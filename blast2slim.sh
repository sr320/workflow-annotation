#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./blast2goslim.sh -i query.fasta -o outdir [--protein] [--dbdir path] [--threads N]
#
# Notes:
#   - Default assumes nucleotide queries -> BLASTX to Swiss-Prot (reviewed).
#   - Use --protein if your FASTA is protein -> BLASTP.
#   - Requires: BLAST+ (makeblastdb, blastx/blastp), curl, gzip, python3 (+ requests, pandas, goatools).

# ---------- args ----------
QUERY=""
OUTDIR="output"  # base directory now fixed; all run dirs will live under top-level output/
DBDIR="blastdb"
THREADS=8
MODE="blastx"   # or blastp with --protein

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input) QUERY="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    --dbdir) DBDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --protein) MODE="blastp"; shift 1;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "${QUERY}" ]]; then
        echo "ERROR: Provide -i query.fasta"; exit 1
fi

# Normalize user-supplied OUTDIR so it always resides under top-level output/
# Examples:
#   (default) output -> output/run_<ts>
#   -o myproj        -> output/myproj/run_<ts>
#   -o output/x      -> output/x/run_<ts>
if [[ "${OUTDIR}" == /* ]]; then
    # Convert absolute path to just its basename components under output/
    # (keeps last path element so user still distinguishes runs) 
    OUTDIR="output/$(basename "${OUTDIR}")"
else
    OUTDIR="${OUTDIR#./}"          # strip leading ./
    OUTDIR="${OUTDIR%/}"           # trim trailing /
    # If user did not start with output/, prefix it
    if [[ "${OUTDIR}" != output* ]]; then
        OUTDIR="output/${OUTDIR}"
    fi
fi

# Create per-run timestamped directory inside normalized OUTDIR
BASE_OUTDIR="${OUTDIR%/}"
TS=$(date +%Y%m%d_%H%M%S)
RUN_DIR="${BASE_OUTDIR}/run_${TS}"
if [[ -d "${RUN_DIR}" ]]; then
    i=1
    while [[ -d "${RUN_DIR}" ]]; do
        RUN_DIR="${BASE_OUTDIR}/run_${TS}_${i}"; ((i++))
    done
fi
OUTDIR="${RUN_DIR}"
echo "[INFO] Output directory: ${OUTDIR}"
mkdir -p "${OUTDIR}" "${DBDIR}"

# If input is an HTTP/HTTPS URL, download it locally into OUTDIR
if [[ "${QUERY}" =~ ^https?:// ]]; then
    echo "[INFO] Downloading input FASTA from URL: ${QUERY}"
    fname=$(basename "${QUERY}")
    # fallback name if URL ends with /
    if [[ -z "${fname}" || "${fname}" == */ ]]; then
        fname="query.fasta"
    fi
    DL_PATH="${OUTDIR}/${fname}"
    curl -L -o "${DL_PATH}" "${QUERY}"
    QUERY="${DL_PATH}"
    echo "[INFO] Saved URL FASTA to ${QUERY}"
fi

# ---------- 1) Download latest Swiss-Prot & build/refresh BLAST DB ----------
# (kept simple; re-download if missing; you can pin a release by renaming file)
SP_FASTA_GZ="${DBDIR}/uniprot_sprot.fasta.gz"
SP_FASTA="${DBDIR}/uniprot_sprot.fasta"
SP_DB_PREFIX="${DBDIR}/uniprot_sprot"

if [[ ! -f "${SP_FASTA_GZ}" && ! -f "${SP_FASTA}" ]]; then
  echo "[INFO] Downloading current Swiss-Prot..."
  curl -L -o "${SP_FASTA_GZ}" "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
fi

if [[ ! -f "${SP_FASTA}" ]]; then
  echo "[INFO] Unzipping Swiss-Prot..."
  gunzip -kc "${SP_FASTA_GZ}" > "${SP_FASTA}"
fi

if [[ ! -f "${SP_DB_PREFIX}.pin" ]]; then
  echo "[INFO] Building BLAST protein DB..."
  makeblastdb -in "${SP_FASTA}" -dbtype prot -out "${SP_DB_PREFIX}"
fi

# ---------- 2) Run BLAST ----------
BASENAME=$(basename "${QUERY%.*}")
BLAST_TSV="${OUTDIR}/${BASENAME}.blast.tsv"

# outfmt: qseqid sacc pident length evalue bitscore stitle
if [[ "${MODE}" == "blastx" ]]; then
  echo "[INFO] Running BLASTX..."
  blastx -query "${QUERY}" -db "${SP_DB_PREFIX}" -evalue 1e-20 -num_threads "${THREADS}" \
         -max_target_seqs 1 -outfmt "6 qseqid sacc pident length evalue bitscore stitle" > "${BLAST_TSV}"
else
  echo "[INFO] Running BLASTP..."
  blastp -query "${QUERY}" -db "${SP_DB_PREFIX}" -evalue 1e-20 -num_threads "${THREADS}" \
         -max_target_seqs 1 -outfmt "6 qseqid sacc pident length evalue bitscore stitle" > "${BLAST_TSV}"
fi

# ---------- 3) Python post-processing: UniProt GO + GO-Slim ----------
PY="${OUTDIR}/postprocess_uniprot_go.py"
cat > "${PY}" << 'PYCODE'
import sys, os, time, json, math
import pandas as pd
import requests
from pathlib import Path

BLAST_TSV = sys.argv[1]
OUTDIR = sys.argv[2]

os.makedirs(OUTDIR, exist_ok=True)
hits = pd.read_csv(BLAST_TSV, sep='\t', header=None,
                   names=['query','accession','pident','length','evalue','bitscore','title'])

# Normalize accession (handle formats like "sp|P69907|HBA_PANTR") for UniProt REST API
def _norm_acc(a: str):
    if not isinstance(a, str):
        return a
    if '|' in a:
        parts = a.split('|')
        if len(parts) >= 3 and parts[0] in {'sp','tr'}:
            return parts[1]
    return a
hits['accession'] = hits['accession'].apply(_norm_acc)
hits['query'] = hits['query'].astype(str)

# Deduplicate accessions for API efficiency
accs = sorted(hits['accession'].dropna().unique().tolist())
# UniProt REST batch: use query with OR list; fields: accession, id, protein_name, organism_name, go_id, go_p, go_c, go_f
# See UniProt REST docs. We will page in chunks to be safe.
# https://www.uniprot.org/help/api_queries ; https://www.uniprot.org/help/return_fields
def fetch_uniprot_batch(acc_batch):
    """Fetch UniProt annotations for a list of accessions, splitting recursively if query causes 400."""
    if not acc_batch:
        return pd.DataFrame()
    # Base case: single accession -> always try directly
    if len(acc_batch) == 1:
        query = f"accession:{acc_batch[0]}"
    else:
        query = " OR ".join([f"accession:{a}" for a in acc_batch])
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,id,reviewed,protein_name,organism_name,go_id,go_p,go_c,go_f"
    }
    try:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
    except requests.HTTPError as e:
        # If too large / bad request and we have multiple accs, split and recurse
        if r.status_code == 400 and len(acc_batch) > 1:
            mid = len(acc_batch)//2
            left = fetch_uniprot_batch(acc_batch[:mid])
            right = fetch_uniprot_batch(acc_batch[mid:])
            return pd.concat([left, right], ignore_index=True)
        raise
    if not r.text.strip():
        return pd.DataFrame(columns=["accession","id","reviewed","protein_name","organism","go_ids","go_bp","go_cc","go_mf"])
    from io import StringIO
    df = pd.read_csv(StringIO(r.text), sep='\t')
    # normalize column names
    rename_map = {
        "Entry":"accession",
        "Entry Name":"id",
        "Status":"reviewed",
        "Protein names":"protein_name",
        "Organism":"organism",
        "Gene Ontology IDs":"go_ids",
        "Gene Ontology (biological process)":"go_bp",
        "Gene Ontology (cellular component)":"go_cc",
        "Gene Ontology (molecular function)":"go_mf"
    }
    df = df.rename(columns={k:v for k,v in rename_map.items() if k in df.columns})
    return df

frames = []
CHUNK = 50  # smaller initial chunk to avoid long URLs
for i in range(0, len(accs), CHUNK):
    batch = accs[i:i+CHUNK]
    frames.append(fetch_uniprot_batch(batch))
    time.sleep(0.1)

if frames:
    u = pd.concat(frames, ignore_index=True)
else:
    # Ensure merge does not fail if there are no hits
    u = pd.DataFrame(columns=['accession','id','reviewed','protein_name','organism','go_ids','go_bp','go_cc','go_mf'])
# Join back to BLAST best hits
merged = hits.merge(u, how='left', on='accession')

# Write intermediate with full GO
full_path = Path(OUTDIR) / "annotation_full_go.tsv"
merged.to_csv(full_path, sep='\t', index=False)

# ---- Map full GO -> GO-Slim (generic) using GOATOOLS ----
# inputs: go-basic.obo + goslim_generic.obo
from goatools.obo_parser import GODag
from collections import defaultdict

obo_basic = Path(OUTDIR) / "go-basic.obo"
obo_slim  = Path(OUTDIR) / "goslim_generic.obo"

# download if missing
import urllib.request
if not obo_basic.exists():
    urllib.request.urlretrieve("http://purl.obolibrary.org/obo/go/go-basic.obo", obo_basic.as_posix())
if not obo_slim.exists():
    urllib.request.urlretrieve("https://go.princeton.edu/GOTermMapper/goSlimFiles/goslim_generic.obo", obo_slim.as_posix())

obodag = GODag(obo_basic.as_posix(), optional_attrs={'relationship'})
slimdag = GODag(obo_slim.as_posix(), optional_attrs={'relationship'})

# Build query: list of GO terms per gene (merge bp/cc/mf) as GO IDs
def split_goids(s):
    if isinstance(s, float) and math.isnan(s): return []
    # UniProt 'go_ids' is semi-colon separated IDs (e.g., "GO:0005524; GO:0004672")
    return [x.strip() for x in str(s).split(';') if x.strip().startswith('GO:')]

gene2gos = defaultdict(set)
for _, row in merged.iterrows():
    gs = set(split_goids(row.get('go_ids', '')))
    if gs:
        gene2gos[row['query']].update(gs)

# Map to slim: returns counts per slim term; we also want per-gene mapping
# We'll perform per-gene slim mapping using GOATOOLS internal search
def _ancestors(goid, dag):
    term = dag.get(goid)
    if term is None:
        return []
    seen = set()
    stack = [term]
    while stack:
        t = stack.pop()
        if t.id in seen:
            continue
        seen.add(t.id)
        for p in t.parents:
            stack.append(p)
    return seen

slim_set = {t.id for t in slimdag.values()}
slim_map_rows = []
for gene, gos in gene2gos.items():
    slim_hits = set()
    for go in gos:
        for anc in _ancestors(go, obodag):
            if anc in slim_set:
                slim_hits.add(anc)
    slim_map_rows.append({
        "query": gene,
        "goslim_ids": "; ".join(sorted(slim_hits))
    })
if slim_map_rows:
    slim_df = pd.DataFrame(slim_map_rows)
    # Add names for slim terms
    def name_of(goid):
        n = slimdag.get(goid)
        return n.name if n is not None else ""
    slim_df['goslim_names'] = slim_df['goslim_ids'].apply(
        lambda s: "; ".join([name_of(x) for x in s.split('; ')]) if s else ""
    )
else:
    # Ensure required columns exist for merge
    slim_df = pd.DataFrame(columns=["query","goslim_ids","goslim_names"])
    slim_df = slim_df.astype({"query":"string","goslim_ids":"string","goslim_names":"string"})

final = merged.merge(slim_df, on='query', how='left')
final['query'] = final['query'].astype(str)

# Tidy final columns
desired_cols = [
    "query","accession","id","reviewed","protein_name","organism",
    "pident","length","evalue","bitscore","title",
    "go_ids","go_bp","go_cc","go_mf",
    "goslim_ids","goslim_names"
]
for c in desired_cols:
    if c not in final.columns:
        final[c] = pd.NA
final = final[desired_cols]

final.to_csv(Path(OUTDIR) / "annotation_with_goslim.tsv", sep='\t', index=False)
print(f"[DONE] Wrote:\n - {full_path}\n - {Path(OUTDIR) / 'annotation_with_goslim.tsv'}")
PYCODE

# Ensure Python deps (use a venv if you like)
python3 - <<'PYCHK'
import sys, subprocess
def ensure(pkg):
    try:
        __import__(pkg)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
for p in ("requests","pandas","goatools"):
    ensure(p)
print("deps ok")
PYCHK

python3 "${PY}" "${BLAST_TSV}" "${OUTDIR}"

echo "[OK] All done. See ${OUTDIR}/annotation_with_goslim.tsv"