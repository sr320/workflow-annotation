#!/usr/bin/env bash
set -euo pipefail

# Record start time for performance tracking
START_TIME=$(date +%s)
START_DATE=$(date '+%Y-%m-%d %H:%M:%S')

# Usage:
#   ./blast2slim.sh -i query.fasta -o outdir [--protein] [--dbdir path] [--threads N] [--diamond]
#
# Input formats supported:
#   - .fasta, .fa, .fas, .cds (nucleotide or protein sequences)
#   - HTTP/HTTPS URLs to FASTA files (auto-downloaded)
#
# Notes:
#   - Default assumes nucleotide queries -> BLASTX to Swiss-Prot (reviewed).
#   - Use --protein if your input contains protein sequences -> BLASTP.
#   - Use --diamond for faster searches (especially for large datasets).
#   - CDS files (.cds) are treated as nucleotide sequences by default.
#   - Requires: BLAST+ (makeblastdb, blastx/blastp) OR DIAMOND, curl, gzip, python3 (+ requests, pandas, goatools).

# ---------- args ----------
QUERY=""
OUTDIR="output"  # base directory now fixed; all run dirs will live under top-level output/
DBDIR="blastdb"
THREADS=40
MODE="blastx"   # or blastp with --protein
USE_DIAMOND=false

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--input) QUERY="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    --dbdir) DBDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --protein) MODE="blastp"; shift 1;;
    --diamond) USE_DIAMOND=true; shift 1;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "${QUERY}" ]]; then
    echo "ERROR: Provide input file with -i (supports .fasta, .fa, .fas, .cds, or URL)"; exit 1
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

# Build appropriate database based on tool choice
if [[ "${USE_DIAMOND}" == "true" ]]; then
  # DIAMOND database
  DIAMOND_DB="${DBDIR}/uniprot_sprot.dmnd"
  if [[ ! -f "${DIAMOND_DB}" ]]; then
    echo "[INFO] Building DIAMOND protein DB..."
    diamond makedb --in "${SP_FASTA}" -d "${DBDIR}/uniprot_sprot"
  fi
else
  # BLAST+ database
  if [[ ! -f "${SP_DB_PREFIX}.pin" ]]; then
    echo "[INFO] Building BLAST protein DB..."
    makeblastdb -in "${SP_FASTA}" -dbtype prot -out "${SP_DB_PREFIX}"
  fi
fi

# ---------- 2) Run BLAST/DIAMOND ----------
BASENAME=$(basename "${QUERY%.*}")
BLAST_TSV="${OUTDIR}/${BASENAME}.blast.tsv"

# Run search with appropriate tool
if [[ "${USE_DIAMOND}" == "true" ]]; then
  DIAMOND_DB="${DBDIR}/uniprot_sprot.dmnd"
  if [[ "${MODE}" == "blastx" ]]; then
    echo "[INFO] Running DIAMOND BLASTX..."
    diamond blastx -d "${DIAMOND_DB}" -q "${QUERY}" -o "${BLAST_TSV}" \
      -e 1e-20 -p "${THREADS}" -k 1 \
      -f 6 qseqid sseqid pident length evalue bitscore stitle
  else
    echo "[INFO] Running DIAMOND BLASTP..."
    diamond blastp -d "${DIAMOND_DB}" -q "${QUERY}" -o "${BLAST_TSV}" \
      -e 1e-20 -p "${THREADS}" -k 1 \
      -f 6 qseqid sseqid pident length evalue bitscore stitle
  fi
else
  # Standard BLAST+
  if [[ "${MODE}" == "blastx" ]]; then
    echo "[INFO] Running BLASTX..."
    blastx -query "${QUERY}" -db "${SP_DB_PREFIX}" -evalue 1e-20 -num_threads "${THREADS}" \
           -max_target_seqs 1 -outfmt "6 qseqid sacc pident length evalue bitscore stitle" > "${BLAST_TSV}"
  else
    echo "[INFO] Running BLASTP..."
    blastp -query "${QUERY}" -db "${SP_DB_PREFIX}" -evalue 1e-20 -num_threads "${THREADS}" \
           -max_target_seqs 1 -outfmt "6 qseqid sacc pident length evalue bitscore stitle" > "${BLAST_TSV}"
  fi
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

# Generate summary report
blast_hits = len(final[final['accession'].notna()])
go_matches = len(final[final['go_ids'].notna() & (final['go_ids'] != '')])
goslim_matches = len(final[final['goslim_ids'].notna() & (final['goslim_ids'] != '')])
total_seqs = len(final)

# GO-Slim category counts for visualization (Biological Process only)
goslim_counts = {}
if goslim_matches > 0:
    for _, row in final[final['goslim_names'].notna() & (final['goslim_names'] != '')].iterrows():
        terms = [x.strip() for x in str(row['goslim_names']).split(';') if x.strip()]
        for term in terms:
            goslim_counts[term] = goslim_counts.get(term, 0) + 1

# Filter GO-Slim counts for Biological Process terms only
bp_goslim_counts = {}
if goslim_counts:
    # Get corresponding GO IDs for filtering by namespace
    for _, row in final[final['goslim_ids'].notna() & (final['goslim_ids'] != '') & 
                       final['goslim_names'].notna() & (final['goslim_names'] != '')].iterrows():
        slim_ids = [x.strip() for x in str(row['goslim_ids']).split(';') if x.strip()]
        slim_names = [x.strip() for x in str(row['goslim_names']).split(';') if x.strip()]
        
        for slim_id, slim_name in zip(slim_ids, slim_names):
            if slim_id in slimdag:
                term_obj = slimdag[slim_id]
                # Filter for biological process namespace
                if hasattr(term_obj, 'namespace') and term_obj.namespace == 'biological_process':
                    bp_goslim_counts[slim_name] = bp_goslim_counts.get(slim_name, 0) + 1

# Use BP counts for visualization, fall back to all counts if no BP terms
chart_counts = bp_goslim_counts if bp_goslim_counts else goslim_counts

# Write summary stats for main script
summary_stats = {
    'total_sequences': total_seqs,
    'blast_hits': blast_hits,
    'go_matches': go_matches,
    'goslim_matches': goslim_matches,
    'goslim_counts': chart_counts  # Use filtered BP counts for chart
}

import json
with open(Path(OUTDIR) / "summary_stats.json", 'w') as f:
    json.dump(summary_stats, f)

print(f"[INFO] Summary stats: {blast_hits}/{total_seqs} BLAST hits, {go_matches} GO annotations, {goslim_matches} GO-Slim mappings")
PYCODE

# Ensure Python deps (use a venv if you like)
python3 - <<'PYCHK'
import sys, subprocess
def ensure(pkg):
    try:
        __import__(pkg)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
for p in ("requests","pandas","goatools","matplotlib"):
    ensure(p)
print("deps ok")
PYCHK

python3 "${PY}" "${BLAST_TSV}" "${OUTDIR}"

# ---------- Generate Summary Report ----------
END_TIME=$(date +%s)
END_DATE=$(date '+%Y-%m-%d %H:%M:%S')
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

# Read summary stats from Python output
if [[ -f "${OUTDIR}/summary_stats.json" ]]; then
    # Extract stats using Python
    python3 - "${OUTDIR}" "$(basename "${QUERY}")" "${START_DATE}" "${END_DATE}" "${HOURS}h ${MINUTES}m ${SECONDS}s" "${THREADS}" "${DURATION}" "${USE_DIAMOND}" "${MODE}" "$(basename "${BLAST_TSV}")" <<'SUMMARY_SCRIPT'
import json
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Get arguments from bash
outdir = sys.argv[1]
input_basename = sys.argv[2]
start_time = sys.argv[3]
end_time = sys.argv[4]
duration_str = sys.argv[5]
cpu_count = int(sys.argv[6])
duration_seconds = int(sys.argv[7])
use_diamond = sys.argv[8]
mode = sys.argv[9]
blast_file = sys.argv[10]

# Load summary statistics
with open(f"{outdir}/summary_stats.json", 'r') as f:
    stats = json.load(f)

# Generate GO-Slim bar chart
goslim_counts = stats.get('goslim_counts', {})
if goslim_counts:
    # Sort by count and take top 15
    sorted_items = sorted(goslim_counts.items(), key=lambda x: x[1], reverse=True)[:15]
    
    if sorted_items:
        terms, counts = zip(*sorted_items)
        
        plt.figure(figsize=(12, 8))
        bars = plt.barh(range(len(terms)), counts, color='steelblue', alpha=0.7)
        plt.yticks(range(len(terms)), [term[:50] + '...' if len(term) > 50 else term for term in terms])
        plt.xlabel('Number of Sequences')
        plt.title('Top GO-Slim Biological Process Categories')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(f"{outdir}/goslim_chart.png", dpi=150, bbox_inches='tight')
        print("[INFO] Generated GO-Slim chart: goslim_chart.png")

# Generate tool info string
if use_diamond == "true":
    tool_prefix = "DIAMOND "
else:
    tool_prefix = "BLAST+ "

if mode == "blastp":
    tool_suffix = "BLASTP (protein)"
else:
    tool_suffix = "BLASTX (nucleotide)"

tool_info = tool_prefix + tool_suffix

md_content = f"""# Annotation Summary Report

## Job Information
- **Input file**: {input_basename}
- **Start time**: {start_time}
- **End time**: {end_time}
- **Duration**: {duration_str}
- **CPUs used**: {cpu_count}
- **Tool**: {tool_info}

## Results Overview
- **Total sequences**: {stats['total_sequences']:,}
- **BLAST hits found**: {stats['blast_hits']:,} ({100*stats['blast_hits']/stats['total_sequences']:.1f}%)
- **GO annotations**: {stats['go_matches']:,} ({100*stats['go_matches']/stats['total_sequences']:.1f}%)
- **GO-Slim mappings**: {stats['goslim_matches']:,} ({100*stats['goslim_matches']/stats['total_sequences']:.1f}%)

## Output Files
- **Main results**: annotation_with_goslim.tsv
- **Full GO data**: annotation_full_go.tsv  
- **Raw BLAST**: {blast_file}
- **Processing script**: postprocess_uniprot_go.py

"""

if goslim_counts:
    md_content += """## Top GO-Slim Categories

| GO-Slim Term | Count |
|--------------|-------|
"""
    for term, count in sorted(goslim_counts.items(), key=lambda x: x[1], reverse=True)[:15]:
        md_content += f"| {term} | {count} |\n"
    
    md_content += f"""
![GO-Slim Categories](goslim_chart.png)

*Top 15 GO-Slim categories by sequence count*
"""
else:
    md_content += "## GO-Slim Categories\n\nNo GO-Slim mappings found.\n"

md_content += f"""
## Performance
- **BLAST throughput**: {stats['total_sequences']/max(1, duration_seconds):.1f} sequences/second
- **Annotation rate**: {stats['blast_hits']/max(1, duration_seconds):.1f} hits/second

---
*Generated by blast2slim.sh on {end_time}*
"""

with open(f"{outdir}/summary.md", 'w') as f:
    f.write(md_content)

print(f"[INFO] Generated summary report: summary.md")
SUMMARY_SCRIPT
fi

echo "[OK] All done. See ${OUTDIR}/annotation_with_goslim.tsv and ${OUTDIR}/summary.md"