# blast2slim.sh - Swiss-Prot Annotation Pipeline

A comprehensive bash script that annotates nucleotide or protein FASTA sequences with Swiss-Prot best hits, UniProt protein names, full Gene Ontology (GO) annotations, and GO-Slim terms.

## Script Overview

The `blast2slim.sh` script is a self-contained pipeline that:   
1. **Downloads & indexes Swiss-Prot database** (one-time setup)   
2. **Runs BLAST searches** (BLASTX for nucleotides, BLASTP for proteins)   
3. **Retrieves UniProt annotations** via REST API   
4. **Maps GO terms to GO-Slim** categories using GOATOOLS  
5. **Outputs structured TSV files** ready for downstream analysis  

## Command Syntax

```bash
bash blast2slim.sh -i <input.fasta> [OPTIONS]
```

## Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `-i, --input FILE` | Input FASTA file (supports .fasta, .fa, .fas, .cds extensions) | `-i transcripts.cds` |

## Optional Arguments

| Argument | Default | Description | Example |
|----------|---------|-------------|---------|
| `-o, --outdir DIR` | `output` | Output subdirectory name (always under `output/`) | `-o my_project` |
| `--dbdir DIR` | `blastdb` | Directory for Swiss-Prot database files | `--dbdir /shared/db` |
| `--threads N` | `40` | Number of CPU threads for BLAST/DIAMOND | `--threads 16` |
| `--protein` | (none) | Use BLASTP mode for protein input (default: BLASTX) | `--protein` |
| `--diamond` | (none) | Use DIAMOND instead of BLAST+ (faster for large datasets) | `--diamond` |

## Usage Examples

### Basic Examples

#### Nucleotide Sequences (Default: BLASTX)
```bash
# FASTA format
bash blast2slim.sh -i transcripts.fasta

# CDS format (coding sequences)
bash blast2slim.sh -i genes.cds -o cds_analysis --threads 16

# Other supported formats
bash blast2slim.sh -i sequences.fa
bash blast2slim.sh -i transcripts.fas
```

#### Protein Sequences (BLASTP Mode)
```bash
# Protein sequences require --protein flag
bash blast2slim.sh -i proteins.faa --protein

# With custom output directory
bash blast2slim.sh -i proteins.faa --protein -o protein_annotation --threads 12
```

#### Fast Searches with DIAMOND
```bash
# DIAMOND BLASTX (much faster for large nucleotide datasets)
bash blast2slim.sh -i large_transcriptome.fasta --diamond --threads 40

# DIAMOND BLASTP (faster protein searches)
bash blast2slim.sh -i proteins.faa --protein --diamond -o fast_protein_search

# Combined with other options
bash blast2slim.sh -i genes.cds --diamond -o diamond_analysis --threads 64
```

#### Direct URL Input
```bash
# Download and annotate FASTA directly from URL
bash blast2slim.sh -i "https://example.org/sequences.fasta" --protein -o url_analysis

# Real example with public data
bash blast2slim.sh -i "https://raw.githubusercontent.com/urol-e5/timeseries_molecular/refs/heads/main/D-Apul/data/Machinery.fasta" --protein -o machinery_proteins
```

### Advanced Examples

#### Custom Database Location
```bash
# Use shared database directory to avoid re-downloading
bash blast2slim.sh -i sequences.fasta --dbdir /shared/blastdb -o shared_db_run

# Useful for multiple projects using same reference
bash blast2slim.sh -i project1.fasta --dbdir /data/uniprot_db -o project1
bash blast2slim.sh -i project2.fasta --dbdir /data/uniprot_db -o project2
```

#### High-Performance Settings
```bash
# Maximum threads for large datasets
bash blast2slim.sh -i large_transcriptome.fasta --threads $(nproc) -o large_analysis

# Monitor progress with logging
bash blast2slim.sh -i sequences.fasta -o logged_run --threads 16 2>&1 | tee analysis.log
```

#### Batch Processing
```bash
# Process multiple FASTA files
for fasta in data/*.fasta; do
    basename=$(basename "$fasta" .fasta)
    bash blast2slim.sh -i "$fasta" -o "batch_${basename}" --threads 4
done
```

## Detailed Parameter Descriptions

### Input Options

**`-i, --input FILE`** (REQUIRED)
- Accepts local FASTA files (.fasta, .fa, .fas, .cds) or HTTP/HTTPS URLs
- For URLs, file is automatically downloaded to output directory
- Supports both nucleotide and protein sequences
- CDS files (.cds) are treated as nucleotide sequences by default
- Sequence headers: text before first space becomes query ID

**`--protein`** (FLAG)
- Changes from BLASTX (nucleotide→protein) to BLASTP (protein→protein)
- Required when input contains protein sequences
- Omit this flag for nucleotide sequences (transcripts, genes, etc.)

**`--diamond`** (FLAG)
- Use DIAMOND instead of BLAST+ for much faster searches
- Recommended for large datasets (>10,000 sequences)
- ~10-20x faster than BLAST+ with similar sensitivity
- Creates `.dmnd` database files instead of BLAST+ `.pin/.psq` files
- Compatible with both `--protein` and nucleotide modes

### Output Options

**`-o, --outdir DIR`** (Default: `output`)
- Creates timestamped run directory: `output/{DIR}/run_YYYYMMDD_HHMMSS/`
- Examples:
  - `-o myproject` → `output/myproject/run_20250826_143022/`
  - Default → `output/run_20250826_143022/`
- All files organized under this directory structure

**`--dbdir DIR`** (Default: `blastdb`)
- Location for Swiss-Prot database files (~90MB compressed, ~260MB uncompressed)
- Reusable across multiple runs
- Files created: `uniprot_sprot.fasta.gz`, `uniprot_sprot.fasta`, `uniprot_sprot.p*`

### Performance Options

**`--threads N`** (Default: `40`)
- Number of CPU threads for BLAST search
- Higher values = faster processing (up to available cores)
- Recommended: match to available CPU cores (`$(nproc)`)
- Memory usage scales with thread count

## Fixed Internal Parameters

The script uses optimized BLAST parameters:

| Parameter | Value | Purpose |
|-----------|-------|---------|
| E-value threshold | `1e-20` | Stringent significance cutoff |
| Max target sequences | `1` | Only best hit per query |
| Output format | Tabular | `qseqid sacc pident length evalue bitscore stitle` |
| Database | Swiss-Prot | Reviewed UniProt entries only |

## Output Structure

```
output/
└── {project_name}/
    └── run_YYYYMMDD_HHMMSS/
        ├── {basename}.blast.tsv           # Raw BLAST results
        ├── annotation_full_go.tsv         # BLAST + UniProt + full GO
        ├── annotation_with_goslim.tsv     # Main output with GO-Slim
        ├── summary.md                     # Concise summary report with stats & visualization
        ├── goslim_chart.png              # GO-Slim category bar chart
        ├── postprocess_uniprot_go.py      # Generated processing script
        ├── go-basic.obo                   # GO ontology
        └── goslim_generic.obo            # GO-Slim ontology
```

### Output File Descriptions

**`summary.md`** - Concise summary report
- Job metadata (duration, CPUs used, input file)
- Key statistics (BLAST hits, GO matches, GO-Slim mappings)
- Top GO-Slim categories table and visualization
- Performance metrics

**`annotation_with_goslim.tsv`** - Primary results file
- BLAST hits with UniProt annotations and GO-Slim terms
- Columns: query, accession, protein_name, organism, BLAST metrics, GO terms, GO-Slim categories

**`annotation_full_go.tsv`** - Intermediate results
- Same as above but with full GO terms instead of GO-Slim

**`{basename}.blast.tsv`** - Raw BLAST output
- Tabular BLAST results before UniProt annotation
- Useful for troubleshooting or custom post-processing

## Dependencies

### Required Tools
- **BLAST+**: `makeblastdb`, `blastx`, `blastp` (NCBI BLAST+ suite) OR
- **DIAMOND**: `diamond` (faster alternative for large datasets)
- **Python 3** (≥3.8) with packages:
  - `requests` (UniProt REST API)
  - `pandas` (data processing)
  - `goatools` (GO term mapping)
  - `matplotlib` (summary visualization)
- **System tools**: `bash`, `curl`, `gzip`

### Installation Examples

**Via Conda (Recommended)**
```bash
conda create -n annotation python=3.11 -y
conda activate annotation
# Install BLAST+ (default option)
conda install -c bioconda blast -y
# OR install DIAMOND (faster option)
conda install -c bioconda diamond -y
# Install Python packages
pip install requests pandas goatools matplotlib
```

**Via Package Managers**
```bash
# Ubuntu/Debian
sudo apt install ncbi-blast+ python3-pip
pip3 install requests pandas goatools

# macOS with Homebrew
brew install blast python3
pip3 install requests pandas goatools
```

## Workflow Details

### Step 1: Database Setup (First Run Only)
1. Downloads Swiss-Prot FASTA from UniProt FTP
2. Uncompresses and builds BLAST protein database
3. Database files reused for subsequent runs

### Step 2: BLAST Search
- **BLASTX**: Nucleotide queries → Swiss-Prot proteins (6-frame translation)
- **BLASTP**: Protein queries → Swiss-Prot proteins (direct comparison)
- Parameters: E-value ≤1e-20, best hit only

### Step 3: UniProt Annotation
- Batch queries UniProt REST API for accessions found
- Retrieves: protein names, organism, GO terms (BP/CC/MF)
- Handles large result sets with automatic query splitting

### Step 4: GO-Slim Mapping
- Downloads current GO ontology (go-basic.obo)
- Downloads GO-Slim generic subset (goslim_generic.obo)
- Maps full GO terms to high-level GO-Slim categories using GOATOOLS

## Performance Characteristics

### BLAST+ vs DIAMOND Comparison

| Metric | BLAST+ | DIAMOND | Notes |
|--------|--------|---------|-------|
| **Speed** | Baseline | 10-20x faster | DIAMOND optimized for large datasets |
| **Sensitivity** | High | Very similar | Minimal difference in annotation quality |
| **Memory** | 1-4GB | 2-6GB | DIAMOND uses slightly more RAM |
| **Database size** | ~350MB | ~120MB | DIAMOND databases are smaller |
| **Best for** | <10K sequences | >10K sequences | Recommendation by dataset size |

### Timing Estimates

**BLAST+ Mode:**
- **Database setup**: 5-15 minutes (first run only)
- **Search speed**: ~0.1-1 second per query sequence
- **UniProt annotation**: ~0.1 seconds per unique hit
- **GO mapping**: ~1-10 seconds total

**DIAMOND Mode:**
- **Database setup**: 2-5 minutes (first run only)
- **Search speed**: ~0.005-0.05 seconds per query sequence
- **UniProt annotation**: ~0.1 seconds per unique hit
- **GO mapping**: ~1-10 seconds total

### Resource Usage
- **Disk space**: ~350MB (BLAST+) or ~120MB (DIAMOND) for database files
- **Memory**: ~1-4GB (BLAST+) or ~2-6GB (DIAMOND) depending on dataset size
- **Network**: ~90MB download (first run), API calls during annotation

### Scalability
- **Small datasets** (1-1000 seqs): 1-10 minutes total (either tool)
- **Medium datasets** (1000-10000 seqs): 10-60 minutes (BLAST+), 5-15 minutes (DIAMOND)
- **Large datasets** (10000+ seqs): 1+ hours (BLAST+), 10-30 minutes (DIAMOND)

## Troubleshooting

### Common Issues

**"command not found: blastx" or "command not found: diamond"**
```bash
# Install BLAST+
conda install -c bioconda blast
# OR install DIAMOND (faster)
conda install -c bioconda diamond
# OR via package manager
sudo apt install ncbi-blast+ diamond-aligner
```

**"ModuleNotFoundError: No module named 'requests'"**
```bash
pip install requests pandas goatools matplotlib
```

**"Permission denied"**
```bash
chmod +x blast2slim.sh
```

**"No space left on device"**
- Check disk space: `df -h`
- Use custom `--dbdir` on larger partition
- Clean up old output directories

**Long UniProt API delays**
- Network-dependent; just wait
- Large numbers of unique hits take longer
- Check firewall/proxy settings if timing out

### Monitoring Long Runs

```bash
# Run with logging
bash blast2slim.sh -i large_file.fasta -o big_run --threads 16 2>&1 | tee run.log

# Monitor progress in another terminal
tail -f run.log

# Check BLAST progress
ps aux | grep blast
```

### Validation

```bash
# Check output file completeness
wc -l output/*/run_*/annotation_with_goslim.tsv
head -1 output/*/run_*/annotation_with_goslim.tsv

# Verify database files
ls -lh blastdb/
```

## Integration Examples

### With R/Bioconductor
```r
library(readr)
results <- read_tsv("output/myproject/run_*/annotation_with_goslim.tsv")
```

### With Python/Pandas
```python
import pandas as pd
results = pd.read_csv("output/myproject/run_*/annotation_with_goslim.tsv", sep='\t')
```

### Pipeline Integration
```bash
# Part of larger workflow
bash preprocess_sequences.sh input.fasta > cleaned.fasta
bash blast2slim.sh -i cleaned.fasta -o analysis_$(date +%Y%m%d)
Rscript analyze_results.R output/analysis_*/run_*/annotation_with_goslim.tsv
```

## Citation

When using this workflow in publications:

> "Sequences were annotated using Swiss-Prot BLAST best hits (E-value ≤1e-20) with UniProt GO term retrieval and GO-Slim mapping via GOATOOLS."

## License

This script is provided as-is for research and educational use.