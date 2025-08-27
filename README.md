# workflow-annotation

A streamlined bioinformatics workflow to annotate nucleotide or protein FASTA sequences with Swiss-Prot best hits, UniProt protein names, full Gene Ontology (GO) annotations, and GO-Slim terms.

## Overview

The `blast2slim.sh` script performs:
1. Downloads Swiss-Prot database and builds BLAST index (first run only)
2. Runs BLASTX (nucleotide) or BLASTP (protein) for best hits
3. Retrieves UniProt annotations via REST API
4. Maps full GO terms to GO-Slim categories
5. Outputs clean TSV files ready for analysis

## Using in Quarto Documents

This workflow is designed to integrate seamlessly with Quarto documents. Below are the key steps to get started.

### Clone the Repository

```{bash}
#| eval: false
git clone https://github.com/sr320/workflow-annotation.git
cd workflow-annotation
```

### Check Requirements

Ensure you have the required tools installed:

```{bash}
#| eval: false
# Check BLAST+ installation
which makeblastdb blastx blastp

# Check Python and install required packages
python3 -c "import requests, pandas; print('Python packages OK')" || \
python3 -m pip install requests pandas goatools
```

### Basic Usage Examples

#### Annotate Nucleotide Sequences (BLASTX)

```{bash}
#| eval: false
# Download example transcriptome
wget -O transcripts.fasta "https://example.org/transcripts.fasta"

# Run annotation workflow
bash blast2slim.sh -i transcripts.fasta -o my_transcripts --threads 8
```

#### Annotate Protein Sequences (BLASTP)

```{bash}
#| eval: false
# Download example protein set
wget -O proteins.faa "https://raw.githubusercontent.com/urol-e5/timeseries_molecular/refs/heads/main/D-Apul/data/Machinery.fasta"

# Run protein annotation
bash blast2slim.sh -i proteins.faa --protein -o protein_annotation --threads 8
```

#### Annotate from URL (Direct Download)

```{bash}
#| eval: false
# Annotate directly from a URL
bash blast2slim.sh -i "https://example.org/sequences.fasta" --protein -o url_analysis
```

### Output Structure

All outputs are organized under the `output/` directory:

```{bash}
#| eval: false
# List output structure
ls -la output/
ls -la output/*/run_*/
```

### Key Output Files

```{bash}
#| eval: false
# View the main annotation results
head output/my_project/run_*/annotation_with_goslim.tsv

# Check raw BLAST results
head output/my_project/run_*/yourfile.blast.tsv
```

### Load Results in R

```{r}
#| eval: false
library(readr)
library(dplyr)

# Read annotation results
annotations <- read_tsv("output/my_project/run_20250826_123456/annotation_with_goslim.tsv")

# Quick summary
annotations %>%
  count(organism, sort = TRUE) %>%
  head(10)

# GO-Slim category counts
annotations %>%
  separate_rows(goslim_names, sep = "; ") %>%
  count(goslim_names, sort = TRUE) %>%
  filter(!is.na(goslim_names), goslim_names != "") %>%
  head(15)
```

### Performance Tips

```{bash}
#| eval: false
# Use more threads for faster processing
bash blast2slim.sh -i sequences.fasta --threads $(nproc)

# For large datasets, monitor progress
bash blast2slim.sh -i large_dataset.fasta -o big_analysis --threads 16 2>&1 | tee analysis.log
```

### Custom Database Location

```{bash}
#| eval: false
# Use custom database directory (shared across projects)
bash blast2slim.sh -i sequences.fasta --dbdir /shared/blastdb -o custom_analysis
```

## Output Files Explained

| File | Description |
|------|-------------|
| `annotation_with_goslim.tsv` | **Main result**: BLAST hits + UniProt info + GO-Slim terms |
| `annotation_full_go.tsv` | BLAST hits + UniProt info + full GO terms |
| `sequences.blast.tsv` | Raw BLAST tabular output |
| `postprocess_uniprot_go.py` | Generated Python script (for reproducibility) |
| `go-basic.obo`, `goslim_generic.obo` | Ontology files |

## Key Columns in Results

- **query**: Your original sequence ID
- **accession**: UniProt accession from best BLAST hit
- **protein_name**: Protein description
- **organism**: Source organism
- **pident**, **evalue**, **bitscore**: BLAST alignment metrics
- **go_ids**: Full GO term IDs (semicolon-separated)
- **goslim_ids**, **goslim_names**: Mapped GO-Slim categories

## Requirements

- **BLAST+**: `makeblastdb`, `blastx`, `blastp`
- **Python 3** (≥3.8) with packages: `requests`, `pandas`, `goatools`
- **Standard tools**: `bash`, `curl`, `gzip`

### Install with Conda

```{bash}
#| eval: false
conda create -n annotation python=3.11 -y
conda activate annotation
conda install -c bioconda blast -y
pip install requests pandas goatools
```

## Troubleshooting

### Common Issues

```{bash}
#| eval: false
# BLAST not found
conda install -c bioconda blast

# Python packages missing
python3 -m pip install requests pandas goatools

# Check script permissions
chmod +x blast2slim.sh

# Monitor large runs
tail -f analysis.log
```

## Integration with Analysis Pipelines

### In Quarto Reports

````markdown
```{bash}
#| label: run-annotation
#| eval: true
#| output: false
bash blast2slim.sh -i data/transcripts.fasta -o transcriptome_analysis --threads 8
```

```{r}
#| label: load-results
library(readr)
results <- read_tsv("output/transcriptome_analysis/run_*/annotation_with_goslim.tsv")
```
````

### Batch Processing

```{bash}
#| eval: false
# Process multiple files
for fasta in data/*.fasta; do
    basename=$(basename "$fasta" .fasta)
    bash blast2slim.sh -i "$fasta" -o "batch_${basename}" --threads 4
done
```

## License

This workflow is provided as-is for research and educational use.

## Citation

When using this workflow, cite: "Sequences were annotated using Swiss-Prot BLAST best hits with UniProt GO term retrieval and GO-Slim mapping via GOATOOLS."

---

**Quick Start**: Clone repo → Run `bash blast2slim.sh -i your_sequences.fasta` → Analyze `output/*/annotation_with_goslim.tsv`