# RNA-seq Snakemake Pipeline — Chromosome 19 | Alzheimer's Disease

Reproducible end-to-end RNA-seq workflow built with Snakemake, analyzing differential gene expression on Chromosome 19 in an Alzheimer's Disease mouse model (Mus musculus GRCm39). Designed to run on low-resource environments without HPC access.

---

## Pipeline Overview

```
Raw FASTQ → FastQC QC → HISAT2 Alignment → SAMtools BAM processing
→ featureCounts Read Counting → DESeq2 Differential Expression
→ clusterProfiler Pathway Analysis → MultiQC Report
```

---

## Tools & Technologies

| Tool | Purpose |
|---|---|
| Snakemake | Workflow management and reproducibility |
| FastQC | Per-sample read quality control |
| HISAT2 | Splice-aware RNA-seq alignment |
| SAMtools | BAM sorting and indexing |
| featureCounts (subread) | Gene-level read counting |
| DESeq2 (R/Bioconductor) | Differential expression analysis |
| clusterProfiler (R) | GO and KEGG pathway enrichment |
| MultiQC | Aggregated QC report across all samples |

---

## Repository Structure

```
rnaseq-chr19-snakemake/
├── Snakefile               # Snakemake workflow — all pipeline rules
├── config.yaml             # Parameters: reference, annotation, threads
├── deseq2.R                # DESeq2 differential expression script
├── pathway_analysis.R      # GO/KEGG pathway enrichment script
├── raw_counts.txt          # featureCounts output — input to DESeq2
├── .gitignore              # Excludes BAMs, FASTQs, index files
└── README.md               # This file
```

---

## How to Run

**1. Clone the repository**
```bash
git clone https://github.com/tejasw2204/rnaseq-chr19-snakemake
cd rnaseq-chr19-snakemake
```

**2. Install dependencies**
```bash
conda install -c conda-forge -c bioconda snakemake hisat2 samtools subread fastqc multiqc -y
```

**3. Add your FASTQ files**

Place paired-end FASTQ files in the working directory following this naming convention:
```
samplename_1.fastq
samplename_2.fastq
```
Sample names are detected **automatically** — no hardcoding needed.

**4. Update config.yaml**
```yaml
reference_index: "chr19_index"   # path to your HISAT2 index prefix
annotation:      "chr19.gtf"     # path to your GTF annotation file
threads:         4
```

**5. Dry run first (always)**
```bash
snakemake --dry-run
```

**6. Run the pipeline**
```bash
snakemake --cores 4
```

**7. Run DESeq2 and pathway analysis**
```bash
Rscript deseq2.R
Rscript pathway_analysis.R
```

---

## Output Files

| File | Description |
|---|---|
| `raw_counts.txt` | Read count matrix across all samples |
| `multiqc_report/multiqc_report.html` | Aggregated QC report |
| `results/deseq2_results.csv` | Differentially expressed genes |
| `results/pathway_results.csv` | Enriched GO/KEGG pathways |

---

## Key Features

- **Auto-detects samples** from FASTQ files in the working directory — works on any dataset without code changes
- **Skips completed steps** — reruns only what changed, saving compute time
- **Low-resource optimised** — tested on 4GB RAM without HPC access
- **Fully parameterised** via `config.yaml` — no hardcoded paths in pipeline code

---

## Data Availability

Raw FASTQ files, BAM alignments, and genome index files are not included in this repository due to file size. Reference genome: *Mus musculus* GRCm39, Chromosome 19 (`chr19.fa`). Annotation: Ensembl GRCm39 release 110 (`Mus_musculus.GRCm39.110.gtf`).

The `raw_counts.txt` count matrix is included and is sufficient to reproduce all downstream DESeq2 and pathway analysis results directly.

---

## Author

**Tejas Manoj Wankhade**  
M.Sc. Bioinformatics — Guru Nanak Khalsa College, Mumbai  
[LinkedIn](https://www.linkedin.com/in/tejas-wankhade-642699246) · [GitHub](https://github.com/tejasw2204)
