# AutoFlow-SNAI1

End-to-end automated RNA-Seq and drug discovery pipeline for transcriptomic analysis of SNAI1 Knockout models in Breast Cancer. Transitions from raw differential expression to clinical-grade drug repurposing candidates in a single command.

---

## Overview

AutoFlow-SNAI1 is a fully automated, reproducible bioinformatics workflow built on Snakemake. It processes raw count matrices through differential expression, strict biomarker filtration, functional pathway enrichment, and automated drug repurposing — querying the ChEMBL API to surface Phase 3 and Phase 4 clinical inhibitors for overexpressed gene targets.

---

## Pipeline

01 — Differential Expression
Processes raw count matrices using DESeq2 to identify statistically significant gene expression changes.

02 — Strict Filtration
Isolates high-confidence biomarkers at padj < 0.05 and |log2FC| >= 2.0 thresholds.

03 — Pathway Enrichment
Functional annotation via GSEApy against KEGG and Gene Ontology databases.

04 — Drug Repurposing
Automated ChEMBL API queries to retrieve clinical-grade inhibitors for overexpressed targets.

---

## Tech Stack

| Category        | Tools                      |
|-----------------|----------------------------|
| Workflow        | Snakemake                  |
| Languages       | Python 3.x, R 4.x          |
| Bioinformatics  | DESeq2, GSEApy, BiomaRt    |
| Data Science    | Pandas, ChEMBL WebClient   |

---

## Project Structure

AutoFlow-SNAI1/
  ├── Scripts/     # Core analysis and API scripts
  ├── data/        # Input count matrices and metadata
  ├── results/     # Automated plots and reports
  ├── Snakefile    # Workflow automation engine
  └── config.yaml  # Analysis parameters

---

## Execution

Run the complete pipeline across 4 parallel threads:

    snakemake --cores 4

---

## Author

Amna Asghar
Bioinformatics and Computational Oncology
https://www.linkedin.com/in/amna-asghar-030771274
