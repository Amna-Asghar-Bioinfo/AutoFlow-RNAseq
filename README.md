#AutoFlow-SNAI1: Automated RNA-Seq & Drug Discovery

An end-to-end, automated bioinformatics workflow designed to analyze transcriptomic changes in **SNAI1 Knockout** models of Breast Cancer. This pipeline transitions from raw differential expression to targeted drug repurposing candidates using Snakemake.

## Key Highlights
* **Full Automation:** Single-command execution via Snakemake.
* **Clinical Relevance:** Filters for Phase 3/4 drug candidates using the ChEMBL API.
* **Reproducible Research:** Configuration-driven analysis for consistent results.

## Pipeline Logic
1. **Differential Expression:** Processes raw count matrices using `DESeq2`.
2. **Strict Filtration:** Identifies biomarkers with `padj < 0.05` and `|log2FC| >= 2.0`.
3. **Pathway Enrichment:** Functional annotation via `GSEApy` (KEGG/GO).
4. **Drug Repurposing:** Automated discovery of clinical-grade inhibitors for overexpressed genes.

## Tech Stack
| Category | Tools |
| :--- | :--- |
| **Workflow** | Snakemake |
| **Languages** | Python 3.x, R 4.x |
| **Bioinformatics** | DESeq2, GSEApy, BiomaRt |
| **Data Science** | Pandas, ChEMBL WebClient |

##Project Structure


├── Scripts/            # Core analysis and API scripts
├── data/               # Input matrices and metadata
├── results/            # Automated Plots & Reports
├── Snakefile           # The automation engine
└── config.yaml         # Analysis parameters


##Execution
To run the entire pipeline across 4 threads:

snakemake --cores 4

**Author:** Amna Asghar  
**Focus:** Bioinformatics & Computational Oncology
**Connect:** [![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/amna-asghar-030771274)



