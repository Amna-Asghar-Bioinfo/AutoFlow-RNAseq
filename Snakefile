docstring: """
SNAI1 Knockout RNA-Seq Analysis Pipeline
Author: Amna Asghar
Date: April 2026

This workflow automates:
1. Differential Expression Analysis (DESeq2)
2. Significant Gene Filtration (Log2FC >= 2.0)
3. Pathway Enrichment Analysis (GSEApy/Enrichr)
4. Drug Repurposing Discovery (ChEMBL API)
"""

# Rest of your Snakefile starts here...
configfile: "config.yaml"

# ------------------------------------------------------------
# FINAL TARGET
# ------------------------------------------------------------
rule all:
    input:
        "results/final_structured_report.csv"


# ------------------------------------------------------------
# 1. DESEQ2 (R)
# ------------------------------------------------------------
rule deseq:
    input:
        counts=config["counts_file"],
        metadata=config["metadata_file"]
    output:
        "results/deseq.csv"
    shell:
        """
        Rscript Scripts/Analyzing_with_DESeq.R {input.counts} {input.metadata} {output}
        """


# ------------------------------------------------------------
# 2. FILTER GENES (R)
# ------------------------------------------------------------
rule filter:
    input:
        "results/deseq.csv"
    output:
        up="results/up.csv",
        down="results/down.csv"
    shell:
        """
        Rscript Scripts/Filteration.R {input} {output.up} {output.down}
        """


# ------------------------------------------------------------
# 3. PATHWAY ANALYSIS (PYTHON)
# ------------------------------------------------------------
rule pathway:
    input:
        up="results/up.csv",
        down="results/down.csv"
    output:
        up="results/up_pathways.csv",
        down="results/down_pathways.csv"
    shell:
        """
        python Scripts/Extracting_Genes_with_Significant_Pathway.py {input.up} {input.down} {output.up} {output.down}
        """



# ------------------------------------------------------------
# 4. DRUG ANALYSIS (PYTHON)
# ------------------------------------------------------------

rule drug:
    input:
        up="results/up_pathways.csv",
        down="results/down_pathways.csv"
    output:
        "results/final_drug_report.csv"
    log: "logs/drug_search.log"
    shell:
        """
        python Scripts/Extracting_Drugs_Fast.py \
        {input.up} {input.down} {output} 2>&1 | tee results/drug_search.log
        """

# Rule 5: Uses that report
rule final:
    input:
        "results/final_drug_report.csv"   
    output:
        "results/final_structured_report.csv"
    shell:
        """
        python Scripts/final.py {input} {output}
        """
