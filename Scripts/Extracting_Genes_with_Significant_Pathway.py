"""
Functional Enrichment Analysis and Gene-Pathway Mapping Pipeline

Description:
This script performs over-representation analysis (ORA) using GSEApy to identify 
statistically significant biological pathways and gene ontology terms associated 
with a given gene set. It maps identified pathways back to the original gene 
expression data to facilitate biological interpretation.

Key Features:
- Utilizes KEGG_2021_Human and GO_Biological_Process_2021 libraries.
- Implements a statistical fallback: Prioritizes Adjusted P-values (< 0.05), 
  defaulting to raw P-values and top 10 terms if significance is not met.
- Performs many-to-one mapping: Concatenates multiple significant pathways for 
  each gene into a single 'Significant_Pathways' column.
- Supports Snakemake/CLI execution via positional system arguments.

Dependencies:
- gseapy: Interface for Enrichr API.
- pandas: Data merging and set-based filtering.
- os/sys: File system validation and command-line interface.

Author: Amna Asghar
"""
import gseapy as gp
import pandas as pd
import os
import sys

# ------------------------------------------------------------
# 1. COMMAND LINE INPUTS
# ------------------------------------------------------------
# usage:
# python script.py up_file down_file up_output down_output

UP_INPUT = sys.argv[1]
DOWN_INPUT = sys.argv[2]
UP_OUTPUT = sys.argv[3]
DOWN_OUTPUT = sys.argv[4]

MY_DISEASE = "SNAI1_KO_Analysis"

# ------------------------------------------------------------
# 2. FUNCTION FOR ENRICHMENT
# ------------------------------------------------------------
def generate_pathway_report(file_path, output_path, direction_label):
"""
    Executes enrichment, maps genes to terms, and exports a consolidated report.

    Args:
        file_path (str): CSV containing gene names (must have 'Gene.Name' column).
        output_path (str): Destination for the enriched data CSV.
        direction_label (str): Metadata label (e.g., 'Upregulated').
    """
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return

    df = pd.read_csv(file_path)
    gene_list = df['Gene.Name'].dropna().unique().tolist()

    if not gene_list:
        print(f"No genes in {direction_label} file.")
        return

    print(f"\n[INFO] Running GSEApy for {direction_label} ({len(gene_list)} genes)")

    try:
        # --------------------------------------------------------
        # 3. ENRICHMENT ANALYSIS
        # --------------------------------------------------------
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
            organism='human'
        )

        sig_results = enr.results[enr.results['Adjusted P-value'] < 0.05]

        # fallback logic
        if sig_results.empty:
            print(f"   -> No significant pathways (Adj P < 0.05). Using fallback.")
            sig_results = enr.results[enr.results['P-value'] < 0.05] \
                .sort_values('P-value') \
                .head(10)

        if sig_results.empty:
            print(f"   -> No enrichment found for {direction_label}")
            return

        # --------------------------------------------------------
        # 4. MAP GENES TO PATHWAYS
        # --------------------------------------------------------
        gene_to_pathway = {}

        for _, row in sig_results.iterrows():
            pathway = row['Term']
            genes_in_pathway = str(row['Genes']).split(';')

            for g in genes_in_pathway:
                g = g.upper()

                if g in gene_to_pathway:
                    if pathway not in gene_to_pathway[g]:
                        gene_to_pathway[g] += "; " + pathway
                else:
                    gene_to_pathway[g] = pathway

        # --------------------------------------------------------
        # 5. MERGE WITH INPUT DATA
        # --------------------------------------------------------
        df['Gene_Upper'] = df['Gene.Name'].str.upper()
        pathway_df = df[df['Gene_Upper'].isin(gene_to_pathway.keys())].copy()

        if pathway_df.empty:
            print(f"   -> No overlap between enrichment and gene list")
            return

        pathway_df['Significant_Pathways'] = pathway_df['Gene_Upper'].map(gene_to_pathway)
        pathway_df['Disease'] = MY_DISEASE
        pathway_df['Analysis_Type'] = f"RNA-Seq ({direction_label})"

        pathway_df.drop(columns=['Gene_Upper'], inplace=True)

        # --------------------------------------------------------
        # 6. SAVE OUTPUT
        # --------------------------------------------------------
        pathway_df.to_csv(output_path, index=False)

        print(f"   -> Saved: {output_path}")
        print(f"   -> Top pathways: {', '.join(sig_results['Term'].head(2).tolist())}")

    except Exception as e:
        print(f"[ERROR] GSEApy failed: {e}")


# ------------------------------------------------------------
# 3. RUN BOTH ANALYSES
# ------------------------------------------------------------
generate_pathway_report(UP_INPUT, UP_OUTPUT, "Upregulated")
generate_pathway_report(DOWN_INPUT, DOWN_OUTPUT, "Downregulated")
