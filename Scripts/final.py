"""
Consolidation and Structural Formatting of Drug Repurposing Reports

Description:
This script serves as the final formatting layer of the pipeline. It ingests 
raw drug-target interaction data, standardizes column nomenclature, and 
applies clinical relevance filters (Phase 3 & 4). The resulting report is 
structured for clinical interpretation and publication-ready documentation.

Key Features:
- Dynamic column cleaning and name mapping to ensure cross-script compatibility.
- Adaptive phase filtering: Prioritizes late-stage clinical candidates (Phases 3 and 4).
- Metadata enrichment: Injects fixed context (Disease, Analysis type) into the final dataset.
- Resilient column selection: Validates existence of fields (e.g., Repurposing_Score) 
  before sorting or exporting to prevent runtime errors.

Dependencies:
- pandas: Data restructuring and formatting.
- sys: Command-line argument management for pipeline integration.

Author: Amna Asghar
"""
import pandas as pd
import sys

def build_final_report():
"""
    Standardizes and filters the drug-target dataset into a structured CSV report.
    
    Reads input from system arguments or default local paths, handles column 
    renaming variations, and exports the final filtered dataframe.
    """
    try:
        INPUT_FILE = sys.argv[1]
        OUTPUT_FILE = sys.argv[2]
    except IndexError:
        INPUT_FILE = "results/final_drug_report.csv"
        OUTPUT_FILE = "results/final_structured_report.csv"

    print(f"Reading from: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)

    # ---- Clean column names (removes spaces/extra characters) ----
    df.columns = [c.strip().replace(" ", "_") for c in df.columns]
    
    # DEBUG: Show what columns were actually found
    print("Found columns:", df.columns.tolist())

    # ---- Ensure numeric phase (Only if the column exists) ----
    # We check for common variations: 'Clinical_Phase' or 'Phase'
    phase_col = None
    for col in ["Clinical_Phase", "Phase", "ClinicalPhase"]:
        if col in df.columns:
            phase_col = col
            break
    
    if phase_col:
        df[phase_col] = pd.to_numeric(df[phase_col], errors="coerce")
        # Keep only Phase 3 & 4
        df = df[df[phase_col].isin([3, 4])]
        # Rename it to our standard for the final report
        df = df.rename(columns={phase_col: "Clinical_Phase"})
    else:
        print("Warning: No Clinical Phase column found. Skipping phase filtering.")

    # ---- Add fixed metadata ----
    df["Disease"] = "Breast Cancer"
    df["Analysis"] = "RNA-Seq Differential Expression"

    # ---- Rename biological columns for final output ----
    # Mapping old names (left) to new names (right)
    rename_map = {
        "Gene": "Gene Name",
        "P-Value": "P-value",
        "pvalue": "P-value",
        "Log2FC": "LogFC",
        "log2FoldChange": "LogFC"
    }
    df = df.rename(columns=rename_map)

    # ---- Final column structure ----
    final_cols = [
        "Disease", "Analysis", "Gene Name", "P-value", "LogFC", 
        "Direction", "Pathway", "Drug_Name", "Drug_ID", 
        "Clinical_Phase", "Repurposing_Score"
    ]

    # Only select columns that actually exist to avoid another KeyError
    existing_cols = [c for c in final_cols if c in df.columns]
    final_df = df[existing_cols]

    if "Repurposing_Score" in final_df.columns:
        final_df = final_df.sort_values(by="Repurposing_Score", ascending=False)

    final_df.to_csv(OUTPUT_FILE, index=False)
    print("[DONE] Final structured report saved:", OUTPUT_FILE)

if __name__ == "__main__":
    build_final_report()
