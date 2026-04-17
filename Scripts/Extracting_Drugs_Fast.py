"""
ChEMBL Drug-Gene Interaction Discovery Pipeline

Description:
This script identifies potential drug candidates for upregulated and downregulated 
gene sets derived from differential expression analysis. It utilizes the ChEMBL 
database API to map gene symbols to protein targets and retrieves associated 
bioactive molecules.

Key Features:
- Implements a signal-based timeout mechanism to handle API latency.
- Filters for 'SINGLE PROTEIN' target types to ensure mapping specificity.
- Consolidates gene expression metrics (Log2FC, p-values) with pharmacological data.
- Handles API rate-limiting via controlled time delays and record capping.

Dependencies:
- pandas: Data manipulation and CSV I/O.
- chembl_webresource_client: Interface for the ChEMBL database.
- signal: Unix-based timeout management.

Author: Amna Asghar
"""


import pandas as pd
import time
from chembl_webresource_client.new_client import new_client
import signal

# -------------------------
# SETTINGS
# -------------------------
UP_FILE = "results/up_pathways.csv"
DOWN_FILE = "results/down_pathways.csv"
OUTPUT_FILE = "results/final_drug_report.csv"

target_api = new_client.target
activity_api = new_client.activity

# -------------------------
# TIMEOUT HANDLER (IMPORTANT)
# -------------------------
class TimeoutError(Exception):
    pass

def handler(signum, frame):
    raise TimeoutError()

signal.signal(signal.SIGALRM, handler)

# -------------------------
# SAFE CHEMBL QUERY WRAPPER
# -------------------------
def safe_target_search(gene):
    signal.alarm(15)  # 15 sec timeout
    try:
        res = list(target_api.search(gene))
        signal.alarm(0)
        return res[:3]  # limit explosion
    except:
        signal.alarm(0)
        return []

def safe_activity_query(tid):
    signal.alarm(15)
    try:
        res = list(activity_api.filter(target_chembl_id=tid))
        signal.alarm(0)
        return res[:20]  # prevent huge dump
    except:
        signal.alarm(0)
        return []

# -------------------------
# MAIN FUNCTION
# -------------------------
def find_drugs(file_path, label):
"""
    Iterates through a list of genes to query ChEMBL for bioactive compounds.
    
    Args:
        file_path (str): Path to the input CSV containing gene names and statistics.
        label (str): Classification label for the gene set (e.g., 'Upregulated').
        
    Returns:
        pd.DataFrame: A table of gene-drug interactions with associated metadata.
    """
    df = pd.read_csv(file_path)

    results = []
    seen = set()

    print(f"\n[INFO] Processing {label}: {len(df)} genes")

    for _, row in df.iterrows():
        gene = str(row["Gene.Name"]).strip()
        print(f"   -> {gene}")

        try:
            targets = safe_target_search(gene)

            for t in targets:
                if t.get("target_type") != "SINGLE PROTEIN":
                    continue

                tid = t.get("target_chembl_id")
                if not tid:
                    continue

                activities = safe_activity_query(tid)

                for a in activities:
                    drug = a.get("molecule_chembl_id")
                    if not drug:
                        continue

                    key = f"{gene}_{drug}"
                    if key in seen:
                        continue

                    results.append({
                        "Gene": gene,
                        "Log2FC": row.get("log2FoldChange"),
                        "P-Value": row.get("padj"),
                        "Direction": label,
                        "Pathway": row.get("Significant_Pathways"),
                        "Drug_ID": drug,
                        "Target_ID": tid
                    })

                    seen.add(key)

            time.sleep(0.2)  # small delay (safe for API)

        except TimeoutError:
            print(f"   ⚠ Timeout: {gene}")
            continue
        except Exception as e:
            print(f"   ⚠ Error: {gene} → {e}")
            continue

    return pd.DataFrame(results)

# -------------------------
# RUN PIPELINE
# -------------------------
up = find_drugs(UP_FILE, "Upregulated")
down = find_drugs(DOWN_FILE, "Downregulated")

final = pd.concat([up, down], ignore_index=True)

final.to_csv(OUTPUT_FILE, index=False)

print("\n DONE!")
print(f"File saved → {OUTPUT_FILE}")
print(f"Total drug hits → {len(final)}")
