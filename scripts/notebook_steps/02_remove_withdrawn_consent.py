# Extracted from WMH_Analysis.ipynb, code cell 1.
# Notebook heading: # Removal of Participants Who Withdrew Consent
# Run this file from the repository root unless a local CONFIG section is edited.

# Removal of Participants Who Withdrew Consent
"""
Publication-Ready Description
-----------------------------
Objective
    Remove UK Biobank participants who have withdrawn consent from a single
    analysis dataset and report post-removal sample sizes by treatment arm.

Inputs
    1) data20250415.csv
       - Contains at minimum:
           • A participant identifier column (auto-detected; see below)
           • A treatment assignment column named exactly: treatment_var
             (e.g., 'Study' vs 'Control' or 1 vs 0)
    2) w48286_20250818.csv
       - Headerless CSV; the FIRST column lists participant IDs who withdrew consent.

Outputs
    • data.csv                             (cleaned analysis dataset)
    • Console report (publication-grade English):
        - N removed in total
        - Post-removal sample sizes by treatment_var

ID Handling
    • The participant ID column is detected heuristically from common names:
      ['Participant_ID', 'participant_id', 'eid', 'EID', 'ID', 'id'].
      If none are present, the first column is used.
    • All IDs are compared as strings to avoid leading-zero mismatches.

Notes
    • This script performs a *strict* row-wise removal: any row whose ID appears
      in the withdrawal list is dropped. No cluster or matched-set logic is applied.
    • The script prints both pre- and post-removal counts by treatment_var.

Usage
    python remove_withdrawals_simple.py
"""

from pathlib import Path
from typing import Dict, Set
import pandas as pd

# -----------------------
# Configuration
# -----------------------
DATA_IN = "data20250415.csv"
WITHDRAWN_FN = "w48286_20250818.csv"   # headerless; first column are IDs
DATA_OUT = "data.csv"

# -----------------------
# Utilities
# -----------------------
def detect_id_col(df: pd.DataFrame) -> str:
    """
    Detect the participant ID column from common names.
    Fallback: use the first column if none of the common names are present.
    """
    candidates = ["Participant_ID", "participant_id", "eid", "EID", "ID", "id"]
    for c in candidates:
        if c in df.columns:
            return c
    return df.columns[0]

def load_withdrawn_ids(path: str | Path) -> Set[str]:
    """
    Load withdrawn IDs from a headerless CSV where the FIRST column contains IDs.
    Return as a set of strings (NaN and empty values discarded).
    """
    wd = pd.read_csv(path, header=None)
    if wd.empty:
        return set()
    col0 = wd.columns[0]
    return set(
        wd[col0]
        .astype(str)
        .str.strip()
        .replace({"nan": ""})
        .dropna()
        .loc[lambda s: s.ne("")]
        .unique()
        .tolist()
    )

def summarize_by_treatment(df: pd.DataFrame, title: str) -> Dict[str, int]:
    """
    Print and return counts by treatment_var.
    """
    print(f"\n=== {title} ===")
    if "treatment_var" not in df.columns:
        raise KeyError("Required column 'treatment_var' not found in the dataset.")
    counts = df["treatment_var"].value_counts(dropna=False).to_dict()
    total = len(df)
    # Pretty print
    print(f"Total sample size: {total}")
    for k, v in counts.items():
        print(f"  treatment_var = {k!r}: n = {v}")
    # Return a flattened dict for optional downstream use
    out = {"total": total}
    out.update({f"treatment_{k}": int(v) for k, v in counts.items()})
    return out

# -----------------------
# Main
# -----------------------
if __name__ == "__main__":
    # Existence checks
    for p in [DATA_IN, WITHDRAWN_FN]:
        if not Path(p).exists():
            raise FileNotFoundError(f"Required input not found: {p}")

    # Load main dataset
    df = pd.read_csv(DATA_IN)
    if df.empty:
        raise ValueError("Input dataset is empty: data20250415.csv")

    # Detect ID column and normalize to string for safe comparison
    id_col = detect_id_col(df)
    df[id_col] = df[id_col].astype(str).str.strip()

    # Sanity: ensure treatment_var is present
    if "treatment_var" not in df.columns:
        raise KeyError(
            "The dataset must contain a 'treatment_var' column "
            "(e.g., 'Study' vs 'Control' or 1 vs 0)."
        )

    # Load withdrawal list
    withdrawn_ids = load_withdrawn_ids(WITHDRAWN_FN)
    print(f"Loaded withdrawal list: {len(withdrawn_ids)} unique IDs.")

    # Pre-removal report
    summarize_by_treatment(df, "Pre-removal sample size by treatment_var")

    # Remove withdrawn IDs
    before_n = len(df)
    mask_keep = ~df[id_col].isin(withdrawn_ids)
    removed_n = int((~mask_keep).sum())
    df_clean = df.loc[mask_keep].copy()

    # Post-removal report
    summarize_by_treatment(df_clean, "Post-removal sample size by treatment_var")

    # High-level removal statement (publication-ready)
    print("\n--- Publication-Grade Summary ---")
    print(f"Participants removed due to withdrawn consent: n = {removed_n} "
          f"(from N = {before_n} to N = {len(df_clean)}).")
    print("Post-removal counts are reported above by 'treatment_var'.")

    # Save cleaned dataset
    df_clean.to_csv(DATA_OUT, index=False, encoding="utf-8-sig")
    print(f"\n✅ Cleaned dataset saved to: {DATA_OUT}")
