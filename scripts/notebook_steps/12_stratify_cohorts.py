# Extracted from WMH_Analysis.ipynb, code cell 11.
# Notebook heading: # Stratify process
# Run this file from the repository root unless a local CONFIG section is edited.

# Stratify process
# ==============================================================================================
# This script:
#   • Reuses cleaning outputs (SA_ascertain_group, Years_since_sleep_disorder) or derives them if missing.
#   • Creates Study-only stratification columns (Controls = NaN):
#       - SA_MRI_interval_stratum   (LE10_years / GT10_years)
#       - Snoring_Stratum
#       - Sleepy_SA_Stratum
#   • Additionally prints summary counts for:
#       - SA_ascertain_group (source of SA diagnosis)
#       - Sex (sex distribution)
#       - Self_Report_Only proportion within Study group
#   • Saves *_stratified.csv for downstream analyses.

import pandas as pd
import numpy as np
import os
import re

# ---------- Input & Output ----------
cohorts = {
    "primary_cohort.csv": "primary_cohort_stratified.csv",
    "sensitivity_cohort.csv": "sensitivity_cohort_stratified.csv"
}

# ---------- Helpers ----------
def _clean_cols(df: pd.DataFrame) -> pd.DataFrame:
    cols = (df.columns
              .str.replace(r"[^0-9a-zA-Z]+", "_", regex=True)
              .str.replace(r"_{2,}", "_", regex=True)
              .str.strip("_"))
    out = df.copy(); out.columns = cols
    return out

def _derive_sa_flags_and_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Derive group, treatment_var, and SA_ascertain_group (same as cleaning)."""
    out = df.copy()
    col_any, col_main, col_sec = "Diagnoses_ICD10", "Diagnoses_main_ICD10", "Diagnoses_secondary_ICD10"

    sr_pat = re.compile(r"^Non_cancer_illness_code_self_reported_Instance_[0-2]_Array_([0-9]|[1-2][0-9]|3[0-3])$")
    self_report_cols = [c for c in out.columns if sr_pat.match(c)]
    if self_report_cols:
        sr_block = out[self_report_cols].astype(str)
        sr_has_1123 = sr_block.apply(lambda s: s.eq("1123")).any(axis=1)
    else:
        sr_has_1123 = pd.Series(False, index=out.index)

    has_any  = out.get(col_any,  pd.Series("", index=out.index)).fillna("").str.contains("G473", regex=False)
    has_main = out.get(col_main, pd.Series("", index=out.index)).fillna("").str.contains("G473", regex=False)
    has_sec  = out.get(col_sec,  pd.Series("", index=out.index)).fillna("").str.contains("G473", regex=False)

    out["treatment_var"] = np.where(has_any | sr_has_1123, 1, 0)
    out["group"] = np.where(out["treatment_var"].eq(1), "Study", "Control")

    out["SA_ascertain_group"] = pd.NA
    out.loc[out["group"] == "Control", "SA_ascertain_group"] = "No_SA"
    out.loc[(out["group"] == "Study") & has_main,              "SA_ascertain_group"] = "Hospital_Primary"
    out.loc[(out["group"] == "Study") & (~has_main) & has_sec, "SA_ascertain_group"] = "Hospital_Secondary"
    out.loc[(out["group"] == "Study") & (~has_any),            "SA_ascertain_group"] = "Self_Report_Only"
    return out

def _ensure_years_since_sa(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure Years_since_sleep_disorder exists (same exclusion rules as cleaning)."""
    out = df.copy()
    mri_col = "Date_of_attending_assessment_centre_Instance_2"
    sa_col  = "Sleep_Disorder_Diagnosis_Date"

    if "Years_since_sleep_disorder" not in out.columns:
        if mri_col in out.columns:
            out[mri_col] = pd.to_datetime(out[mri_col], errors="coerce")
        if sa_col in out.columns:
            out[sa_col] = pd.to_datetime(out[sa_col], errors="coerce")

        if (mri_col in out.columns) and (sa_col in out.columns):
            delta_days = (out[mri_col] - out[sa_col]).dt.days
            out["Years_since_sleep_disorder"] = np.where(
                (~out[mri_col].isna()) & (~out[sa_col].isna()),
                delta_days / 365.25,
                np.nan
            )
        else:
            out["Years_since_sleep_disorder"] = np.nan

    if "group" not in out.columns or "treatment_var" not in out.columns:
        out = _derive_sa_flags_and_groups(out)

    neg_mask     = out["Years_since_sleep_disorder"] < 0
    study_mask   = out.get("treatment_var", pd.Series(dtype=int)) == 1
    control_mask = out.get("treatment_var", pd.Series(dtype=int)) == 0

    out.loc[control_mask & neg_mask, "Years_since_sleep_disorder"] = np.nan
    excl_n = int((study_mask & neg_mask).sum())
    if excl_n > 0:
        out = out[~(study_mask & neg_mask)]
        print(f"[info] Dropped {excl_n} Study participants with negative Years_since_sleep_disorder (<0).")

    return out

# ---------- Main ----------
for infile, outfile in cohorts.items():
    print(f"\n=== Processing: {infile} ===")
    df = pd.read_csv(infile)
    df = _clean_cols(df)

    # Reuse or derive SA grouping
    if ("SA_ascertain_group" not in df.columns) or ("group" not in df.columns) or ("treatment_var" not in df.columns):
        df = _derive_sa_flags_and_groups(df)

    # Ensure Years_since_sleep_disorder consistent with cleaning
    df = _ensure_years_since_sa(df)

    is_study = df["group"].eq("Study")

    # --- 1) SA_MRI_interval_stratum (Study only) ---
    df["SA_MRI_interval_stratum"] = pd.Series(pd.array([pd.NA] * len(df), dtype="object"))
    if "Years_since_sleep_disorder" in df.columns:
        bins = [-float("inf"), 10, float("inf")]
        labels = ["LE10_years", "GT10_years"]
        df.loc[is_study, "SA_MRI_interval_stratum"] = pd.cut(
            df.loc[is_study, "Years_since_sleep_disorder"],
            bins=bins, labels=labels
        )

    # --- 2) Snoring_Stratum (Study only) ---
    df["Snoring_Stratum"] = pd.Series(pd.array([pd.NA] * len(df), dtype="Int64"))
    if "Snoring_Group" in df.columns:
        df.loc[is_study, "Snoring_Stratum"] = df.loc[is_study, "Snoring_Group"].astype("Int64")
    elif "Snoring_Instance_0" in df.columns:
        df.loc[is_study, "Snoring_Stratum"] = (df.loc[is_study, "Snoring_Instance_0"] == 1).astype("Int64")

    # --- 3) Sleepy_SA_Stratum (Study only) ---
    df["Sleepy_SA_Stratum"] = pd.Series(pd.array([pd.NA] * len(df), dtype="Int64"))
    if "Daytime_dozing_sleeping_Instance_0" in df.columns:
        df.loc[is_study, "Sleepy_SA_Stratum"] = df.loc[is_study, "Daytime_dozing_sleeping_Instance_0"].isin([1, 2]).astype("Int64")

    # --- Save ---
    df.to_csv(outfile, index=False, encoding="utf-8-sig")
    print(f"[OK] Saved Study-only–stratified file: {outfile}")

    # --- Summaries (Study only) ---
    print("\n[Study-only summaries]")
    for col in ["SA_ascertain_group", "Sex", "SA_MRI_interval_stratum", "Snoring_Stratum", "Sleepy_SA_Stratum"]:
        if col in df.columns:
            print(f"{col}:\n", df.loc[is_study, col].value_counts(dropna=False))
        else:
            print(f"{col}: (not present)")

    # --- Extra: Self-report proportion within Study group ---
    if "SA_ascertain_group" in df.columns:
        total_study = is_study.sum()
        self_report_n = (df.loc[is_study, "SA_ascertain_group"] == "Self_Report_Only").sum()
        prop = (self_report_n / total_study * 100) if total_study > 0 else np.nan
        print(f"\n[Info] Self_Report_Only participants: {self_report_n} / {total_study} ({prop:.1f}%) within Study group")
