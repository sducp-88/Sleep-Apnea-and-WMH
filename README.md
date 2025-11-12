## **WMH Analysis (UK Biobank) — Main Cross-Sectional & Short-Interval Change**

Public, reproducible code for the primary analyses of the UK Biobank project on sleep-apnea–related differences in white-matter hyperintensities (WMH).  
This notebook consolidates the main cross-sectional models and the exploratory short-interval longitudinal analyses.  
Extra models (e.g., extended adjustment sets, mediation, effect-modification) have been removed from this public version.

---

### **1) What this repository contains**

**WMH_Analysis_pub.ipynb** — single, publication-ready notebook that:

- Removes UK Biobank withdrawn-consent participants from the working dataset.  
- Performs cleaning & final dataset preparation (variable derivation, date harmonization, sample flow logging).  
- Provides a matching / ATO weighting scaffold (logistic PS with sklearn + balance reports) to support the main matched cohorts.  
- Runs WMH transformation diagnostics (raw vs. log1p(raw × volumetric_scaling_factor), normality and heteroscedasticity checks, AIC comparisons).  
- Fits primary cross-sectional OLS models (head-size–normalized WMH, log1p scale) using a prespecified covariate set; exports publication-grade tables and figures.  
- Implements the short-interval change analysis (I3–I2): annualized slopes and increase/decrease proportions, with supporting visualizations.  
- Generates figure legends, eMethods snippets, and package version prints for reproducibility.  
- **[NEW IN THIS UPDATE]** Adds automated FDR correction (Benjamini–Hochberg) for PWMH/DWMH families; fixes logging of sample counts within each stratum; updates visualization scripts for forest and bar-type plots; and improves compatibility with head-size normalization pipelines.  
- **[NEW FILE]** `run_apoe.sh` — Bash script used on **UK Biobank RAP** (Research Analysis Platform) to extract APOE genotypes (rs429358, rs7412) from BGEN/VCF files.  
  The script supports parallelized `bgenix`/`qctool` workflows and outputs a clean CSV file with `eid`, `APOE_e2e3e4` status, and genotype counts for downstream linkage to phenotype tables.

---

### **2) Data & privacy**

UK Biobank data are **not distributed** in this repository.  
You must have an approved UKB application and access to the corresponding data fields.  
The notebook expects a working CSV with required variables and uses a UKB withdrawn-consent list (e.g., `w48286_YYYYMMDD.csv`) to remove participants.  
Any sharing must comply with **UKB’s Material Transfer and Data Access Policies**.

---

### **3) Required inputs (typical)**

- **Working analysis dataset:** e.g., `data_processed.csv`, containing at least:
  participant identifier, `treatment_var`, `match_id`, WMH raw volumes, volumetric scaling factor, covariates, and dates for longitudinal subset.  
- **Withdrawn list:** e.g., `w48286_YYYYMMDD.csv`.  
- If variable names differ, adapt the **CONFIG** section in the notebook.

---

### **4) Environment**

- **Python ≥ 3.9**  
- **Core packages:** `pandas`, `numpy`, `statsmodels`, `scipy`, `matplotlib`, `seaborn`, `scikit-learn`, `patsy`  
- **Docs/figures:** `python-docx`, `openpyxl`, `graphviz`, `pillow`, `tqdm`

**Install:**
```bash
pip install pandas numpy statsmodels scipy matplotlib seaborn scikit-learn patsy python-docx openpyxl graphviz pillow tqdm
```

---

### **5) How to run**

1. Place your working CSV and withdrawn list in the project directory.  
2. Launch Jupyter:
   ```bash
   jupyter lab
   ```
3. Open **WMH_Analysis_pub.ipynb** and execute cells top-to-bottom.  

To run as a script:
```bash
jupyter nbconvert --to script WMH_Analysis_pub.ipynb
python WMH_Analysis_pub.py
```

To extract APOE genotypes on UKB-RAP:
```bash
bash run_apoe.sh
```

---

### **6) Statistical design**

- **Cross-sectional:** log1p(head-size–normalized WMH); OLS with cluster-robust SEs.  
- **Short-interval change:** annualized slope & proportion increase; LOESS and bar-type plots.  
- **BH-FDR** applied to PWMH/DWMH families.

---

### **7) Key outputs**

`OLS_Results_Manuscript_Table.csv/.docx`,  
`Forest_HSN_LogPct.png`,  
`Diagnostics_raw_vs_log1p.csv`,  
`Hist_*.png`, `QQ_*.png`,  
`FigureX_Longitudinal_WMH.png`,  
`Derived_Variables_UKB.csv`,  
`Sample_Flow.txt`,  
`PSM_Summary.txt`,  
and `APOE_Genotypes.csv` (from `run_apoe.sh`).

---

### **8) Reproducibility notes**

Withdrawn-consent removal is mandatory.  
Logs all derived and dropped variables.  
Final cell prints package versions for reproducibility.

---

### **9) Customization**

Edit the **CONFIG** block for variable names and file paths.  
Modify inclusion/exclusion criteria as needed.  
Figures are exported in PNG/PDF/SVG/TIFF (600 dpi, LZW compression).

---

### **10) Citation**

Cheng P, et al. *Sleep Apnea and White Matter Hyperintensities in the UK Biobank: Cross-Sectional Association and Short-Interval Change.* 2025 (manuscript in submission).  
**UK Biobank acknowledgment:** This research was conducted using the UK Biobank resource (Application XXXX).

---

### **11) License**

Use an open-source license such as **MIT**, **Apache-2.0**, or **CC-BY-4.0**.

---

### **12) Contact**

**Peng Cheng, MD**  
Visiting Research Fellow, Mayo Clinic  
Neurologist, Second Hospital of Shandong University  
📧 Cheng.Peng@Mayo.edu | sducp@email.sdu.edu.cn  
