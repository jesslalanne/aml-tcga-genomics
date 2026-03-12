# Data

The dataset used in this project is publicly available through **cBioPortal** and does not require registration.

> ⚠️ Raw data files are **not committed to this repository** (file size). Follow the steps below to reproduce the analysis locally.

---

## Dataset: AML TCGA Pan-Cancer Atlas 2018

**Study ID:** `laml_tcga_pan_can_atlas_2018`  
**Source:** [https://www.cbioportal.org](https://www.cbioportal.org/study/summary?id=laml_tcga_pan_can_atlas_2018)  
**Cohort size:** ~200 AML patients  
**Data types:** Clinical, somatic mutations (MAF), RNA-seq (RSEM)

---

## Download Instructions

### Option A — cBioPortal Web Interface (recommended)

1. Go to: https://www.cbioportal.org/study/summary?id=laml_tcga_pan_can_atlas_2018
2. Click **"Download"** in the top right
3. Select **"All"** data types and download the ZIP archive
4. Extract the archive into this `data/` directory

### Option B — cBioPortal API

```python
import requests, zipfile, io, os

study_id = "laml_tcga_pan_can_atlas_2018"
url = f"https://www.cbioportal.org/api/studies/{study_id}/downloadAll"
response = requests.get(url)
z = zipfile.ZipFile(io.BytesIO(response.content))
z.extractall("data/")
print("Done.")
```

---

## Expected File Structure

After download, your `data/` folder should contain:

```
data/
├── data_clinical_patient.txt        # Clinical metadata: age, OS, FAB subtype, cytogenetics
├── data_clinical_sample.txt         # Sample-level metadata
├── data_mutations.txt               # Somatic mutations in MAF format (Hugo_Symbol, VAF, etc.)
├── data_mrna_seq_v2_rsem.txt        # RNA-seq gene expression (RSEM log2 normalized)
├── data_cna.txt                     # Copy number alterations (optional)
└── meta_study.txt                   # Study metadata
```

---

## Key Variables Used

| Variable | File | Description |
|---|---|---|
| `OS_MONTHS` | clinical_patient | Overall survival in months |
| `OS_STATUS` | clinical_patient | 0 = alive, 1 = deceased |
| `AGE` | clinical_patient | Patient age at diagnosis |
| `FAB_TYPE` | clinical_patient | FAB morphological classification (M0–M7) |
| `Hugo_Symbol` | mutations | Gene name |
| `t_alt_count`, `t_ref_count` | mutations | Tumor allele counts → VAF |
| `Variant_Classification` | mutations | Missense, Nonsense, Frameshift, etc. |

---

## Preprocessing Notes

- **Duplicate patients** in mutations file are expected (one row per mutation per patient) — use `patient_id` as groupby key
- **RSEM values** are log2(x+1) transformed before PCA/clustering
- **Survival analysis:** patients with `OS_MONTHS = 0` or missing are excluded
- **VAF** is computed as: `t_alt_count / (t_alt_count + t_ref_count)`
