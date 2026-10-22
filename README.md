# IBD and Hematological Comorbidities

R code for *Hematological Comorbidities and Blood Cell and Biochemical Profiles in Inflammatory Bowel Disease: Evidence from a Prospective Cohort Study in the UK and a Multicenter Cross-Sectional Study in China* (Jiang, Yang, Chen et al.).

## Scripts (run in order, from the repository root)

| Script | Purpose | Outputs |
|---|---|---|
| `01_cohort_and_outcomes.R` | Build UK Biobank cohort, IBD exposure, outcomes, covariates, blood traits | `data/*.rds` |
| `02_cox_ibd_hematological.R` | Time-varying Cox models, IBD → hematological outcomes; sensitivity analyses | Tables S10–S12 |
| `03_ldsc_commands.sh` | LDSC genetic correlation (`ldsc.py`) | `data/ldsc/logs/` |
| `03_mr_ldsc.R` | Two-sample MR and LDSC summary | Tables S6–S8, S13, S14 |
| `04_blood_trait_profile.R` | 65 blood traits by IBD status, UKB and Chinese cohort | Tables S15–S19 |
| `05_traits_incident_disease.R` | Baseline blood traits → incident hematological disease | Tables S20, S21 |
| `06_mediation.R` | Mediation by baseline blood traits | Table S22 |

## Requirements

R ≥ 4.2: `survival`, `openxlsx`, `data.table`, `TwoSampleMR`, `MRPRESSO`, `MASS`, `R.utils`.
LDSC v1.0.1 (Python 2.7) with European LD scores and the HapMap3 SNP list.

## Data

Individual-level data aren't shared. UK Biobank data are available to approved researchers (Application 157997). FinnGen R10 and the GWAS Catalog are public. ONE-IBD data are available from the corresponding author on reasonable request.

## Contact

Xue Li: xueli157@zju.edu.cn
