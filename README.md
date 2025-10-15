# Integrated Spatial Analysis of Ovarian Precancerous Lesions

This repository contains the R scripts and analysis code used in the study published at  
**[BioRxiv (2025.06.24.661327v1)](https://www.biorxiv.org/content/10.1101/2025.06.24.661327v1)**

---
### `scripts/`
#### `1_STIC_Nanostring_preprocessing.Rmd`
- Converts raw **WTA Nanostring data** into processed count data.  
- **Note:** You do **not** need to run this step if starting from the dataset `GSM9005069`.

#### `2_STIC_Nanostring_NMF.R`
- Performs **Non-negative Matrix Factorization (NMF)** for molecular subtype classification.  
- The resulting subtype information is stored in the `anno_data` column: Molecular_Subtype

#### `3_STIC_Nanostring_analysis.Rmd`
- Generates all main analyses and figures presented in the publication.  
- Before running, please modify the root path in the script:  

------------------------------------------------------------

## 📂Folder Structure and Script Descriptions

### scripts/
#### 1_STIC_Nanostring_preprocessing.Rmd
- Converts raw WTA Nanostring data into processed count data.
- Note: You do not need to run this step if starting from GSM9005069.

#### 2_STIC_Nanostring_NMF.R
- Performs Non-negative Matrix Factorization (NMF) for molecular subtype classification.
- The resulting subtype information is stored in the anno_data column "Molecular_Subtype".

#### 3_STIC_Nanostring_analysis.Rmd
- Performs the main analyses and generates figures included in the publication.
- Before running, please modify the path:
  root <- "Location"
  Replace "Location" with your actual project path.

#### R/STIC_Nanostring_Custom_Functions.R
- Contains custom plotting functions used in 3_STIC_Nanostring_analysis.Rmd,
  including the volcano plot function.

------------------------------------------------------------

### Notes
- The provided scripts assume the same data structure as in the original study.
- Analyses were performed using R version 4.3 or later with standard RNA expression and visualization packages.
- Each R Markdown file contains session information for reproducibility.

------------------------------------------------------------

© 2025 Tu-Yung Chang and collaborators.
This code is provided for academic and reproducibility purposes.