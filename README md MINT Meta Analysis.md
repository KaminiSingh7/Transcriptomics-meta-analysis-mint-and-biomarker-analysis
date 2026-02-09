README.md

##########################Main_structure_model################################

mint-meta-analysis/
│
├── data/
│   ├── example_PD_DATA/
│   │   └── GSE_SAMPLE/
│   │       ├── GSE_SAMPLE_counts.csv
│   │       └── GSE_SAMPLE_metadata.csv
│
├── results/                # auto generated
│
├── scripts/
│   └── mint_analysis_pipeline_enhanced.R
│
├── README.md
├── requirements.R
├── .gitignore
├── config.R                # optional improvement
└── LICENSE


# MINT Meta-Analysis Pipeline for Parkinson's Disease

This repository contains an R pipeline for performing MINT-based meta-analysis of multi-study gene expression datasets.

## 📊 Features
- Integrates multiple GEO datasets
- Performs MINT PLS-DA analysis
- Automatic feature selection
- Generates visualization reports
- Identifies potential biomarkers

---

## 📂 Input Data Structure

Place datasets in:

data/PD_DATA/

Each study folder must follow:

GSE_ID/
   GSE_ID_counts.csv
   GSE_ID_metadata.csv

---

## 📄 Counts File Requirements
- Rows: Genes
- Columns: Samples
- Numeric count data only

---

## 📄 Metadata Requirements
- Must contain sample IDs
- Must contain group column (Control/Disease)

---

## ⚙️ Installation

Install dependencies in R:


