# Input Files Requirements

This pipeline requires gene expression datasets organized in a specific folder structure.

## Folder Structure

data/
└── PD_DATA/
    └── GSE_ID/
        ├── GSE_ID_counts.csv
        └── GSE_ID_metadata.csv

---

## 1. Counts File

### File Naming
GSE_ID_counts.csv

### Description
This file contains gene expression count data.

### Required Format

| Gene | Sample1 | Sample2 | Sample3 |
|-------|----------|-----------|-----------|
| GeneA | 50 | 60 | 45 |
| GeneB | 12 | 18 | 22 |

### Requirements
- First column must contain gene names
- Remaining columns must contain sample expression values
- All values must be numeric
- Sample names must match metadata file sample IDs

---

## 2. Metadata File

### File Naming
GSE_ID_metadata.csv

### Description
Contains sample classification and grouping information.

### Required Format

| SampleID | Group |
|-----------|------------|
| Sample1 | Control |
| Sample2 | Disease |

### Requirements
- Must contain sample ID column
- Must contain group/diagnosis/disease classification column
- Group values must include Control and Disease labels
- Sample IDs must match counts file column names

---

## Notes
- Multiple GSE datasets can be included
- Each dataset must be placed in its own folder
- File naming must remain consistent

