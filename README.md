# Thyroid Cancer Dedifferentiation-State Analysis

Code and processed source-data tables for:

> Yin S, Li S, Shi C. *Integrated transcriptomic and single-cell mapping
> of a one-carbon-metabolism-enriched dedifferentiation state in thyroid
> cancer prioritizes SLC7A5, MTHFD2 and AXL.* Journal of Translational
> Medicine. 2026. (Under review.)

[![DOI](https://zenodo.org/badge/1253988198.svg)](https://doi.org/10.5281/zenodo.20453860)

> ⚠️ Zenodo badge above will activate after the first release (see "Citation" below).

---

## Overview

This repository contains the R analysis scripts and processed source-data
tables that reproduce all figures and tables in the manuscript. The composite
dedifferentiation-state score is defined as:

> z(One-carbon) + z(MAPK) + z(dedifferentiation/EMT) - z(TDS)

Six public cohorts are integrated:

| Cohort       | Modality     | Role                      | N        |
|--------------|--------------|---------------------------|----------|
| TCGA-THCA    | Bulk RNA-seq | Discovery                 | 505      |
| GSE151179    | Bulk array   | RAI exploratory           | 49       |
| GSE33630     | Bulk array   | PTC/ATC validation        | 105      |
| GSE76039     | Bulk array   | PDTC/ATC validation       | 37       |
| GSE184362    | scRNA-seq    | PTC ecology               | 197,955  |
| GSE232237    | scRNA-seq    | PTC/ATC ecology           | 84,803   |

## Repository structure
## Reproducibility

### Environment

- R version: 4.5.3
- OS: Windows 11 / Ubuntu 22.04
- Key packages (versions in `sessionInfo.txt`):
  - TCGAbiolinks 2.38.0
  - GEOquery 2.78.0
  - limma 3.66.0
  - fgsea 1.36.2
  - Seurat 5.5.0
  - copykat 1.1.0

### Data access

Raw transcriptomic data are NOT included in this repository (too large).
They can be obtained from the original sources:

- **TCGA-THCA**: Genomic Data Commons (`TCGAbiolinks::GDCquery`)
- **GSE151179 / GSE33630 / GSE76039 / GSE184362 / GSE232237**: NCBI GEO
- **Human Protein Atlas**: https://www.proteinatlas.org (XML/TSV download)
- **ChEMBL**: https://www.ebi.ac.uk/chembl (REST API)

Run scripts in order (`01_*.R` through `09_*.R`). Each script reads from
the previous step's output and produces intermediate tables under `data/`.

### Random seeds

All stochastic steps (copyKAT sampling, train/test splits) use `set.seed(123)`
unless noted otherwise.

## Citation

If you use this code, please cite both the paper and this repository:
## Contact

- Corresponding authors: Chang Shi (changshi@csu.edu.cn), Shi Yin (shiyin910515@csu.edu.cn)
- Issues / questions: please open a GitHub Issue

## License

MIT License — see [LICENSE](LICENSE) file.
