# Survival Analysis of Wilms' Tumor

## Introduction
This repository contains the R scripts and data analysis for the research project titled "Survival Analysis of Wilms Tumor." The study aims to understand the factors influencing relapse in patients with Wilms' tumor, a type of kidney cancer primarily affecting children.

## Problem Statement
Wilms' tumor remains a significant health concern despite advances in treatment. This study focuses on analyzing the time until relapse and the effectiveness of new treatments, considering factors such as age, stage of the disease, and histology.

## Data
The study was conducted on a cohort of 4028 patients diagnosed with Wilms' tumor. Data were collected on various factors, including age, stage of the disease, and histology.

## Methodology
The research utilized a two-phase study design:
1. **Phase 1**: Initial screening of a large cohort.
2. **Phase 2**: Detailed analysis of a subset of patients.

### Statistical Analysis
- **Survival Analysis**: Kaplan-Meier estimates and Cox proportional hazards models.
- **Handling Censoring**: Appropriate methods for censored data.
- **Software**: R, with key packages including `survival` and `survminer`.

## Results
Key findings include:
- Age at diagnosis, stage of the disease, and histological classification significantly influence the time until tumor relapse.
- The Kaplan-Meier survival curves and Cox model results are provided in the report.

## Conclusion
The study highlights the importance of age, tumor stage, and histological prognosis in predicting survival outcomes for Wilms' tumor patients. Future research should focus on improving measurement accuracy and exploring additional prognostic factors.

## Repository Structure
- `dataset/`: Contains the dataset used for analysis.
- `survival.r`: R scripts for data preparation, analysis, and visualization.
- `results/`: Output files, including survival curves and model diagnostics.


## Authors
- Saliou CISSE
- Ethan ADA
- Mohamed AL JALANJI

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
We would like to thank the National Wilms Tumor Study Group and the Children's Oncology Group for providing the data and resources for this research.
