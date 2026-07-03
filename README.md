# Quartet Reference Materials for MS-based Proteomics

This repository contains the R scripts used for the generation, evaluation, and application of property values for Quartet Reference Materials in MS-based proteomics.

The current version of this repository includes analysis and figure-generation scripts only. Large-scale raw data, processed data, and intermediate results are not included at this stage because of file size limitations. These data will be deposited in an appropriate public proteomics data repository and linked here upon manuscript submission or publication.

## Repository Status

This repository is currently under active development for manuscript submission.

At this stage, the repository provides:

- R scripts for qualitative and quantitative property value generation;
- R scripts for normality assessment, homogeneity assessment, and stability assessment;
- R scripts for certified property value generation;
- R scripts for laboratory performance evaluation;
- R scripts used to generate the main figures.

Large data files are not included in this GitHub repository.

## File Description

| File | Description |
|---|---|
| `1_qualitative_properties.R` | Analysis of qualitative property values for Quartet Reference Materials. |
| `2_quantitative_properties.R` | Analysis of quantitative property values based on relative abundance ratios. |
| `3_lab_performance.R` | Evaluation of laboratory performance using Quartet Reference Materials. |
| `4_shapiro_test.R` | Normality assessment of candidate quantitative property values. |
| `5_homogeneity.R` | Homogeneity assessment of candidate property values. |
| `6_stability.R` | Short-term and long-term stability assessment. |
| `7_certified_values.R` | Generation of certified property values after quality filtering. |
| `8_application.R` | Application analysis using certified property values of Quartet Reference Materials. |
| `DEP.R` | Differential expression/protein analysis helper script. |
| `PCA.R` | Principal component analysis script. |
| `figure2_qualitative_property_values.R` | Script for generating Figure 2. |
| `figure3_ratio_validation.R` | Script for generating Figure 3. |
| `figure4_homogeneity_stability.R` | Script for generating Figure 4. |
| `figure5_application.R` | Script for generating Figure 5. |

## Analysis Overview

The analysis workflow includes the following major steps:

1. Generation of qualitative property values for Quartet Reference Materials;
2. Generation of quantitative property values based on donor-to-donor abundance ratios;
3. Filtering of candidate property values by statistical and quality-control criteria;
4. Normality assessment and homogeneity assessment;
5. Short-term and long-term stability assessment;
6. Generation of certified property values;
7. Application of the certified property values to evaluate MS-based proteomics data quality across laboratories.

## Data Availability

The large-scale data files used in this study are not included in this GitHub repository.

The following data will be made available through a public repository upon manuscript submission or publication:

- Raw mass spectrometry data;
- Processed identification and quantification results;
- Intermediate result tables;
- Source data used for figure generation;
- Search database files and metadata.

Repository accession numbers or DOIs will be added here when available.

```text
Raw and processed proteomics data: [to be added]
Source data: [to be added]
