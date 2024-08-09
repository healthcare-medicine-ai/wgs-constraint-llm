# Unified Genome Constraint, Pathogenicity, and pLoF Model

This repository contains the code and notebooks used for the analysis presented in the paper titled [A Unified Genome Constraint, Pathogenicity, and pLoF Model Identifies New Genes Associated with Epilepsy](https://www.medrxiv.org/content/10.1101/2024.06.27.24309590v2.full-text). The repository is organized into several Jupyter notebooks that correspond to different aspects of the analysis, as well as a results folder containing key outputs.

## Repository Structure
- **AoU Mutation Predictions.ipynb**:
This notebook generates mutation predictions for the All of Us (AoU) dataset. The analysis is part of the unified model described in the paper, focusing on the prediction of pathogenicity across the genome.

- **Constraint vs gnomAD Analysis.ipynb**:
This notebook performs a comparative analysis between the constraint metrics derived from the model and the gnomAD dataset. The analysis aims to validate the model's constraint predictions by comparing them with observed mutation frequencies in gnomAD.

- **Schizophrenia Analysis.ipynb**:
Experimental: This notebook contains preliminary analysis related to schizophrenia. While not included in the paper, it explores potential associations between the unified model's predictions and schizophrenia-related genes.

- **Autism Analysis.ipynb**:
Experimental: Similar to the schizophrenia notebook, this one explores associations between the model's predictions and autism-related genes. This analysis is also experimental and not part of the published paper.

- **Epilepsy Analysis.ipynb**:
This notebook provides the detailed analysis of epilepsy-related genes, as described in the paper. It includes the final methodology and results that were validated and presented in the publication.

- **Epilepsy Analysis - Old.ipynb**:
Deprecated: This is an older version of the epilepsy analysis that is no longer relevant. It is retained for reference but should not be used.

- **RGC + AoU Predictions.ipynb**:
This notebook combines predictions from the RGC (Rare Genomics Consortium) dataset with those from the AoU dataset. The combined analysis provides insights into rare variants' pathogenicity, integrating data across multiple cohorts.

- **Constraint + AM Analysis.ipynb**:
This notebook provides a detailed analysis of constraint metrics in combination with additional annotations (AM). It explores how these metrics correlate with known pathogenicity markers and contributes to the broader analysis framework described in the paper.

- **Bipolar Analysis.ipynb**:
Experimental: This notebook contains an exploratory analysis related to bipolar disorder. Like the schizophrenia and autism analyses, it is experimental and not included in the paper.

- **results/**:
This folder contains select results from the analysis that will be published. It includes key outputs, figures, and tables that support the findings presented in the paper.

## Usage
To run the notebooks, you will need Python and Jupyter Notebook installed. You can install the required Python packages by running the setup cells included in the beginning of each notebook.

Each notebook is self-contained and includes comments and explanations for the code and analysis. Start with the relevant notebook based on the analysis you are interested in, and follow the steps described to reproduce the results.

The data used by the notebook is stored locally (not included in this repository). After obtaining the necessary files, you should redirect the paths to the data defined at the beginning of the notebook.

## Contact
For questions or issues related to this repository, please contact Oscar Aguilar at osthoag[at]stanford[dot]edu.
