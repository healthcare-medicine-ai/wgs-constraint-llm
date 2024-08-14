# Unified Genome Constraint, Pathogenicity, and pLoF Model

This repository contains the code and notebooks used for the analysis presented in the paper titled [A Unified Genome Constraint, Pathogenicity, and pLoF Model Identifies New Genes Associated with Epilepsy](https://www.medrxiv.org/content/10.1101/2024.06.27.24309590v2.full-text). The repository is organized into several Jupyter notebooks that correspond to different aspects of the analysis, as well as a results folder containing key outputs.

## Repository Structure
- **RGC Mutation Predictions.ipynb**:
This notebook generates constraint predictions for the Regeneron Genetics Center Millon Exome (RGC-ME) dataset using a Hidden Markov Model (HMM). The analysis is part of the unified model described in the paper, focusing on the prediction of pathogenicity across the genome.

- **AoU Mutation Predictions.ipynb**:
This notebook generates constraint predictions for the All of Us (AoU) dataset. The predictions are used as an independent WES dataset.

- **RGC + AoU Predictions.ipynb**:
This notebook combines predictions from the RGC dataset with those from the AoU dataset. The combined analysis uses these independent datasets to validate the HMM predictions as a measure of constraint.

- **Constraint vs gnomAD Analysis.ipynb**:
This notebook performs a comparative analysis between the HMM constraint predictions and constraint metrics provided by gnomAD 4.0.

- **Constraint + AM Analysis.ipynb**:
This notebook provides an analysis between the HMM constraint probabilities and pathogenicity predictions from AlphaMissense (AM).

- **Epilepsy Analysis.ipynb**:
This notebook provides the detailed analysis of epilepsy-related genes using HMM constraint probabilities, AM pathogenicity predictions, and information on predicted loss-of-function (pLoF) variants. It includes the final methodology and results that were validated and presented in the publication.

- **results/**:
This folder contains select results from the analysis that will be published. It includes key outputs, figures, and tables that support the findings presented in the paper.

## Usage
To run the notebooks, you will need Python and Jupyter Notebook installed. You can install the required Python packages by running the setup cells included in the beginning of each notebook.

Each notebook is self-contained and includes comments and explanations for the code and analysis. Start with the relevant notebook based on the analysis you are interested in, and follow the steps described to reproduce the results.

The data used by the notebook is stored locally (not included in this repository). After obtaining the necessary files, you should redirect the paths to the data defined at the beginning of the notebook.

## Contact
For questions or issues related to this repository, please contact Oscar Aguilar at osthoag[at]stanford[dot]edu.
