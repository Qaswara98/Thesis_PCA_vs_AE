# Comparative Analysis of Autoencoder and PCA Based Dimensionality Reduction Techniques for Gene Expression Data

**Welcome to the comprehensive README Documentation and User Manual.** This repository contains the code necessary to replicate the results of a thesis project, which compares Autoencoder-Based and PCA-Based Dimensionality Reduction Techniques for Gene Expression Data. It also includes the script for the AutoGeneReducer tool, designed for performing dimensionality reduction on gene expression datasets.

## Contents of this repository 

- `DataPreprocessing.R`: A script for preprocessing the data, performed prior to training the model and implementing PCA.

- `Thesis_Intro_AE_vs_PCA.ipynb`: A comprehensive guide illustrating the comparison between AE and PCA methods. It provides a detailed walkthrough of their implementation, with a special focus on the extraction of significant features from both methods.

- `GeneOverlapAnalysis.R`: A script designed to identify overlapping genes and subsequently perform a statistical analysis. The purpose of this analysis is to evaluate the significance of the overlap between genes.
  
- `GeneEnrichmentAnalysis.R`: This script performs Gene Ontology (GO) and Reactome pathway enrichment analysis on both the top 500 genes from PCA and the AE Model  

- `AutoGeneReducer.py`: A script encapsulating the functionality of the AutoGeneReducer tool. It is designed to train models for dimensionality reduction on gene expression data and can be conveniently operated via the command line.

## Getting Started

For a quick start, it is recommended to first visit the ‘INTRODUCTION’ section to gain an understanding of the project. Following this, proceed to the ‘INSTALL’ section, where detailed instructions are provided for a smooth installation of the Tool on your system. Upon successful installation, navigate to the ‘USE’ section. This section offers a succinct tutorial on the basic operations and functionalities of the Tool, enabling immediate and effective usage.

## Get Data

The dataset used for this project is sourced from the Gene Expression Omnibus (GEO) database. It is identifiable under the accession number `GSE216738`. The specific file used for the analysis is `GSE216738_RAW_COUNTS_Abnormal-AML-50bp.txt.gz`.

For direct access to the dataset, please visit the [GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216738).

# Introduction
