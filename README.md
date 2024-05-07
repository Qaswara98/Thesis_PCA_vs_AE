# Comparative Analysis of Autoencoder and PCA Based Dimensionality Reduction Techniques for Gene Expression Data

**Welcome to the comprehensive README Documentation and User Manual.** This repository contains the code necessary to replicate the results of a thesis project, which compares Autoencoder-Based and PCA-Based Dimensionality Reduction Techniques for Gene Expression Data. It also includes the script for the AutoGeneReducer tool, designed for performing dimensionality reduction on gene expression datasets.

## Contents of this repository 

- `DataPreprocessing.R`: A script for preprocessing the data, performed prior to training the model and implementing PCA.

- `Thesis_Intro_AE_vs_PCA.ipynb`: A comprehensive guide illustrating the comparison between AE and PCA methods. It provides a detailed walkthrough of their implementation, with a special focus on the extraction of significant features from both methods.

- `GeneOverlapAnalysis.R`: A script designed to identify overlapping genes and subsequently perform a statistical analysis. The purpose of this analysis is to evaluate the significance of the overlap between genes.
  
- `Annotation Enrichment Analysis.R`: This script performs Gene Ontology (GO) and Reactome pathway enrichment analysis on both the top 3000 genes from PCA and the AE Model  

- `AutoGeneReducer.py`: A script encapsulating the functionality of the AutoGeneReducer tool. It is designed to train models for dimensionality reduction on gene expression data and can be conveniently operated via the command line.

## Getting Started

For a quick start, it is recommended to first visit the ‘INTRODUCTION’ section to gain an understanding of the project. Following this, proceed to the ‘Installation Instructions’ section, where detailed instructions are provided for a smooth installation of the Tool on your system. Upon successful installation, navigate to the ‘USE’ section. This section offers a succinct tutorial on the basic operations and functionalities of the Tool, enabling immediate and effective usage.
# Introduction
Soon to be updated ....

## Get Data

The dataset used for this project is sourced from the Gene Expression Omnibus (GEO) database. It is identifiable under the accession number `GSE216738`. The specific file used for the analysis is `GSE216738_RAW_COUNTS_Abnormal-AML-50bp.txt.gz`.

For direct access to the dataset, please visit the [GEO database GSE216738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216738).
# Reproducibility of the Results
Soon to be updated ....

# AutoGeneReducer tool


## Features
- **Support for Various File Formats**: The tool can process data in CSV, XLS, or XLSX format.
- **Customizable Model Architecture**: Users can adjust the depth, width, and dropout rate of the neural network.
- **Visualization Tools**: The tool provides the ability to plot training loss and reconstruction error histograms, helping users evaluate model performance.
- **Flexible Data Handling**: The tool can work with separate training and validation datasets.
- **Optional Pre-trained Model Loading**: The tool supports loading a pre-trained model from a specified file path. This allows for model reusability and fine-tuning on new data without starting from scratch.

## Dependencies
- Python 3.6 or later
- Pandas
- Numpy
- Matplotlib
- scikit-learn
- tensorflow
- Keras

## Installation Instructions
You have several options to get started with the AutoGeneReducer tool:

1. **Clone the Repository**: You can clone the repository to your local machine. This gives you access to all the necessary files. Use the following command to clone the repository:

```bash
git clone git@github.com:Qaswara98/Thesis_PCA_vs_AE.git
```

2. **Download the Script**: Alternatively, you can specifically download the `AutoGeneReducer.py` script from the repository. Please ensure that all dependencies are installed before running the script. The required dependencies are listed in the section above.

3. **Use the Docker Image**: For a hassle-free setup and execution, we provide a Docker image. This image comes with a pre-configured environment and all the dependencies required by the AutoGeneReducer tool, ensuring a consistent runtime across different platforms. You can download the Docker image using the following command:

```sh
docker pull qaswara98/ubuntu:AutoGeneReducer
```
## Usage
The script is executed from the command line with the following syntax:
```bash
python AutoGeneReducer.py [arguments]
```
## Arguments
- `data_filepath`: This is the path to the data file you want to use. The data file should be in CSV, XLS, or XLSX format.
- `val_data_filepath`: This is an optional argument. If you have a separate validation data file, you can specify the path here. If not provided, the script will split the training data into training and validation sets.
- `model_path`: This is an optional argument. If you have a pre-trained model that you want to load, you can specify the path here.
- `latent_dim`: This is the dimension of the latent space of the autoencoder. The default value is 20.
- `depth`: This is the depth of the neural network, i.e., the number of layers. The default value is 2.
- `first_layer_size`: This is the size of the first layer of the neural network. The default value is 500.
- `dropout_rate`: This is the dropout rate for training the neural network. The default value is 0.1.
- `epochs`: This is the number of epochs for training the model. The default value is 200.
- `batch_size`: This is the batch size for training the model. The default value is 80.
- `activation`: This is the activation function for the neural network layers. The default value is 'relu'.
- `test_size`: This is the proportion of the dataset to include in the test split. The default value is 0.25.
- `shuffle`: This determines whether to shuffle the training data before each epoch. The default value is True.
- `plot_loss`: If this argument is included, the script will plot the loss curves after training.
- `plot_reconstruction_error`: If this argument is included, the script will plot the reconstruction error histogram after training.

### Running the Script
In this section, you'll find examples demonstrating how to execute the script with various arguments.
### Basic Command
Run the script with minimal setup using the command below. This will run the script with all default parameters:
```bash
python AutoGeneReducer.py path/to/your/data.csv
```
### Advanced Options
Customize the script's behavior using various options. For example:
```bash
python AutoGeneReducer.py data.csv --latent_dim 32 --epochs 150 --plot_loss
```
You can execute the script with all arguments as shown below:
```bash
python AutoGeneReducer.py --data_filepath /path/to/data.csv --val_data_filepath /path/to/validation_data.csv --model_path /path/to/model.keras --latent_dim 20 --depth 2 --first_layer_size 500 --dropout_rate 0.1 --epochs 200 --batch_size 80 --activation relu --test_size 0.25 --shuffle True --plot_loss --plot_reconstruction_error
```
> # Note
Please replace `/path/to/data.csv`, `/path/to/validation_data.csv`, and `/path/to/model.keras` with the actual paths to your files. If you want to use the default values for any of the optional arguments, you can simply omit them from the command.

## Running the AutoGeneReducer Tool Using Docker
### Data Preparation
Please make sure that your training files are either in the working directory or a specified directory on your host system. To make these files accessible to AutoGeneReducer, you'll need to mount this directory to the Docker container.

```sh
docker run -v /path/to/your/host/data:/data/in/container qaswara98/ubuntu:AutoGeneReducer python3 AutoGeneReducer.py --data_filepath /data/in/container/data.csv --val_data_filepath /data/in/container/validation_data.csv --model_path /data/in/container/model.keras --latent_dim 20 --depth 2 --first_layer_size 500 --dropout_rate 0.1 --epochs 200 --batch_size 80 --activation relu --test_size 0.25 --shuffle True --plot_loss --plot_reconstruction_error
```
> # Note
Please replace /path/to/your/host/data and /data/in/container with your actual paths. Also, ensure that the Docker image qaswara98/ubuntu:AutoGeneReducer is available on your system or Docker Hub. 

## Troubleshooting

- **File Format Error**: Ensure your file is correctly formatted. The script only accepts CSV, XLS, or XLSX files.
- **Dependency Errors**: Verify all required Python libraries are installed.
- **Model Performance**: If the model performance is not satisfactory, consider adjusting the hyperparameters such as  `--latent_dim`, `--epochs`, or `--batch_size`.

## License

The project is released under the GPL-3.0 license. For more details on the licensing, please See the [LICENSE](LICENSE) file. 

