#!/usr/bin/env nextflow

params.input = "" // default input file
params.outdir = "Thesis_PCA_vs_AE"
inputFile = Channel.fromPath(params.input)

process DataPreprocessing {
    publishDir "${params.outdir}/DataPreprocessing", mode: 'copy'

    input:
    path inputFile

    output:
    path 'normCounts_res.CSV'

    script:
    """
    #!/usr/bin/env Rscript
    # Load necessary R packages for data manipulation, analysis, and visualization 
    library(dplyr)
    library(tidyr)
    library(GEOquery) 
    library(DESeq2) 
    library(edgeR)
    library(openxlsx) 
    library(biomaRt)
    # Set timeout to 8000 seconds for the BioMart database
    options(timeout=8000) 

    # Load  count data
    counts <- read.table("$inputFile", header = TRUE, as.is = TRUE, row.names = 1")
    dim(counts)
    # Set up the BioMart database
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    # Define the attributes and filters for the query
    attributes <- c("ensembl_gene_id", "gene_biotype")
    filters <- "ensembl_gene_id"

    # Get the gene IDs from  count table
    gene_ids <- rownames(counts)
    head(gene_ids)
    # Query the database to get gene biotypes
    gene_data <- getBM(attributes=attributes, filters=filters, values=gene_ids, mart=mart)


    head(gene_data)

    # Count the number of each gene biotype
    gene_counts <- table(gene_data$gene_biotype)
    print(gene_counts)


    # Define the gene biotypes to be filtered out
    filter_biotypes <- c("artifact","TEC","pseudogene", "processed_pseudogene", "unprocessed_pseudogene", 
                         "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", 
                         "unitary_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", 
                         "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "translated_processed_pseudogene",
                          "rRNA_pseudogene", "rRNA", "Mt_rRNA", "Mt_tRNA", "snRNA", "snoRNA", "misc_RNA","transcribed_unitary_pseudogene","IG_pseudogene")

    # Get the gene biotypes from  gene data
    gene_biotypes <- gene_data$gene_biotype

    # Filter out the specified gene biotypes
    filtered_gene_data <- gene_data[!gene_biotypes %in% filter_biotypes, ]
    dim(filtered_gene_data)
    # Filter out the counts of the filtered genes
    counts_filtered <- counts[rownames(counts) %in% filtered_gene_data$ensembl_gene_id, ]
    dim(counts_filtered )
    # Get the unique biotypes in the filtered data
    remaining_biotypes <- unique(filtered_gene_data$gene_biotype)
    print(remaining_biotypes)
    # Print the remaining biotypes
    #print(paste("The remaining biotypes after filtering are: ", 
            #paste(remaining_biotypes, collapse = ", ")))

    # Fetching associated metadata from GEO
    gse <- getGEO("GSE216738", GSEMatrix = TRUE, getGPL = FALSE)

    # Access the phenotype data
    metadata <- pData(gse[[1]])

    # Manipulate the metadata dataframe to keep only the columns of interest and rename them
    metadata.subset <- metadata[, c(1, 46, 47, 48)]
    colnames(metadata.subset) <- c("samples", "eln_group", "Lncrna_score", "tissue")


    # Check if the dimensions of counts_filtered and metadata.subset are identical 
    length(rownames(metadata.subset)) == length(colnames(counts_filtered))


    # Excluding 'gene_name' and 'gene_biotype' columns
    counts_filtered <- counts_filtered[, !(names(counts_filtered) %in% c("gene_name", "gene_biotype"))]

    # Setting the row names of the metadata to the matching sample names
    rownames(metadata.subset) <- metadata.subset$samples


    # Ensuring the order of metadata matches the order of samples in counts_filtered
    # Extracting the sample names from the column names of the counts_filtered
    counts_filtered_sample_names <- colnames(counts_filtered)

    # Subset the metadata to retain only rows that match the counts_filtered sample names, in the order of counts_filtered
    metadata.subset_aligned <- metadata.subset[counts_filtered_sample_names, ]

    # Find the most frequent category for eln_group
    most_frequent_eln_group <- names(sort(table(metadata.subset_aligned$eln_group), decreasing = TRUE))[1]

    # Impute missing values with the most frequent category
    metadata.subset_aligned$eln_group[is.na(metadata.subset_aligned$eln_group)] <- most_frequent_eln_group

    # Find the most frequent category for Lncrna_score
    most_frequent_Lncrna_score <- names(sort(table(metadata.subset_aligned$Lncrna_score), decreasing = TRUE))[1]

    # Impute missing values with the most frequent category
    metadata.subset_aligned$Lncrna_score[is.na(metadata.subset_aligned$Lncrna_score)] <- most_frequent_Lncrna_score
    # Filter out lowly expressed genes
    meanLog2CPM <- rowMeans(log2(cpm(counts_filtered) + 1))

    hist(meanLog2CPM)
    sum(meanLog2CPM <= 1)
    counts_filtered <- counts_filtered[meanLog2CPM > 1, ]
    nrow(counts_filtered)
    dim(counts_filtered)

    # Create a DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata.subset_aligned, design = ~ eln_group)  


    # Normalize the data
    normCounts <- vst(dds)

    # Plot a histogram of the normalized counts
    hist(assay(normCounts))

    # Export normalized counts for training the AE  model and downstream analysis
    countData <- assay(normCounts)
    dim(countData)
    write.csv(countData, file = "normCounts_res.CSV")
    """
}

process AutoencoderTraining {
    publishDir "${params.outdir}/AutoencoderTraining", mode: 'copy'

    input:
    path preprocessedData from DataPreprocessing.out

    output:
    path 'autoencoder_model.keras'
    path 'combined_transformation_matrix_with_genes.csv'

    script:
    """
    #!/usr/bin/env python

    # Importing necessary libraries
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split
    from tensorflow.keras.models import Model
    from tensorflow.keras.layers import Input, Dense, Dropout
    import tensorflow as tf
    from tensorflow.keras.models import load_model
    import scipy.stats
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import seaborn as sns

    # Define the hyperparameters for the autoencoder model
    params = {
        'latent_dim':85,  # The dimensionality of the latent space
        'depth': 2,  # The number of hidden layers in the encoder and decoder
        'first_layer_size': 500,  # The size of the first hidden layer
        'dropout_rate': 0,  # The dropout rate used after each layer
        'epochs': 200,  # The number of epochs to train for
        'batch_size': 80,  # The batch size used in training
        'activation': 'relu',  # The activation function used in each layer
    }

    def load_data(filepath):
    """Load the data from a CSV file."""
    # Read the CSV file at the given filepath and return it as a pandas DataFrame
    data = pd.read_csv(filepath, index_col=0)
    print(data.head())  # Print the first few rows of the DataFrame
    return data
    
    def build_ae(original_dim, params):
    """Build the autoencoder model."""
    # Define the input layer
    input_layer = Input(shape=(original_dim,))
    
    # Define the first hidden layer and apply dropout
    encoded = Dense(params['first_layer_size'], activation=params['activation'])(input_layer)
    encoded = Dropout(params['dropout_rate'])(encoded)
    
    # Add additional hidden layers with dropout
    for _ in range(params['depth'] - 1):
        encoded = Dense(params['first_layer_size'], activation=params['activation'])(encoded)
        encoded = Dropout(params['dropout_rate'])(encoded)
    
    # Define the latent space layer and the output layer
    encoded = Dense(params['latent_dim'], activation=params['activation'])(encoded)
    decoded = Dense(original_dim, activation=params['activation'])(encoded)
    
    # Define the autoencoder model and compile it
    autoencoder = Model(input_layer, decoded)
    autoencoder.compile(optimizer='adam', loss='mse')
    
    # Print the summary of the model
    autoencoder.summary()
    return autoencoder

    def train_ae(autoencoder, X_train, X_val, params):
    """Train the autoencoder model."""
    # Train the model using the training data (75% of the data) and validate it using the validation data(25% of the data)
    history = autoencoder.fit(X_train, X_train,
                              epochs=params['epochs'],
                              batch_size=params['batch_size'],
                              shuffle=True,
                              validation_data=(X_val, X_val))
    return autoencoder, history

    def plot_loss(history):
    """Plot the training and validation loss."""
    # Create a new figure
    plt.figure()
    
    # Plot the training loss and validation loss
    plt.plot(history.history['loss'], label='Train Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    
    # Add a title and labels to the plot
    plt.title('Loss over epochs')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend()
    #save plot
    plt.savefig('training loss and validation loss plot.png')
    # Display the plot
    plt.show()
   
   def plot_reconstruction_error(autoencoder, X_test):
    """Plot the reconstruction error on the test set."""
    # Use the autoencoder to reconstruct the test data
    reconstructions = autoencoder.predict(X_test)
    
    # Calculate the reconstruction error
    reconstruction_error = np.mean(np.abs(reconstructions - X_test), axis=1)
    
    # Create a new figure
    plt.figure()
    
    # Plot a histogram of the reconstruction error
    plt.hist(reconstruction_error, bins=50)
    
    # Add a title and labels to the plot
    plt.title('Reconstruction error histogram')
    plt.xlabel('Reconstruction error')
    plt.ylabel('Number of examples')
    
    # Display the plot
    plt.show()
    
    def extract_and_multiply_weights(autoencoder):
    """Extracts weights from the trained autoencoder and computes the combined transformation matrix."""
    # Extract weights
    weight_1 = autoencoder.layers[1].get_weights()[0]  # Input to First Hidden Layer
    weight_2 = autoencoder.layers[3].get_weights()[0]  # First Hidden Layer to Second Hidden Layer
    weight_3 = autoencoder.layers[5].get_weights()[0]  # Second Hidden Layer to Latent Layer
    
    # Perform matrix multiplication
    combined_transformation = np.dot(np.dot(weight_1, weight_2), weight_3)
    return combined_transformation

    def main():
        \"\"\"Main function to run the pipeline.\"\"\"
        # Define the filepath of the data
        filepath = os.path.join(DATA_DIR, 'normCounts_res.CSV')
        
        # Load the data
        data = load_data(filepath)
        
        # Get the original dimension of the data (number of genes)
        original_dim = data.shape[0]  # Corrected to shape[1] as we need the number of features, not samples
        
        # Transpose the data because the neural network expects samples as rows
        data_transposed = data.transpose()
        
        # Split the data into training and validation sets (75%, 25%)
        X_train, X_val = train_test_split(data_transposed, test_size=0.25, random_state=42)
        
        # Build the autoencoder model
        autoencoder = build_ae(original_dim, params)
        
        # Train the autoencoder model
        trained_autoencoder, history = train_ae(autoencoder, X_train, X_val, params)

        # Create a separate encoder model from the trained autoencoder
        encoder = Model(inputs=trained_autoencoder.input, outputs=trained_autoencoder.layers[-2].output)

        # Transform the input data into the latent space to get the latent features
        latent_features = encoder.predict(data_transposed)

        # Save the latent features to a CSV file
        pd.DataFrame(latent_features, index=data_transposed.index).to_csv('latent_features.csv')

        # Plot the loss curves and reconstruction error
        plot_loss(history)
        plot_reconstruction_error(trained_autoencoder, X_val)
        
        # Save the trained autoencoder and encoder models
        trained_autoencoder.save('autoencoder_model.keras')
        encoder.save('encoder_model.keras')
        
        combined_transformation_matrix = extract_and_multiply_weights(trained_autoencoder)
        
        print("Combined Transformation Matrix Shape:", combined_transformation_matrix.shape)
        # Convert the numpy array to a pandas DataFrame
        df = pd.DataFrame(combined_transformation_matrix)
        # Extract the combined transformation matrix
        combined_transformation_matrix = extract_and_multiply_weights(trained_autoencoder)
        
        # Convert the numpy array to a pandas DataFrame, using gene names as the index
        combined_transformation_df = pd.DataFrame(combined_transformation_matrix, index=data.index)

        # Save the DataFrame to a CSV file
        combined_transformation_df.to_csv('combined_transformation_matrix_with_genes.csv')

    if __name__ == "__main__":
        main()
    """
}

process IdentifyTopGenes {
    publishDir "${params.outdir}/IdentifyTopGenes", mode: 'copy'

    input:
    path weightsMatrix from AutoencoderTraining.out[1]

    output:
    path 'top_2000_genes_from_autoencoder.csv'

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Load the combined weights matrix
    file_path = '$weightsMatrix'
    weights_df = pd.read_csv(file_path, index_col=0)  # The first column contains gene names

    # Calculate an overall importance score for each gene
    # Here, I'm using the sum of weights across different features.
    gene_importance = weights_df.abs().sum(axis=1)

    # Sort genes by their importance score in descending order
    sorted_genes = gene_importance.sort_values(ascending=False)

    # Select the top 2000 genes
    top_2000_genes = sorted_genes.head(2000)

    # Save the top 2000 genes to a new CSV file for further analysis
    top_2000_genes.to_csv('top_2000_genes_from_autoencoder.csv')
    """
}

process PCAAnalysis {
    publishDir "${params.outdir}/PCAAnalysis", mode: 'copy'

    input:
    path normCounts from DataPreprocessing.out

    output:
    path 'PCA_results' into PCA_results
    output:
    path 'PCA_results' into PCA_results_channel

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Load the data
    data_path = '$normCounts'
    data = pd.read_csv(data_path, index_col=0)

    # Transpose the data so that genes are columns and samples are rows
    data_transposed = data.T

    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data_transposed)

    # Perform PCA
    pca = PCA(n_components=10)
    principalComponents = pca.fit_transform(scaled_data)
    loadings = pca.components_.T

    # Convert to DataFrame for easier handling
    principalDf = pd.DataFrame(data=principalComponents,
                               columns=['PC' + str(i) for i in range(1, 11)])

    # Calculate eigenvalues
    eigenvalues = pca.explained_variance_

    # Create a Scree plot
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(eigenvalues) + 1), eigenvalues, 'ro-', linewidth=2)
    plt.title('Scree Plot')
    plt.xlabel('Principal Component')
    plt.ylabel('Eigenvalue')
    plt.axhline(y=1, color='r', linestyle='--')
    plt.savefig('Scree plot.png')
    plt.close()

    # Calculate and print the explained variance ratio for each principal component
    explained_variance_ratio = pca.explained_variance_ratio_
    with open('explained_variance_ratio.txt', 'w') as f:
        for i, ratio in enumerate(explained_variance_ratio):
            f.write(f"PC{i + 1} explains {ratio * 100:.2f}% of the variance.\\n")

    # Combine the loadings from PCA1, PCA2, PCA3 and PCA4
    combined_loadings = np.abs(loadings[:, 0]) + np.abs(loadings[:, 1]) + np.abs(loadings[:, 2] + np.abs(loadings[:, 3]))

    # Create a DataFrame for easier handling
    loadings_df = pd.DataFrame(combined_loadings, index=data_transposed.columns, columns=['Combined Loadings'])

    # Sort the DataFrame based on the combined loadings
    sorted_loadings_df = loadings_df.sort_values(by='Combined Loadings', ascending=False)

    # Get the top 2000 most influential genes
    top_2000_genes = sorted_loadings_df.head(2000)
    top_2000_genes.to_csv('top_2000_genes_from_PCA.csv')

    # Generate a correlation matrix plot for loadings
    correlation_matrix = np.corrcoef(loadings.T)
    plt.figure(figsize=(12, 10))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f')
    plt.title('Correlation Matrix of Loadings')
    plt.xlabel('Principal Components')
    plt.ylabel('Principal Components')
    plt.savefig('Correlation Matrix of Loadings.png')
    plt.close()
    """

    # Create a directory to store all PCA results
    mkdir -p PCA_results
    mv Scree plot.png PCA_results/
    mv explained_variance_ratio.txt PCA_results/
    mv top_2000_genes_from_PCA.csv PCA_results/
    mv Correlation Matrix of Loadings.png PCA_results/
}

process GeneOverlapAnalysis {
    publishDir "${params.outdir}/GeneOverlapAnalysis", mode: 'copy'

    input:
    path top2000ae from IdentifyTopGenes.out
    path top2000pca from PCAAnalysis.out
    path normCounts from DataPreprocessing.out

    output:
    path 'gene_overlap_results'

    script:
    """
    #!/usr/bin/env Rscript

    # Load necessary libraries 
    library(stats) # For the hypergeometric test
    library(VennDiagram) # For creating Venn diagrams

    # Read in normalized counts, top 2000 genes from autoencoder and PCA
    normCounts <- read.csv('$normCounts')
    top2000ae <- read.csv('$top2000ae')
    top2000pca <- read.csv('$top2000pca')

    # Extract gene names from your datasets
    ae_genes <- as.character(top2000ae$X)
    pca_genes <- as.character(top2000pca$X)

    # The universe is the list of all genes in normCounts dataset
    universe_genes <- as.character(normCounts$X)

    # Calculate the actual overlap
    actual_overlap <- length(intersect(ae_genes, pca_genes))

    # Calculate expected overlap
    total_genes <- length(universe_genes)  # Total number of genes in  universe
    expected_overlap <- (length(ae_genes) * length(pca_genes)) / total_genes

    # Perform the hypergeometric test
    p_value <- phyper(actual_overlap - 1, length(ae_genes), total_genes - length(ae_genes), length(pca_genes), lower.tail=FALSE)

    # Print results
    cat("Actual Overlap:", actual_overlap, "\\n")
    cat("Expected Overlap:", expected_overlap, "\\n")
    cat("P-value:", p_value, "\\n")

    # Interpret the results
    if (p_value < 0.05) {
      cat("The overlap is statistically significant.\\n")
    } else {
      cat("The overlap is not statistically significant.\\n")
    }

    # Create the Venn diagram for the overlap genes 
    venn.plot <- venn.diagram(
      x = list( AE= top2000ae$X, PCA = top2000pca$X),
      filename = NULL,
      fill = c("red", "blue"),
      alpha = 0.5,
      cat.cex = 2,
      cat.col = c("red", "blue")
    )

    # Save the Venn diagram to a file
    png('venn_diagram.png')
    grid.draw(venn.plot)
    dev.off()
    """

    # Create a directory to store all gene overlap results
    mkdir -p gene_overlap_results
    mv venn_diagram.png gene_overlap_results/
}
process GOPathwayAnalysis {
    publishDir "${params.outdir}/GOPathwayAnalysis", mode: 'copy'

    input:
    path top2000ae from IdentifyTopGenes.out
    path top2000pca from PCAAnalysis.out

    output:
    path 'go_pathway_analysis_results'

    script:
    """
    #!/usr/bin/env Rscript

    # Load the necessary library
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(graphite)
    library(DOSE)
    library(enrichplot) 
    library(ReactomePA)
    library(ggplot2)
    library(pathview)

    # Read the files containing the top 2000 genes from autoencoder and PCA
    genes_autoencoder <- read.csv("$top2000ae")
    genes_pca <- read.csv("$top2000pca")

    # Prepare the gene lists
    gene_list_autoencoder <- as.character(genes_autoencoder$X) 
    gene_list_pca <- as.character(genes_pca$X) 

    # Convert ENSEMBL IDs to ENTREZ IDs for both gene lists
    gene_ids_autoencoder <- bitr(gene_list_autoencoder, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    gene_ids_pca <- bitr(gene_list_pca, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    # Check for genes that couldn't be converted in both gene lists
    missing_genes_autoencoder <- gene_list_autoencoder[is.na(gene_ids_autoencoder$ENTREZID)]
    missing_genes_pca <- gene_list_pca[is.na(gene_ids_pca$ENTREZID)]

    if (length(missing_genes_autoencoder) > 0) {
      message(length(missing_genes_autoencoder), " genes from autoencoder data could not be converted to ENTREZ IDs.")
    }

    if (length(missing_genes_pca) > 0) {
      message(length(missing_genes_pca), " genes from PCA data could not be converted to ENTREZ IDs.")
    }

    # Continue with only the genes that could be converted in both gene lists
    converted_gene_ids_autoencoder <- gene_ids_autoencoder$ENTREZID[!is.na(gene_ids_autoencoder$ENTREZID)]
    converted_gene_ids_pca <- gene_ids_pca$ENTREZID[!is.na(gene_ids_pca$ENTREZID)]

    # Function to perform GO analysis
    perform_go_analysis <- function(converted_gene_ids) {
      go_enrichment_results <- enrichGO(gene         = converted_gene_ids,
                                        OrgDb        = org.Hs.eg.db, 
                                        keyType      = "ENTREZID",
                                        ont          = "BP", # for Biological Process
                                        pAdjustMethod = "BH",
                                        qvalueCutoff = 0.05,
                                        readable     = TRUE)
      return(go_enrichment_results)
    }

    # Function to visualize GO analysis results
    visualize_go_analysis <- function(go_enrichment_results, converted_gene_ids) {
      # View the results
      head(go_enrichment_results)
      (head(go_enrichment_results@result))
      
      # Visualize the results
      print(barplot(go_enrichment_results, showCategory=20, font.size = 10) + ggtitle("Biological Process (GO) barplot") +
              theme(plot.title = element_text(hjust = 0.5)))
      
      print(dotplot(go_enrichment_results, showCategory=20, font.size = 10) + ggtitle("Biological Process (GO) dotplot") +
              theme(plot.title = element_text(hjust = 0.5)))
      
      print(cnetplot(go_enrichment_results, gene_ids= converted_gene_ids, font.size = 10) + ggtitle("cnetplot") +
              theme(plot.title = element_text(hjust = 0.5)))
      
      
      pwt<-pairwise_termsim(go_enrichment_results) 
      print(emapplot(pwt,showCategory = 10)+ # Overlapping gene sets will have thicker edges
              ggtitle("Enrichment map")+
              theme(plot.title = element_text(color="black", size=10, face="bold.italic")))
   """

    }

    process NetworkAnalysis {
    publishDir "${params.outdir}/NetworkAnalysis", mode: 'copy'

    input:
    path top2000ae from IdentifyTopGenes.out
    path top2000pca from PCAAnalysis.out

    output:
    path 'network_analysis_results'

    script:
    """
    #!/usr/bin/env Rscript

    # Load necessary libraries
    library(STRINGdb)
    library(igraph)
    library(BioNAR)
    library(randomcoloR)

    # Function to perform operations
    perform_operations <- function(gene_list, string_db) {
      # Create data frame for mapping
      gene_df <- data.frame(query = gene_list)
      
      # Map the genes using the STRINGdb map function
      top_genes_mapped <- string_db$map(gene_df, my_data_frame_id_col_names = "query")
      
      # Retrieve interactions and the networks are already igraph objects
      g <- string_db$get_subnetwork(top_genes_mapped$STRING_id)
      
      # Calculate network metrics (betweenness centrality)
      btwn <- betweenness(g)
      
      # Identify hub genes for network
      degree <- degree(g)
      hub_genes <- names(degree[degree > quantile(degree, 0.9)])
      
      # Community detection
      community <- cluster_louvain(g)
      
      # Calculate degree for each node
      degree <- degree(g)
      
      # Return list of results
      return(list(g = g, btwn = btwn, hub_genes = hub_genes, community = community, degree = degree))
    }

    # Initialize STRINGdb object for human (species = 9606)
    string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400)

    # Read gene lists
    genes_autoencoder <- read.csv("$top2000ae")
    genes_pca <- read.csv("$top2000pca")

    # gene symbols are in the first column
    gene_list_autoencoder <- genes_autoencoder[,1]
    gene_list_pca <- genes_pca[,1]

    # Perform operations for autoencoder and PCA gene lists
    results_autoencoder <- perform_operations(gene_list_autoencoder, string_db)
    results_pca <- perform_operations(gene_list_pca, string_db)

    results_autoencoder$g 
    results_autoencoder$btwn 
    results_autoencoder$hub_genes 
    results_autoencoder$degree 

    # Print additional network metrics
    print(paste("Closeness centrality (Autoencoder):", mean(closeness(results_autoencoder$g), na.rm = TRUE)))
    print(paste("Closeness centrality (PCA):", mean(closeness(results_pca$g), na.rm = TRUE)))
    print(paste("Eigenvector centrality (Autoencoder):", mean(eigen_centrality(results_autoencoder$g)$vector)))
    print(paste("Eigenvector centrality (PCA):", mean(eigen_centrality(results_pca$g)$vector)))
    print(paste("Number of communities (Autoencoder):", max(membership(results_autoencoder$community))))
    print(paste("Number of communities (PCA):", max(membership(results_pca$community))))
    print(paste("Average Path Length (Autoencoder):", mean_distance(results_autoencoder$g, directed = FALSE)))
    print(paste("Average Path Length (PCA):", mean_distance(results_pca$g, directed = FALSE)))
    print(paste("Clustering Coefficient (Autoencoder):", transitivity(results_autoencoder$g)))
    print(paste("Clustering Coefficient (PCA):", transitivity(results_pca$g)))
    print(paste("Modularity (Autoencoder):", modularity(results_autoencoder$community)))
    print(paste("Modularity (PCA):", modularity(results_pca$community)))
    print(paste("Degree Distribution Summary (Autoencoder):"))
    print(summary(degree_distribution(results_autoencoder$g)))
    print(paste("Degree Distribution Summary (PCA):"))
    print(summary(degree_distribution(results_pca$g)))

    # Additional analysis using BioNAR package for autoencoder
    clusters <- calcAllClustering(results_autoencoder$g)
    pFit <- fitDegree(as.vector(igraph::degree(graph=clusters)), threads=1, Nsim=5, plot=TRUE)
    alg = "louvain"
    clusters <- calcCentrality(results_autoencoder$g)
    getCentralityMatrix(clusters)
    clusters <- calcClustering(clusters, alg)
    summary(clusters)
    V(clusters)$louvain

    mem_df <- data.frame(names=V(clusters)$name, membership=as.numeric(V(clusters)$louvain))
    palette <- distinctColorPalette(max(as.numeric(mem_df$membership)))
    lay <- layoutByCluster(clusters, mem_df, layout = layout_nicely)
    idx <- base::match(V(clusters)$name, mem_df$names)
    cgg <- getCommunityGraph(clusters, mem_df$membership[idx])
    D0 = unname(degree(cgg))
    plot(cgg, vertex.size=sqrt(V(cgg)$size), vertex.cex = 0.8, vertex.color=round(log(D0)) + 1, layout=layout_with_kk, margin=0)

    plot(clusters, vertex.size=3, layout=lay, vertex.label=NA, vertex.color=palette[as.numeric(mem_df$membership)], edge.color='grey95')
    legend('topright', legend=names(table(mem_df$membership)), col=palette, pch=19, ncol = 2,cex=0.2)


    # Additional analysis using BioNAR package for PCA
    clusters <- calcAllClustering(results_pca$g)
    pFit <- fitDegree(as.vector(igraph::degree(graph=clusters)), threads=1, Nsim=5, plot=TRUE)
    alg = "louvain"
    clusters <- calcCentrality(results_pca$g)
    getCentralityMatrix(clusters)
    clusters <- calcClustering(clusters, alg)
    summary(clusters)
    V(clusters)$louvain

    mem_df <- data.frame(names=V(clusters)$name, membership=as.numeric(V(clusters)$louvain))
    palette <- distinctColorPalette(max(as.numeric(mem_df$membership)))
    lay <- layoutByCluster(clusters, mem_df, layout = layout_nicely)
    idx <- base::match(V(clusters)$name, mem_df$names)
    cgg <- getCommunityGraph(clusters, mem_df$membership[idx])
    D0 = unname(degree(cgg))
    plot(cgg, vertex.size=sqrt(V(cgg)$size), vertex.cex = 0.8, vertex.color=round(log(D0)) + 1, layout=layout_with_kk, margin=0)

    plot(clusters, vertex.size=3, layout=lay, vertex.label=NA, vertex.color=palette[as.numeric(mem_df$membership)], edge.color='grey95')
    legend('topright', legend=names(table(mem_df$membership)), col=palette, pch=19, ncol = 2)

    # Compare hub genes
    common_hub_genes <- intersect(results_autoencoder$hub_genes, results_pca$hub_genes)
    print(paste("Number of common hub genes:", length(common_hub_genes)))
    if(length(common_hub_genes) > 0) {
      print("Common hub genes:")
      print(common_hub_genes)
    }

    # Compare community structures
    # Note: This is a simple comparison based on the number of communities and their sizes.
    # More sophisticated comparisons could consider the actual composition of the communities.
    autoencoder_comm_sizes <- table(membership(results_autoencoder$community))
    pca_comm_sizes <- table(membership(results_pca$community))
    print(paste("Number of communities in Autoencoder network:", length(autoencoder_comm_sizes)))
    print(paste("Number of communities in PCA network:", length(pca_comm_sizes)))
    print("Sizes of communities in Autoencoder network:")
    print(autoencoder_comm_sizes)
    print("Sizes of communities in PCA network:")
    print(pca_comm_sizes)


    # Load necessary libraries
    library(clusterProfiler)
    library(org.Hs.eg.db)

    # Extract protein identifiers from STRING IDs
    string_ids <- results_autoencoder$hub_genes
    string_ids <-results_pca$hub_genes
    protein_ids <- sapply(strsplit(string_ids, "[.]"), "[", 2)

    # Convert protein identifiers to gene symbols
    gene_symbols <- mapIds(org.Hs.eg.db, keys=protein_ids, column="SYMBOL", keytype="ENSEMBLPROT", multiVals="first")

    # Perform GO enrichment analysis
    # Perform GO enrichment analysis
    go_enrich_results <- enrichGO(gene         = gene_symbols,
                                  OrgDb        = org.Hs.eg.db,
                                  keyType      = "SYMBOL",
                                  ont          = "BP", # Change this to "CC" for Cellular Component or "MF" for Molecular Function
                                  pAdjustMethod = "BH", # Method for p-value adjustment
                                  pvalueCutoff = 0.05, # Adjusted p-value cutoff
                                  qvalueCutoff = 0.2, # q-value cutoff
                                  readable     = TRUE) # Convert gene IDs to readable gene symbols

    # View the results
    head(go_enrich_results)
    dotplot(go_enrich_results)
    
    """  
   
}


workflow {
    inputFile = Channel.fromPath(params.input)
    DataPreprocessing(inputFile)
    AutoencoderTraining(DataPreprocessing.out)
    IdentifyTopGenes(AutoencoderTraining.out[1])
    PCAAnalysis(DataPreprocessing.out)
    GeneOverlapAnalysis(IdentifyTopGenes.out, PCAAnalysis.out, DataPreprocessing.out)
    GOPathwayAnalysis(IdentifyTopGenes.out, PCAAnalysis.out)
    NetworkAnalysis(IdentifyTopGenes.out, PCAAnalysis.out)
}

