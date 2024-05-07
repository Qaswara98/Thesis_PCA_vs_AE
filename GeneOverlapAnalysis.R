# Load necessary libraries 
library(stats) # For the hypergeometric test
library(VennDiagram) # For creating Venn diagrams

# Set the working directory 
 setwd("C:/Users/ABDULLAHI HAJI/OneDrive/Documents/Thesis_poject/GSE216738_RAW_COUNTS_Abnormal-AML-50bp")


 
 # Read in normalized counts, top 500 genes from autoencoder and PCA
normCounts <- read.csv("normCounts_res.CSV")
top3000ae <- read.csv("top_3000_genes_from_autoencoder.csv")
top3000pca <- read.csv("top_3000_genes_from_PCA.csv")

# Extract gene names from your datasets
ae_genes <- as.character(top3000ae$X)
pca_genes <- as.character(top3000pca$X)

# The universe is the list of all genes in normCounts dataset
universe_genes <- as.character(normCounts$X)

# Calculate the actual overlap
actual_overlap <- length(intersect(ae_genes, pca_genes))

# Calculate expected overlap
total_genes <- length(universe_genes)  # Total number of genes in  universe
expected_overlap <- (length(ae_genes) * length(pca_genes)) / total_genes

# Perform the hypergeometric test
# The parameters are: number of shared genes, size of the first gene list,
# the difference between the universe size and the size of the first gene list,
# and the size of the second gene list.
p_value <- phyper(actual_overlap - 1, length(ae_genes), total_genes - length(ae_genes), length(pca_genes), lower.tail=FALSE)

# Print results
cat("Actual Overlap:", actual_overlap, "\n")
cat("Expected Overlap:", expected_overlap, "\n")
cat("P-value:", p_value, "\n")


# Interpret the results
if (p_value < 0.05) {
  cat("The overlap is statistically significant.\n")
} else {
  cat("The overlap is not statistically significant.\n")
}

# Create the Venn diagram for the overlap genes 
venn.plot <- venn.diagram(
  x = list( AE= top3000ae$X, PCA = top3000pca$X),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.cex = 2,
  cat.col = c("red", "blue")
)

# Display the plot
grid.draw(venn.plot)
