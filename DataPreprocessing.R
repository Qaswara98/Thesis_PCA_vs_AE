# Load necessary R packages for data manipulation, analysis, and visualization 
library(dplyr)
library(tidyr)
library(GEOquery) # For fetching gene expression data from GEO
library(DESeq2) # For normalization and differential expression analysis
library(edgeR)  # For normalization and additional RNA-seq data handling
library(openxlsx) # For reading and writing Excel files

# Set the working directory 
setwd("C:/Users/ABDULLAHI HAJI/OneDrive/Documents/Thesis_poject/GSE216738_RAW_COUNTS_Abnormal-AML-50bp")

# Load gene expression count data
countTable <- read.table("GSE216738_RAW_COUNTS_Abnormal-AML-50bp.txt", header = TRUE, as.is = TRUE, row.names = 1, sep = "\t")

# Fetching associated metadata from GEO
gse <- getGEO("GSE216738", GSEMatrix = TRUE, getGPL = FALSE)

# Access the phenotype data
metadata <- pData(gse[[1]])

# Manipulate the metadata dataframe to keep only the columns of interest and rename them
metadata.subset <- metadata[, c(1, 46, 47, 48)]
colnames(metadata.subset) <- c("samples", "eln_group", "Lncrna_score", "tissue")

# Check if the dimensions of countTable and metadata.subset are identical 
length(rownames(metadata.subset)) == length(colnames(countTable))

# Excluding 'gene_name' and 'gene_biotype' columns
countTable <- countTable[, !(names(countTable) %in% c("gene_name", "gene_biotype"))]

# Setting the row names of the metadata to the matching sample names
rownames(metadata.subset) <- metadata.subset$samples

# Ensuring the order of metadata matches the order of samples in countTable
# Extracting the sample names from the column names of the countTable
countTable_sample_names <- colnames(countTable)

# Subset the metadata to retain only rows that match the countTable sample names, in the order of countTable
metadata.subset_aligned <- metadata.subset[countTable_sample_names, ]

# Find the most frequent category for eln_group
most_frequent_eln_group <- names(sort(table(metadata.subset_aligned$eln_group), decreasing = TRUE))[1]

# Impute missing values with the most frequent category
metadata.subset_aligned$eln_group[is.na(metadata.subset_aligned$eln_group)] <- most_frequent_eln_group

# Find the most frequent category for Lncrna_score
most_frequent_Lncrna_score <- names(sort(table(metadata.subset_aligned$Lncrna_score), decreasing = TRUE))[1]

# Impute missing values with the most frequent category
metadata.subset_aligned$Lncrna_score[is.na(metadata.subset_aligned$Lncrna_score)] <- most_frequent_Lncrna_score

# Filter out lowly expressed genes
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1))
hist(meanLog2CPM)
sum(meanLog2CPM <= 1)
countTable <- countTable[meanLog2CPM > 1, ]
dim(countTable)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countTable, colData = metadata.subset_aligned, design = ~ eln_group)  

# Normalize the data
normCounts <- vst(dds)

# Plot a histogram of the normalized counts
hist(assay(normCounts))
# Export normalized counts for training the AE  model and downstream analysis
countData <- assay(normCounts)

write.csv(countData, file = "normCounts_res.CSV")

