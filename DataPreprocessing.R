# Load necessary libraries
library(dplyr)
library(tidyr)
library(GEOquery) 
library(DESeq2) 
library(edgeR)
library(openxlsx) 
library(biomaRt)
# Set timeout to 8000 seconds for the BioMart database
options(timeout=8000) 
# Set the working directory
setwd("C:/Users/ABDULLAHI HAJI/OneDrive/Documents/Thesis_poject/GSE216738_RAW_COUNTS_Abnormal-AML-50bp")

# Load  count data
counts <- read.table("GSE216738_RAW_COUNTS_Abnormal-AML-50bp.txt", header=TRUE, row.names=1)
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
filter_biotypes <- c("artifact","TEC","pseudogene", "processed_pseudogene", "unprocessed_pseudogene", "TR_C_gene","sRNA","ribozyme","TR_D_gene",
                     "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "vault_RNA","miRNA","scaRNA","IG_D_gene",
                     "unitary_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "lncRNA","IG_V_gene","vault_RNA","TR_V_gene",
                     "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "translated_processed_pseudogene","scRNA","IG_C_gene","IG_J_gene","TR_J_gene",
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

# Export normalized counts for training the AE  model , PCA and downstream analysis
countData <- assay(normCounts)
dim(countData)
write.csv(countData, file = "normCounts_res.CSV")
