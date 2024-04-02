#Change the working directory
setwd("C:/Users/ABDULLAHI HAJI/OneDrive/Documents/Thesis_poject/GSE216738_RAW_COUNTS_Abnormal-AML-50bp")
# Load libraries needed in this task 
library(dplyr)
library(GEOquery)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
# Import the counts into a data frame
countTable<- read.table("GSE216738_RAW_COUNTS_Abnormal-AML-50bp.txt", header = TRUE,as.is = TRUE,row.names = 1,sep="\t")

# Inspecting the class of the data
class(countTable)
# Inspecting the dimensions of the countTable
dim(countTable)
# Display the first 10 rows and 5 columns of the countTable
countTable[1:10, 1:5]
# Get metadata 
gse<- getGEO("GSE216738", GSEMatrix =TRUE, getGPL=FALSE)

# Access the phenotype data
metadata<- pData(gse[[1]])
# View the metadata
View(metadata)
# Manipulate the metadata dataframe to keep only the columns of interest and rename them
metadata.subset <- metadata %>%
  dplyr::select(1, 8, 46, 47, 48) %>%
  rename(
    samples = title,
    source_name = source_name_ch1,
    eln_group = `eln_group:ch1`,
    Lncrna_score = `lncrna_score:ch1`,
    tissue = `tissue:ch1`
  )
# Inspecting the metadata.subset 
head(metadata.subset)

# Check if the dimensions of countTable and metadata.subset are identical 
length(rownames(metadata.subset)) == length(colnames(countTable))
# Checking the number of samples in metadata.subset
print(paste("Number of samples in metadata.subset: ", nrow(metadata.subset)))

# Checking  the number of samples in count data
print(paste("Number of samples in count data: ", ncol(countTable)))
# Check the first few sample names in count data
print(head(colnames(countTable)))
# Inspecting  the first few sample names in metadata.subset
print(head(rownames(metadata.subset)))
# Excluding 'gene_name' and 'gene_biotype' from the countTable column names
sample_ids<- setdiff(colnames(countTable), c("gene_name", "gene_biotype"))
# Checking if the lengths are the same now
length(rownames(metadata.subset)) == length(sample_ids)
# Creating a named vector where names are 'sample_ids' and values are 'rownames(metadata.subset)'
id_mapping <- setNames(rownames(metadata.subset), sample_ids)
# Renaming the sample columns in 'countTable'
colnames(countTable)[match(sample_ids, colnames(countTable))] <- id_mapping
# Inspecting the first few rows and columns of 'countTable'
head(countTable)
# Check if all column names in 'countTable' are present in 'metadata.subset'
all(colnames(countTable) %in% rownames(metadata.subset))
# Check the first few column names in countTable
print(head(colnames(countTable)))

# Exclude 'gene_name' and 'gene_biotype' from the check
all(setdiff(colnames(countTable), c("gene_name", "gene_biotype")) %in% rownames(metadata.subset))
# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countTable[, -c(1,2)], colData = metadata.subset, design = ~ 1)# i used this beacuse i could not use any other parameter for the desing  
dds <- estimateSizeFactors(dds)
dds <- vst(dds)


