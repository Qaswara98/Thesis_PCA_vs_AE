#!/usr/bin/env nextflow

params.input = "" //  input file
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
    counts <- read.table("$inputFile", header = TRUE, as.is = TRUE, row.names = 1)
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
    """
}

workflow {
    inputFile = Channel.fromPath(params.input)
    DataPreprocessing(inputFile)
}
