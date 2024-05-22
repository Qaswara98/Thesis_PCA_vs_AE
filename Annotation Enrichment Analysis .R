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

# Set the working directory 
setwd("C:/Users/ABDULLAHI HAJI")

# Read the files containing the top 3000 genes from autoencoder and PCA
genes_autoencoder <- read.csv("top_3000_genes_from_autoencoder.csv")

genes_pca <- read.csv("top_3000_genes_from_pca.csv")

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
  print(head(g_enrichment_results@result,n=100))
}

# Function to visualize GO analysis results
visualize_go_analysis <- function(go_enrichment_results, converted_gene_ids) {
  # View the results
  head(go_enrichment_results)
  print(head(go_enrichment_results@result, n=100))
  
  # Visualize the results
  print(barplot(go_enrichment_results, showCategory=20, font.size = 10) + ggtitle("Barplot for the Biological Process (GO) Using an PCA Gene List") +
          theme(plot.title = element_text(hjust = 0.5)))
  
  print(dotplot(go_enrichment_results, showCategory=20, font.size = 10) + ggtitle("Dotplot for the Biological Process (GO) Using an PCA Gene List") +
          theme(plot.title = element_text(hjust = 0.5)))
  print(goplot(go_enrichment_results, showCategory=20, font.size = 10) + ggtitle("goplot for the Biological Process (GO) Using an PCA Gene List") +
          theme(plot.title = element_text(hjust = 0.5)))
  
  
  print(cnetplot(go_enrichment_results, gene_ids= converted_gene_ids, font.size = 10) + ggtitle("cnetplot") +
          theme(plot.title = element_text(hjust = 0.5)))
  
  
  pwt<-pairwise_termsim(go_enrichment_results) 
  print(emapplot(pwt,showCategory = 10)+ # Overlapping gene sets will have thicker edges
          ggtitle("Enrichment map")+
          theme(plot.title = element_text(color="black", size=10, face="bold.italic")))
}

# Perform GO analysis for both gene lists
go_enrichment_results_autoencoder <- perform_go_analysis(converted_gene_ids_autoencoder)
go_enrichment_results_pca <- perform_go_analysis(converted_gene_ids_pca)

# Visualize GO analysis results for both gene lists
visualize_go_analysis(go_enrichment_results_autoencoder, converted_gene_ids_autoencoder)
visualize_go_analysis(go_enrichment_results_pca, converted_gene_ids_pca)

# Function to perform Reactome pathway enrichment analysis
perform_reactome_analysis <- function(converted_gene_ids) {
  reactome_res <- enrichPathway(gene =  converted_gene_ids,
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",  # for Benjamini & Hochberg (BH) adjustment
                                organism = "human")
  return(reactome_res)
  print(head(reactome_res@result))
}

# Function to visualize Reactome analysis results
visualize_reactome_analysis <- function(reactome_res, converted_gene_ids) {
  # View the results
  head(reactome_res)
  
  # For Reactome results
  print(barplot(reactome_res, showCategory=20))
  
  # For Reactome results
  print(dotplot(reactome_res, showCategory=20))
  
  # For Reactome results
  #print(cnetplot(reactome_res, gene_ids=converted_gene_ids))
  
  # For Reactome results
  reactome_pwt <- pairwise_termsim(reactome_res)
  print(emapplot(reactome_pwt, showCategory=20))
}

# Perform Reactome pathway enrichment analysis for both gene lists
reactome_res_autoencoder <- perform_reactome_analysis(converted_gene_ids_autoencoder)
reactome_res_pca <- perform_reactome_analysis(converted_gene_ids_pca)

# Visualize Reactome analysis results for both gene lists
visualize_reactome_analysis(reactome_res_autoencoder, gene_ids_autoencoder)
visualize_reactome_analysis(reactome_res_pca, gene_ids_pca)

# Function to perform KEGG pathway analysis
perform_kegg_analysis <- function(converted_gene_ids) {
  kegg <- enrichKEGG(gene         = converted_gene_ids,
                     organism     = 'hsa', # for Homo sapiens
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05)
  return(kegg)
}

# Function to visualize KEGG analysis results
visualize_kegg_analysis <- function(kegg, converted_gene_ids) {
  # View the results
  print(head(kegg))
  
  # Visualize the results
  print(barplot(kegg, showCategory=20))
  
  print(dotplot(kegg, showCategory=20))
  
  #print(cnetplot(kegg, gene_ids=converted_gene_ids))
  
  # Use color.params instead of gene_ids
  print(cnetplot(kegg, color.params = list(gene_ids = converted_gene_ids)))
  
  # pathview of KEGG pathway
  #pathview(gene.data = kegg, pathway.id = "hsa04520")  
  #pathview(gene.data = kegg, pathway.id = "hsa04914")
}

# Perform KEGG pathway analysis for both gene lists
kegg_autoencoder <- perform_kegg_analysis(converted_gene_ids_autoencoder)
kegg_pca <- perform_kegg_analysis(converted_gene_ids_pca)

# Visualize KEGG analysis results for both gene lists
visualize_kegg_analysis(kegg_autoencoder, gene_ids_autoencoder)
visualize_kegg_analysis(kegg_pca, gene_ids_pca)

# check the the intersect between GO terms of PCA and AE 
intersect_GO_terms<-intersect(go_enrichment_results_autoencoder@result$ID[1:1000],go_enrichment_results_pca@result$ID[1:1000])
length(intersect_GO_terms)
print(intersect_GO_terms)
View(go_enrichment_results_autoencoder)
View(go_enrichment_results_pca)


# check the the intersect between GO terms of PCA and AE 
intersect_kegg_terms<-intersect(kegg_autoencoder@result$ID[1:100],kegg_pca@result$ID[1:100])
length(intersect_kegg_terms)
print(intersect_kegg_terms)

View(kegg_autoencoder)
View(kegg_pca)

