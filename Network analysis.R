# Load necessary libraries
library(STRINGdb)
library(igraph)
library(BioNAR)
library(randomcoloR)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
# sets timeout limit to 300 seconds
options(timeout = 2000)  


# Set the working directory
setwd("C:/Users/ABDULLAHI HAJI")

# Function to perform operations
perform_operations <- function(gene_list, string_db) {
  # Create data frame for mapping
  gene_df <- data.frame(query = gene_list)
  
  # Map the genes using the STRINGdb map function
  top_genes_mapped <- string_db$map(gene_df, my_data_frame_id_col_names = "query")
  
  # Retrieve interactions and the networks are already igraph objects
  g <- string_db$get_subnetwork(top_genes_mapped$STRING_id)
  
  # Extract Largest Connected Component
  components <- components(g)
  largest_component <- which.max(components$csize)
  g <- induced_subgraph(g, which(components$membership == largest_component))
  
  # Calculate network metrics (betweenness centrality, degree, and hub genes within LCC)
  btwn <- betweenness(g)
  degree <- degree(g)
  hub_genes <- names(degree[degree > quantile(degree, 0.9)])
  community <- cluster_louvain(g)
  
  # Return list of results
  return(list(g = g, btwn = btwn, hub_genes = hub_genes, community = community, degree = degree))
}

# Initialize STRINGdb object for human (species = 9606)
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400)


# Read gene lists files 
genes_autoencoder <- read.csv("top_3000_genes_from_autoencoder.csv")
genes_pca <- read.csv("top_3000_genes_from_pca.csv")

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

# Plots a network of hub genes identified by the autoencoder
string_db$plot_network(results_autoencoder$hub_genes)
# Plots a network of hub genes identified by the PCA
string_db$plot_network(results_pca$hub_genes)

# Identify common hub genes between PCA and AE 
common_hub_genes <- intersect(results_autoencoder$hub_genes, results_pca$hub_genes)
print(paste("Number of common hub genes:", length(common_hub_genes)))
print(common_hub_genes)

# Plots a network of common hub genes identified by the autoencoder
string_db$plot_network(common_hub_genes)

# Compare community structures
#  This is a simple comparison based on the number of communities and their sizes.
autoencoder_comm_sizes <- table(membership(results_autoencoder$community))
pca_comm_sizes <- table(membership(results_pca$community))
print(paste("Number of communities in Autoencoder network:", length(autoencoder_comm_sizes)))
print(paste("Number of communities in PCA network:", length(pca_comm_sizes)))
print("Sizes of communities in Autoencoder network:")
print(autoencoder_comm_sizes)
print("Sizes of communities in PCA network:")
print(pca_comm_sizes)

# Extract protein identifiers from STRING IDs for the hub genes
string_ids <- results_autoencoder$hub_genes
string_ids <-results_pca$hub_genes

# Extract protein identifiers from STRING IDs for the community
string_ids <- results_autoencoder$community
string_ids <- results_pca$community

# Extract protein identifiers from STRING IDs FOR the common hub genes
string_ids_common_hub_genes <- common_hub_genes

string_ids <- common_hub_genes

protein_ids <- sapply(strsplit(string_ids, "[.]"), "[", 2)


# Convert protein identifiers to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys=protein_ids, column="SYMBOL", keytype="ENSEMBLPROT", multiVals="first")

# Perform GO enrichment analysis
go_enrich_results <- enrichGO(gene         = gene_symbols,
                              OrgDb        = org.Hs.eg.db,
                              keyType      = "SYMBOL",
                              ont          = "BP",
                              pAdjustMethod = "BH", # Method for p-value adjustment
                              pvalueCutoff = 0.05, # Adjusted p-value cutoff
                              qvalueCutoff = 0.2, # q-value cutoff
                              readable     = TRUE) # Convert gene IDs to readable gene symbols

# View the results
head(go_enrich_results)
dotplot(go_enrich_results)


# Convert protein identifiers to Entrez IDs for Pathway analysis 
entrez_ids <- mapIds(org.Hs.eg.db, keys=protein_ids, column="ENTREZID", keytype="ENSEMBLPROT", multiVals="first")

# Remove NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# Now use the converted Entrez IDs in enrichKEGG
kegg <- enrichKEGG(gene         = entrez_ids,
                   organism     = 'hsa', # for Homo sapiens
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)


head(kegg)

barplot(kegg)

# Perform Reactome pathway enrichment analysis
reactome_res <- enrichPathway(gene =  entrez_ids,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",  # for Benjamini & Hochberg (BH) adjustment
                              organism = "human")

# View the results
head(reactome_res)


# For Reactome results
barplot(reactome_res, showCategory=20)


# Compare community structures
# This is a simple comparison based on the number of communities and their sizes.
autoencoder_comm_sizes <- table(membership(results_autoencoder$community))
pca_comm_sizes <- table(membership(results_pca$community))
print(paste("Number of communities in Autoencoder network:", length(autoencoder_comm_sizes)))
print(paste("Number of communities in PCA network:", length(pca_comm_sizes)))
print("Sizes of communities in Autoencoder network:")
print(autoencoder_comm_sizes)
print("Sizes of communities in PCA network:")
print((pca_comm_sizes))

# More sophisticated comparison considering the actual composition of the communities
autoencoder_communities <- split(V(results_autoencoder$g)$name, membership(results_autoencoder$community))
pca_communities <- split(V(results_pca$g)$name, membership(results_pca$community))


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



string_ids_autoencoder <- results_autoencoder$hub_genes
protein_ids_autoencoder <- sapply(strsplit(string_ids_autoencoder, "[.]"), "[", 2)
write.csv(protein_ids_autoencoder, "hub_genes_autoencoder.csv")

string_ids__pca<-results_pca$hub_genes
protein_ids_pca <- sapply(strsplit(string_ids__pca, "[.]"), "[", 2)
write.csv(protein_ids_pca, "hub_genes_pca.csv")

string_ids__common<-common_hub_genes
protein_ids_common <- sapply(strsplit(string_ids__common, "[.]"), "[", 2)
write.csv(protein_ids_common, "common_hub_genes.csv")

