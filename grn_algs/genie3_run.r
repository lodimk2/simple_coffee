library(GENIE3)
library(dplyr)
library(Seurat)
library(doRNG)

args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
#out_dir <- args[2]
setwd(input_dir)
#print(cat("Working Directory:", input_dir))
print(getwd())
data_mat <- read.csv(paste0(getwd(), "/cluster_top_exp.csv"), row.names = 1)

seurat_obj <- CreateSeuratObject(counts = data_mat)

# Normalize gene expression
seurat_obj <- NormalizeData(seurat_obj)


log_transformed_matrix <- as.data.frame(GetAssayData(seurat_obj, slot = "data"))
geneNames <- rownames(log_transformed_matrix)
rownames(log_transformed_matrix) <- c(geneNames)


weightMat <- GENIE3(as.matrix(log_transformed_matrix), nCores=8, verbose=TRUE)

linkList <- as.data.frame(getLinkList(weightMat))


# Check if the directory exists
if (!file.exists("grn_networks")) {
  # If not, create the directory
  dir.create("grn_networks")
} else {
  cat(paste("Directory '", "/grn_networks", "' already exists.\n", sep = ""))
}


##inkList$Gene1 <- gsub("\"", "", linkList$Gene1)
#linkList$Gene2 <- gsub("\"", "", linkList$Gene2)

write.table(linkList, paste0(getwd(), "/grn_networks/genie3_network.csv"), sep = "\t", row.names = FALSE, quote = FALSE)

