# Use PPCOR package to calculate GRN
#setwd('~/Documents/grn_clustering')
library(ppcor)
library(dplyr)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
setwd(input_dir)

data_mat <- read.csv(paste0(getwd(), "/cluster_top_exp.csv"), row.names = 1)

seurat_obj <- CreateSeuratObject(counts = data_mat)

# Normalize gene expression
seurat_obj <- NormalizeData(seurat_obj)


log_transformed_matrix <- as.data.frame(GetAssayData(seurat_obj, slot = "data"))
geneNames <- rownames(log_transformed_matrix)
rownames(log_transformed_matrix) <- c(geneNames)

pcorResults=  pcor(x= t(as.matrix(log_transformed_matrix)), method = "spearman")

# From BEELINE Script 

DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
outDF <- DF[order(DF$corVal, decreasing=TRUE), ]

# Adjust the threshold as needed
#Filter based on significant pvalue of 0.05 


# Check if the directory exists
# Check if the directory exists
if (!file.exists("grn_networks")) {
  # If not, create the directory
  dir.create("grn_networks")
} else {
  cat(paste("Directory '", "/grn_networks", "' already exists.\n", sep = ""))
}
colnames(outDF) <- c("Gene1", "Gene2", "EdgeWeight", "PValue")

outDF <- outDF

write.table(outDF[, c(1,2,3)], paste0(getwd(), "/grn_networks/ppcor_network.csv"), row.names = FALSE, sep = "\t", quote = FALSE)

