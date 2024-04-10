# Used to determine 500 genes varying across pseudotime points for gene filtering using Slingshot
if(!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager") 
}
BiocManager::install("tradeSeq")

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(mclust, quietly = TRUE)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
input_file <- args[2]

setwd(input_dir)

# Read in raw count matrix

data <- read.csv(input_file, row.names = 1)


sce <- SingleCellExperiment(assays = list(counts = as.matrix(data)))
rownames(sce) <- rownames(data)
counts <- as.matrix(counts(sce))
# Filter Genes
geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

# Normalize Data
FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

#plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)


reducedDims(sce) <- SimpleList(PCA = rd1)

sce <- slingshot(sce,reducedDim = 'PCA')

summary(sce$slingPseudotime_1)

sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)

topgenes <- c(rownames(ATres[order(ATres$pvalue), ])[1:100])

## Subset original count matrix 

subset_df <- data %>% filter(rownames(data) %in% topgenes)

write.csv(subset_df, paste0(input_dir, "/cluster_top_exp.csv"), row.names = TRUE)
