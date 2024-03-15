# 14-03-2024
# Differential gene expression analysis using DESeq2 package
# setwd("C:/Users/salon/OneDrive/Desktop/DESeq2")

# Load libraries
library(BiocManager)
BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
BiocManager::install("airway", install_args = "-t 900")  # In case of an error, increase timeout to 15 minutes from 5 mins , i.e, 300s (default)
library(airway)

# Obtain data from airway package
data(airway)
airway
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('CellLine', 'Dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = TRUE, row.names = T, quote = F)
countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# Step 1: Preparing count data -----------------------

# read in sample info
colData <- read.csv('sample_info.csv')

# read in counts data
counts_data <- read.csv('counts_data.csv')

# making sure the row names in colData matches to column names in counts_data
all(rownames(colData) %in% colnames(counts_data))

# are they in the same order?
all(rownames(colData) == colnames(counts_data))

# If not in the same order
# counts_data <-counts_data[, rownames(colData)]


# Step 2: Construct a DESeqDataSet object ------------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                        colData = colData,
                        design = ~Dexamethasone)
dds

# NOTE: Collapse technical replicates, if any. 
# dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design ~ condition + replicate)
# dds <- collapseReplicates(dds, groupBy = "replicate")

# Removing rows (genes) with low gene counts
# Keeping rows that have at least 10 reads in total

dds_10 <- rowSums(counts(dds)) >= 10
dds_filtered <- dds[dds_10,]

# Set the factor level
dds_filtered$Dexamethasone <- relevel(dds_filtered$Dexamethasone, ref = 'untreated')
dds_filtered


# Step 3: Run DESeq ------------------------------

dds_filtered <- DESeq(dds_filtered)

# Explore results --------------------------------
res <- results(dds_filtered)
res
summary(res)

res0.01 <- results(dds_filtered, alpha = 0.01)
summary(res0.01)

# Saving the result
write.csv(as.data.frame(res), file="results.csv")

# contrasts
resultsNames(dds_filtered)

# In case we have multiple levels, for example: treated_4hrs, treated_8hrs and untreated
# results(dds_filtered, contrast = c('Dexamethasone', 'treated_4hrs', 'untreated'))
# results(dds_filtered, contrast = c('Dexamethasone', 'treated_8hrs', 'untreated'))


# MA plot: To visualize genes that are differentially expressed

plotMA(res, ylim = c(-5,5))

# Identify the gene with the minimum adjusted p-value (most significant) 
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="red", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=1, cex=0.8, col="black")
})


# plotDispEsts: To visualize dispersion estimates
plotDispEsts(dds_filtered, ylim = c(1e-6, 1e1))

# Normalized counts for the most significant gene over treatment group
install.packages('ggbeeswarm')
library(ggbeeswarm)
geneCounts <- plotCounts(dds_filtered, gene = topGene, intgroup = c("Dexamethasone","CellLine"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Dexamethasone, y = count, color = CellLine, group = CellLine)) + geom_point(size = 3) + geom_line()


# Gene Clustering
# transform the raw count data dds_filtered
# vst function will perform variance stabilizing transformation
BiocManager::install("genefilter")
library(genefilter)
install.packages("pheatmap")
library(pheatmap)
vsd <- vst(dds_filtered, blind = F)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("CellLine","Dexamethasone")])
pheatmap(mat, annotation_col = anno)


# plotPCA: To visualize how our samples group by treatment
plotPCA(vsd, intgroup="Dexamethasone") 


#reset par
par(mfrow=c(1,1))

# Basic volcano plot: To identify genes with statistically significant differential expression
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


