# Differential Gene Expression Analysis
This repository contains one script that performs differential gene expressional analysis.

# Study Design
RNA-sequencing performed on 4 primary human airway smooth muscle cell lines treated with 1 micromolar dexamethasone for 18 hours.
For each of the 4 cell lines, study has a treated and untreated sample.

Goal: To understand the transcriptional changes occurring due to treatment with dexamethasone - a glucocorticosteroid.


## Requirements
- R version 4.2.1
- Package:
BiocManager
DESeq2 
tidyverse
airway
ggbeeswarm
genefilter
pheatmap

## Data used
airway package data

## DESeq2 workflow
- Estimate size factors
- Estimate dispersions
- Gene-wise dispersion estimates
- Mean-dispersion relationships
- Final dispersion estimates
- Fitting model
- Testing

The null hypothesis for DESeq2 is that there is no differential expression between 2 sample groups (treated and untreated). 
So log-fold change is computed from data using generalised linear models and a statistical test like Wald test that provides us sufficient evidence to conclude that the log-fold change is not zero. 
It is different from zero and what we observe is highly significant and is greater than what we would have observed due to random variation.


## Data visualization 
Plots used:

   a. MA plot: A scatter plot between mean of normalized count vs log two-fold change that shows which genes are differentially expressed.
      The gene with Ensembl ID ‘ENSG00000152583’ has the lowest p-adjusted value and hence is the most significant.
   
   b. Heatmap: The blocks represents genes that covary across patients.
      A set of genes in the heatmap are separating the N061011 cell line from the others, and there is another set of genes for which the dexamethasone treated samples have higher gene expression.
   
   c. PCA plot: Visualizes the variation in gene expression profiles across samples that received Dexamethasone treatment compared to those that didn't.
      Since, the "treated" and "untreated" samples form distinct clusters in the plot, it suggests significant gene expression changes upon Dexamethasone treatment. 

   d. Volcano plot: Visually identify genes with statistically significant differential expression, prioritize genes based on the magnitude of their expression change and
      focus on genes that are both statistically significant and have a biologically relevant fold change.
      Each data point in the plot represents a gene.
      The overall shape of the plot resembles a volcano, with many genes clustered around the center (low fold change, non-significant) and
      fewer genes scattered towards the extremes (large fold change, highly significant).
   
   
## Discussion
Thus, few of the important glucocorticoid-responsive genes which could be used for further studies are:
ANGPTL7 (Angiopoietin-like 7), GPX3 (Glutathione peroxidase 3), SPARCL1 (SPARC-like 1), FKBP5 (FK506-binding protein 5), KLF15 (Kruppel-like factor 15), SAMHD1 (SAM domain and HD domain 1)	and 
ZBTB16 (Zinc finger and BTB domain containing 16).
