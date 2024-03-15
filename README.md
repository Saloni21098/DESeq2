# Differential Gene Expression Analysis
This repository contains one R script (deseq2.R) that performs differential gene expressional analysis.

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

## Result
After running DESeq, a dataframe of 22368 genes were obtained with 6 columns described below:

1. baseMean: The average of the normalised count values dividing by size factors taken over all the samples.

2. log2 fold change (Maximum likelihood estimation): The effect size estimate and is calculated in the design factor Dexamethasone between treated and untreated levels. 
   It tells us how much the gene’s expression seems to have changed due to treatment with Dexamethasone in comparison to control.
   The positive values are the upregulated genes in the treated conditions and the negative values are the downregulated genes in the treated conditions.

3. lfcSE: The uncertainty associated with the effect size estimate is provided in the column lfcSE which tells us about the standard error estimate for the log two-fold change.

4. stat: The uncertainty of a particular effect size estimate can also be provided as the result of a statistical test – Wald test.
   The purpose of a test for differential expression is to test whether the data provides sufficient evidence to conclude that this value is really different from zero.
   DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is no effect of the treatment on the gene and that the observed          
   difference between treatment and control was merely caused by experimental variability (i.e., the type of variability that we can just as well expect between different samples in the same treatment group).
   
5. pvalue: The result of the statistical test is reported as a p-value (alpha significance is 0.05) for the genes.
   A p value indicates the probability that a fold change as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.

6. padj: The p adjusted value is the corrected p-values for multiple testing.
   We need to correct p-values for multiple testing because whenever we perform a statistical test we use a p-value of 0.05.
   So, 5% of differentially expressed genes are not really differentially expressed but they are only due to random chance and the drug has no effect on them.

In the dataset that we used, there are around 22,000 genes and 5% of 22,368 is 1118 genes. 
So, in our list of the genes being differentially expressed, 1118 of these genes are false positive which is really a huge number. 
In order to account for this problem, there are various methods to adjust the p-values and DESeq2 employs these methods to adjust p-values to avoid detection of false positive genes.


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
   
   
## Conclusion
Dexamethasone, a glucocoticosteroid can target the primary tissue (airway smooth muscle) during asthma treatment. Few of the important glucocorticoid-responsive genes which could be used for further studies are:
ANGPTL7 (Angiopoietin-like 7), GPX3 (Glutathione peroxidase 3), SPARCL1 (SPARC-like 1), FKBP5 (FK506-binding protein 5), KLF15 (Kruppel-like factor 15), SAMHD1 (SAM domain and HD domain 1)	and 
ZBTB16 (Zinc finger and BTB domain containing 16).
