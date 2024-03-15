# Differential Gene Expression Analysis
This repository contains one script that performs differential gene expressional analysis.

# Study Design
RNA-sequencing performed on 4 primary human airway smooth muscle cell lines treated with 1 micromolar dexamethasone for 18 hours.
For each of the 4 cell lines, study has a treated and untreated sample.

Goal: To understand the transcriptional changes occurring due to treatment with dexamethasone - a glucocorticosteroid.


## Requirements
- R version 4.2.1
- Package
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
Thus, few of the important glucocorticoid-responsive genes observed are listed below and these genes could be used for further studies.   

Ensembl ID	    GeneSymbol	GeneName	                                GenomicCoordinates	      baseMean	  log2FoldChange	lfcSE	      stat	      pvalue	  padj
ENSG00000171819	ANGPTL7	    Angiopoietin-like 7	                      chr1:11189341-11195981	  254.8825987	5.796857869	    1.397218194	4.148856558	3.34E-05	0.000600311
ENSG00000211445	GPX3	      Glutathione peroxidase 3	                chr5:151020438-151028993	12285.70005	3.742803103	    0.376835773	9.932186301	3.02E-23	1.21E-20
ENSG00000152583	SPARCL1	    SPARC-like 1	                            chr4:87473335-87531061	  997.4447326	4.602554541	    0.211741943	21.73662182	9.25E-105	1.70E-100
ENSG00000096060	FKBP5	      FK506-binding protein 5	                  chr6:35573585-35728583	  2564.384375	3.938431015	    0.296906744	13.26487559	3.70E-40	5.24E-37
ENSG00000163884	KLF15	      Kruppel-like factor 15	                  chr3:126342635-126357442	561.1110483	4.58218143	    0.412315586	11.11328696	1.08E-28	8.66E-26
ENSG00000101347	SAMHD1	    SAM domain and HD domain 1	              chr20:36890229-36951843	  12703.41283	3.857750212	    0.279324769	13.81098504	2.19E-43	4.03E-40
ENSG00000109906	ZBTB16	    Zinc finger and BTB domain containing 16	chr11:114059593-114250676	385.0713277	7.176149396	    0.4963593	  14.45757015	2.25E-47	5.17E-44





