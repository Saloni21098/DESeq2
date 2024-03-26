# 20-03-2024
# script to convert Ensembl IDs to gene IDs
# setwd("C:/Users/salon/OneDrive/Desktop/DeSeq2")

# Load the libraries
library(biomaRt)
library(tidyverse)

library(devtools)
devtools::install_github("stephenturner/annotables")
library(annotables)

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)


# creating text file for list of Ensembl ID's
ensembl_id <- c("ENSG00000171819", "ENSG00000211445", "ENSG00000152583", "ENSG00000096060", 
                 "ENSG00000163884", "ENSG00000101347", "ENSG00000109906")
df_ensembl_id <- data.frame(ensembl_id)
ensembl_ids <- write_delim(df_ensembl_id, "ensembl_ids.txt", delim = "\t")


devtools::install_version("dbplyr", version = "2.3.4")
BiocManager::install("Bioconductor/BiocFileCache")


# method1: biomaRt

listEnsembl() 
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl_con <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl_con)
fil <- listFilters(ensembl_con)

getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
      filters = 'ensembl_gene_id',
      values = ensembl_ids$ensembl_id,
      mart = ensembl_con)


# method 2: annotables
grch38 %>% 
  filter(ensgene %in% ensembl_ids$ensembl_id)


# method 3: annotation DBs

keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

mapIds(org.Hs.eg.db,
       keys = ensembl_ids$ensembl_id,
       keytype = 'ENSEMBL',
       column = "SYMBOL")

columns(EnsDb.Hsapiens.v86)
keytypes(EnsDb.Hsapiens.v86)

mapIds(EnsDb.Hsapiens.v86,
       keys = ensembl_ids$ensembl_id,
       column = 'SYMBOL',
       keytype = 'GENEID')