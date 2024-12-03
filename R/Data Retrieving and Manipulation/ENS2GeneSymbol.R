###############################
# Enemble IDs to Gene Symbols #
###############################

# Methods Available:
#  1. BioMart web and biomaRt R Package
#  2. Annotables R package
#  3. AnnotationDbi

# Load libraries

bio_lib = c('biomaRt', 'org.Hs.eg.db', 'EnsDb.Hsapiens.v86')
library('biomaRt')
library('org.Hs.eg.db')
library('EnsDb.Hsapiens.v86')
library('tidyverse')

# Input list of ensemble ids
gene_ids = c('ENSG00000139618', 'ENSG00000171094', 'ENSG00000109132', 'ENSG00000152256', 'ENSG00000105369')
gene_df = as.data.frame(gene_ids)

# method 1:
listEnsembl()
ensembl <- useEnsembl(biomart = 'genes')
datasets <- listDatasets(ensembl)

ensembl_con <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl_con)
filt <- listFilters(ensembl_con)

getBM(attributes = c('ensembl_gene_id','external_gene_name'),
      filters = 'ensembl_gene_id',
      values = gene_ids,
      mart = ensembl_con)

# Method 2:
library('annotables')

grch38 |> 
  filter(ensgene %in% gene_ids)

# Method 3:
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

mapIds(org.Hs.eg.db, 
       keys = gene_ids,
       keytype = 'ENSEMBL',
       column = 'SYMBOL')

mapIds(EnsDb.Hsapiens.v86, 
       keys = gene_ids,
       keytype = 'GENEID',
       column = 'SYMBOL')