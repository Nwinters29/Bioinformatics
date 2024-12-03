# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

data(airway)

sample_info <- as.data.frame(colData(airway))

head(sample_info)

sample_info <- sample_info[, c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = 'sample_info.csv', sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# Step 1: preparing count data

# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# read in sample info data 
col_data <- read.csv('sample_info.csv')
head(col_data)

# Ensure data is the same between counts data and sample data
all(colnames(counts_data) %in% rownames(col_data))

# Are they in the same order?
all(colnames(counts_data) == rownames(col_data))

# Step 2: Construct a DESeqDataSet Object 

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = col_data,
                       design = ~ dexamethasone)
dds

# Recommended Steps
# pre-filter and remove rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# set a factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = 'untreated')

# NOTE: collapse technical replicates with the collapse function (never collapse biological replicates)

# Step 3: Run DESeq
dds <- DESeq(dds)
res <- results(dds)

res

# Explore Results
summary(res)

res0.01 <- results(dds, alpha=0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# ex for multiple level contrasts (treated4hrs, treated8hrs, untreated)
#### results(dds, contrast = c('dexamethasone', 'treated4hrs', 'untreated'))
#### allows for comparison between levels

# MA plot
plotMA(res)
