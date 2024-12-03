"Load required libraries"
library(tidyverse)
library(GEOquery)
library(ggplot2)

##################################
#          Get the Data          #
##################################

# Read in gene dataset
data = read_csv("GSE183947_fpkm.csv")

# Get gene meta data
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)

##################################
#      Manipulate the Data       #
##################################

# Assign first element in gse to new object
metaData <- pData(phenoData(gse[[1]]))
head(metaData)

# Select only the necessary variables from the metaData
cleanMeta <- metaData |> 
  select(c(1, 10, 11, 17)) |> 
  rename(tissue = characteristics_ch1,
         metastasis = characteristics_ch1.1,
         samples = description) |> 
  mutate(tissue = gsub("tissue: ", "", tissue),
         metastasis = gsub("metastasis: ", "", metastasis))

head(cleanMeta)

dim(data)

# reshape the gene data into long format
data_long <- data |> 
  rename(gene = ...1) |> 
  gather(key = "samples", value = "FPKM", -gene)

head(data_long)
dim(data_long)

# Join the gene and metadata DFs
data_comb <- data_long |> 
  left_join(cleanMeta)

##################################
# Explore and Visualize the Data #
##################################

# Explore the Data
data_comb |> 
  filter(gene == 'BRCA1' | gene == 'BRCA2') |> 
  group_by(gene, tissue) |> 
  summarize(mean_FPKM = mean(FPKM)) |> 
  head()

# Barplot
data_comb |> 
  filter(gene == "BRCA1") |> 
  ggplot(
    aes(x = samples, y = FPKM, fill = tissue)) + 
  geom_col()

# Density plot
data_comb |> 
  filter(gene == "BRCA1") |> 
  ggplot(
    aes(x = FPKM, fill = tissue)) + 
  geom_density(alpha = 0.5)

# Boxplot
data_comb |> 
  filter(gene == "BRCA1") |> 
  ggplot(
    aes(x = metastasis, y = FPKM)) + 
  geom_boxplot()

# Scatterplot
data_comb |> 
  filter(gene == "BRCA1" | gene == "BRCA2") |> 
  spread(key = gene, value = FPKM) |> 
  ggplot(
    aes(x = BRCA1, y = BRCA2)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)

data_comb |> 
  filter(gene == "BRCA1" | gene == "BRCA2") |> 
  spread(key = gene, value = FPKM) |> 
  ggplot(
    aes(x = BRCA1, y = BRCA2, color = tissue)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F)

# Heatmap
genes_of_interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

data_long |> 
  filter(gene %in% genes.of.interest) |> 
  ggplot(aes(x = samples, y = gene, fill = FPKM)) + 
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')







