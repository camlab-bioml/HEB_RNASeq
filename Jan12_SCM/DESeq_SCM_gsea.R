# Script to perform differential gene expression analysis using DESeq2 package and GSEA with fgsea


# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("org.Mm.eg.db")
# devtools::install_github("slowkow/ggrepel")
# BiocManager::install("GO.db")
# BiocManager::install("GOstats")
# BiocManager::install("pathview")
# BiocManager::install("gage")
# BiocManager::install("gageData")
# install.packages("devtools")
# devtools::install_github("stephenturner/annotables")
# install.packages("RColorBrewer")
# BiocManager::install("RColorBrewer")
# BiocManager::install("DOSE")
# BiocManager::install("genekitr")


# load libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(kableExtra)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(conflicted)
library(DT)
library(ggpubr)
library(conflicted)
library(DOSE)
conflict_prefer("PlotMA", "DESeq2")

# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv('salmon.merged.gene_counts_length_scaled_removeID.tsv', header = T, row.names = 1,
                        sep = "\t")
head(counts_data)
dim(counts_data)

# read in sample info
colData <- read.csv('Sample_info.txt', header = T, row.names = 1, sep = "\t")


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))


# are they in the same order?
all(colnames(counts_data) == rownames(colData))



# Step 2: construct a DESeqDataSet object ----------

# Create & subset SCM dataset only
dds <- DESeqDataSetFromMatrix(countData = round(counts_data[,colData$CellType=="SCM"]),
                              colData = colData[colData$CellType=="SCM",],
                              design = ~ Status)

dds 

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 100 reads total
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]

dds

# Check reference levels
dds$Status

# Set the factor level
dds$Status <- relevel(dds$Status, ref = "WT")
dds$Status

# data tranfromation
vsd <- vst(dds, blind=FALSE)



# Step 3: Run DESeq ----------------------

dds <- DESeq(dds)
res <- results(dds)

res

# Check contrasts
resultsNames(dds)



# Explore Results ----------------

summary(res)

# Summary of results based on p=0.05
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

# MA plot
plotMA(res)
plotMA(res0.05)

# write.csv(res0.05,"Results_DGE_p0.05")



# GSEA (https://www.biostars.org/p/467197/) ----------------

library(tibble)
library(dplyr)
library(tidyr)
library(fgsea)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(DT)
library(AnnotationDbi)
library(org.Mm.eg.db)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")

# already we performed DESeq2 analysis and have statistics for working on it
res0.05$row <- rownames(res0.05)

# Map Ensembl gene IDs to the symbol. First, create a mapping table.
columns(org.Mm.eg.db)

ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
                                    key=res0.05$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
names(ens2symbol)[1] <- "row"

ens2symbol <- as_tibble(ens2symbol)
ens2symbol
# joining
res0.05 <- merge(data.frame(res0.05), ens2symbol, by=c("row"))

# remove the NAs, averaging statistics for a multi-hit symboln ###Do this before DESEQ step? And ADD COUNTS instead of averaging?
res4 <- res0.05 %>% 
  select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

res4

# creating  a named vector [ranked genes]
ranks <- res4$stat
names(ranks) <- res4$SYMBOL



# Hallmark Pathways ----------------

# Load the pathway (gene set) into a named list (downloaded mysigdb were located in my "~" directory):
pathways.hallmark <- gmtPathways("msigdb_v2022.1.Mm_files_to_download_locally/mh.all.v2022.1.Mm.symbols.gmt")

# show a few lines from the pathways file
head(pathways.hallmark)

# Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, minSize = 15)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

### BAR PLOT ###

# Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

fgseaResTidy %>% 
  mutate(pathway=gsub("HALLMARK_","", pathway)) %>% 
  filter(NES>1|NES<0) %>% 
  ggplot(aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA (SCM)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic() +
  theme(text = element_text(size = 35))



# Biocarta Pathways ----------------

# Load the pathway (gene set) into a named list (downloaded mysigdb were located in my "~" directory):
pathways.biocarta <- gmtPathways("MSIGDB genesets/m2.cp.biocarta.v2023.1.Mm.symbols.gmt")

# show a few lines from the pathways file
head(pathways.biocarta)

# Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways.biocarta, stats=ranks, minSize = 20)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

### BAR PLOT ###

# Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

fgseaResTidy %>% 
  mutate(pathway=gsub("BIOCARTA_","", pathway)) %>% 
  filter(NES>1|NES<0) %>% 
  ggplot(aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Biocarta pathways Enrichment Score from GSEA (SCM)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic() +
  theme(text = element_text(size = 35))


