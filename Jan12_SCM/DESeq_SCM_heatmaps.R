# Script to perform differential gene expression analysis using DESeq2 package - Visualization of data


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
# BiocManager::install("enrichplot")


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
library(pheatmap)
library(enrichplot)
conflict_prefer("pheatmap", "pheatmap")


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

# Match ENSEMBLE ID to Gene ID & remove duplicate gene names
counts_data2 <- as.data.frame(counts_data)
counts_data2 <- rownames_to_column(counts_data2, 'ensgene')
counts_data3 <- inner_join(counts_data2, annotables::grcm38)

head(counts_data3)

counts_data3 <- select(counts_data3, -c(ensgene, entrez, chr:description)) |>
  as.data.frame()

counts_data4 <- counts_data3 |> 
  pivot_longer(-symbol, names_to = "condition", values_to = "expression") |>
  group_by(symbol, condition) |>
  summarize(expression = sum(expression)) |>
  ungroup() |>
  pivot_wider(names_from = "condition", values_from = "expression")

counts_data5 <- as.data.frame(counts_data4) |>
  column_to_rownames("symbol") |>
  as.matrix()

head(counts_data5)



# Step 2: construct a DESeqDataSet object ----------

# Create & subset SCM dataset only
dds <- DESeqDataSetFromMatrix(countData = round(counts_data5[,colData$CellType=="SCM"]),
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

write.csv(res0.05,"Results_DGE_p0.05")



# Visualizing Results ----------------

# Heatmap of NFAT pathway:
gene.list <- read.table("~/SRI work/RNAseq/Genesets/Pathway genesets/SCM/NFAT.txt")
genes.to.plot <- as.character(gene.list$V1)

sampleInfo <- as.data.frame(colData(dds)[,c("Status","CellType")])

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

sample_order <- c("WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2")

expr_df <- assay(vsd)[genes.to.plot, ] |> t() |>
  as.data.frame() |>
  rownames_to_column("sample")

expr_df$sample <- str_split_fixed(expr_df$sample, "_", 2)[,1]

mat <- pivot_longer(expr_df, -sample, names_to = "gene", values_to = "expression") |>
  group_by(sample, gene) |>
  summarize(expression = mean(expression)) |>
  pivot_wider(names_from = "gene", values_from = "expression") |>
  as.data.frame() |>
  column_to_rownames("sample") |>
  as.matrix() |>
  t()

pheatmap(assay(vsd)[genes.to.plot, sample_order], cluster_cols = F, cluster_rows = T, 
         cellwidth = 20, cellheight = 20, angle_col = "45", color = colfunc(10), fontsize = 25, 
         scale = 'row', main = "NFAT Pathway")


# Heatmap of PI3K/MTOR/AKT pathway:
gene.list <- read.table("~/SRI work/RNAseq/Genesets/Pathway genesets/SCM/PI3K_MTOR_AKT.txt")
genes.to.plot <- as.character(gene.list$V1)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)
sample_order <- c("WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2")

expr_df <- assay(vsd)[genes.to.plot, ] |> t() |>
  as.data.frame() |>
  rownames_to_column("sample")

expr_df$sample <- str_split_fixed(expr_df$sample, "_", 2)[,1]

mat <- pivot_longer(expr_df, -sample, names_to = "gene", values_to = "expression") |>
  group_by(sample, gene) |>
  summarize(expression = mean(expression)) |>
  pivot_wider(names_from = "gene", values_from = "expression") |>
  as.data.frame() |>
  column_to_rownames("sample") |>
  as.matrix() |>
  t()

pheatmap(assay(vsd)[genes.to.plot, sample_order], cluster_cols = F, cluster_rows = T, 
         cellwidth = 15, cellheight = 15, angle_col = "90", color = colfunc(10), fontsize = 20, 
         scale = 'row', main = "PI3K/MTOR/AKT Pathway")

