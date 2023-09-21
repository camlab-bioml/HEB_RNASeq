# Script to perform differential gene expression analysis using DESeq2 package


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
# if (!requireNamespace('remotes', quietly = TRUE))
# install.packages('remotes')
# remotes::install_github('kevinblighe/PCAtools')
# BiocManager::install("dittoSeq")
# install.packages("binovisualfields")
# BiocManager::install("DEGreport")


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
library(PCAtools)
library(dittoSeq)
library(binovisualfields)
library(pheatmap)
library(dplyr)
conflict_prefer("pheatmap", "pheatmap")
library(AnnotationDbi)
library(org.Mm.eg.db)

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

# Create dataset 
dds <- DESeqDataSetFromMatrix(countData = round(counts_data5),
                              colData = colData,
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



# PCATools ----------------------

vst <- assay(vst(dds))
p <- pca(vst, metadata = colData, removeVar = 0.1)

screeplot(p, axisLabSize = 25, titleLabSize = 30)

biplot(p, labSize = 5)
biplot(p, x = "PC1", y = "PC4", labSize = 4, shape = 'Status', 
       legendPosition = 'right', colby = 'CellType', pointSize = 5)
biplot(p, x = "PC1", y = "PC4", labSize = 0, shape = 'Status', 
       legendPosition = 'right', colby = 'CellType', pointSize = 7, drawConnectors = F,
       legendLabSize = 30, legendTitleSize = 30, axisLabSize = 30)

biplot(p, showLoadings = TRUE,
       labSize = 4, pointSize = 4, sizeLoadingsNames = 5)
biplot(p, x = "PC1", y = "PC4", labSize = 4, shape = 'Status', 
       legendPosition = 'right', colby = 'CellType', pointSize = 5, 
       showLoadings = TRUE, sizeLoadingsNames = 4)

pairsplot(p)

plotloadings(p, labSize = 3)

plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

eigencorplot(p, metavars = c('CellType', 'Status'))

eigencorplot(p,
             components = getComponents(p, 1:10),
             metavars = c('CellType', 'Status'),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PCA correlations',
             colFrame = 'white',
             plotRsquared = FALSE) 

p$rotated[1:16, 1:5]

p$loadings[1:100, 1:5]



# Explore Results ----------------

summary(res)

# Summary of results based on p=0.05
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

# MA plot
plotMA(res)
plotMA(res0.05)

write.csv(res0.05,"Results_DGE_p0.05")

save.image("Data_setup.RData")
load("Data_setup.RData")



# Visualizing Results (Box plots of individual genes) ----------------

# Box plot (WT vs KO):
dds$Status <- factor(dds$Status, levels=c("WT", "KO"))

plot_my_counts <- function(dds, gene_name) {
  plotCounts(dds, gene_name, intgroup = c("Status"), transform = TRUE, returnData = TRUE) %>%   
    ggplot(aes(x = Status, y = count,col=Status)) + 
    geom_jitter(width=0.1) +
    labs(title = gene_name, y = "log-normalized count") + 
    scale_y_log10() +
    geom_boxplot(aes(fill=Status)) +
    theme_bw() +
    stat_compare_means(method = "t.test", aes(label = ..p.signif..), hide.ns = TRUE, label.x = 1.5, label.y. = 1, size = 8)
}

# Call function:
plot_my_counts(dds, "Gzmb")


# Box plot (WT vs KO by cell type:)
dds$CellType <- factor(dds$CellType, levels=c("Naive", "SCM", "CM", "EM"))

plot_my_counts <- function(dds, gene_name) {
  plotCounts(dds, gene_name, intgroup = c("Status","CellType"), transform = TRUE, returnData = TRUE) %>%   
    ggplot(aes(x = Status, y = count,col=Status)) + 
    geom_jitter(width=0.1) + 
    facet_wrap(~CellType, scales = "free_y") +
    labs(title = gene_name, y = "log-normalized count") + 
    scale_y_log10(expand = expansion(mult = c(0.05,0.3))) +
    geom_boxplot(aes(fill=Status)) +
    theme_bw() +
    theme(strip.background=element_rect(fill="white"), text = element_text(size = 30))
}

# Call function:
plot_my_counts(dds, "Sell")


# Facet-wrap box plot with multiple genes (WT vs KO by cell type):
genes.to.plot <- c("Sell", "Il7r", "Tcf7", "Eomes", "Bcl6", 
                   "Ly6a", "Il6ra", "Bcl2")

genes.to.plot <- c("Ifng", "Zeb2", "Gzma", 
                   "Cd44", "Prdm1", "Tbx21", "Cxcr3", "Il2rb", "Cd28", "Ccl5")

dds$CellType <- factor(dds$CellType, levels=c("Naive", "SCM", "CM", "EM"))

# To check if a gene is missing or has no count:
genes.to.plot[!(genes.to.plot %in% rownames(assay(vsd)))]

comparisons <- list(c("Naive", "CM"), c("Naive", "EM"), c("SCM", "CM"), c("SCM", "EM"))

lapply(genes.to.plot, function(x){
  plotCounts(dds, x, intgroup = c("CellType"), transform = TRUE, returnData = TRUE) |>
    mutate(gene = x)
}) |> bind_rows() %>%
  ggplot(aes(x = CellType, y = count, color = CellType)) +
  geom_boxplot(aes(fill=CellType)) +
  scale_y_log10(expand = expansion(mult = c(0.05,0.1))) +
  ylab("log10 Normalized Counts") +
  ggtitle("Naive-like T cell gene markers") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  theme(strip.background=element_rect(fill="white"), 
        text = element_text(size = 30),  axis.text = element_text(size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# Visualizing Results (Heatmap) ----------------

# Heatmap of DEG from each cell type:
geneslist <- read.table("Celltype DEG genelist/WT_v_KO_genelist_dup_TCR_removed.txt")
genes.to.plot <- as.character(geneslist$V1)

# To check if a gene is missing or has no count:
genes.to.plot[!(genes.to.plot %in% rownames(assay(vsd)))]

vsd$CellType <- factor(vsd$CellType, levels=c("Naive", "SCM", "CM", "EM"))
vsd$Status <-  relevel(vsd$Status, ref = "WT")

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

sample_order <- c("WT_Naive_1", "WT_Naive_2", "KO_Naive_1", "KO_Naive_2",
                  "WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2",
                  "WT_CM_1", "WT_CM_2", "KO_CM_1", "KO_CM_2",
                  "WT_EM_1", "WT_EM_2", "KO_EM_1", "KO_EM_2")

pheatmap(assay(vsd)[genes.to.plot, sample_order], annotation_col = sampleInfo, cluster_cols = F, 
         cluster_rows = T, fontsize_row = 2, scale = "row", treeheight_row = 10, 
         clustering_distance_rows = "correlation", 
         annotation_colors = list("Status" = c("WT" = "#FF595E", "KO" = "#ffca3a"),
                                  "CellType" = c("Naive" = "#2297E6", "SCM" = "#DF536B", 
                                                 "CM" = "#61D04F", "EM" = "#CD0BBC")))


# Heatmap of SCM gene signature from Pace et al:
geneslist <- read.table("~/SRI work/RNAseq/DESeq2/July4_WT_v_KO_DEG_Heatmap/SCM geneset_pace_remove_no_count.txt")
genes.to.plot <- as.character(geneslist$V1)

# To check if a gene is missing or has no count:
genes.to.plot[!(genes.to.plot %in% rownames(assay(vsd)))]

vsd$CellType <- factor(vsd$CellType, levels=c("Naive", "SCM", "CM", "EM"))
vsd$Status <-  relevel(vsd$Status, ref = "WT")

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

sample_order <- c("WT_Naive_1", "WT_Naive_2", "KO_Naive_1", "KO_Naive_2",
                  "WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2",
                  "WT_CM_1", "WT_CM_2", "KO_CM_1", "KO_CM_2",
                  "WT_EM_1", "WT_EM_2", "KO_EM_1", "KO_EM_2")

breaksList = seq(-3, 2, by = 0.5)

pheatmap(assay(vsd)[genes.to.plot, sample_order], annotation_col = sampleInfo, cluster_cols = F, 
         cluster_rows = T, fontsize_row = 10, scale = "row", 
         color = colorRampPalette(rev(brewer.pal(n = 100, name = "RdYlBu")))(length(breaksList)),
         cellwidth = 15, 
         cellheight = 8,
         fontsize = 15,
         breaks = breaksList, annotation_colors = list("Status" = c("WT" = "#FF595E", "KO" = "#ffca3a"),
                                                       "CellType" = c("Naive" = "#2297E6", "SCM" = "#DF536B", 
                                                                      "CM" = "#61D04F", "EM" = "#CD0BBC")))


# Heatmap of GO BP geneset for T cell activation:
geneslist <- read.table("~/SRI work/RNAseq/Genesets/MSIGDB genesets/GOBP_ALPHA_BETA_T_CELL_ACTIVATION.v2023.txt")
genes.to.plot <- as.character(geneslist$V1)

# To check if a gene is missing or has no count:
genes.to.plot[!(genes.to.plot %in% rownames(assay(vsd)))]

vsd$CellType <- factor(vsd$CellType, levels=c("Naive", "SCM", "CM", "EM"))
vsd$Status <-  relevel(vsd$Status, ref = "WT")

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

sample_order <- c("WT_Naive_1", "WT_Naive_2", "KO_Naive_1", "KO_Naive_2",
                  "WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2",
                  "WT_CM_1", "WT_CM_2", "KO_CM_1", "KO_CM_2",
                  "WT_EM_1", "WT_EM_2", "KO_EM_1", "KO_EM_2")

breaksList = seq(-2, 2, by = 0.1)

pheatmap(assay(vsd)[genes.to.plot, sample_order], annotation_col = sampleInfo, cluster_cols = F, 
         cluster_rows = T, fontsize_row = 5, scale = "row", color = colorRampPalette(rev(brewer.pal(n = 100, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList, cellwidth = 15, cellheight = 4, fontsize = 15,
         annotation_colors = list("Status" = c("WT" = "#FF595E", "KO" = "#ffca3a"),
                                  "CellType" = c("Naive" = "#2297E6", "SCM" = "#DF536B", 
                                                 "CM" = "#61D04F", "EM" = "#CD0BBC")))

# Heatmap of T cell markers:
genes.to.plot <- c("Ccr7", "Lef1", "Sell", "Foxo1", "Ltb", "Il7r", "Tcf7", "Eomes", "Foxp1", "Id3", "Bcl6", "Stat3", 
                   "Klf7", "Gzmm", "Ifng", "Ccl4", "Ccl5", "Gzmb", "Nkg7", "Zeb2", "Gzma", "Prf1", "Il12rb1", "Id2",
                   "Cd44", "Prdm1", "Tbx21", "Stat4", "Il2rb", "Cxcr3", "Bcl2", "Ccr5", "Fas", "Cd27", "Fasl", "Klf2", 
                   "Cxcr5", "Slamf6", "Tnf", "Cx3cr1", "Runx1")

# To check if a gene is missing or has no count:
genes.to.plot[!(genes.to.plot %in% rownames(assay(vsd)))]

vsd$CellType <- factor(vsd$CellType, levels=c("Naive", "SCM", "CM", "EM"))
vsd$Status <-  relevel(vsd$Status, ref = "WT")

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

sample_order <- c("WT_Naive_1", "WT_Naive_2", "KO_Naive_1", "KO_Naive_2",
                  "WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2",
                  "WT_CM_1", "WT_CM_2", "KO_CM_1", "KO_CM_2",
                  "WT_EM_1", "WT_EM_2", "KO_EM_1", "KO_EM_2")

colfunc <- colorRampPalette(c("blue", "black", "red"))  

breaksList = seq(-2, 2, by = 0.1)

pheatmap(assay(vsd)[genes.to.plot, sample_order], 
         annotation_col = sampleInfo, 
         cluster_cols = F, 
         cluster_rows = T, 
         fontsize_row = 10, 
         scale = "row", 
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList, 
         cellwidth = 10, 
         cellheight = 10,
         annotation_colors = list("Status" = c("WT" = "#FF595E", "KO" = "#ffca3a"),
                                  "CellType" = c("Naive" = "#2297E6", "SCM" = "#DF536B", "CM" = "#61D04F", "EM" = "#CD0BBC")))


# Heatmap of GO.Shiny gene clusters:
#Cluster 3
genes.to.plot <- c("Stat1", "Gbp2", "Gbp3", "Ifit1", "Ifi209", "Iigp1", "Ifi214", "Ifit3", "Igtp", "Gbp6", "Bst2", "Igtp", "Isg15", "Oas2", "Oas3", 
                   "Usp18", "Samhd1", "Zbp1")

#Cluster 5 (IFNY + Lymphocyte activation)
genes.to.plot <- c("Il12rb1", "Axl", "Il12rb2", "Xcl1", "Klrk1", "Cd160", "Irf8", "Ccr2", "Tbx21",
                   "Dlg5", "Cdkn1a", "Adam8", "Lag3", "Klrc1", "Itgal", "Fgl2", "Ctla2a", "Nedd4", "Rbpj", "Il2rb")

# To check if a gene is missing or has no count:
genes.to.plot[!(genes.to.plot %in% rownames(assay(vsd)))]

vsd$CellType <- factor(vsd$CellType, levels=c("Naive", "SCM", "CM", "EM"))
vsd$Status <-  relevel(vsd$Status, ref = "WT")

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

sample_order <- c("WT_Naive_1", "WT_Naive_2", "KO_Naive_1", "KO_Naive_2",
                  "WT_SCM_1", "WT_SCM_2", "KO_SCM_1", "KO_SCM_2",
                  "WT_CM_1", "WT_CM_2", "KO_CM_1", "KO_CM_2",
                  "WT_EM_1", "WT_EM_2", "KO_EM_1", "KO_EM_2")

colfunc <- colorRampPalette(c("blue", "black", "red"))  

pheatmap(assay(vsd)[genes.to.plot, sample_order], annotation_col = sampleInfo, cluster_cols = F, 
         cluster_rows = F, fontsize = 20, scale = "row", cellwidth = 20, cellheight = 20, 
         angle_col = "45", color = colfunc(10), 
         annotation_colors = list("Status" = c("WT" = "#FF595E", "KO" = "#ffca3a"),
                                  "CellType" = c("Naive" = "#2297E6", "SCM" = "#DF536B", 
                                                 "CM" = "#61D04F", "EM" = "#CD0BBC")))



# Volcano Plot ----------------

res1 = as.data.frame(res0.05)
res1 <- na.omit(res1) # check if correct code?

# Plot data based on p<0.05, log2Foldchange >0.5
# add a new column using the mutate function in dplyr
res1 = mutate(res1, sig=ifelse(res1$padj<0.05, "FDR<0.05", "Not Sig"))
res1[which(abs(res1$log2FoldChange)<0.5), 'sig'] = "Not Sig"

library(ggplot2)
p = ggplot(res1, aes(log2FoldChange, -log10(padj))) + 
  geom_point(aes(col=sig)) + 
  scale_color_manual(values = c("red", "black"))

p

# Annotate volcano plot
columns(org.Mm.eg.db)

res$ensembl <- gsub("\\..*","", row.names(res))
res$entrez <- mapIds(org.Mm.eg.db,
                     keys= res$ensembl,
                     column="ENTREZID",
                     keytype="SYMBOL", #Out ID is ENSMBL
                     multiVals="first")

res$symbol <- mapIds(org.Mm.eg.db,
                     keys= res$ensembl,
                     column="SYMBOL",
                     keytype="SYMBOL", #Out ID is ENSMBL
                     multiVals="first")

res1 = as.data.frame(res)
res1 <- na.omit(res1) # check if correct code?

# add a new column using the mutate function in dplyr
res1 = mutate(res1, sig=ifelse(res1$padj<0.05, "FDR<0.05", "Not Sig"))
res1[which(abs(res1$log2FoldChange)<0.5),'sig'] <- "Not Sig"

labels <- res1 %>% 
  filter(!grepl("Trav|Trdv|Trdc|Trbv", symbol),
         padj < 1e-5)

p = ggplot(res1, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig), size = 2) +
  scale_color_manual(values=c("red", "black")) +
  geom_vline(xintercept = -0.5, col = "black", linetype="dashed") +
  geom_vline(xintercept = 0.5,  col = "black", linetype="dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 30))

p+geom_text_repel(data=labels, aes(label=symbol), size=8,
                  max.overlaps = Inf) 

