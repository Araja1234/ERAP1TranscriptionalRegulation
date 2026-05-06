# =============================================================================
# CD19+ B Cell Transcriptome Analysis 
# =============================================================================
setwd("/Users/araja2/Desktop/CD19_DISC")
library(DESeq2); library(edgeR); library(ggplot2); library(ggrepel)
library(dplyr); library(ggsci); library(WGCNA); library(flashClust)
library(ReactomePA); library(clusterProfiler); library(org.Hs.eg.db)
library(AnnotationDbi); library(enrichplot); library(readr)
library(biomaRt); library(openxlsx)

# =============================================================================
# 1. LOAD DATA FROM RDS
# =============================================================================

dat         <- readRDS("CD19_data.rds")
Counts      <- dat$counts       
metadatanew <- dat$metadata     

stopifnot(all(colnames(Counts) == rownames(metadatanew)))
cat("Data Dimesnsions:", nrow(Counts), "genes x", ncol(Counts), "samples\n")

# =============================================================================
# 2. DESeq2 — check ERAP1/ERAP2 expression per condition (Wald test)
# =============================================================================

cond    <- factor(metadatanew$Condition)
colData <- data.frame(row.names = colnames(Counts), cond)
dds     <- DESeqDataSetFromMatrix(Counts, colData = colData, design = ~ cond)
dds$Condition <- relevel(dds$cond, "HC")
dds1    <- DESeq(dds, test = "Wald", modelMatrixType = "standard")

plotCounts(dds1, gene = "ENSG00000164307", intgroup = "cond", main = "ERAP1",
           normalized = FALSE, transform = TRUE)
plotCounts(dds1, gene = "ENSG00000164308", intgroup = "cond", main = "ERAP2",
           normalized = FALSE, transform = TRUE)

# =============================================================================
# 3. DESeq2 — pan-UV vs HC 
# =============================================================================

metadata_combined <- metadatanew
metadata_combined$Condition <- ifelse(
  metadata_combined$Condition %in% c("AU", "BS", "IU"), "pan_UV",
  metadata_combined$Condition)

dds_combined  <- DESeqDataSetFromMatrix(Counts, colData = metadata_combined,
                                        design = ~ Condition)
dds_combined$Condition <- relevel(dds_combined$Condition, "HC")
dds_combined1 <- DESeq(dds_combined, parallel = TRUE, test = "LRT",
                        reduced = ~ 1, modelMatrixType = "standard")
res_combined  <- results(dds_combined1,
                         contrast = c("Condition", "pan_UV", "HC"))

# Significant genes (p < 0.05)
res_combined_subset       <- subset(res_combined,
                                    pvalue <= 0.05 & !is.na(pvalue) & !is.na(padj))
significant_genes_combined <- sum(res_combined_subset$pvalue <= 0.05)
cat("Significant genes (p<0.05):", significant_genes_combined, "\n")

# ERAP1 and ERAP2 p-values
ERAP1 <- "ENSG00000164307"
ERAP2 <- "ENSG00000164308"
cat("ERAP1:\n"); print(subset(res_combined, rownames(res_combined) == ERAP1))
cat("ERAP2:\n"); print(subset(res_combined, rownames(res_combined) == ERAP2))

# Annotate with gene symbols
convertedGenes         <- bitr(rownames(res_combined), fromType = "ENSEMBL",
                               toType = "SYMBOL", OrgDb = org.Hs.eg.db)
res_combinedvsd_scaled <- as.data.frame(res_combined[order(res_combined$padj), ])
rescombined_anno       <- merge(res_combinedvsd_scaled, convertedGenes,
                                by.x = "row.names", by.y = "ENSEMBL")

# Write supplementary table
wb <- createWorkbook()
addWorksheet(wb, "Table2")
writeData(wb, sheet = "Table2", rescombined_anno)
saveWorkbook(wb, "SupplementaryFile1.xlsx", overwrite = TRUE)

# --- Supplementary Figure: Volcano -------------------------------------------
volcano_data <- data.frame(
  gene           = rescombined_anno$SYMBOL,
  log2FoldChange = replace(rescombined_anno$log2FoldChange,
                           is.na(rescombined_anno$log2FoldChange), 0),
  pvalue         = replace(rescombined_anno$pvalue,
                           is.na(rescombined_anno$pvalue), 1))
volcano_data$Significance <- with(volcano_data,
  ifelse(pvalue < 0.05 & log2FoldChange >  0.5, "Upregulated",
  ifelse(pvalue < 0.05 & log2FoldChange < -0.5, "Downregulated",
         "Non-significant")))

supp_volcano <- ggplot(data = volcano_data,
                       aes(x = log2FoldChange, y = -log10(pvalue),
                           color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated"     = "#D55E00",
                                "Downregulated"   = "#0072B2",
                                "Non-significant" = "grey")) +
  geom_point(data = subset(volcano_data, gene %in% c("ERAP1", "ERAP2")),
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = "black", size = 3) +
  geom_text_repel(data = subset(volcano_data, gene %in% c("ERAP1", "ERAP2")),
                  aes(label = gene), size = 5, color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  labs(title = "Autoimmune Uveitis Patients vs Controls in CD19+ B cells",
       x = "Log2 Fold Change", y = "-Log10 P-Value") +
  theme_classic() +
  theme(text       = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 12))
print(supp_volcano)
ggsave("SuppFig_CD19_volcano.svg", plot = supp_volcano, width = 8, height = 6)

# =============================================================================
# 4. WGCNA — VST -> network -> module detection
# =============================================================================

dds2      <- DESeqDataSetFromMatrix(countData = Counts,
                                    colData   = metadatanew, design = ~ 1)
dds75c    <- dds2[rowSums(counts(dds2) >= 1) >= 30, ]  # 75% of 42 samples
dds_normc <- vst(dds75c)
datExpr01 <- t(assay(dds_normc))   # rows = samples, cols = genes

# QC
gsg <- goodSamplesGenes(datExpr01, verbose = 3)
gsg$allOK
if (!gsg$allOK) datExpr01 <- datExpr01[gsg$goodSamples, gsg$goodGenes]

# Outlier detection
sampleTree  <- hclust(dist(datExpr01), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
clust       <- cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
keepSamples <- (clust == 1)
datExpr01   <- datExpr01[keepSamples, ]
dim(datExpr01)
cat("WGCNA input:", nrow(datExpr01), "samples x", ncol(datExpr01), "genes\n")

# Soft-thresholding power
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft    <- pickSoftThreshold(datExpr01, dataIsExpr = TRUE,
                            powerVector = powers, corFnc = cor,
                            corOptions  = list(use = "p"),
                            networkType = "unsigned")

tiff("powerAnalysis.tif", width = 10, height = 10,
     units = "in", res = 300, compression = "lzw")
par(mfrow = c(1, 2)); cex1 <- 0.9
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "blue")
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
dev.off()

# Network construction (softPower = 9)
softPower <- 9
adj       <- adjacency(datExpr01, type = "unsigned", power = softPower)
TOM       <- TOMsimilarity(adj)
colnames(TOM) <- colnames(datExpr01)
rownames(TOM) <- colnames(datExpr01)
dissTOM   <- 1 - TOM
geneTree  <- flashClust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab = "", sub = "", cex = 0.3)

# Dynamic tree cut
minModuleSize  <- 20
dynamicModules <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                method = "hybrid", deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
dynamicColors  <- labels2colors(dynamicModules)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Merge modules (cutHeight = 0.1 i.e. r > 0.9)
MEList       <- moduleEigengenes(datExpr01, colors = dynamicColors)
MEs          <- MEList$eigengenes
merge        <- mergeCloseModules(datExpr01, dynamicColors,
                                  cutHeight = 0.1, verbose = 3)
mergedColors <- merge$colors
mergedMEs    <- merge$newMEs
mergedColors[dynamicColors == "grey"] <- "grey"
moduleColors <- mergedColors

# Assign gene names
names(datExpr01) <- rownames(dds_normc)

# Check ERAP1/ERAP2 module
cat("ERAP1 module:", unique(moduleColors[names(datExpr01) == ERAP1]), "\n")
cat("ERAP2 module:", unique(moduleColors[names(datExpr01) == ERAP2]), "\n")
cat("Genes in red module:", sum(moduleColors == "red"), "\n")

# Gene-module assignment table
annotated_genes    <- bitr(colnames(datExpr01), fromType = "ENSEMBL",
                           toType = c("SYMBOL", "ENTREZID"),
                           OrgDb = org.Hs.eg.db)
moduleLabels       <- as.numeric(factor(moduleColors))
gene_MEs           <- as.data.frame(t(rbind(colnames(datExpr01),
                                            moduleColors, moduleLabels)))
colnames(gene_MEs)[1:3] <- c("ensembl_gene_id", "mod_color", "mod_number")
gene_MEs_annotated <- merge(gene_MEs, annotated_genes,
                            by.x = "ensembl_gene_id", by.y = "ENSEMBL")
write.table(gene_MEs_annotated, "CD19_allModules.csv",
            sep = ",", row.names = FALSE)

tiff("NewgeneTree_merged_dynamic.tif", width = 6, height = 5,
     units = "in", res = 300, compression = "lzw")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# =============================================================================
# 5. kME + ERAP1 CORRELATION -> merge_data
# =============================================================================

MEs         <- moduleEigengenes(datExpr01, colors = moduleColors)$eigengenes
datKME      <- signedKME(datExpr01, MEs, outputColumnName = "kME")
red_genes   <- moduleColors == "red"
ensembl_ids <- colnames(datExpr01)[red_genes]
kME_values  <- datKME[red_genes, "kMEred"]

# Annotate (AnnotationDbi)
annotationsKME  <- select(org.Hs.eg.db, keys = ensembl_ids,
                          keytype = "ENSEMBL", columns = "SYMBOL")
annotationsKME  <- annotationsKME[!duplicated(annotationsKME$ENSEMBL), ]
matched_indices <- match(annotationsKME$ENSEMBL, ensembl_ids)
valid_indices   <- !is.na(matched_indices)

result_df <- data.frame(
  Module        = rep("red", sum(valid_indices)),
  Ensembl_ID    = ensembl_ids[matched_indices[valid_indices]],
  Gene_Name     = annotationsKME$SYMBOL[valid_indices],
  Module_Number = moduleLabels[red_genes][matched_indices[valid_indices]],
  kME           = kME_values[matched_indices[valid_indices]])

# Pearson correlation with ERAP1 expression
ERAP1_expression  <- datExpr01[, ERAP1]
red_gene_ids      <- colnames(datExpr01)[moduleColors == "red"]
gene_correlations <- sapply(red_gene_ids, function(gene) {
  cor(datExpr01[, gene], ERAP1_expression, method = "pearson")
})
correlation_df <- data.frame(Gene = red_gene_ids,
                             Correlation = gene_correlations)

merge_data <- merge(result_df, correlation_df,
                    by.x = "Ensembl_ID", by.y = "Gene", all.x = TRUE)
write.table(merge_data, "May2025Red_Gene_ERAP1_CorrelationsKME.txt",
            sep = "\t", row.names = FALSE)

# Add to supplementary file
wb <- loadWorkbook("SupplementaryFile1.xlsx")
addWorksheet(wb, "Table 3")
writeData(wb, "Table 3", merge_data)
saveWorkbook(wb, "SupplementaryFile1.xlsx", overwrite = TRUE)

save(merge_data, file = "merge_data.RData")

# =============================================================================
# 6. TF FILTERING — ENCODE ChIP-seq >= 2 experiments at ERAP1 promoter
# =============================================================================

file_paths  <- c("~/Downloads/EH38E2393417.tsv", "~/Downloads/EH38E2393464.tsv", "~/Downloads/EH38E2393465.tsv",
                 "~/Downloads/EH38E2393466.tsv", "~/Downloads/EH38E2393467.tsv")
df_combined <- bind_rows(lapply(file_paths, read_delim, delim = "\t",
                                escape_double = FALSE, trim_ws = TRUE))
TFs_with_support <- df_combined %>%
  filter(`# of experiments that support TF binding` >= 2) %>%
  pull(Factor) %>% unique()
save(TFs_with_support, file = "TFs_with_support.RData")

# =============================================================================
# 7. FIGURE 2A — kME vs ERAP1 correlation scatter
# =============================================================================

highlight_genes <- c("XBP1", "IRF4")

redM <- ggplot() +
  geom_point(data = merge_data, aes(x = kME, y = Correlation),
             color = "red", alpha = 0.3) +
  geom_point(data = merge_data[merge_data$Gene_Name %in% TFs_with_support, ],
             aes(x = kME, y = Correlation),
             size = 4, shape = 16, color = "darkred") +
  geom_text_repel(data = merge_data[merge_data$Gene_Name %in% TFs_with_support, ],
                  aes(x = kME, y = Correlation, label = Gene_Name),
                  size = 6, color = "black") +
  geom_point(data = merge_data[merge_data$Gene_Name == "ERAP1", ],
             aes(x = kME, y = Correlation),
             size = 4, shape = 16, color = "black") +
  geom_text_repel(data = merge_data[merge_data$Gene_Name == "ERAP1", ],
                  aes(x = kME, y = Correlation, label = Gene_Name),
                  size = 6, color = "black") +
  geom_point(data = merge_data[merge_data$Gene_Name %in% highlight_genes, ],
             aes(x = kME, y = Correlation),
             size = 4, shape = 16, color = "black") +
  geom_text_repel(data = merge_data[merge_data$Gene_Name %in% highlight_genes, ],
                  aes(x = kME, y = Correlation, label = Gene_Name),
                  size = 6, color = "black") +
  xlim(-1, 1) + ylim(-1, 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Module Membership (kME)", y = "Correlation with ERAP1")

print(redM)
ggsave("FigureACD19.svg", plot = redM, width = 8, height = 6, dpi = 300)

# =============================================================================
# 8. FIGURE 2B — Reactome enrichment full ERAP1 module
# =============================================================================

entrez_idsred       <- bitr(merge_data$Gene_Name, fromType = "SYMBOL",
                            toType = "ENTREZID", OrgDb = org.Hs.eg.db)
reactome_resultsRed <- enrichPathway(gene = entrez_idsred$ENTREZID,
                                     organism = "human",
                                     pvalueCutoff = 0.1, readable = TRUE)
reactome_resultsReddf <- as.data.frame(reactome_resultsRed)
write.csv(reactome_resultsReddf, "reactome_resultsReddf.csv")
save(reactome_resultsReddf, file = "reactome_resultsReddf.RData")

selected_pathways <- c("R-HSA-983169","R-HSA-983170","R-HSA-381070",
                       "R-HSA-381038","R-HSA-5633008","R-HSA-381119",
                       "R-HSA-8939236","R-HSA-532668","R-HSA-5362768")
top_dfred <- reactome_resultsReddf %>%
  filter(ID %in% selected_pathways) %>%
  arrange(p.adjust) %>% head(10)

top_dfred$Description <- gsub(
  "Antigen Presentation: Folding, assembly and peptide loading of class I MHC",
  "Antigen Presentation", top_dfred$Description)
top_dfred$Description <- gsub(
  "N-glycan trimming in the ER and Calnexin/Calreticulin cycle",
  "N-glycan trimming in ER", top_dfred$Description)
top_dfred$Description <- gsub(
  "TP53 Regulates Transcription of Cell Death Genes",
  "TP53 Transcription regulation", top_dfred$Description)
top_dfred$Description <- gsub(
  "Class I MHC mediated antigen processing & presentation",
  "MHC-I antigen processing", top_dfred$Description)

nejm_colors   <- pal_nejm("default")(8)
color_palette <- colorRampPalette(nejm_colors[c(1, 5, 7, 2)])
p_values_red  <- top_dfred$p.adjust

pathway_plot_red <- ggplot(top_dfred,
                           aes(x = reorder(Description, Count), y = Count,
                               fill = p_values_red)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.8,
           color = "black", linewidth = 0.7, position = "dodge") +
  scale_fill_gradient(low = "#AA1D1D", high = "#FFE5B4", name = "p-adjusted") +
  labs(x = NULL, y = "Gene Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title  = element_text(face = "bold"),
        panel.grid  = element_blank())
print(pathway_plot_red)
ggsave("Plot-RedGenesPathways.svg", plot = pathway_plot_red,
       width = 8, height = 6, dpi = 300)

# =============================================================================
# 9. FIGURE 2C — Reactome enrichment TF subset (JASPAR filtered)
# =============================================================================

genes        <- read.table("TFcheckpoint_CD19.txt", header = TRUE,
                           sep = "\t", check.names = FALSE)
genes        <- genes[!duplicated(genes$Gene_symbol), ]
gene_symbols <- as.character(genes$Gene_symbol)
common_TFs   <- intersect(gene_symbols, merge_data$Gene_Name)
tf_output_df <- merge_data[merge_data$Gene_Name %in% common_TFs, ]
write.table(tf_output_df, "CD19ERAP1-tf.txt", sep = "\t", row.names = TRUE)
save(tf_output_df, file = "tf_output_df.RData")

entrez_idstf          <- bitr(tf_output_df$Gene_Name, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Hs.eg.db)
reactome_resultstffiltered <- enrichPathway(gene = entrez_idstf$ENTREZID,
                                            organism = "human",
                                            pvalueCutoff = 0.1,
                                            readable = TRUE)
reactome_resultstffiltereddf <- as.data.frame(reactome_resultstffiltered) %>%
  mutate(GeneCount = as.numeric(sub("/.*", "", GeneRatio)))
save(reactome_resultstffiltereddf, file = "reactome_resultstffiltereddf.RData")

top_df        <- reactome_resultstffiltereddf %>%
  arrange(desc(GeneCount)) %>% head(10)
color_palette <- colorRampPalette(nejm_colors[c(1, 3, 7, 2)])

plotCD19 <- ggplot(top_df, aes(x = reorder(Description, GeneCount),
                               y = GeneCount,
                               fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradientn(colors = c("red","royalblue"),
                       name = "-log10(Adjusted P-Value)") +
  labs(x = NULL, y = "Gene Count") +
  theme_minimal() +
  theme(axis.text.y      = element_text(size = 14, hjust = 1),
        axis.text.x      = element_text(size = 12),
        legend.title     = element_text(size = 12, face = "bold"),
        legend.text      = element_text(size = 12),
        plot.margin      = margin(1, 1, 1, 1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio     = 1.2) +
  coord_flip()
print(plotCD19)
ggsave("Plot-RedTFPathways.svg", plot = plotCD19, width = 12, height = 8, dpi = 600)

# =============================================================================
# SAVE WORKSPACE 
# =============================================================================

save(rescombined_anno,
     merge_data,
     TFs_with_support,
     tf_output_df,
     reactome_resultsReddf,
     reactome_resultstffiltereddf,
     file = "CD19_workspace.RData")

