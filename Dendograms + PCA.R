library(DESeq2)
library(data.table)
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(gridExtra)
library(pheatmap)
library(viridis)
library(MetaboAnalystR)
library(WGCNA)
library(ggalluvial)
library(EnrichmentBrowser)
library(KEGGREST)
library(FELLA)
library(CorLevelPlot)
library(methylKit)
library(simplifyEnrichment)
library(ggh4x)
library(ggtext)

##### THEME   ----          
Style_format_theme <- theme(
  ## Title of the axis
  axis.title.x       = element_text(color="black", size=30),
  axis.title.y       = element_text(color="black", size=30),
  ## Text in the axis
  axis.text          = element_text(color="black", size=26),
  axis.text.x        = element_text(color="black", size=26),
  axis.text.y        = element_text(color="black", size=26),
  ## Line of the axis
  axis.line          = element_line(color="black", size= 1),
  axis.line.x        = element_line(color="black", size= 1),
  axis.line.y        = element_line(color="black", size= 1),
  ## Ticks of axis
  axis.ticks         = element_line(color="black", size= 1),
  axis.ticks.x       = element_line(color="black", size= 1),
  axis.ticks.y       = element_line(color="black", size= 1),
  axis.ticks.length  = unit(0.2,"cm"),
  ## Legend
  legend.title       = element_text(size=10),
  legend.text        = element_text(size=10),
  legend.background  = element_rect( colour ="black", size=0.8, linetype="solid"),
  legend.position    ="none",
  #legend.position="right",
  #legend.position = c(0.6, 0.9),
  ## Panel
  panel.grid.major   = element_line(color="white"),
  panel.grid.minor   = element_line(color="white"),
  panel.background   = element_blank() ,
  ## Margins
  plot.margin       = unit(c(1,1,1,1), "cm"))
####


# EPIGENETICS ----
# EPIGENETICS F14R ----
# load("D:/Decicomp/R/Melthilkit projects/coverage8X/meth_F14R.RData")

file.list <- list(
  # C3
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C3_T0_S1_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C3_T0_S2_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C3_T0_S3_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C3_T0_S4_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C3_T0_S5_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C3_T0_S6_DNA.myCpG.txt.gz",
  
  # C7
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C7_T0_S1_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C7_T0_S2_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C7_T0_S3_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C7_T0_S4_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C7_T0_S5_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C7_T0_S6_DNA.myCpG.txt.gz",
  
  #C8
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C8_T0_S1_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C8_T0_S2_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C8_T0_S3_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C8_T0_S4_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C8_T0_S5_DNA.myCpG.txt.gz",
  "E:/Decicomp/R/Melthilkit projects/coverage8X/F14R_C8_T0_S6_DNA.myCpG.txt.gz")

myobj_F14R <- methRead(file.list,
                       sample.id=list("F14R_C3_T0_S1_DNA","F14R_C3_T0_S2_DNA", "F14R_C3_T0_S3_DNA","F14R_C3_T0_S4_DNA", "F14R_C3_T0_S5_DNA", "F14R_C3_T0_S6_DNA",
                                      "F14R_C7_T0_S1_DNA","F14R_C7_T0_S2_DNA", "F14R_C7_T0_S3_DNA","F14R_C7_T0_S4_DNA", "F14R_C7_T0_S5_DNA", "F14R_C7_T0_S6_DNA",
                                      "F14R_C8_T0_S1_DNA","F14R_C8_T0_S2_DNA", "F14R_C8_T0_S3_DNA","F14R_C8_T0_S4_DNA", "F14R_C8_T0_S5_DNA", "F14R_C8_T0_S6_DNA"),
                       pipeline = "amp",
                       assembly="crg",
                       treatment=c(1,1,1,1,1,1, 
                                   1,1,1,1,1,1, 
                                   1,1,1,1,1,1),
                       context="CpG",
                       mincov = 8)

filtered.myobj_F14R          <- filterByCoverage(myobj_F14R, lo.count=8, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
normalize.myobj_F14R         <- normalizeCoverage(filtered.myobj_F14R, method = "median")
meth_F14R                    <- methylKit::unite(normalize.myobj_F14R, destrand=T)

sample_pca_epigenetics            <- PCASamples(meth_F14R, scale=F, obj.return = TRUE)
summary_sample_pca_epigenetics    <- summary(sample_pca_epigenetics)

coldata_epigenetics <- rownames(sample_pca_epigenetics$x) %>% as.data.frame() %>% 
  dplyr::rename("sample"=".") %>% 
  dplyr::mutate(group = case_when(
    grepl("F14R_C3_T0", sample) ~ "F14R_C3",
    grepl("F14R_C7_T0", sample) ~ "F14R_C7",
    grepl("F14R_C8_T0", sample) ~ "F14R_C8", TRUE ~ NA_character_)) %>% 
  tibble::column_to_rownames(var = "sample")

pca_scores_epigenetics          <- data.frame(PC1 = sample_pca_epigenetics$x[, 1], 
                                              PC2 = sample_pca_epigenetics$x[, 2], 
                                              PC3 = sample_pca_epigenetics$x[, 3],
                                              PC4 = sample_pca_epigenetics$x[, 4])

pca_scores_epigenetic_coldata_F14R <- cbind(pca_scores_epigenetics, coldata_epigenetics) %>% 
  tibble::rownames_to_column(var = "sample")  %>%   
  dplyr::mutate(sample = gsub("_DNA", "", sample))

PCA_F14R_epigenetics  <- ggplot(data = pca_scores_epigenetic_coldata_F14R, aes(x = PC1 , y = PC4)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = group, color=group), linetype = 1)  +
  geom_point(aes(fill = group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (14.0%)", y = "PC4 (7.8%)") +
  #scale_y_continuous(limits=c(-700,700)) +
  #scale_x_continuous(limits=c(-500,500)) +
  #scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "F14R_C7"="#8faadc", "F14R_C8"="#2f5597", "F14R_C3"="#dae3f3" )) +
  scale_color_manual(values=c(  "F14R_C7"="#8faadc", "F14R_C8"="#2f5597", "F14R_C3"="#dae3f3" )) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_F14R_epigenetics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_F14R_epigenetics_PC1_PC2.tiff", PCA_F14R_epigenetics, width = 10, height = 8, dpi = 600)
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/pca/PCA_F14R_epigenetics_2_vs_3.tiff", PCA_F14R_epigenetics, width = 10, height = 8, dpi = 600)


Coordenates_all_cpg_meth_F14R <- data.frame(
  chr   = meth_F14R$chr, 
  start = meth_F14R$start, 
  end   = meth_F14R$end)

Methylation_meth_F14R <- percMethylation(meth_F14R)

Betas_all_cpg_meth_F14R <- cbind(Coordenates_all_cpg_meth_F14R, Methylation_meth_F14R) %>%
  dplyr::select(-end) %>%
  mutate(pos = paste(chr, start, sep = "_")) %>%
  dplyr::select(-chr, -start) %>%
  column_to_rownames(var = "pos")

dist_mat <- dist(t(Betas_all_cpg_meth_F14R), method = "euclidean")
hc       <- hclust(dist_mat, method = "ward.D2")

sample_colors <- c(
  "F14R_C3_T0_S1_DNA" = "#dae3f3",
  "F14R_C3_T0_S2_DNA" = "#dae3f3",
  "F14R_C3_T0_S3_DNA" = "#dae3f3",
  "F14R_C3_T0_S4_DNA" = "#dae3f3",
  "F14R_C3_T0_S5_DNA" = "#dae3f3",
  "F14R_C3_T0_S6_DNA" = "#dae3f3",
  
  "F14R_C7_T0_S1_DNA" = "#8faadc",
  "F14R_C7_T0_S2_DNA" = "#8faadc",
  "F14R_C7_T0_S3_DNA" = "#8faadc",
  "F14R_C7_T0_S4_DNA" = "#8faadc",
  "F14R_C7_T0_S5_DNA" = "#8faadc",
  "F14R_C7_T0_S6_DNA" = "#8faadc",
  
  "F14R_C8_T0_S1_DNA" = "#2f5597",
  "F14R_C8_T0_S2_DNA" = "#2f5597",
  "F14R_C8_T0_S3_DNA" = "#2f5597",
  "F14R_C8_T0_S4_DNA" = "#2f5597",
  "F14R_C8_T0_S5_DNA" = "#2f5597",
  "F14R_C8_T0_S6_DNA" = "#2f5597")


dend <- as.dendrogram(hc)

# Orden real del dendrograma
labels_dend <- labels(dend)

# Reordenar colores según ese orden
label_colors <- sample_colors[labels_dend]

# Aplicar colores a labels
labels_colors(dend) <- label_colors

# 🔥 NEGRITA (forma correcta en dendextend)
dend <- set(dend, "labels_cex", 1.5)

par(mar = c(14, 4, 4, 2)) 
plot(dend, main = " ", ylab = "Height", cex = 1)


# Plot
tiff("C:/Users/avaldivi/Desktop/BMC_Biology_revision/pca/F14R_dendrogram.tiff",
     width = 3000,
     height = 2500,
     res = 300)

par(mar = c(14, 4, 4, 2))

plot(dend,
     main = " ",
     ylab = "Height",
     cex = 1)

dev.off()




# EPIGENETICS H2D ----
# load("D:/Decicomp/R/Melthilkit projects/coverage8X/meth_H2D.RData")

file.list <- list(
  # C3
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C3_T0_S1_DNA.myCpG.txt.gz",
  #"D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C3_T0_S2_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C3_T0_S3_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C3_T0_S4_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C3_T0_S5_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C3_T0_S6_DNA.myCpG.txt.gz",
  
  # C7
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C7_T0_S1_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C7_T0_S2_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C7_T0_S3_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C7_T0_S4_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C7_T0_S5_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C7_T0_S6_DNA.myCpG.txt.gz",
  
  #C8
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C8_T0_S1_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C8_T0_S2_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C8_T0_S3_DNA.myCpG.txt.gz",
  #"D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C8_T0_S4_DNA.myCpG.txt.gz",
  #"D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C8_T0_S5_DNA.myCpG.txt.gz",
  "D:/Decicomp/R/Melthilkit projects/coverage8X/H2D_C8_T0_S6_DNA.myCpG.txt.gz")

myobj_H2D <- methRead(file.list,
                      sample.id=list("H2D_C3_T0_S1_DNA",                    "H2D_C3_T0_S3_DNA","H2D_C3_T0_S4_DNA", "H2D_C3_T0_S5_DNA", "H2D_C3_T0_S6_DNA",
                                     "H2D_C7_T0_S1_DNA","H2D_C7_T0_S2_DNA", "H2D_C7_T0_S3_DNA","H2D_C7_T0_S4_DNA", "H2D_C7_T0_S5_DNA", "H2D_C7_T0_S6_DNA",
                                     "H2D_C8_T0_S1_DNA","H2D_C8_T0_S2_DNA", "H2D_C8_T0_S3_DNA",                                        "H2D_C8_T0_S6_DNA"),
                      pipeline = "amp",
                      assembly="crg",
                      treatment=c(1,1,1,1,1, 
                                  1,1,1,1,1,1, 
                                  1,1,1,1),
                      context="CpG",
                      mincov = 8)

filtered.myobj_H2D          <- filterByCoverage(myobj_H2D, lo.count=8, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
normalize.myobj_H2D         <- normalizeCoverage(filtered.myobj_H2D, method = "median")
meth_H2D                    <- methylKit::unite(normalize.myobj_H2D, destrand=T)


sample_pca_epigenetics            <- PCASamples(meth_H2D, scale=F, obj.return = TRUE)
summary_sample_pca_epigenetics    <- summary(sample_pca_epigenetics)

coldata_epigenetics <- rownames(sample_pca_epigenetics$x) %>% as.data.frame() %>% 
  dplyr::rename("sample"=".") %>% 
  dplyr::mutate(group = case_when(
    grepl("H2D_C3_T0", sample) ~ "H2D_C3",
    grepl("H2D_C7_T0", sample) ~ "H2D_C7",
    grepl("H2D_C8_T0", sample) ~ "H2D_C8", TRUE ~ NA_character_)) %>% 
  tibble::column_to_rownames(var = "sample")

pca_scores_epigenetics          <- data.frame(PC1 = sample_pca_epigenetics$x[, 1], 
                                              PC2 = sample_pca_epigenetics$x[, 2], 
                                              PC3 = sample_pca_epigenetics$x[, 3],
                                              PC4 = sample_pca_epigenetics$x[, 4])
                                          

pca_scores_epigenetic_coldata_H2D <- cbind(pca_scores_epigenetics, coldata_epigenetics) %>% 
  tibble::rownames_to_column(var = "sample")  %>%   
  dplyr::mutate(sample = gsub("_DNA", "", sample))

PCA_H2D_epigenetics  <- ggplot(data = pca_scores_epigenetic_coldata_H2D, aes(x = PC1 , y = PC4)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = group, color=group), linetype = 1)  +
  geom_point(aes(fill = group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (10.8%)", y = "PC4 (8.6%)") +
  #scale_y_continuous(limits=c(-700,700)) +
  #scale_x_continuous(limits=c(-500,500)) +
  # scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "H2D_C7"="#A9D18E", "H2D_C8"="#385700", "H2D_C3"="#E2F0D9")) +
  scale_color_manual(values=c(  "H2D_C7"="#A9D18E", "H2D_C8"="#385700", "H2D_C3"="#E2F0D9")) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_H2D_epigenetics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_H2D_epigenetics.tiff", PCA_H2D_epigenetics, width = 10, height = 8, dpi = 600)
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/pca/PCA_H2D_epigenetics_2_vs_3.tiff", PCA_H2D_epigenetics, width = 10, height = 8, dpi = 600)


Coordenates_all_cpg_meth_H2D <- data.frame(
  chr   = meth_H2D$chr, 
  start = meth_H2D$start, 
  end   = meth_H2D$end)

Methylation_meth_H2D <- percMethylation(meth_H2D)

Betas_all_cpg_meth_H2D <- cbind(Coordenates_all_cpg_meth_H2D, Methylation_meth_H2D) %>%
  dplyr::select(-end) %>%
  mutate(pos = paste(chr, start, sep = "_")) %>%
  dplyr::select(-chr, -start) %>%
  column_to_rownames(var = "pos")

dist_mat <- dist(t(Betas_all_cpg_meth_H2D), method = "euclidean")
hc       <- hclust(dist_mat, method = "ward.D2")

sample_colors <- c(
  "H2D_C3_T0_S1_DNA" = "#E2F0D9",
  "H2D_C3_T0_S2_DNA" = "#E2F0D9",
  "H2D_C3_T0_S3_DNA" = "#E2F0D9",
  "H2D_C3_T0_S4_DNA" = "#E2F0D9",
  "H2D_C3_T0_S5_DNA" = "#E2F0D9",
  "H2D_C3_T0_S6_DNA" = "#E2F0D9",
  
  "H2D_C7_T0_S1_DNA" = "#A9D18E",
  "H2D_C7_T0_S2_DNA" = "#A9D18E",
  "H2D_C7_T0_S3_DNA" = "#A9D18E",
  "H2D_C7_T0_S4_DNA" = "#A9D18E",
  "H2D_C7_T0_S5_DNA" = "#A9D18E",
  "H2D_C7_T0_S6_DNA" = "#A9D18E",
  
  "H2D_C8_T0_S1_DNA" ="#385700",
  "H2D_C8_T0_S2_DNA" ="#385700",
  "H2D_C8_T0_S3_DNA" ="#385700",
  "H2D_C8_T0_S4_DNA" ="#385700",
  "H2D_C8_T0_S5_DNA" ="#385700",
  "H2D_C8_T0_S6_DNA" ="#385700")


dend <- as.dendrogram(hc)

# Orden real del dendrograma
labels_dend <- labels(dend)

# Reordenar colores según ese orden
label_colors <- sample_colors[labels_dend]

# Aplicar colores a labels
labels_colors(dend) <- label_colors

# 🔥 NEGRITA (forma correcta en dendextend)
dend <- set(dend, "labels_cex", 1.5)

par(mar = c(14, 4, 4, 2)) 
plot(dend, main = " ", ylab = "Height", cex = 1)


# Plot
tiff("C:/Users/avaldivi/Desktop/BMC_Biology_revision/pca/H2D_dendrogram.tiff",
     width = 3000,
     height = 2500,
     res = 300)

par(mar = c(14, 4, 4, 2))

plot(dend,
     main = " ",
     ylab = "Height",
     cex = 1)

dev.off()

