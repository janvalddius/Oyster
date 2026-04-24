# load("D:/Decicomp/R/projects/Immunome_good_analysis_transcriptome.RData")

set.seed(123)
random_numbers <- runif(5)
print(random_numbers)

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


##### DATA BASES   ----  
# ROSETA  
# Biomart list
all_genes_enseble_from_biomart <- read.delim("D:/Gestinnov/1) Oyster proyect/R in datarmor/Roseta_stone/all_genes_enseble_from_biomart.txt")
all_genes_enseble_from_biomart <- as.data.frame(all_genes_enseble_from_biomart[,c(1)])
colnames(all_genes_enseble_from_biomart) <- c("Gene.stable.ID")
dim(all_genes_enseble_from_biomart) # 36012     

# Gene list IHPE
Annot_Cg_Ros_IHPEV1           <- read_excel("D:/Gestinnov/1) Oyster proyect/R in datarmor/Roseta_stone/Annot_Cg_Ros_IHPEV1.xlsx",  sheet="Annot_B2GO")
Annot_Cg_Ros_IHPEV1           <- Annot_Cg_Ros_IHPEV1[,c(2,4)]
colnames(Annot_Cg_Ros_IHPEV1) <- c("Gene.stable.ID", "description")
Roseta                        <- merge(x = all_genes_enseble_from_biomart, y = Annot_Cg_Ros_IHPEV1, by = "Gene.stable.ID", all.x = TRUE )
colnames(Roseta)              <-  c("Gene", "Description")
Annot_Cg_Ros_IHPEV1           <- Annot_Cg_Ros_IHPEV1 %>% rename(Gene= Gene.stable.ID)

# MWU_GO      
List_genes           <- read.delim("D:/Decicomp/R/GO_terms/all_go.tab", header=F) %>% as.data.frame() %>% dplyr::select(1) %>% dplyr::rename("Gene"= "V1")
GO_terms_big         <- read_excel("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_correlates/GO_terms_big.xlsx", sheet= "Feuil3")

# KEGGS M. gigas
KEGGS_Magallana_gigas <- read_excel("D:/Decicomp/R/keegs_extraction/KEGGS_Magallas_gigas.xlsx") %>%
                         mutate(code = as.numeric(code)) # %>% mutate(code = paste0("crg", code))

KEGG_pathways_oyster_all_ld <- read.delim("D:/Decicomp/R/keegs_extraction/KEGG_pathways_oyster_all_ld.txt") %>%
                               mutate(code = as.numeric(code)) %>% rename("KEGG.name"="Pathway_term")

# Genes Conversion GXX to LOC 
Genes_conversion <- read_excel("D:/Decicomp/R/Genes_conversion/Genes_conversion_latest.xlsx", sheet = "Converted") %>%  
                    dplyr::select(1,2,3) %>% mutate(Gene_LOC = as.character(Gene_LOC)) %>%  as.data.frame()  

Gene_conversion_Manu <- read.delim2("D:/Decicomp/R/Genes_conversion/annot_roslin_cds_ok_short__v__GCF_963853765.1_xbMagGiga1.1_short.tsv")

### DATA CLEAN OF OUTLIERS  ----
countData            <- read_excel("D:/Decicomp/Decicomp_matrix.xlsx", sheet = "matrix_outliers")
countData            <- as.data.frame(countData)
row.names(countData) <- countData[, 1]
countData            <- countData[,-1]
countData            <- countData %>% #dplyr::select(matches(grep("T0|T3|T6|T12", colnames(countData), value = TRUE))) %>% 
                        dplyr::select(-c("F14R_C3_T24_S5", "F14R_C3_T24_S6",
                                         "F14R_C7_T0_S3",  "F14R_C7_T24_S4",
                                         "F14R_C8_T6_S3",  "F14R_C8_T24_S1",
                                         "H2D_C3_T24_S4")) 
coldata <- colnames(countData) %>% 
  as.data.frame() %>%
  dplyr::rename("sample"=".") %>%
  dplyr::mutate(group = case_when(
    grepl("H2D_C3_T0",  sample) ~  "H2D_C3_T0",
    grepl("H2D_C3_T3",  sample) ~  "H2D_C3_T3",
    grepl("H2D_C3_T6",  sample) ~  "H2D_C3_T6",
    grepl("H2D_C3_T12", sample) ~  "H2D_C3_T12",
    grepl("H2D_C3_T24", sample) ~  "H2D_C3_T24",
    
    grepl("H2D_C7_T0",  sample) ~  "H2D_C7_T0",
    grepl("H2D_C7_T3",  sample) ~  "H2D_C7_T3",
    grepl("H2D_C7_T6",  sample) ~  "H2D_C7_T6",
    grepl("H2D_C7_T12", sample) ~  "H2D_C7_T12",
    grepl("H2D_C7_T24", sample) ~  "H2D_C7_T24",
    
    grepl("H2D_C8_T0",  sample) ~  "H2D_C8_T0", 
    grepl("H2D_C8_T3",  sample) ~  "H2D_C8_T3",
    grepl("H2D_C8_T6",  sample) ~  "H2D_C8_T6",
    grepl("H2D_C8_T12", sample) ~  "H2D_C8_T12",
    grepl("H2D_C8_T24", sample) ~  "H2D_C8_T24",
    
    grepl("F14R_C3_T0", sample) ~  "F14R_C3_T0",
    grepl("F14R_C3_T3", sample) ~  "F14R_C3_T3",
    grepl("F14R_C3_T6", sample) ~  "F14R_C3_T6",
    grepl("F14R_C3_T12",sample) ~  "F14R_C3_T12",
    grepl("F14R_C3_T24",sample) ~  "F14R_C3_T24",
    
    grepl("F14R_C7_T0", sample) ~  "F14R_C7_T0",
    grepl("F14R_C7_T3", sample) ~  "F14R_C7_T3",
    grepl("F14R_C7_T6", sample) ~  "F14R_C7_T6",
    grepl("F14R_C7_T12",sample) ~  "F14R_C7_T12",
    grepl("F14R_C7_T24",sample) ~  "F14R_C7_T24",
    
    grepl("F14R_C8_T0", sample) ~  "F14R_C8_T0", 
    grepl("F14R_C8_T3", sample) ~  "F14R_C8_T3",
    grepl("F14R_C8_T6", sample) ~  "F14R_C8_T6",
    grepl("F14R_C8_T12",sample) ~  "F14R_C8_T12", 
    grepl("F14R_C8_T24",sample) ~  "F14R_C8_T24", TRUE ~ NA_character_)) %>% 
    tibble::column_to_rownames(var = "sample")

coldata$group        <- as.factor(coldata$group)

dds                  <- DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ group)
keep                 <- rowSums(counts(dds)) >= 15 # dim
dds                  <- dds[keep,]
dds                  <- DESeq(dds)
sizeFactors(dds)
resultsNames(dds) # importante
normalizedCounts    <- counts(dds, normalized = TRUE) # dim(normalizedCounts)
sampleData          <- colData(dds)
dds_norm            <- vst(dds, blind = F) # Blind =F, to take into account the experimental design

# Correlation matrix to detect outliers
counts_norm_F14R           <- assay(dds_norm) %>% as.data.frame() %>% dplyr::select(grep("F14R", colnames(normalizedCounts), value = TRUE))
counts_norm_F14R_outliers  <- cor(counts_norm_F14R)

counts_norm_F14R_outliers_heatmap     <- Heatmap(counts_norm_F14R_outliers, col = magma(10), name = "Correlation", 
                                                        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6)) 
counts_norm_F14R_outliers_heatmap


counts_norm_H2D           <- assay(dds_norm) %>% as.data.frame() %>% dplyr::select(grep("H2D", colnames(normalizedCounts), value = TRUE))
counts_norm_H2D_outliers  <- cor(counts_norm_H2D)

counts_norm_H2D_outliers_heatmap     <- Heatmap(counts_norm_H2D_outliers, col = magma(10), name = "Correlation", 
                                                 row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6)) 
counts_norm_H2D_outliers_heatmap


counts_norm               <- assay(dds_norm) %>% as.data.frame()
counts_norm_outliers      <- cor(counts_norm)

counts_norm_outliers_heatmap     <- Heatmap(counts_norm_outliers, col = magma(10), name = "Correlation", 
                                                row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6)) 
counts_norm_outliers_heatmap


# F14R Family ----
# C3  ----
### F14R_C3_T3 vs F14R_C3_T0 GO_MWU  ----
F14R_C3_T3_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C3_T3_vs_T0_FC    <- left_join(List_genes, F14R_C3_T3_vs_T0_FC)
#write(write.table(List_genes_F14R_C3_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T3_vs_T0/List_genes_F14R_C3_T3_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
#write(write.table(List_genes_F14R_C3_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T3_vs_T0/Immunome/reference/List_genes_F14R_C3_T3_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T3_vs_T0/Immunome/reference/")
input="List_genes_F14R_C3_T3_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C3_T3 vs F14R_C3_T0 GO_MWU
F14R_C3_T3_vs_T0          <- results(dds, contrast=c("group", "F14R_C3_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C3_T3_vs_T0_sig_DEGS <- F14R_C3_T3_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C3_T3_vs_T0_sig_DEGS)

F14R_C3_T3_vs_T0_sig_DEGS_result <- F14R_C3_T3_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                    Down  = sum(log2FoldChange < 0),
                                    Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    dplyr::rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                    mutate(Family = "F14R", Age = "4", Time = "3") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C3_T3_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T3_vs_T0/Immunome/reference/BP_List_genes_F14R_C3_T3_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C3_T3_vs_T0, by = "Gene") 

gene_counts_F14R_C3_T3_vs_T0_immunome    <- BP_List_genes_F14R_C3_T3_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T3_vs_T0_immunome  <- BP_List_genes_F14R_C3_T3_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C3_T3_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C3_T3_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C3_T3_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C3_T3_vs_T0_sig_DEGS_result        <- merge(F14R_C3_T3_vs_T0_sig_DEGS_result,F14R_C3_T3_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C3_T3_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T3_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C3_T3_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C3_T3_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



### F14R_C3_T6 vs F14R_C3_T0 GO_MWU  ----
F14R_C3_T6_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C3_T6_vs_T0_FC    <- left_join(List_genes, F14R_C3_T6_vs_T0_FC)
write(write.table(List_genes_F14R_C3_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T6_vs_T0/List_genes_F14R_C3_T6_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C3_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T6_vs_T0/Immunome/reference/List_genes_F14R_C3_T6_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T6_vs_T0/Immunome/reference")
input="List_genes_F14R_C3_T6_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C3_T6 vs F14R_C3_T0 GO_MWU
F14R_C3_T6_vs_T0          <- results(dds, contrast=c("group", "F14R_C3_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C3_T6_vs_T0_sig_DEGS <- F14R_C3_T6_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C3_T6_vs_T0_sig_DEGS)

F14R_C3_T6_vs_T0_sig_DEGS_result <- F14R_C3_T6_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                    Down  = sum(log2FoldChange < 0),
                                    Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>%  rownames_to_column(var="DEGS") %>%
                                    mutate(Family = "F14R", Age = "4", Time = "6") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C3_T6 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C3_T6_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T6_vs_T0/Immunome/reference/BP_List_genes_F14R_C3_T6_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C3_T6_vs_T0, by = "Gene") 

gene_counts_F14R_C3_T6_vs_T0_immunome    <- BP_List_genes_F14R_C3_T6_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T6_vs_T0_immunome  <- BP_List_genes_F14R_C3_T6_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C3_T6_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C3_T6_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C3_T6_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C3_T6_vs_T0_sig_DEGS_result         <- merge(F14R_C3_T6_vs_T0_sig_DEGS_result, F14R_C3_T6_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C3_T6 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C3_T6_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T6_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C3_T6_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C3_T6_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C3_T12 vs F14R_C3_T0 GO_MWU  ----
F14R_C3_T12_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                           as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                           filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C3_T12_vs_T0_FC    <- left_join(List_genes, F14R_C3_T12_vs_T0_FC)
write(write.table(List_genes_F14R_C3_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T12_vs_T0/List_genes_F14R_C3_T12_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C3_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T12_vs_T0/Immunome/reference/List_genes_F14R_C3_T12_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T12_vs_T0/Immunome/reference/")
input="List_genes_F14R_C3_T12_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C3_T12 vs F14R_C3_T0 GO_MWU
F14R_C3_T12_vs_T0          <- results(dds, contrast=c("group", "F14R_C3_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C3_T12_vs_T0_sig_DEGS <- F14R_C3_T12_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C3_T12_vs_T0_sig_DEGS)

F14R_C3_T12_vs_T0_sig_DEGS_result <- F14R_C3_T12_vs_T0_sig_DEGS %>%
                                     summarise(Up    = sum(log2FoldChange > 0),
                                     Down  = sum(log2FoldChange < 0),
                                     Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                     rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                     mutate(Family = "F14R", Age = "4", Time = "12") %>%
                                     dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C3_T12 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C3_T12_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T12_vs_T0/Immunome/reference/BP_List_genes_F14R_C3_T12_vs_T0.txt") %>% 
                                             dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C3_T12_vs_T0, by = "Gene") 

gene_counts_F14R_C3_T12_vs_T0_immunome    <- BP_List_genes_F14R_C3_T12_vs_T0_immunome %>%
                                             dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T12_vs_T0_immunome  <- BP_List_genes_F14R_C3_T12_vs_T0_immunome %>% 
                                             left_join(gene_counts_F14R_C3_T12_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C3_T12_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C3_T12_vs_T0_immunome %>% 
                                             filter(Gene_significant == "Yes") %>% 
                                             summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                             mutate(Family = "F14R")

F14R_C3_T12_vs_T0_sig_DEGS_result         <- merge(F14R_C3_T12_vs_T0_sig_DEGS_result, F14R_C3_T12_vs_T0_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T12 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C3_T12_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T12_vs_T0.txt", sep="") %>% 
                                             mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                             mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                             left_join(BP_List_genes_F14R_C3_T12_vs_T0_immunome, by = c("term", "name")) %>% 
                                             mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                             mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                             mutate(Comparison = "F14R_C3_T12_vs_T0") %>% 
                                             left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C3_T24 vs F14R_C3_T0 GO_MWU  ----
F14R_C3_T24_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                           as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                           filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C3_T24_vs_T0_FC    <- left_join(List_genes, F14R_C3_T24_vs_T0_FC)
write(write.table(List_genes_F14R_C3_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T24_vs_T0/List_genes_F14R_C3_T24_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C3_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T24_vs_T0/Immunome/reference/List_genes_F14R_C3_T24_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T24_vs_T0/Immunome/reference/")
input="List_genes_F14R_C3_T24_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C3_T24 vs F14R_C3_T0 GO_MWU
F14R_C3_T24_vs_T0          <- results(dds, contrast=c("group", "F14R_C3_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C3_T24_vs_T0_sig_DEGS <- F14R_C3_T24_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C3_T24_vs_T0_sig_DEGS)

F14R_C3_T24_vs_T0_sig_DEGS_result <- F14R_C3_T24_vs_T0_sig_DEGS %>%
                                     summarise(Up    = sum(log2FoldChange > 0),
                                               Down  = sum(log2FoldChange < 0),
                                               Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                     rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                     mutate(Family = "F14R", Age = "4", Time = "24") %>%
                                     dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C3_T24 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C3_T24_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T24_vs_T0/Immunome/reference/BP_List_genes_F14R_C3_T24_vs_T0.txt") %>% 
                                             dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C3_T24_vs_T0, by = "Gene") 

gene_counts_F14R_C3_T24_vs_T0_immunome    <- BP_List_genes_F14R_C3_T24_vs_T0_immunome %>%
                                             dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T24_vs_T0_immunome  <- BP_List_genes_F14R_C3_T24_vs_T0_immunome %>% 
                                             left_join(gene_counts_F14R_C3_T24_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C3_T24_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C3_T24_vs_T0_immunome %>% 
                                             filter(Gene_significant == "Yes") %>% 
                                             summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                              mutate(Family = "F14R")

F14R_C3_T24_vs_T0_sig_DEGS_result         <- merge(F14R_C3_T24_vs_T0_sig_DEGS_result, F14R_C3_T24_vs_T0_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T24 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C3_T24_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T24_vs_T0.txt", sep="") %>% 
                                             mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                             mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                             left_join(BP_List_genes_F14R_C3_T24_vs_T0_immunome, by = c("term", "name")) %>% 
                                             mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                             mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                             mutate(Comparison = "F14R_C3_T24_vs_T0") %>% 
                                             left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


# C7 ----

### F14R_C7_T0 vs F14R_C3_T0 GO_MWU  ----
F14R_C7_T0_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C7_T0_vs_T0_FC    <- left_join(List_genes, F14R_C7_T0_vs_T0_FC)
write(write.table(List_genes_F14R_C7_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T0_vs_T0/List_genes_F14R_C7_T0_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C7_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T0_vs_T0/Immunome/reference/List_genes_F14R_C7_T0_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T0_vs_T0/Immunome/reference/")
input="List_genes_F14R_C7_T0_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C7_T0 vs F14R_C3_T0 GO_MWU
F14R_C7_T0_vs_T0          <- results(dds, contrast=c("group", "F14R_C7_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C7_T0_vs_T0_sig_DEGS <- F14R_C7_T0_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C7_T0_vs_T0_sig_DEGS)

F14R_C7_T0_vs_T0_sig_DEGS_result <- F14R_C7_T0_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                              Down  = sum(log2FoldChange < 0),
                                              Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                    mutate(Family = "F14R", Age = "16", Time = "0") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C7_T0 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C7_T0_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T0_vs_T0/Immunome/reference/BP_List_genes_F14R_C7_T0_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C7_T0_vs_T0, by = "Gene") 

gene_counts_F14R_C7_T0_vs_T0_immunome    <- BP_List_genes_F14R_C7_T0_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T0_vs_T0_immunome  <- BP_List_genes_F14R_C7_T0_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C7_T0_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C7_T0_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C7_T0_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C7_T0_vs_T0_sig_DEGS_result         <- merge(F14R_C7_T0_vs_T0_sig_DEGS_result, F14R_C7_T0_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C7_T0 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C7_T0_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T0_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C7_T0_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C7_T0_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C7_T3 vs F14R_C3_T0 GO_MWU  ----
F14R_C7_T3_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C7_T3_vs_T0_FC    <- left_join(List_genes, F14R_C7_T3_vs_T0_FC)
write(write.table(List_genes_F14R_C7_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T3_vs_T0/List_genes_F14R_C7_T3_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C7_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T3_vs_T0/Immunome/reference/List_genes_F14R_C7_T3_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T3_vs_T0/Immunome/reference/")
input="List_genes_F14R_C7_T3_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C7_T3 vs F14R_C3_T0 GO_MWU
F14R_C7_T3_vs_T0          <- results(dds, contrast=c("group", "F14R_C7_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C7_T3_vs_T0_sig_DEGS <- F14R_C7_T3_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C7_T3_vs_T0_sig_DEGS)

F14R_C7_T3_vs_T0_sig_DEGS_result <- F14R_C7_T3_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                              Down  = sum(log2FoldChange < 0),
                                              Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                    mutate(Family = "F14R", Age = "16", Time = "3") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C7_T3 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C7_T3_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T3_vs_T0/Immunome/reference/BP_List_genes_F14R_C7_T3_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C7_T3_vs_T0, by = "Gene") 

gene_counts_F14R_C7_T3_vs_T0_immunome    <- BP_List_genes_F14R_C7_T3_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T3_vs_T0_immunome  <- BP_List_genes_F14R_C7_T3_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C7_T3_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C7_T3_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C7_T3_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C7_T3_vs_T0_sig_DEGS_result         <- merge(F14R_C7_T3_vs_T0_sig_DEGS_result, F14R_C7_T3_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C7_T3 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C7_T3_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T3_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C7_T3_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C7_T3_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



### F14R_C7_T6 vs F14R_C3_T0 GO_MWU  ----
F14R_C7_T6_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C7_T6_vs_T0_FC    <- left_join(List_genes, F14R_C7_T6_vs_T0_FC)
write(write.table(List_genes_F14R_C7_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T6_vs_T0/List_genes_F14R_C7_T6_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C7_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T6_vs_T0/Immunome/reference/List_genes_F14R_C7_T6_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T6_vs_T0/Immunome/reference/")
input="List_genes_F14R_C7_T6_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C7_T6 vs F14R_C3_T0 GO_MWU
F14R_C7_T6_vs_T0          <- results(dds, contrast=c("group", "F14R_C7_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C7_T6_vs_T0_sig_DEGS <- F14R_C7_T6_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C7_T6_vs_T0_sig_DEGS)

F14R_C7_T6_vs_T0_sig_DEGS_result <- F14R_C7_T6_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                    Down  = sum(log2FoldChange < 0),
                                    Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                    mutate(Family = "F14R", Age = "16", Time = "6") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C7_T6 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C7_T6_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T6_vs_T0/Immunome/reference/BP_List_genes_F14R_C7_T6_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C7_T6_vs_T0, by = "Gene") 

gene_counts_F14R_C7_T6_vs_T0_immunome    <- BP_List_genes_F14R_C7_T6_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T6_vs_T0_immunome  <- BP_List_genes_F14R_C7_T6_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C7_T6_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C7_T6_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C7_T6_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C7_T6_vs_T0_sig_DEGS_result         <- merge(F14R_C7_T6_vs_T0_sig_DEGS_result, F14R_C7_T6_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C7_T6 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C7_T6_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T6_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C7_T6_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C7_T6_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C7_T12 vs F14R_C3_T0 GO_MWU  ----
F14R_C7_T12_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                           as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                           filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C7_T12_vs_T0_FC    <- left_join(List_genes, F14R_C7_T12_vs_T0_FC)
write(write.table(List_genes_F14R_C7_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T12_vs_T0/List_genes_F14R_C7_T12_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C7_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T12_vs_T0/Immunome/reference/List_genes_F14R_C7_T12_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T12_vs_T0/Immunome/reference")
input="List_genes_F14R_C7_T12_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C7_T12 vs F14R_C3_T0 GO_MWU
F14R_C7_T12_vs_T0          <- results(dds, contrast=c("group", "F14R_C7_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C7_T12_vs_T0_sig_DEGS <- F14R_C7_T12_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C7_T12_vs_T0_sig_DEGS)

F14R_C7_T12_vs_T0_sig_DEGS_result <- F14R_C7_T12_vs_T0_sig_DEGS %>%
                                     summarise(Up    = sum(log2FoldChange > 0),
                                               Down  = sum(log2FoldChange < 0),
                                               Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                     rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                     mutate(Family = "F14R", Age = "16", Time = "12") %>%
                                     dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C7_T12 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C7_T12_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T12_vs_T0/Immunome/reference/BP_List_genes_F14R_C7_T12_vs_T0.txt") %>% 
                                             dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C7_T12_vs_T0, by = "Gene") 

gene_counts_F14R_C7_T12_vs_T0_immunome    <- BP_List_genes_F14R_C7_T12_vs_T0_immunome %>%
                                             dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T12_vs_T0_immunome  <- BP_List_genes_F14R_C7_T12_vs_T0_immunome %>% 
                                             left_join(gene_counts_F14R_C7_T12_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C7_T12_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C7_T12_vs_T0_immunome %>% 
                                             filter(Gene_significant == "Yes") %>% 
                                             summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                             mutate(Family = "F14R")

F14R_C7_T12_vs_T0_sig_DEGS_result         <- merge(F14R_C7_T12_vs_T0_sig_DEGS_result, F14R_C7_T12_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C7_T12 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C7_T12_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T12_vs_T0.txt", sep="") %>% 
                                             mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                             mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                             left_join(BP_List_genes_F14R_C7_T12_vs_T0_immunome, by = c("term", "name")) %>% 
                                             mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                             mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                             mutate(Comparison = "F14R_C7_T12_vs_T0") %>% 
                                             left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C7_T24 vs F14R_C3_T0 GO_MWU  ----
F14R_C7_T24_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                           as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                           filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C7_T24_vs_T0_FC    <- left_join(List_genes, F14R_C7_T24_vs_T0_FC)
write(write.table(List_genes_F14R_C7_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T24_vs_T0/List_genes_F14R_C7_T24_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C7_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T24_vs_T0/Immunome/reference/List_genes_F14R_C7_T24_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T24_vs_T0/Immunome/reference/")
input="List_genes_F14R_C7_T24_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C7_T24 vs F14R_C3_T0 GO_MWU
F14R_C7_T24_vs_T0          <- results(dds, contrast=c("group", "F14R_C7_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C7_T24_vs_T0_sig_DEGS <- F14R_C7_T24_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C7_T24_vs_T0_sig_DEGS)

F14R_C7_T24_vs_T0_sig_DEGS_result <- F14R_C7_T24_vs_T0_sig_DEGS %>%
                                     summarise(Up    = sum(log2FoldChange > 0),
                                               Down  = sum(log2FoldChange < 0),
                                               Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                     rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                     mutate(Family = "F14R", Age = "16", Time = "24") %>%
                                     dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C7_T24 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C7_T24_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T24_vs_T0/Immunome/reference/BP_List_genes_F14R_C7_T24_vs_T0.txt") %>% 
                                             dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C7_T24_vs_T0, by = "Gene") 

gene_counts_F14R_C7_T24_vs_T0_immunome    <- BP_List_genes_F14R_C7_T24_vs_T0_immunome %>%
                                             dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T24_vs_T0_immunome  <- BP_List_genes_F14R_C7_T24_vs_T0_immunome %>% 
                                             left_join(gene_counts_F14R_C7_T24_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C7_T24_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C7_T24_vs_T0_immunome %>% 
                                             filter(Gene_significant == "Yes") %>% 
                                             summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                             mutate(Family = "F14R")

F14R_C7_T24_vs_T0_sig_DEGS_result         <- merge(F14R_C7_T24_vs_T0_sig_DEGS_result, F14R_C7_T24_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C7_T24 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C7_T24_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T24_vs_T0.txt", sep="") %>% 
                                             mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                             mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                             left_join(BP_List_genes_F14R_C7_T24_vs_T0_immunome, by = c("term", "name")) %>% 
                                             mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                             mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                             mutate(Comparison = "F14R_C7_T24_vs_T0") %>% 
                                             left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


# C8  ----

### F14R_C8_T0 vs F14R_C3_T0 GO_MWU  ----
F14R_C8_T0_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
  filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T0_vs_T0_FC    <- left_join(List_genes, F14R_C8_T0_vs_T0_FC)
write(write.table(List_genes_F14R_C8_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T0_vs_T0/List_genes_F14R_C8_T0_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C8_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T0_vs_T0/Immunome/reference/List_genes_F14R_C8_T0_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T0_vs_T0/Immunome/reference/")
input="List_genes_F14R_C8_T0_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C8_T0 vs F14R_C3_T0 GO_MWU
F14R_C8_T0_vs_T0          <- results(dds, contrast=c("group", "F14R_C8_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
                             F14R_C8_T0_vs_T0_sig_DEGS <- F14R_C8_T0_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T0_vs_T0_sig_DEGS)

F14R_C8_T0_vs_T0_sig_DEGS_result <- F14R_C8_T0_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                              Down  = sum(log2FoldChange < 0),
                                              Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                    mutate(Family = "F14R", Age = "28", Time = "0") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C8_T0 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C8_T0_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T0_vs_T0/Immunome/reference/BP_List_genes_F14R_C8_T0_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T0_vs_T0, by = "Gene") 

gene_counts_F14R_C8_T0_vs_T0_immunome    <- BP_List_genes_F14R_C8_T0_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T0_vs_T0_immunome  <- BP_List_genes_F14R_C8_T0_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C8_T0_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T0_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_T0_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C8_T0_vs_T0_sig_DEGS_result         <- merge(F14R_C8_T0_vs_T0_sig_DEGS_result, F14R_C8_T0_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C8_T0 vs F14R_C8_T0 GO_MWU
GO_terms_F14R_C8_T0_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T0_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C8_T0_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C8_T0_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 

### F14R_C8_T3 vs F14R_C3_T0 GO_MWU  ----
F14R_C8_T3_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T3_vs_T0_FC    <- left_join(List_genes, F14R_C8_T3_vs_T0_FC)
write(write.table(List_genes_F14R_C8_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T3_vs_T0/List_genes_F14R_C8_T3_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C8_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T3_vs_T0/Immunome/reference/List_genes_F14R_C8_T3_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T3_vs_T0/Immunome/reference/")
input="List_genes_F14R_C8_T3_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C8_T3 vs F14R_C3_T0 GO_MWU
F14R_C8_T3_vs_T0          <- results(dds, contrast=c("group", "F14R_C8_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
                             F14R_C8_T3_vs_T0_sig_DEGS <- F14R_C8_T3_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T3_vs_T0_sig_DEGS)

F14R_C8_T3_vs_T0_sig_DEGS_result <- F14R_C8_T3_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                              Down  = sum(log2FoldChange < 0),
                                              Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                    mutate(Family = "F14R", Age = "28", Time = "3") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C8_T3 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C8_T3_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T3_vs_T0/Immunome/reference/BP_List_genes_F14R_C8_T3_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T3_vs_T0, by = "Gene") 

gene_counts_F14R_C8_T3_vs_T0_immunome    <- BP_List_genes_F14R_C8_T3_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T3_vs_T0_immunome  <- BP_List_genes_F14R_C8_T3_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C8_T3_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T3_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_T3_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C8_T3_vs_T0_sig_DEGS_result         <- merge(F14R_C8_T3_vs_T0_sig_DEGS_result, F14R_C8_T3_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C8_T3 vs F14R_C8_T0 GO_MWU
GO_terms_F14R_C8_T3_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T3_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C8_T3_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C8_T3_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



### F14R_C8_T6 vs F14R_C3_T0 GO_MWU  ----
F14R_C8_T6_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T6_vs_T0_FC    <- left_join(List_genes, F14R_C8_T6_vs_T0_FC)
write(write.table(List_genes_F14R_C8_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T6_vs_T0/List_genes_F14R_C8_T6_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C8_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T6_vs_T0/Immunome/reference/List_genes_F14R_C8_T6_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T6_vs_T0/Immunome/reference/")
input="List_genes_F14R_C8_T6_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C8_T6 vs F14R_C3_T0 GO_MWU
F14R_C8_T6_vs_T0          <- results(dds, contrast=c("group", "F14R_C8_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T6_vs_T0_sig_DEGS <- F14R_C8_T6_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T6_vs_T0_sig_DEGS)

F14R_C8_T6_vs_T0_sig_DEGS_result <- F14R_C8_T6_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                              Down  = sum(log2FoldChange < 0),
                                              Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                    mutate(Family = "F14R", Age = "28", Time = "6") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C8_T6 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C8_T6_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T6_vs_T0/Immunome/reference/BP_List_genes_F14R_C8_T6_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T6_vs_T0, by = "Gene") 

gene_counts_F14R_C8_T6_vs_T0_immunome    <- BP_List_genes_F14R_C8_T6_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T6_vs_T0_immunome  <- BP_List_genes_F14R_C8_T6_vs_T0_immunome %>% 
                                            left_join(gene_counts_F14R_C8_T6_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T6_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_T6_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "F14R")

F14R_C8_T6_vs_T0_sig_DEGS_result         <- merge(F14R_C8_T6_vs_T0_sig_DEGS_result, F14R_C8_T6_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C8_T6 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C8_T6_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T6_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_F14R_C8_T6_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "F14R_C8_T6_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C8_T12 vs F14R_C3_T0 GO_MWU  ----
F14R_C8_T12_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                           as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                           filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T12_vs_T0_FC    <- left_join(List_genes, F14R_C8_T12_vs_T0_FC)
write(write.table(List_genes_F14R_C8_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T12_vs_T0/List_genes_F14R_C8_T12_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C8_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T12_vs_T0/Immunome/reference/List_genes_F14R_C8_T12_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T12_vs_T0/Immunome/reference/")
input="List_genes_F14R_C8_T12_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C8_T12 vs F14R_C3_T0 GO_MWU
F14R_C8_T12_vs_T0          <- results(dds, contrast=c("group", "F14R_C8_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T12_vs_T0_sig_DEGS <- F14R_C8_T12_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T12_vs_T0_sig_DEGS)

F14R_C8_T12_vs_T0_sig_DEGS_result <- F14R_C8_T12_vs_T0_sig_DEGS %>%
                                     summarise(Up    = sum(log2FoldChange > 0),
                                               Down  = sum(log2FoldChange < 0),
                                               Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                     rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                     mutate(Family = "F14R", Age = "28", Time = "12") %>%
                                     dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C8_T12 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C8_T12_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T12_vs_T0/Immunome/reference/BP_List_genes_F14R_C8_T12_vs_T0.txt") %>% 
                                             dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T12_vs_T0, by = "Gene") 

gene_counts_F14R_C8_T12_vs_T0_immunome    <- BP_List_genes_F14R_C8_T12_vs_T0_immunome %>%
                                             dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T12_vs_T0_immunome  <- BP_List_genes_F14R_C8_T12_vs_T0_immunome %>% 
                                             left_join(gene_counts_F14R_C8_T12_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T12_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_T12_vs_T0_immunome %>% 
                                             filter(Gene_significant == "Yes") %>% 
                                             summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                             mutate(Family = "F14R")

F14R_C8_T12_vs_T0_sig_DEGS_result         <- merge(F14R_C8_T12_vs_T0_sig_DEGS_result, F14R_C8_T12_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C8_T12 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C8_T12_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T12_vs_T0.txt", sep="") %>% 
                                             mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                             mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                             left_join(BP_List_genes_F14R_C8_T12_vs_T0_immunome, by = c("term", "name")) %>% 
                                             mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                             mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                             mutate(Comparison = "F14R_C8_T12_vs_T0") %>% 
                                             left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### F14R_C8_T24 vs F14R_C3_T0 GO_MWU  ----
F14R_C8_T24_vs_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                           as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                           filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T24_vs_T0_FC    <- left_join(List_genes, F14R_C8_T24_vs_T0_FC)
write(write.table(List_genes_F14R_C8_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T24_vs_T0/List_genes_F14R_C8_T24_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_F14R_C8_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T24_vs_T0/Immunome/reference/List_genes_F14R_C8_T24_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T24_vs_T0/Immunome/reference/")
input="List_genes_F14R_C8_T24_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use F14R_C8_T24 vs F14R_C3_T0 GO_MWU
F14R_C8_T24_vs_T0          <- results(dds, contrast=c("group", "F14R_C8_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T24_vs_T0_sig_DEGS <- F14R_C8_T24_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T24_vs_T0_sig_DEGS)

F14R_C8_T24_vs_T0_sig_DEGS_result <- F14R_C8_T24_vs_T0_sig_DEGS %>%
                                     summarise(Up    = sum(log2FoldChange > 0),
                                               Down  = sum(log2FoldChange < 0),
                                               Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                     rename("Number"="V1") %>% rownames_to_column(var ="DEGS") %>%
                                     mutate(Family = "F14R", Age = "28", Time = "24") %>%
                                     dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome F14R_C8_T24 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C8_T24_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T24_vs_T0/Immunome/reference/BP_List_genes_F14R_C8_T24_vs_T0.txt") %>% 
                                             dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T24_vs_T0, by = "Gene") 

gene_counts_F14R_C8_T24_vs_T0_immunome    <- BP_List_genes_F14R_C8_T24_vs_T0_immunome %>%
                                             dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T24_vs_T0_immunome  <- BP_List_genes_F14R_C8_T24_vs_T0_immunome %>% 
                                             left_join(gene_counts_F14R_C8_T24_vs_T0_immunome, by = "term") %>% 
                                             mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T24_vs_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_T24_vs_T0_immunome %>% 
                                             filter(Gene_significant == "Yes") %>% 
                                             summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                             mutate(Family = "F14R")

F14R_C8_T24_vs_T0_sig_DEGS_result         <- merge(F14R_C8_T24_vs_T0_sig_DEGS_result, F14R_C8_T24_vs_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C8_T24 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C8_T24_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T24_vs_T0.txt", sep="") %>% 
                                             mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                             mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                             left_join(BP_List_genes_F14R_C8_T24_vs_T0_immunome, by = c("term", "name")) %>% 
                                             mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                             mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                             mutate(Comparison = "F14R_C8_T24_vs_T0") %>% 
                                             left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 

# Heatmap_F14R -----

GO_terms_F14R_C3_T3_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T3_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C3_T6_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T6_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C3_T12_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T12_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_F14R_C3_T24_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C3_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C3_T24_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_F14R_C7_T0_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T0_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C7_T3_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T3_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C7_T6_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T6_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C7_T12_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T12_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_F14R_C7_T24_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C7_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C7_T24_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_F14R_C8_T0_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T0_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C8_T3_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T3_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C8_T6_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T6_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_F14R_C8_T12_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T12_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_F14R_C8_T24_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_F14R_C8_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_F14R_C8_T24_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)

dim(GO_terms_F14R_C3_T3_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C3_T6_vs_T0_immunome_sig)
dim(GO_terms_F14R_C3_T12_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C3_T24_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C7_T0_vs_T0_immunome_sig)
dim(GO_terms_F14R_C7_T3_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C7_T6_vs_T0_immunome_sig)
dim(GO_terms_F14R_C7_T12_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C7_T24_vs_T0_immunome_sig)
dim(GO_terms_F14R_C8_T0_vs_T0_immunome_sig)
dim(GO_terms_F14R_C8_T3_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C8_T6_vs_T0_immunome_sig)
dim(GO_terms_F14R_C8_T12_vs_T0_immunome_sig) 
dim(GO_terms_F14R_C8_T24_vs_T0_immunome_sig) 

# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


GO_F14R_infection_reference <- rbind(GO_terms_F14R_C3_T3_vs_T0_immunome, 
                                     GO_terms_F14R_C3_T6_vs_T0_immunome, 
                                     GO_terms_F14R_C3_T12_vs_T0_immunome, 
                                     GO_terms_F14R_C3_T24_vs_T0_immunome,
                                     GO_terms_F14R_C7_T0_vs_T0_immunome, 
                                     GO_terms_F14R_C7_T3_vs_T0_immunome, 
                                     GO_terms_F14R_C7_T6_vs_T0_immunome,  
                                     GO_terms_F14R_C7_T12_vs_T0_immunome, 
                                     GO_terms_F14R_C7_T24_vs_T0_immunome,
                                     GO_terms_F14R_C8_T0_vs_T0_immunome, 
                                     GO_terms_F14R_C8_T3_vs_T0_immunome, 
                                     GO_terms_F14R_C8_T6_vs_T0_immunome,  
                                     GO_terms_F14R_C8_T12_vs_T0_immunome, 
                                     GO_terms_F14R_C8_T24_vs_T0_immunome) %>%
                               filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% distinct(name, .keep_all = TRUE) %>%
                               dplyr::select(Comparison, name, Enrichment) %>% 
                               pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
                               as.data.frame() %>% mutate_all(~replace(., is.na(.), 0)) %>%
                               mutate(F14R_C3_T0_vs_T0 = 0, F14R_C3_T12_vs_T0 = 0, F14R_C7_T6_vs_T0  = 0) %>%
                               column_to_rownames(var = "name") %>%
                               dplyr::select(13,1,2,14,3,4,5,15,6,7,8,9,10,11,12)
         
# Colum Colors heatmap   
annotation_col <- data.frame(
  Time = c( "T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24"),
  Age  = c(rep("4 months", 5), rep("16 months", 5), rep("28 months", 5)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
                             Age    = c("4 months"   = "#dae3f3", "16 months"  = "#8faadc", "28 months"  ="#2f5597"),
                             Time   = c(T0 = "white",  T3           = "#CD5C5C", T6 = "#B22222", T12="darkred",  T24 ="red"))

rownames(annotation_col) <- colnames(GO_F14R_infection_reference) 


# Gradient enrichment                        
myBreaks <- c(seq(min(GO_F14R_infection_reference), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_F14R_infection_reference)/100, max(GO_F14R_infection_reference), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

GO_F14R_infection_reference_heatmap <- pheatmap(GO_F14R_infection_reference, 
                           cluster_cols = F, 
                           scale = "none",
                           cluster_rows = T, 
                           fontsize_row = 14, 
                           color = mycolor, 
                           breaks = myBreaks,
                           border_color = "black",
                           clustering_distance_rows = "euclidean",
                           show_colnames = F, 
                           show_rownames = T,
                           annotation_col = annotation_col, 
                           #annotation_row = annotation_row,
                           annotation_colors = annotation_colors,
                           #cutree_rows = 4,
                           gaps_col =  c(5,10,15))
GO_F14R_infection_reference_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/GO_F14R_infection_reference_heatmap.png", GO_F14R_infection_reference_heatmap, width = 13, height = 8, dpi = 300)


# Filtering T0
GO_F14R_infection_reference_T0 <- rbind(GO_terms_F14R_C7_T0_vs_T0_immunome, GO_terms_F14R_C8_T0_vs_T0_immunome) %>%
  filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% distinct(name, .keep_all = TRUE) %>%
  dplyr::select(Comparison, name, Enrichment) %>% 
  pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
  as.data.frame() %>% mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(F14R_C3_T0_vs_T0 = 0) %>%
  column_to_rownames(var = "name") %>%
  dplyr::select(3,1,2)


# Colum Colors heatmap   
annotation_col <- data.frame(
 #Time = c("T0", "T0", "T0"),
  Age  = c(rep("4 months", 1), rep("16 months", 1), rep("28 months", 1)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
  Age    = c("4 months"   = "#dae3f3", "16 months"  = "#8faadc", "28 months"  ="#2f5597"))
  #Time   = c(T0 = "white"))

rownames(annotation_col) <- colnames(GO_F14R_infection_reference_T0) 


# Gradient enrichment                        
myBreaks <- c(seq(min(GO_F14R_infection_reference_T0), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_F14R_infection_reference_T0)/100, max(GO_F14R_infection_reference_T0), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

GO_F14R_infection_reference_T0_heatmap <- pheatmap(GO_F14R_infection_reference_T0, 
                                                   cluster_cols = F, 
                                                   scale = "none",
                                                   cluster_rows = T, 
                                                   fontsize_row = 14, 
                                                   color = mycolor, 
                                                   breaks = myBreaks,
                                                   border_color = "black",
                                                   clustering_distance_rows = "euclidean",
                                                   show_colnames = F, 
                                                   show_rownames = T,
                                                   annotation_col = annotation_col, 
                                                   #annotation_row = annotation_row,
                                                   annotation_names_col = FALSE,
                                                   annotation_legend = FALSE,
                                                   legend = F,
                                                   legend_title = NULL,
                                                   annotation_colors = annotation_colors)
#cutree_rows = 4,
#gaps_col =  c(4,9,14))
GO_F14R_infection_reference_T0_heatmap
 # ggsave("D:/Decicomp/R/MOFA_omics/GO_F14R_infection_reference_T0_heatmap.tiff", GO_F14R_infection_reference_T0_heatmap, width = 9, height = 8, dpi = 300)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



# H2D Family ----
# C3  ----
### H2D_C3_T3 vs H2D_C3_T0 GO_MWU  ----
H2D_C3_T3_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                        filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C3_T3_vs_T0_FC    <- left_join(List_genes, H2D_C3_T3_vs_T0_FC)
write(write.table(List_genes_H2D_C3_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T3_vs_T0/List_genes_H2D_C3_T3_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C3_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T3_vs_T0/Immunome/reference/List_genes_H2D_C3_T3_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T3_vs_T0/Immunome/reference/")
input="List_genes_H2D_C3_T3_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C3_T3 vs H2D_C3_T0 GO_MWU
H2D_C3_T3_vs_T0          <- results(dds, contrast=c("group", "H2D_C3_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C3_T3_vs_T0_sig_DEGS <- H2D_C3_T3_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C3_T3_vs_T0_sig_DEGS)

H2D_C3_T3_vs_T0_sig_DEGS_result <- H2D_C3_T3_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                             Down  = sum(log2FoldChange < 0),
                                             Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                             rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                             mutate(Family = "H2D", Age = "4", Time = "3") %>%
                                            dplyr::select(3,4,5,1,2)


## GO terms and genes involved immunome H2D_C3_T3 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C3_T3_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T3_vs_T0/Immunome/reference/BP_List_genes_H2D_C3_T3_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C3_T3_vs_T0, by = "Gene") 

gene_counts_H2D_C3_T3_vs_T0_immunome    <- BP_List_genes_H2D_C3_T3_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T3_vs_T0_immunome  <- BP_List_genes_H2D_C3_T3_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C3_T3_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C3_T3_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C3_T3_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "H2D")

H2D_C3_T3_vs_T0_sig_DEGS_result        <- merge(H2D_C3_T3_vs_T0_sig_DEGS_result,H2D_C3_T3_vs_T0_immunome_DEGS)


## Significant GO terms immunome H2D_C3_T3 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C3_T3_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T3_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C3_T3_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C3_T3_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



### H2D_C3_T6 vs H2D_C3_T0 GO_MWU  ----
H2D_C3_T6_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C3_T6_vs_T0_FC    <- left_join(List_genes, H2D_C3_T6_vs_T0_FC)
write(write.table(List_genes_H2D_C3_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T6_vs_T0/List_genes_H2D_C3_T6_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C3_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T6_vs_T0/Immunome/reference/List_genes_H2D_C3_T6_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T6_vs_T0/Immunome/reference/")
input="List_genes_H2D_C3_T6_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C3_T6 vs H2D_C3_T0 GO_MWU
H2D_C3_T6_vs_T0          <- results(dds, contrast=c("group", "H2D_C3_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C3_T6_vs_T0_sig_DEGS <- H2D_C3_T6_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C3_T6_vs_T0_sig_DEGS)

H2D_C3_T6_vs_T0_sig_DEGS_result <- H2D_C3_T6_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                             Down  = sum(log2FoldChange < 0),
                                             Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                   rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                   mutate(Family = "H2D", Age = "4", Time = "6") %>%
                                   dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C3_T6 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C3_T6_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T6_vs_T0/Immunome/reference/BP_List_genes_H2D_C3_T6_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C3_T6_vs_T0, by = "Gene") 

gene_counts_H2D_C3_T6_vs_T0_immunome    <- BP_List_genes_H2D_C3_T6_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T6_vs_T0_immunome  <- BP_List_genes_H2D_C3_T6_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C3_T6_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C3_T6_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C3_T6_vs_T0_immunome %>% 
                                           filter(Gene_significant == "Yes") %>% 
                                           summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                           mutate(Family = "H2D")

H2D_C3_T6_vs_T0_sig_DEGS_result        <- merge(H2D_C3_T6_vs_T0_sig_DEGS_result,H2D_C3_T6_vs_T0_immunome_DEGS)


## Significant GO terms immunome H2D_C3_T6 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C3_T6_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T6_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C3_T6_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C3_T6_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C3_T12 vs H2D_C3_T0 GO_MWU  ----
H2D_C3_T12_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C3_T12_vs_T0_FC    <- left_join(List_genes, H2D_C3_T12_vs_T0_FC)
write(write.table(List_genes_H2D_C3_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T12_vs_T0/List_genes_H2D_C3_T12_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C3_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T12_vs_T0/Immunome/reference/List_genes_H2D_C3_T12_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T12_vs_T0/Immunome/reference/")
input="List_genes_H2D_C3_T12_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C3_T12 vs H2D_C3_T0 GO_MWU
H2D_C3_T12_vs_T0          <- results(dds, contrast=c("group", "H2D_C3_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C3_T12_vs_T0_sig_DEGS <- H2D_C3_T12_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C3_T12_vs_T0_sig_DEGS)

H2D_C3_T12_vs_T0_sig_DEGS_result <- H2D_C3_T12_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                   Down  = sum(log2FoldChange < 0),
                                   Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                   rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                   mutate(Family = "H2D", Age = "4", Time = "12") %>%
                                   dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C3_T12 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C3_T12_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T12_vs_T0/Immunome/reference/BP_List_genes_H2D_C3_T12_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C3_T12_vs_T0, by = "Gene") 

gene_counts_H2D_C3_T12_vs_T0_immunome    <- BP_List_genes_H2D_C3_T12_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T12_vs_T0_immunome  <- BP_List_genes_H2D_C3_T12_vs_T0_immunome %>% 
                                            left_join(gene_counts_H2D_C3_T12_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C3_T12_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C3_T12_vs_T0_immunome %>% 
                                           filter(Gene_significant == "Yes") %>% 
                                           summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                           mutate(Family = "H2D")

H2D_C3_T12_vs_T0_sig_DEGS_result        <- merge(H2D_C3_T12_vs_T0_sig_DEGS_result,H2D_C3_T12_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C3_T12 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C3_T12_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T12_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_H2D_C3_T12_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "H2D_C3_T12_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C3_T24 vs H2D_C3_T0 GO_MWU  ----
H2D_C3_T24_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C3_T24_vs_T0_FC    <- left_join(List_genes, H2D_C3_T24_vs_T0_FC)
write(write.table(List_genes_H2D_C3_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T24_vs_T0/List_genes_H2D_C3_T24_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C3_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T24_vs_T0/Immunome/reference/List_genes_H2D_C3_T24_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T24_vs_T0/Immunome/reference/")
input="List_genes_H2D_C3_T24_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C3_T24 vs H2D_C3_T0 GO_MWU
H2D_C3_T24_vs_T0          <- results(dds, contrast=c("group", "H2D_C3_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C3_T24_vs_T0_sig_DEGS <- H2D_C3_T24_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C3_T24_vs_T0_sig_DEGS)

H2D_C3_T24_vs_T0_sig_DEGS_result <- H2D_C3_T24_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                   Down  = sum(log2FoldChange < 0),
                                   Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                   rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                   mutate(Family = "H2D", Age = "4", Time = "24") %>%
                                   dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C3_T24 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C3_T24_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T24_vs_T0/Immunome/reference/BP_List_genes_H2D_C3_T24_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C3_T24_vs_T0, by = "Gene") 

gene_counts_H2D_C3_T24_vs_T0_immunome    <- BP_List_genes_H2D_C3_T24_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T24_vs_T0_immunome  <- BP_List_genes_H2D_C3_T24_vs_T0_immunome %>% 
                                            left_join(gene_counts_H2D_C3_T24_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C3_T24_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C3_T24_vs_T0_immunome %>% 
                                            filter(Gene_significant == "Yes") %>% 
                                            summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                            mutate(Family = "H2D")

H2D_C3_T24_vs_T0_sig_DEGS_result        <- merge(H2D_C3_T24_vs_T0_sig_DEGS_result,H2D_C3_T24_vs_T0_immunome_DEGS)


## Significant GO terms immunome H2D_C3_T24 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C3_T24_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T24_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_H2D_C3_T24_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "H2D_C3_T24_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


# C7 ----


### H2D_C7_T0 vs H2D_C3_T0 GO_MWU  ----
H2D_C7_T0_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T0", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C7_T0_vs_T0_FC    <- left_join(List_genes, H2D_C7_T0_vs_T0_FC)
write(write.table(List_genes_H2D_C7_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T0_vs_T0/List_genes_H2D_C7_T0_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C7_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T0_vs_T0/Immunome/reference/List_genes_H2D_C7_T0_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T0_vs_T0/Immunome/reference/")
input="List_genes_H2D_C7_T0_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C7_T0 vs H2D_C3_T0 GO_MWU
H2D_C7_T0_vs_T0          <- results(dds, contrast=c("group", "H2D_C7_T0", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C7_T0_vs_T0_sig_DEGS <- H2D_C7_T0_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C7_T0_vs_T0_sig_DEGS)

H2D_C7_T0_vs_T0_sig_DEGS_result <- H2D_C7_T0_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                    Down  = sum(log2FoldChange < 0),
                                    Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                    mutate(Family = "H2D", Age = "16", Time = "0") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C7_T0 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C7_T0_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T0_vs_T0/Immunome/reference/BP_List_genes_H2D_C7_T0_vs_T0.txt") %>% 
  dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C7_T0_vs_T0, by = "Gene") 

gene_counts_H2D_C7_T0_vs_T0_immunome    <- BP_List_genes_H2D_C7_T0_vs_T0_immunome %>%
  dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T0_vs_T0_immunome  <- BP_List_genes_H2D_C7_T0_vs_T0_immunome %>% 
  left_join(gene_counts_H2D_C7_T0_vs_T0_immunome, by = "term") %>% 
  mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C7_T0_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C7_T0_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C7_T0_vs_T0_sig_DEGS_result        <- merge(H2D_C7_T0_vs_T0_sig_DEGS_result, H2D_C7_T0_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C7_T0 vs H2D_C7_T0 GO_MWU
GO_terms_H2D_C7_T0_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T0_vs_T0.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_H2D_C7_T0_vs_T0_immunome, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C7_T0_vs_T0") %>% 
  left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C7_T3 vs H2D_C3_T0 GO_MWU  ----
H2D_C7_T3_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C7_T3_vs_T0_FC    <- left_join(List_genes, H2D_C7_T3_vs_T0_FC)
write(write.table(List_genes_H2D_C7_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T3_vs_T0/List_genes_H2D_C7_T3_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C7_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T3_vs_T0/Immunome/reference/List_genes_H2D_C7_T3_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T3_vs_T0/Immunome/reference/")
input="List_genes_H2D_C7_T3_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C7_T3 vs H2D_C7_T0 GO_MWU
H2D_C7_T3_vs_T0          <- results(dds, contrast=c("group", "H2D_C7_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C7_T3_vs_T0_sig_DEGS <- H2D_C7_T3_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C7_T3_vs_T0_sig_DEGS)

H2D_C7_T3_vs_T0_sig_DEGS_result <- H2D_C7_T3_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                   Down  = sum(log2FoldChange < 0),
                                   Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                   rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                   mutate(Family = "H2D", Age = "16", Time = "3") %>%
                                   dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C7_T3 vs H2D_C7_T0 GO_MWU
BP_List_genes_H2D_C7_T3_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T3_vs_T0/Immunome/reference/BP_List_genes_H2D_C7_T3_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C7_T3_vs_T0, by = "Gene") 

gene_counts_H2D_C7_T3_vs_T0_immunome    <- BP_List_genes_H2D_C7_T3_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T3_vs_T0_immunome  <- BP_List_genes_H2D_C7_T3_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C7_T3_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C7_T3_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C7_T3_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C7_T3_vs_T0_sig_DEGS_result        <- merge(H2D_C7_T3_vs_T0_sig_DEGS_result, H2D_C7_T3_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C7_T3 vs H2D_C7_T0 GO_MWU
GO_terms_H2D_C7_T3_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T3_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C7_T3_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C7_T3_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



### H2D_C7_T6 vs H2D_C3_T0 GO_MWU  ----
H2D_C7_T6_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C7_T6_vs_T0_FC    <- left_join(List_genes, H2D_C7_T6_vs_T0_FC)
write(write.table(List_genes_H2D_C7_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T6_vs_T0/List_genes_H2D_C7_T6_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C7_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T6_vs_T0/Immunome/reference/List_genes_H2D_C7_T6_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T6_vs_T0/Immunome/reference/")
input="List_genes_H2D_C7_T6_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C7_T6 vs H2D_C3_T0 GO_MWU
H2D_C7_T6_vs_T0          <- results(dds, contrast=c("group", "H2D_C7_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C7_T6_vs_T0_sig_DEGS <- H2D_C7_T6_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C7_T6_vs_T0_sig_DEGS)

H2D_C7_T6_vs_T0_sig_DEGS_result <- H2D_C7_T6_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                   Down  = sum(log2FoldChange < 0),
                                   Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                   rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                   mutate(Family = "H2D", Age = "16", Time = "6") %>%
                                   dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C7_T6 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C7_T6_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T6_vs_T0/Immunome/reference/BP_List_genes_H2D_C7_T6_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C7_T6_vs_T0, by = "Gene") 

gene_counts_H2D_C7_T6_vs_T0_immunome    <- BP_List_genes_H2D_C7_T6_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T6_vs_T0_immunome  <- BP_List_genes_H2D_C7_T6_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C7_T6_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C7_T6_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C7_T6_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C7_T6_vs_T0_sig_DEGS_result        <- merge(H2D_C7_T6_vs_T0_sig_DEGS_result, H2D_C7_T6_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C7_T6 vs H2D_C7_T0 GO_MWU
GO_terms_H2D_C7_T6_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T6_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C7_T6_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C7_T6_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C7_T12 vs H2D_C3_T0 GO_MWU  ----
H2D_C7_T12_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C7_T12_vs_T0_FC    <- left_join(List_genes, H2D_C7_T12_vs_T0_FC)
write(write.table(List_genes_H2D_C7_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T12_vs_T0/List_genes_H2D_C7_T12_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C7_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T12_vs_T0/Immunome/reference/List_genes_H2D_C7_T12_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T12_vs_T0/Immunome/reference/")
input="List_genes_H2D_C7_T12_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C7_T12 vs H2D_C3_T0 GO_MWU
H2D_C7_T12_vs_T0          <- results(dds, contrast=c("group", "H2D_C7_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C7_T12_vs_T0_sig_DEGS <- H2D_C7_T12_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C7_T12_vs_T0_sig_DEGS)

H2D_C7_T12_vs_T0_sig_DEGS_result <- H2D_C7_T12_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                    Down  = sum(log2FoldChange < 0),
                                    Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                    mutate(Family = "H2D", Age = "16", Time = "12") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C7_T12 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C7_T12_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T12_vs_T0/Immunome/reference/BP_List_genes_H2D_C7_T12_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C7_T12_vs_T0, by = "Gene") 

gene_counts_H2D_C7_T12_vs_T0_immunome    <- BP_List_genes_H2D_C7_T12_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T12_vs_T0_immunome  <- BP_List_genes_H2D_C7_T12_vs_T0_immunome %>% 
                                            left_join(gene_counts_H2D_C7_T12_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C7_T12_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C7_T12_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C7_T12_vs_T0_sig_DEGS_result        <- merge(H2D_C7_T12_vs_T0_sig_DEGS_result, H2D_C7_T12_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C7_T12 vs H2D_C7_T0 GO_MWU
GO_terms_H2D_C7_T12_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T12_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_H2D_C7_T12_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "H2D_C7_T12_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C7_T24 vs H2D_C3_T0 GO_MWU  ----
H2D_C7_T24_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C7_T24_vs_T0_FC    <- left_join(List_genes, H2D_C7_T24_vs_T0_FC)
write(write.table(List_genes_H2D_C7_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T24_vs_T0/List_genes_H2D_C7_T24_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C7_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T24_vs_T0/Immunome/reference/List_genes_H2D_C7_T24_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T24_vs_T0/Immunome/reference/")
input="List_genes_H2D_C7_T24_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C7_T24 vs H2D_C3_T0 GO_MWU
H2D_C7_T24_vs_T0          <- results(dds, contrast=c("group", "H2D_C7_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C7_T24_vs_T0_sig_DEGS <- H2D_C7_T24_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C7_T24_vs_T0_sig_DEGS)

H2D_C7_T24_vs_T0_sig_DEGS_result <- H2D_C7_T24_vs_T0_sig_DEGS %>%
                                    summarise(Up    = sum(log2FoldChange > 0),
                                              Down  = sum(log2FoldChange < 0),
                                              Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                    mutate(Family = "H2D", Age = "16", Time = "24") %>%
                                    dplyr::select(3,4,5,1,2)


## GO terms and genes involved immunome H2D_C7_T24 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C7_T24_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T24_vs_T0/Immunome/reference/BP_List_genes_H2D_C7_T24_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C7_T24_vs_T0, by = "Gene") 

gene_counts_H2D_C7_T24_vs_T0_immunome    <- BP_List_genes_H2D_C7_T24_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T24_vs_T0_immunome  <- BP_List_genes_H2D_C7_T24_vs_T0_immunome %>% 
                                            left_join(gene_counts_H2D_C7_T24_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C7_T24_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C7_T24_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C7_T24_vs_T0_sig_DEGS_result        <- merge(H2D_C7_T24_vs_T0_sig_DEGS_result, H2D_C7_T24_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C7_T24 vs H2D_C7_T0 GO_MWU
GO_terms_H2D_C7_T24_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T24_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_H2D_C7_T24_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "H2D_C7_T24_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


# C8  ----


### H2D_C8_T0 vs H2D_C3_T0 GO_MWU  ----
H2D_C8_T0_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T0", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C8_T0_vs_T0_FC    <- left_join(List_genes, H2D_C8_T0_vs_T0_FC)
write(write.table(List_genes_H2D_C8_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T0_vs_T0/List_genes_H2D_C8_T0_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C8_T0_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T0_vs_T0/Immunome/reference/List_genes_H2D_C8_T0_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T0_vs_T0/Immunome/reference/")
input="List_genes_H2D_C8_T0_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C8_T0 vs H2D_C3_T0 GO_MWU
H2D_C8_T0_vs_T0          <- results(dds, contrast=c("group", "H2D_C8_T0", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C8_T0_vs_T0_sig_DEGS <- H2D_C8_T0_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C8_T0_vs_T0_sig_DEGS)

H2D_C8_T0_vs_T0_sig_DEGS_result <- H2D_C8_T0_vs_T0_sig_DEGS %>%
                                   summarise(Up    = sum(log2FoldChange > 0),
                                             Down  = sum(log2FoldChange < 0),
                                             Total = Up +Down) %>% t() %>% as.data.frame() %>% 
                                    rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                    mutate(Family = "H2D", Age = "28", Time = "0") %>%
                                    dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C8_T0 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C8_T0_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T0_vs_T0/Immunome/reference/BP_List_genes_H2D_C8_T0_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C8_T0_vs_T0, by = "Gene") 

gene_counts_H2D_C8_T0_vs_T0_immunome    <- BP_List_genes_H2D_C8_T0_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T0_vs_T0_immunome  <- BP_List_genes_H2D_C8_T0_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C8_T0_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C8_T0_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C8_T0_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C8_T0_vs_T0_sig_DEGS_result        <- merge(H2D_C8_T0_vs_T0_sig_DEGS_result, H2D_C8_T0_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C8_T0 vs H2D_C8_T0 GO_MWU
GO_terms_H2D_C8_T0_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T0_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C8_T0_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C8_T0_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C8_T3 vs H2D_C3_T0 GO_MWU  ----
H2D_C8_T3_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C8_T3_vs_T0_FC    <- left_join(List_genes, H2D_C8_T3_vs_T0_FC)
write(write.table(List_genes_H2D_C8_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T3_vs_T0/List_genes_H2D_C8_T3_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C8_T3_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T3_vs_T0/Immunome/reference/List_genes_H2D_C8_T3_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T3_vs_T0/Immunome/reference/")
input="List_genes_H2D_C8_T3_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C8_T3 vs H2D_C8_T0 GO_MWU
H2D_C8_T3_vs_T0          <- results(dds, contrast=c("group", "H2D_C8_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C8_T3_vs_T0_sig_DEGS <- H2D_C8_T3_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C8_T3_vs_T0_sig_DEGS)

H2D_C8_T3_vs_T0_sig_DEGS_result <- H2D_C8_T3_vs_T0_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up +Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Family = "H2D", Age = "28", Time = "3") %>%
  dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C8_T3 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C8_T3_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T3_vs_T0/Immunome/reference/BP_List_genes_H2D_C8_T3_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C8_T3_vs_T0, by = "Gene") 

gene_counts_H2D_C8_T3_vs_T0_immunome    <- BP_List_genes_H2D_C8_T3_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T3_vs_T0_immunome  <- BP_List_genes_H2D_C8_T3_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C8_T3_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C8_T3_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C8_T3_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C8_T3_vs_T0_sig_DEGS_result        <- merge(H2D_C8_T3_vs_T0_sig_DEGS_result, H2D_C8_T3_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C8_T3 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C8_T3_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T3_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C8_T3_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C8_T3_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



### H2D_C8_T6 vs H2D_C3_T0 GO_MWU  ----
H2D_C8_T6_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                         filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C8_T6_vs_T0_FC    <- left_join(List_genes, H2D_C8_T6_vs_T0_FC)
write(write.table(List_genes_H2D_C8_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T6_vs_T0/List_genes_H2D_C8_T6_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C8_T6_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T6_vs_T0/Immunome/reference/List_genes_H2D_C8_T6_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T6_vs_T0/Immunome/reference/")
input="List_genes_H2D_C8_T6_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C8_T6 vs H2D_C3_T0 GO_MWU
H2D_C8_T6_vs_T0          <- results(dds, contrast=c("group", "H2D_C8_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C8_T6_vs_T0_sig_DEGS <- H2D_C8_T6_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C8_T6_vs_T0_sig_DEGS)

H2D_C8_T6_vs_T0_sig_DEGS_result <- H2D_C8_T6_vs_T0_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up +Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Family = "H2D", Age = "28", Time = "6") %>%
  dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C8_T6 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C8_T6_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T6_vs_T0/Immunome/reference/BP_List_genes_H2D_C8_T6_vs_T0.txt") %>% 
                                           dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C8_T6_vs_T0, by = "Gene") 

gene_counts_H2D_C8_T6_vs_T0_immunome    <- BP_List_genes_H2D_C8_T6_vs_T0_immunome %>%
                                           dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                           nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T6_vs_T0_immunome  <- BP_List_genes_H2D_C8_T6_vs_T0_immunome %>% 
                                           left_join(gene_counts_H2D_C8_T6_vs_T0_immunome, by = "term") %>% 
                                           mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C8_T6_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C8_T6_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C8_T6_vs_T0_sig_DEGS_result        <- merge(H2D_C8_T6_vs_T0_sig_DEGS_result, H2D_C8_T6_vs_T0_immunome_DEGS)


## Significant GO terms immunome H2D_C8_T6 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C8_T6_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T6_vs_T0.txt", sep="") %>% 
                                           mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                           mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                           left_join(BP_List_genes_H2D_C8_T6_vs_T0_immunome, by = c("term", "name")) %>% 
                                           mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                           mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                           mutate(Comparison = "H2D_C8_T6_vs_T0") %>% 
                                           left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C8_T12 vs H2D_C3_T0 GO_MWU  ----
H2D_C8_T12_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C8_T12_vs_T0_FC    <- left_join(List_genes, H2D_C8_T12_vs_T0_FC)
write(write.table(List_genes_H2D_C8_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T12_vs_T0/List_genes_H2D_C8_T12_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C8_T12_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T12_vs_T0/Immunome/reference/List_genes_H2D_C8_T12_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T12_vs_T0/Immunome/reference/")
input="List_genes_H2D_C8_T12_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C8_T12 vs H2D_C3_T0 GO_MWU
H2D_C8_T12_vs_T0          <- results(dds, contrast=c("group", "H2D_C8_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C8_T12_vs_T0_sig_DEGS <- H2D_C8_T12_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C8_T12_vs_T0_sig_DEGS)

H2D_C8_T12_vs_T0_sig_DEGS_result <- H2D_C8_T12_vs_T0_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up +Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Family = "H2D", Age = "28", Time = "12") %>%
  dplyr::select(3,4,5,1,2)


## GO terms and genes involved immunome H2D_C8_T12 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C8_T12_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T12_vs_T0/Immunome/reference/BP_List_genes_H2D_C8_T12_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C8_T12_vs_T0, by = "Gene") 

gene_counts_H2D_C8_T12_vs_T0_immunome    <- BP_List_genes_H2D_C8_T12_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T12_vs_T0_immunome  <- BP_List_genes_H2D_C8_T12_vs_T0_immunome %>% 
                                            left_join(gene_counts_H2D_C8_T12_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C8_T12_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C8_T12_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C8_T12_vs_T0_sig_DEGS_result        <- merge(H2D_C8_T12_vs_T0_sig_DEGS_result, H2D_C8_T12_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C8_T12 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C8_T12_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T12_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_H2D_C8_T12_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "H2D_C8_T12_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


### H2D_C8_T24 vs H2D_C3_T0 GO_MWU  ----
H2D_C8_T24_vs_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                          as.data.frame()      %>% rownames_to_column(var = "gene")        %>% 
                          filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_H2D_C8_T24_vs_T0_FC    <- left_join(List_genes, H2D_C8_T24_vs_T0_FC)
write(write.table(List_genes_H2D_C8_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T24_vs_T0/List_genes_H2D_C8_T24_vs_T0.txt",          row.names = F, quote = F, col.names=F, sep = ","))
write(write.table(List_genes_H2D_C8_T24_vs_T0_FC, "D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T24_vs_T0/Immunome/reference/List_genes_H2D_C8_T24_vs_T0.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome
setwd("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T24_vs_T0/Immunome/reference/")
input="List_genes_H2D_C8_T24_vs_T0.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

## Genes list to use H2D_C8_T24 vs H2D_C3_T0 GO_MWU
H2D_C8_T24_vs_T0          <- results(dds, contrast=c("group", "H2D_C8_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
H2D_C8_T24_vs_T0_sig_DEGS <- H2D_C8_T24_vs_T0 %>% dplyr::filter(padj <= 0.05) 
dim(H2D_C8_T24_vs_T0_sig_DEGS)

H2D_C8_T24_vs_T0_sig_DEGS_result <- H2D_C8_T24_vs_T0_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up +Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Family = "H2D", Age = "28", Time = "24") %>%
  dplyr::select(3,4,5,1,2)

## GO terms and genes involved immunome H2D_C8_T24 vs H2D_C3_T0 GO_MWU
BP_List_genes_H2D_C8_T24_vs_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T24_vs_T0/Immunome/reference/BP_List_genes_H2D_C8_T24_vs_T0.txt") %>% 
                                            dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(H2D_C8_T24_vs_T0, by = "Gene") 

gene_counts_H2D_C8_T24_vs_T0_immunome    <- BP_List_genes_H2D_C8_T24_vs_T0_immunome %>%
                                            dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                            nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T24_vs_T0_immunome  <- BP_List_genes_H2D_C8_T24_vs_T0_immunome %>% 
                                            left_join(gene_counts_H2D_C8_T24_vs_T0_immunome, by = "term") %>% 
                                            mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

H2D_C8_T24_vs_T0_immunome_DEGS           <- BP_List_genes_H2D_C8_T24_vs_T0_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Family = "H2D")

H2D_C8_T24_vs_T0_sig_DEGS_result        <- merge(H2D_C8_T24_vs_T0_sig_DEGS_result, H2D_C8_T24_vs_T0_immunome_DEGS)

## Significant GO terms immunome H2D_C8_T24 vs H2D_C3_T0 GO_MWU
GO_terms_H2D_C8_T24_vs_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T24_vs_T0.txt", sep="") %>% 
                                            mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                            mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                            left_join(BP_List_genes_H2D_C8_T24_vs_T0_immunome, by = c("term", "name")) %>% 
                                            mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                            mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                            mutate(Comparison = "H2D_C8_T24_vs_T0") %>% 
                                            left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 

# Heatmap_H2D -----

GO_terms_H2D_C3_T3_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T3_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C3_T6_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T6_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C3_T12_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T12_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_H2D_C3_T24_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C3_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C3_T24_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_H2D_C7_T0_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T0_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C7_T3_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T3_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C7_T6_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T6_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C7_T12_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T12_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_H2D_C7_T24_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C7_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C7_T24_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_H2D_C8_T0_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T0_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T0_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C8_T3_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T3_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T3_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C8_T6_vs_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T6_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T6_vs_T0.txt", sep="")   %>% filter(p.adj <= 0.05)
GO_terms_H2D_C8_T12_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T12_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T12_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)
GO_terms_H2D_C8_T24_vs_T0_immunome_sig <- read.csv("D:/Decicomp/R/GO_terms/RNA/Infection/GO_RNA_H2D_C8_T24_vs_T0/Immunome/reference/MWU_BP_List_genes_H2D_C8_T24_vs_T0.txt", sep="") %>% filter(p.adj <= 0.05)

dim(GO_terms_H2D_C3_T3_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C3_T6_vs_T0_immunome_sig)
dim(GO_terms_H2D_C3_T12_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C3_T24_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C7_T3_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C7_T6_vs_T0_immunome_sig)
dim(GO_terms_H2D_C7_T12_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C7_T24_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C8_T3_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C8_T6_vs_T0_immunome_sig)
dim(GO_terms_H2D_C8_T12_vs_T0_immunome_sig) 
dim(GO_terms_H2D_C8_T24_vs_T0_immunome_sig) 

# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


GO_H2D_infection_reference <- rbind(GO_terms_H2D_C3_T3_vs_T0_immunome, GO_terms_H2D_C3_T6_vs_T0_immunome, GO_terms_H2D_C3_T12_vs_T0_immunome, GO_terms_H2D_C3_T24_vs_T0_immunome,
                              GO_terms_H2D_C7_T0_vs_T0_immunome, GO_terms_H2D_C7_T3_vs_T0_immunome, GO_terms_H2D_C7_T6_vs_T0_immunome, GO_terms_H2D_C7_T12_vs_T0_immunome, GO_terms_H2D_C7_T24_vs_T0_immunome,
                              GO_terms_H2D_C8_T0_vs_T0_immunome, GO_terms_H2D_C8_T3_vs_T0_immunome, GO_terms_H2D_C8_T6_vs_T0_immunome, GO_terms_H2D_C8_T12_vs_T0_immunome, GO_terms_H2D_C8_T24_vs_T0_immunome) %>%
                              filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% distinct(name, .keep_all = TRUE) %>%
                              dplyr::select(Comparison, name, Enrichment) %>% 
                              pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
                              as.data.frame() %>% mutate_all(~ replace(., is.na(.), 0)) %>%
                              mutate(H2D_C3_T0_vs_T0 = 0) %>%
                              column_to_rownames(var = "name")  %>%
                              dplyr::select(15,1:14)


# Colum Colors heatmap   
annotation_col <- data.frame(
  Time = c("T0","T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24"),
  Age  = c(rep("4 months", 5), rep("16 months", 5), rep("28 months", 5)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"),
  Time   = c(T0 = "white",  T3           = "#CD5C5C", T6 = "#B22222", T12="darkred",  T24 ="red"))

rownames(annotation_col) <- colnames(GO_H2D_infection_reference) 


# Gradient enrichment                        
myBreaks <- c(seq(min(GO_H2D_infection_reference), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_H2D_infection_reference)/100, max(GO_H2D_infection_reference), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

GO_H2D_infection_reference_heatmap <- pheatmap(GO_H2D_infection_reference, 
                                      cluster_cols = F, 
                                      scale = "none",
                                      cluster_rows = T, 
                                      fontsize_row = 14, 
                                      color = mycolor, 
                                      breaks = myBreaks,
                                      border_color = "black",
                                      clustering_distance_rows = "euclidean",
                                      show_colnames = F, 
                                      show_rownames = T,
                                      annotation_col = annotation_col, 
                                      #annotation_row = annotation_row,
                                      annotation_colors = annotation_colors,
                                      #cutree_rows = 4,
                                      gaps_col =  c(5,10,15))
GO_H2D_infection_reference_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/GO_H2D_infection_reference_heatmap.png", GO_H2D_infection_reference_heatmap, width = 13, height = 8, dpi = 300)



# Filtering at T0
GO_H2D_infection_reference_T0 <- rbind(GO_terms_H2D_C7_T0_vs_T0_immunome, GO_terms_H2D_C8_T0_vs_T0_immunome) %>%
  filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% distinct(name, .keep_all = TRUE) %>%
  dplyr::select(Comparison, name, Enrichment) %>% 
  pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
  as.data.frame() %>% mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(H2D_C3_T0_vs_T0 = 0) %>%
  column_to_rownames(var = "name") %>%
  dplyr::select(3,1,2)


# Colum Colors heatmap   
annotation_col <- data.frame(
  #Time = c("T0", "T0", "T0"),
  Age  = c(rep("4 months", 1), rep("16 months", 1), rep("28 months", 1)))

annotation_colors  <-  list( #Family = c(H2D         = "#2f5597"),
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"))
 # Time   = c(T0 = "white"))

rownames(annotation_col) <- colnames(GO_H2D_infection_reference_T0) 


# Gradient enrichment                        
myBreaks <- c(seq(min(GO_H2D_infection_reference_T0), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_H2D_infection_reference_T0)/100, max(GO_H2D_infection_reference_T0), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

GO_H2D_infection_reference_T0_heatmap <- pheatmap(GO_H2D_infection_reference_T0, 
                                                  cluster_cols = F, 
                                                  scale = "none",
                                                  cluster_rows = T, 
                                                  fontsize_row = 14, 
                                                  color = mycolor, 
                                                  breaks = myBreaks,
                                                  border_color = T,
                                                  clustering_distance_rows = "euclidean",
                                                  show_colnames = F, 
                                                  show_rownames = T,
                                                  annotation_col = annotation_col, 
                                                  #annotation_row = annotation_row,
                                                  annotation_names_col = FALSE,
                                                  annotation_legend = FALSE,
                                                  legend = F,
                                                  legend_title = NULL,
                                                  annotation_colors = annotation_colors)
#cutree_rows = 4,
#gaps_col =  c(4,9,14))
GO_H2D_infection_reference_T0_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/GO_H2D_infection_reference_T0_heatmap.tiff", GO_H2D_infection_reference_T0_heatmap, width = 9, height = 8, dpi = 300)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


# Both families heatmap  wierd plots every time changes ----

GO_F14R_infection_frame <- GO_F14R_infection_reference  %>% rownames_to_column(var = "name")  # dim(GO_F14R_infection_frame)
GO_H2D_infection_frame  <- GO_H2D_infection_reference   %>% rownames_to_column(var = "name")  # dim(GO_H2D_infection_frame)

GO_F14R_H2D_infection_reference  <- merge(GO_F14R_infection_frame, GO_H2D_infection_frame, by = "name", all = TRUE) %>% 
                                    mutate_all(~replace(., is.na(.), 0)) %>%
                                    column_to_rownames(var = "name") 

# Define annotation for the columns
annotation_col <- data.frame( 
  Time = rep(c("T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24"), times = 2),
  Age  = c(rep("4 months_F14R", 5), rep("16 months_F14R", 5), rep("28 months_F14R", 5),
           rep("4 months_H2D", 5), rep("16 months_H2D", 5), rep("28 months_H2D", 5)))
# Define colors for annotations
annotation_colors <- list(
  Age = c(
    "4 months_F14R" = "#dae3f3", "16 months_F14R" = "#8faadc", "28 months_F14R" = "#2f5597",
    "4 months_H2D"  = "#E2F0D9", "16 months_H2D"  = "#A9D18E", "28 months_H2D"  = "#385700"),
  Time = c(T0 ="white", T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred", T24 = "red"))

# Assign rownames to annotation_col based on the colnames of the data
rownames(annotation_col) <- colnames(GO_F14R_H2D_infection_reference)



# Gradient enrichment                        
myBreaks <- c(seq(min(GO_F14R_H2D_infection_reference), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_F14R_H2D_infection_reference)/100, max(GO_F14R_H2D_infection_reference), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

GO_F14R_H2D_infection_reference_heatmap <- pheatmap(GO_F14R_H2D_infection_reference, 
                                          cluster_cols = F, 
                                          scale = "none",
                                          cluster_rows = T, 
                                          fontsize_row = 14, 
                                          color = mycolor, 
                                          breaks = myBreaks,
                                          border_color = "black",
                                          clustering_distance_rows = "euclidean",
                                          show_colnames = F, 
                                          show_rownames = T,
                                          annotation_col = annotation_col, 
                                          #annotation_row = annotation_row,
                                          annotation_colors = annotation_colors,
                                          cutree_rows = 1,
                                          gaps_col =  c(5,10,15,20,25))
GO_F14R_H2D_infection_reference_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/GO_F14R_H2D_infection_reference_heatmap.png", GO_F14R_H2D_infection_reference_heatmap, width = 15, height = 8, dpi = 300)




# DEGS ----

F14R_C3_T0_vs_T0_sig_DEGS_result <- data.frame(
  Family = c("F14R", "F14R", "F14R"),
  Age    = c(4, 4, 4),
  Time   = c(0, 0, 0),
  DEGS   = c("Up", "Down", "Total"),
  Number = c(0, 0, 0),
  DEGs_immuno = c(0, 0, 0))

# DEGS F14R ----
F14R_DEGS_context_reference <- rbind( 
  F14R_C3_T0_vs_T0_sig_DEGS_result,
  F14R_C3_T3_vs_T0_sig_DEGS_result,
  F14R_C3_T6_vs_T0_sig_DEGS_result,
  F14R_C3_T12_vs_T0_sig_DEGS_result,
  F14R_C3_T24_vs_T0_sig_DEGS_result,
  F14R_C7_T0_vs_T0_sig_DEGS_result,
  F14R_C7_T3_vs_T0_sig_DEGS_result,
  F14R_C7_T6_vs_T0_sig_DEGS_result,
  F14R_C7_T12_vs_T0_sig_DEGS_result,
  F14R_C7_T24_vs_T0_sig_DEGS_result,
  F14R_C8_T0_vs_T0_sig_DEGS_result,
  F14R_C8_T3_vs_T0_sig_DEGS_result,
  F14R_C8_T6_vs_T0_sig_DEGS_result,
  F14R_C8_T12_vs_T0_sig_DEGS_result,
  F14R_C8_T24_vs_T0_sig_DEGS_result) %>% 
  mutate(ratio = DEGs_immuno / Number) %>% 
  filter(DEGS != "Total")

F14R_DEGS_context_reference$Age  <- factor(F14R_DEGS_context_reference$Age, levels = c(4, 16, 28))
F14R_DEGS_context_reference$Time <- factor(F14R_DEGS_context_reference$Time, levels = c(0, 3, 6, 12, 24))
F14R_DEGS_context_reference$DEGS <- factor(F14R_DEGS_context_reference$DEGS, levels = c("Up", "Down"))


F14R_DEGS_context_reference_plot  <- ggplot(F14R_DEGS_context_reference, aes(x = Time, y = Number, fill = DEGS)) +
  geom_hline(yintercept = c(2500, 5000, 7500), linetype = "dashed", color = "red", linewidth=1) +  
   geom_bar(stat = "identity", position = "stack", color="black") +
  geom_text(aes(label = DEGs_immuno), position = position_stack(vjust = 1), size = 5, color = "black", vjust = -0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10000)) +
  facet_wrap(~ Age, strip.position = "bottom", scales = "free_x") +
  labs(x = " ", y = "Number of DEGs") +
  scale_fill_manual(values = c('Up' = scales::alpha('yellow', 0.7), 'Down' = scales::alpha('blue', 0.7))) +
  Style_format_theme +
  theme(
    strip.placement = "outside",  # Place facet labels outside the panel
    strip.background = element_blank(),  # Remove the background color of facet labels
    strip.text.x = element_text(size = 32),  # Adjust size of facet labels
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
    panel.spacing.x = unit(1, "lines"))  # Rotate x-axis labels for readability

F14R_DEGS_context_reference_plot
# ggsave("D:/Decicomp/R/MOFA_omics/F14R_DEGS_context_reference_plot.png", F14R_DEGS_context_reference_plot, width = 12, height = 8, dpi = 600)





# DEGS ----

H2D_C3_T0_vs_T0_sig_DEGS_result <- data.frame(
  Family = c("H2D", "H2D", "H2D"),
  Age    = c(4, 4, 4),
  Time   = c(0, 0, 0),
  DEGS   = c("Up", "Down", "Total"),
  Number = c(0, 0, 0),
  DEGs_immuno = c(0, 0, 0))

# DEGS H2D ----
H2D_DEGS_context_reference <- rbind( 
  H2D_C3_T0_vs_T0_sig_DEGS_result,
  H2D_C3_T3_vs_T0_sig_DEGS_result,
  H2D_C3_T6_vs_T0_sig_DEGS_result,
  H2D_C3_T12_vs_T0_sig_DEGS_result,
  H2D_C3_T24_vs_T0_sig_DEGS_result,
  H2D_C7_T0_vs_T0_sig_DEGS_result,
  H2D_C7_T3_vs_T0_sig_DEGS_result,
  H2D_C7_T6_vs_T0_sig_DEGS_result,
  H2D_C7_T12_vs_T0_sig_DEGS_result,
  H2D_C7_T24_vs_T0_sig_DEGS_result,
  H2D_C8_T0_vs_T0_sig_DEGS_result,
  H2D_C8_T3_vs_T0_sig_DEGS_result,
  H2D_C8_T6_vs_T0_sig_DEGS_result,
  H2D_C8_T12_vs_T0_sig_DEGS_result,
  H2D_C8_T24_vs_T0_sig_DEGS_result) %>% 
  mutate(ratio = DEGs_immuno / Number) %>% 
  filter(DEGS != "Total")

H2D_DEGS_context_reference$Age  <- factor(H2D_DEGS_context_reference$Age, levels = c(4, 16, 28))
H2D_DEGS_context_reference$Time <- factor(H2D_DEGS_context_reference$Time, levels = c(0, 3, 6, 12, 24))
H2D_DEGS_context_reference$DEGS <- factor(H2D_DEGS_context_reference$DEGS, levels = c("Up", "Down"))


H2D_DEGS_context_reference_plot  <- ggplot(H2D_DEGS_context_reference, aes(x = Time, y = Number, fill = DEGS)) +
  geom_hline(yintercept = c(2500, 5000, 7500), linetype = "dashed", color = "red", linewidth=1) +
  geom_bar(stat = "identity", position = "stack", color="black") +
  geom_text(aes(label = DEGs_immuno), position = position_stack(vjust = 1), size = 5, color = "black", vjust = -0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10000)) +
  facet_wrap(~ Age, strip.position = "bottom", scales = "free_x") +
  labs(x = " ", y = "Number of DEGs") +
  scale_fill_manual(values = c('Up' = scales::alpha('yellow', 0.7), 'Down' = scales::alpha('blue', 0.7))) +
  Style_format_theme +
  theme(
    strip.placement = "outside",  # Place facet labels outside the panel
    strip.background = element_blank(),  # Remove the background color of facet labels
    strip.text.x = element_text(size = 32),  # Adjust size of facet labels
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
    panel.spacing.x = unit(1, "lines"))  # Rotate x-axis labels for readability

H2D_DEGS_context_reference_plot
# ggsave("D:/Decicomp/R/MOFA_omics/H2D_DEGS_context_reference_plot.png", H2D_DEGS_context_reference_plot, width = 12, height = 8, dpi = 600)

# Intersect T3 and T6 DEGS 
F14R_C3_T3_vs_T0_immunome_DEGS_list           <- BP_List_genes_F14R_C3_T3_vs_T0_immunome %>% 
                                                 filter(Gene_significant == "Yes") %>% 
                                                 mutate(Family = "F14R") %>% select(Gene) %>% unique()

F14R_C3_T6_vs_T0_immunome_DEGS_list           <- BP_List_genes_F14R_C3_T6_vs_T0_immunome %>% 
                                                 filter(Gene_significant == "Yes") %>% 
                                                 mutate(Family = "F14R") %>% select(Gene) %>% unique()

F14R_C3_T3_T6_vs_T0_immunome_DEGS_list <- rbind(F14R_C3_T3_vs_T0_immunome_DEGS_list, F14R_C3_T6_vs_T0_immunome_DEGS_list) %>% unique()

H2D_C3_T3_vs_T0_immunome_DEGS_list           <- BP_List_genes_H2D_C3_T3_vs_T0_immunome %>% 
                                                filter(Gene_significant == "Yes") %>% 
                                                mutate(Family = "H2D") %>% select(Gene) %>% unique()

H2D_C3_T6_vs_T0_immunome_DEGS_list           <- BP_List_genes_H2D_C3_T6_vs_T0_immunome %>% 
                                                filter(Gene_significant == "Yes") %>% 
                                                mutate(Family = "H2D") %>% select(Gene) %>% unique()

H2D_C3_T3_T6_vs_T0_immunome_DEGS_list <- rbind(H2D_C3_T3_vs_T0_immunome_DEGS_list, H2D_C3_T6_vs_T0_immunome_DEGS_list) %>% unique()
Intersect_F14R_H2D_T3_T6              <- intersect(F14R_C3_T3_T6_vs_T0_immunome_DEGS_list, H2D_C3_T3_T6_vs_T0_immunome_DEGS_list)





# COMPARING FAMILIES  Same AGE at T0 ----

# Now that we have check the immunome, we want to see when we compare same age but different families at T0

# 4 months between families (F14R C3  H2D C3 at T0 ) ----
F14R_C3_vs_H2D_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame() %>% rownames_to_column(var = "gene") %>% 
                              filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C3_vs_H2D_C3_T0_FC    <- left_join(List_genes, F14R_C3_vs_H2D_C3_T0_FC)
write(write.table(List_genes_F14R_C3_vs_H2D_C3_T0_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C3_vs_H2D_C3_T0/Immunome/List_genes_F14R_C3_vs_H2D_C3_T0_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C3_vs_H2D_C3_T0/Immunome/")
input="List_genes_F14R_C3_vs_H2D_C3_T0_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 


## Genes list to use F14R_C3_T3 vs F14R_C3_T0 GO_MWU
F14R_C3_vs_H2D_C3_T0          <- results(dds, contrast=c("group", "H2D_C3_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                                 as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C3_vs_H2D_C3_T0_sig_DEGS <- F14R_C3_vs_H2D_C3_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C3_vs_H2D_C3_T0_sig_DEGS)

F14R_C3_vs_H2D_C3_T0_sig_DEGS_result <- F14R_C3_vs_H2D_C3_T0_sig_DEGS %>%
                                        summarise(Up    = sum(log2FoldChange > 0),
                                                  Down  = sum(log2FoldChange < 0),
                                                  Total = Up + Down) %>% t() %>% as.data.frame() %>% 
                                        rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                        mutate(Contrast = "F14R_H2D", Age = "4") %>%
                                        dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C3_vs_H2D_C3_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C3_vs_H2D_C3_T0/Immunome/BP_List_genes_H2D_C3_F14R_C3_T0.txt") %>% 
                                                dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C3_vs_H2D_C3_T0, by = "Gene") 

gene_counts_F14R_C3_vs_H2D_C3_T0_immunome    <- BP_List_genes_F14R_C3_vs_H2D_C3_T0_immunome %>%
                                                dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                                nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_vs_H2D_C3_T0_immunome  <- BP_List_genes_F14R_C3_vs_H2D_C3_T0_immunome %>% 
                                                left_join(gene_counts_F14R_C3_vs_H2D_C3_T0_immunome, by = "term") %>% 
                                                mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C3_vs_H2D_C3_T0_immunome_DEGS           <- BP_List_genes_F14R_C3_vs_H2D_C3_T0_immunome %>% 
                                                filter(Gene_significant == "Yes") %>% 
                                                summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                                mutate(Contrast = "F14R_H2D")

F14R_C3_vs_H2D_C3_T0_sig_DEGS_result          <- merge(F14R_C3_vs_H2D_C3_T0_sig_DEGS_result, F14R_C3_vs_H2D_C3_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C3_vs_H2D_C3_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C3_vs_H2D_C3_T0/Immunome/MWU_BP_List_genes_F14R_C3_vs_H2D_C3_T0_FC.txt", sep="") %>% 
                                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                                left_join(BP_List_genes_F14R_C3_vs_H2D_C3_T0_immunome, by = c("term", "name")) %>% 
                                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                                mutate(Comparison = "F14R_C3_vs_H2D_C3_T0") %>% 
                                                left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


# 16 months between families (F14R C7  H2D C7 at T0 ) ----
F14R_C7_vs_H2D_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T0", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame() %>% rownames_to_column(var = "gene") %>% 
                              filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C7_vs_H2D_C7_T0_FC    <- left_join(List_genes, F14R_C7_vs_H2D_C7_T0_FC)
write(write.table(List_genes_F14R_C7_vs_H2D_C7_T0_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C7_vs_H2D_C7_T0/Immunome/List_genes_F14R_C7_vs_H2D_C7_T0_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C7_vs_H2D_C7_T0/Immunome/")
input="List_genes_F14R_C7_vs_H2D_C7_T0_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 


## Genes list to use F14R_C7_T3 vs F14R_C7_T0 GO_MWU
F14R_C7_vs_H2D_C7_T0          <- results(dds, contrast=c("group", "H2D_C7_T0", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                                 as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C7_vs_H2D_C7_T0_sig_DEGS <- F14R_C7_vs_H2D_C7_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C7_vs_H2D_C7_T0_sig_DEGS)

F14R_C7_vs_H2D_C7_T0_sig_DEGS_result <- F14R_C7_vs_H2D_C7_T0_sig_DEGS %>%
                                        summarise(Up    = sum(log2FoldChange > 0),
                                                  Down  = sum(log2FoldChange < 0),
                                                  Total = Up + Down) %>% t() %>% as.data.frame() %>% 
                                        rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                        mutate(Contrast = "F14R_H2D", Age = "16") %>%
                                        dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C7_T3 vs F14R_C7_T0 GO_MWU
BP_List_genes_F14R_C7_vs_H2D_C7_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C7_vs_H2D_C7_T0/Immunome/BP_List_genes_F14R_C7_vs_H2D_C7_T0_FC.txt") %>% 
                                                dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C7_vs_H2D_C7_T0, by = "Gene") 

gene_counts_F14R_C7_vs_H2D_C7_T0_immunome    <- BP_List_genes_F14R_C7_vs_H2D_C7_T0_immunome %>%
                                                dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                                nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_vs_H2D_C7_T0_immunome  <- BP_List_genes_F14R_C7_vs_H2D_C7_T0_immunome %>% 
                                                left_join(gene_counts_F14R_C7_vs_H2D_C7_T0_immunome, by = "term") %>% 
                                                mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C7_vs_H2D_C7_T0_immunome_DEGS           <- BP_List_genes_F14R_C7_vs_H2D_C7_T0_immunome %>% 
                                                filter(Gene_significant == "Yes") %>% 
                                                summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                                mutate(Contrast = "F14R_H2D")

F14R_C7_vs_H2D_C7_T0_sig_DEGS_result          <- merge(F14R_C7_vs_H2D_C7_T0_sig_DEGS_result, F14R_C7_vs_H2D_C7_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C7_T3 vs F14R_C7_T0 GO_MWU
GO_terms_F14R_C7_vs_H2D_C7_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C7_vs_H2D_C7_T0/Immunome/MWU_BP_List_genes_F14R_C7_vs_H2D_C7_T0_FC.txt", sep="") %>% 
                                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                                left_join(BP_List_genes_F14R_C7_vs_H2D_C7_T0_immunome, by = c("term", "name")) %>% 
                                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                                mutate(Comparison = "F14R_C7_vs_H2D_C7_T0") %>% 
                                                left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 



# 28 months between families (F14R C8  H2D C8 at T0 ) ----
F14R_C8_vs_H2D_C8_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T0", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame() %>% rownames_to_column(var = "gene") %>% 
                              filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_vs_H2D_C8_T0_FC    <- left_join(List_genes, F14R_C8_vs_H2D_C8_T0_FC)
write(write.table(List_genes_F14R_C8_vs_H2D_C8_T0_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C8_vs_H2D_C8_T0/Immunome/List_genes_F14R_C8_vs_H2D_C8_T0_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C8_vs_H2D_C8_T0/Immunome/")
input="List_genes_F14R_C8_vs_H2D_C8_T0_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 


## Genes list to use F14R_C8_T3 vs F14R_C8_T0 GO_MWU
F14R_C8_vs_H2D_C8_T0          <- results(dds, contrast=c("group", "H2D_C8_T0", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                                 as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_vs_H2D_C8_T0_sig_DEGS <- F14R_C8_vs_H2D_C8_T0 %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_vs_H2D_C8_T0_sig_DEGS)

F14R_C8_vs_H2D_C8_T0_sig_DEGS_result <- F14R_C8_vs_H2D_C8_T0_sig_DEGS %>%
                                        summarise(Up    = sum(log2FoldChange > 0),
                                        Down  = sum(log2FoldChange < 0),
                                        Total = Up + Down) %>% t() %>% as.data.frame() %>% 
                                        rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                        mutate(Contrast = "F14R_H2D", Age = "16") %>%
                                        dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C8_T3 vs F14R_C8_T0 GO_MWU
BP_List_genes_F14R_C8_vs_H2D_C8_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C8_vs_H2D_C8_T0/Immunome/BP_List_genes_F14R_C8_vs_H2D_C8_T0_FC.txt") %>% 
                                                dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_vs_H2D_C8_T0, by = "Gene") 

gene_counts_F14R_C8_vs_H2D_C8_T0_immunome    <- BP_List_genes_F14R_C8_vs_H2D_C8_T0_immunome %>%
                                                dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                                nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_vs_H2D_C8_T0_immunome  <- BP_List_genes_F14R_C8_vs_H2D_C8_T0_immunome %>% 
                                                left_join(gene_counts_F14R_C8_vs_H2D_C8_T0_immunome, by = "term") %>% 
                                                mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_vs_H2D_C8_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_vs_H2D_C8_T0_immunome %>% 
                                                filter(Gene_significant == "Yes") %>% 
                                                summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                                mutate(Contrast = "F14R_H2D")

F14R_C8_vs_H2D_C8_T0_sig_DEGS_result          <- merge(F14R_C8_vs_H2D_C8_T0_sig_DEGS_result, F14R_C8_vs_H2D_C8_T0_immunome_DEGS)


## Significant GO terms immunome F14R_C8_T3 vs F14R_C8_T0 GO_MWU
GO_terms_F14R_C8_vs_H2D_C8_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C8_vs_H2D_C8_T0/Immunome/MWU_BP_List_genes_F14R_C8_vs_H2D_C8_T0_FC.txt", sep="") %>% 
                                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                                left_join(BP_List_genes_F14R_C8_vs_H2D_C8_T0_immunome, by = c("term", "name")) %>% 
                                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                                mutate(Comparison = "F14R_C8_vs_H2D_C8_T0") %>% 
                                                left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 


# Heatmap comparing families  ---- 
GO_terms_F14R_C3_vs_H2D_C3_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C3_vs_H2D_C3_T0/Immunome/MWU_BP_List_genes_F14R_C3_vs_H2D_C3_T0_FC.txt", sep="")  %>% filter(p.adj <= 0.05)
GO_terms_F14R_C7_vs_H2D_C7_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C7_vs_H2D_C7_T0/Immunome/MWU_BP_List_genes_F14R_C7_vs_H2D_C7_T0_FC.txt", sep="")  %>% filter(p.adj <= 0.05)
GO_terms_F14R_C8_vs_H2D_C8_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_F14R_C8_vs_H2D_C8_T0/Immunome/MWU_BP_List_genes_F14R_C8_vs_H2D_C8_T0_FC.txt", sep="")  %>% filter(p.adj <= 0.05)

dim(GO_terms_F14R_C3_vs_H2D_C3_T0_immunome_sig) 
dim(GO_terms_F14R_C7_vs_H2D_C7_T0_immunome_sig)
dim(GO_terms_F14R_C8_vs_H2D_C8_T0_immunome_sig) 


GO_ages_F14R_vs_H2D <- rbind(       GO_terms_F14R_C3_vs_H2D_C3_T0_immunome, 
                                    GO_terms_F14R_C3_vs_H2D_C3_T0_immunome, 
                                    GO_terms_F14R_C3_vs_H2D_C3_T0_immunome) %>%
                                    filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% 
                                    distinct(name, .keep_all = TRUE) %>%
                                    dplyr::select(Comparison, name, Enrichment) %>%
                                   # rbind(new_rows)  %>%
                                    pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
                                    as.data.frame() %>% mutate_all(~ replace(., is.na(.), 0)) %>%
                                    column_to_rownames(var = "name")  %>% 
                                    mutate(H2D_C7_F14R_C7_T0 =0,
                                           H2D_C8_F14R_C8_T0 =0) 
                                     
# Column Colors heatmap
annotation_col           <- data.frame( Age = rep(c("4 months", "16 months", "28 months"), each = 1))
annotation_colors        <-  list(      Age = c(    "4 months" = "#FDC299",    "16 months" = "#EE820D", "28 months"  ="#FD6A00"))
rownames(annotation_col) <- colnames(GO_ages_F14R_vs_H2D) 

# Gradient enrichment                        
myBreaks <- c(seq(min(GO_ages_F14R_vs_H2D), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_ages_F14R_vs_H2D)/100, max(GO_ages_F14R_vs_H2D), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)


immuno_heatmap_families <- pheatmap(GO_ages_F14R_vs_H2D, 
                           cluster_cols = F, 
                           scale = "none",
                           cluster_rows = T, 
                           fontsize_row = 18, 
                           color = mycolor, 
                           breaks = myBreaks,
                           border_color = "black",
                           clustering_distance_rows = "euclidean",
                           show_colnames = F, 
                           show_rownames = T,
                           annotation_col = annotation_col, 
                           #annotation_row = annotation_row,
                           annotation_legend = FALSE,
                           legend = F,
                           legend_title = F,
                           annotation_names_col = FALSE,
                           annotation_colors = annotation_colors)
                           #cutree_rows = 4,
                           #gaps_col =  c(1,2))
immuno_heatmap_families
# ggsave("D:/Decicomp/R/MOFA_omics/immuno_heatmap_families.tiff", immuno_heatmap_families, width = 9, height = 7, dpi = 600)


# COMPARING H2D C3 vs FR14R C8 different times ----

# At T0 ----
F14R_C8_T0_vs_H2D_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T0", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                                 as.data.frame() %>% rownames_to_column(var = "gene") %>% 
                                 filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC    <- left_join(List_genes, F14R_C8_T0_vs_H2D_C3_T0_FC)
write(write.table(List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T0_vs_F14R_C8_T0/Immunome/List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T0_vs_F14R_C8_T0/Immunome/")
input="List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 



## Genes list to use F14R_C3_T3 vs F14R_C3_T0 GO_MWU
F14R_C8_T0_vs_H2D_C3_T0          <- results(dds, contrast=c("group", "H2D_C3_T0", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                                    as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T0_vs_H2D_C3_T0_sig_DEGS <- F14R_C8_T0_vs_H2D_C3_T0  %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T0_vs_H2D_C3_T0_sig_DEGS)

F14R_C8_T0_vs_H2D_C3_T0_sig_DEGS_result <- F14R_C8_T0_vs_H2D_C3_T0_sig_DEGS %>%
                                           summarise(Up    = sum(log2FoldChange > 0),
                                                     Down  = sum(log2FoldChange < 0),
                                                     Total = Up + Down) %>% t() %>% as.data.frame() %>% 
                                           rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
                                           mutate(Contrast = "H2D_C3_T0_F14R_C8_T0)", Time = "0") %>%
                                           dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T0 GO_MWU
BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T0_vs_F14R_C8_T0/Immunome/BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC.txt") %>% 
                                                   dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T0_vs_H2D_C3_T0, by = "Gene") 

gene_counts_F14R_C8_T0_vs_H2D_C3_T0_immunome    <- BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_immunome %>%
                                                   dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                                   nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_immunome  <- BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_immunome %>% 
                                                   left_join(gene_counts_F14R_C8_T0_vs_H2D_C3_T0_immunome, by = "term") %>% 
                                                   mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T0_vs_H2D_C3_T0_immunome_DEGS           <- BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_immunome %>% 
                                                   filter(Gene_significant == "Yes") %>% 
                                                   summarise(DEGs_immuno = n_distinct(Gene)) %>%
                                                   mutate(Contrast = "F14R_H2D")

F14R_C8_T0_vs_H2D_C3_T0_sig_DEGS_result          <- merge(F14R_C8_T0_vs_H2D_C3_T0_sig_DEGS_result, F14R_C8_T0_vs_H2D_C3_T0_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T0 GO_MWU
GO_terms_F14R_C8_T0_vs_H2D_C3_T0_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T0_vs_F14R_C8_T0/Immunome/MWU_BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC.txt", sep="") %>% 
                                                   mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                                   mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                                   left_join(BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_immunome, by = c("term", "name")) %>% 
                                                   mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                                   mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                                   mutate(Comparison = "H2D_C3_T0_F14R_C8_T0") %>% 
                                                   left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 
GO_terms_F14R_C8_T0_vs_H2D_C3_T0_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T0_vs_F14R_C8_T0/Immunome/MWU_BP_List_genes_F14R_C8_T0_vs_H2D_C3_T0_FC.txt", sep="") %>% filter(p.adj <= 0.05)


# At T3 ----
F14R_C8_T3_vs_H2D_C3_T3_FC    <- results(dds, contrast=c("group", "H2D_C3_T3", "F14R_C8_T3"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC    <- left_join(List_genes, F14R_C8_T3_vs_H2D_C3_T3_FC)
write(write.table(List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T3_vs_F14R_C8_T3/Immunome/List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T3_vs_F14R_C8_T3/Immunome/")
input="List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 


## Genes list to use F14R_C3_T3 vs F14R_C3_T3 GO_MWU
F14R_C8_T3_vs_H2D_C3_T3          <- results(dds, contrast=c("group", "H2D_C3_T3", "F14R_C8_T3"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T3_vs_H2D_C3_T3_sig_DEGS <- F14R_C8_T3_vs_H2D_C3_T3  %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T3_vs_H2D_C3_T3_sig_DEGS)

F14R_C8_T3_vs_H2D_C3_T3_sig_DEGS_result <- F14R_C8_T3_vs_H2D_C3_T3_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up + Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Contrast = "H2D_C3_T3_F14R_C8_T3)", Time = "3") %>%
  dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T3 GO_MWU
BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T3_vs_F14R_C8_T3/Immunome/BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC.txt") %>% 
  dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T3_vs_H2D_C3_T3, by = "Gene") 

gene_counts_F14R_C8_T3_vs_H2D_C3_T3_immunome    <- BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_immunome %>%
  dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_immunome  <- BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_immunome %>% 
  left_join(gene_counts_F14R_C8_T3_vs_H2D_C3_T3_immunome, by = "term") %>% 
  mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T3_vs_H2D_C3_T3_immunome_DEGS           <- BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Contrast = "F14R_H2D")

F14R_C8_T3_vs_H2D_C3_T3_sig_DEGS_result          <- merge(F14R_C8_T3_vs_H2D_C3_T3_sig_DEGS_result, F14R_C8_T3_vs_H2D_C3_T3_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T3 GO_MWU
GO_terms_F14R_C8_T3_vs_H2D_C3_T3_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T3_vs_F14R_C8_T3/Immunome/MWU_BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_immunome, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C3_T3_F14R_C8_T3") %>% 
  left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 
GO_terms_F14R_C8_T3_vs_H2D_C3_T3_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T3_vs_F14R_C8_T3/Immunome/MWU_BP_List_genes_F14R_C8_T3_vs_H2D_C3_T3_FC.txt", sep="") %>% filter(p.adj <= 0.05)



# At T6 ----
F14R_C8_T6_vs_H2D_C3_T6_FC    <- results(dds, contrast=c("group", "H2D_C3_T6", "F14R_C8_T6"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC    <- left_join(List_genes, F14R_C8_T6_vs_H2D_C3_T6_FC)
write(write.table(List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T6_vs_F14R_C8_T6/Immunome/List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T6_vs_F14R_C8_T6/Immunome/")
input="List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 



## Genes list to use F14R_C3_T3 vs F14R_C3_T6 GO_MWU
F14R_C8_T6_vs_H2D_C3_T6          <- results(dds, contrast=c("group", "H2D_C3_T6", "F14R_C8_T6"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T6_vs_H2D_C3_T6_sig_DEGS <- F14R_C8_T6_vs_H2D_C3_T6  %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T6_vs_H2D_C3_T6_sig_DEGS)

F14R_C8_T6_vs_H2D_C3_T6_sig_DEGS_result <- F14R_C8_T6_vs_H2D_C3_T6_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up + Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Contrast = "H2D_C3_T6_F14R_C8_T6)", Time = "6") %>%
  dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T6 GO_MWU
BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T6_vs_F14R_C8_T6/Immunome/BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC.txt") %>% 
  dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T6_vs_H2D_C3_T6, by = "Gene") 

gene_counts_F14R_C8_T6_vs_H2D_C3_T6_immunome    <- BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_immunome %>%
  dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_immunome  <- BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_immunome %>% 
  left_join(gene_counts_F14R_C8_T6_vs_H2D_C3_T6_immunome, by = "term") %>% 
  mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T6_vs_H2D_C3_T6_immunome_DEGS           <- BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Contrast = "F14R_H2D")

F14R_C8_T6_vs_H2D_C3_T6_sig_DEGS_result          <- merge(F14R_C8_T6_vs_H2D_C3_T6_sig_DEGS_result, F14R_C8_T6_vs_H2D_C3_T6_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T6 GO_MWU
GO_terms_F14R_C8_T6_vs_H2D_C3_T6_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T6_vs_F14R_C8_T6/Immunome/MWU_BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_immunome, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C3_T6_F14R_C8_T6") %>% 
  left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 
GO_terms_F14R_C8_T6_vs_H2D_C3_T6_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T6_vs_F14R_C8_T6/Immunome/MWU_BP_List_genes_F14R_C8_T6_vs_H2D_C3_T6_FC.txt", sep="") %>% filter(p.adj <= 0.05)


# At T12 ----
F14R_C8_T12_vs_H2D_C3_T12_FC    <- results(dds, contrast=c("group", "H2D_C3_T12", "F14R_C8_T12"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC    <- left_join(List_genes, F14R_C8_T12_vs_H2D_C3_T12_FC)
write(write.table(List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T12_vs_F14R_C8_T12/Immunome/List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T12_vs_F14R_C8_T12/Immunome/")
input="List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 



## Genes list to use F14R_C3_T3 vs F14R_C3_T12 GO_MWU
F14R_C8_T12_vs_H2D_C3_T12          <- results(dds, contrast=c("group", "H2D_C3_T12", "F14R_C8_T12"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T12_vs_H2D_C3_T12_sig_DEGS <- F14R_C8_T12_vs_H2D_C3_T12  %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T12_vs_H2D_C3_T12_sig_DEGS)

F14R_C8_T12_vs_H2D_C3_T12_sig_DEGS_result <- F14R_C8_T12_vs_H2D_C3_T12_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up + Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Contrast = "H2D_C3_T12_F14R_C8_T12)", Time = "12") %>%
  dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T12 GO_MWU
BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T12_vs_F14R_C8_T12/Immunome/BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC.txt") %>% 
  dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T12_vs_H2D_C3_T12, by = "Gene") 

gene_counts_F14R_C8_T12_vs_H2D_C3_T12_immunome    <- BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_immunome %>%
  dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_immunome  <- BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_immunome %>% 
  left_join(gene_counts_F14R_C8_T12_vs_H2D_C3_T12_immunome, by = "term") %>% 
  mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T12_vs_H2D_C3_T12_immunome_DEGS           <- BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Contrast = "F14R_H2D")

F14R_C8_T12_vs_H2D_C3_T12_sig_DEGS_result          <- merge(F14R_C8_T12_vs_H2D_C3_T12_sig_DEGS_result, F14R_C8_T12_vs_H2D_C3_T12_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T12 GO_MWU
GO_terms_F14R_C8_T12_vs_H2D_C3_T12_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T12_vs_F14R_C8_T12/Immunome/MWU_BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_immunome, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C3_T12_F14R_C8_T12") %>% 
  left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 
GO_terms_F14R_C8_T12_vs_H2D_C3_T12_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T12_vs_F14R_C8_T12/Immunome/MWU_BP_List_genes_F14R_C8_T12_vs_H2D_C3_T12_FC.txt", sep="") %>% filter(p.adj <= 0.05)



# At T24 ----
F14R_C8_T24_vs_H2D_C3_T24_FC    <- results(dds, contrast=c("group", "H2D_C3_T24", "F14R_C8_T24"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  filter(!is.na(padj)) %>% dplyr::select("gene", "log2FoldChange") %>% dplyr::rename(Gene=gene)

List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC    <- left_join(List_genes, F14R_C8_T24_vs_H2D_C3_T24_FC)
write(write.table(List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC, "D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T24_vs_F14R_C8_T24/Immunome/List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC.txt",      row.names = F, quote = F, col.names=F, sep = ","))

### Find GO Immunome 
setwd("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T24_vs_F14R_C8_T24/Immunome/")
input="List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC.txt" 
goAnnotations="immunome.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 



## Genes list to use F14R_C3_T3 vs F14R_C3_T24 GO_MWU
F14R_C8_T24_vs_H2D_C3_T24          <- results(dds, contrast=c("group", "H2D_C3_T24", "F14R_C8_T24"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 
F14R_C8_T24_vs_H2D_C3_T24_sig_DEGS <- F14R_C8_T24_vs_H2D_C3_T24  %>% dplyr::filter(padj <= 0.05) 
dim(F14R_C8_T24_vs_H2D_C3_T24_sig_DEGS)

F14R_C8_T24_vs_H2D_C3_T24_sig_DEGS_result <- F14R_C8_T24_vs_H2D_C3_T24_sig_DEGS %>%
  summarise(Up    = sum(log2FoldChange > 0),
            Down  = sum(log2FoldChange < 0),
            Total = Up + Down) %>% t() %>% as.data.frame() %>% 
  rename("Number"="V1") %>% rownames_to_column(var="DEGS") %>% 
  mutate(Contrast = "H2D_C3_T24_F14R_C8_T24)", Time = "24") %>%
  dplyr::select(3,4,1,2)

## GO terms and genes involved immunome F14R_C3_T3 vs F14R_C3_T24 GO_MWU
BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_immunome  <- read.delim("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T24_vs_F14R_C8_T24/Immunome/BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC.txt") %>% 
  dplyr::rename("Gene" = "seq") %>% dplyr:: filter(lev != -1) %>% left_join(F14R_C8_T24_vs_H2D_C3_T24, by = "Gene") 

gene_counts_F14R_C8_T24_vs_H2D_C3_T24_immunome    <- BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_immunome %>%
  dplyr::group_by(term) %>% dplyr::summarize(nseqs_relative = dplyr::n(),
                                             nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_immunome  <- BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_immunome %>% 
  left_join(gene_counts_F14R_C8_T24_vs_H2D_C3_T24_immunome, by = "term") %>% 
  mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

F14R_C8_T24_vs_H2D_C3_T24_immunome_DEGS           <- BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_immunome %>% 
  filter(Gene_significant == "Yes") %>% 
  summarise(DEGs_immuno = n_distinct(Gene)) %>%
  mutate(Contrast = "F14R_H2D")

F14R_C8_T24_vs_H2D_C3_T24_sig_DEGS_result          <- merge(F14R_C8_T24_vs_H2D_C3_T24_sig_DEGS_result, F14R_C8_T24_vs_H2D_C3_T24_immunome_DEGS)

## Significant GO terms immunome F14R_C3_T3 vs F14R_C3_T24 GO_MWU
GO_terms_F14R_C8_T24_vs_H2D_C3_T24_immunome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T24_vs_F14R_C8_T24/Immunome/MWU_BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_immunome, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C3_T24_F14R_C8_T24") %>% 
  left_join(Roseta, by = c("Gene"))
# # # # --- # # # # --- # # # # --- # # # # --- # # # # --- # # # # --- 
GO_terms_F14R_C8_T24_vs_H2D_C3_T24_immunome_sig  <- read.csv("D:/Decicomp/R/GO_terms/RNA/T0/GO_RNA_H2D_C3_T24_vs_F14R_C8_T24/Immunome/MWU_BP_List_genes_F14R_C8_T24_vs_H2D_C3_T24_FC.txt", sep="") %>% filter(p.adj <= 0.05)




# Heatmap comparing families different ages at same time infection time  ---- 

GO_ages_F14R_vs_H2D_times_infection <- rbind(       
                                    GO_terms_F14R_C8_T0_vs_H2D_C3_T0_immunome,
                                    GO_terms_F14R_C8_T3_vs_H2D_C3_T3_immunome, 
                                    GO_terms_F14R_C8_T6_vs_H2D_C3_T6_immunome, 
                                    GO_terms_F14R_C8_T12_vs_H2D_C3_T12_immunome,
                                    GO_terms_F14R_C8_T24_vs_H2D_C3_T24_immunome) %>%
  filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% 
  distinct(name, .keep_all = TRUE) %>%
  dplyr::select(Comparison, name, Enrichment) %>%
  # rbind(new_rows)  %>%
  pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
  as.data.frame() %>% mutate_all(~ replace(., is.na(.), 0)) %>%
  column_to_rownames(var = "name")  %>% 
  mutate(H2D_C3_T0_F14R_C8_T0 = 0,
         H2D_C3_T3_F14R_C8_T3 = 0,
         H2D_C3_T6_F14R_C8_T6 = 0,
         H2D_C3_T12_F14R_C8_T12=0) %>% 
  dplyr::select(2:5,1)
         

# Colum Colors heatmap   
annotation_col <- data.frame(
  Time = c("T0","T3", "T6", "T12", "T24"))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
  Time   = c(T0 = "white",  T3  = "#CD5C5C", T6 = "#B22222", T12="darkred",  T24 ="red"))

rownames(annotation_col) <- colnames(GO_ages_F14R_vs_H2D_times_infection) 


# Gradient enrichment                        
myBreaks <- c(seq(min(GO_ages_F14R_vs_H2D_times_infection), 0, length.out=ceiling(100/2) + 1),
              seq(max(GO_ages_F14R_vs_H2D_times_infection)/100, max(GO_ages_F14R_vs_H2D_times_infection), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("black", "yellow"))(100)

GO_ages_F14R_vs_H2D_times_infection_heatmap <- pheatmap(GO_ages_F14R_vs_H2D_times_infection, 
                                               cluster_cols = F, 
                                               scale = "none",
                                               cluster_rows = T, 
                                               fontsize_row = 20, 
                                               color = mycolor, 
                                               #breaks = myBreaks,
                                               border_color = "black",
                                               clustering_distance_rows = "euclidean",
                                               show_colnames = F, 
                                               show_rownames = T,
                                               annotation_col = annotation_col, 
                                               #annotation_row = annotation_row,
                                               annotation_colors = annotation_colors,
                                               annotation_names_col = FALSE,
                                               annotation_legend = FALSE,
                                               legend = F,
                                               legend_title = NULL)
                                               #cutree_rows = 4,
                                               #gaps_col =  c(5,10,15))
GO_ages_F14R_vs_H2D_times_infection_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/GO_ages_F14R_vs_H2D_times_infection_heatmap.tiff", GO_ages_F14R_vs_H2D_times_infection_heatmap, width = 13, height = 8, dpi = 300)



### EFFECTS OF AGEING in GENE EXPRESSION PER FAMILY ----

F14R_aging  <- GO_terms_F14R_C8_T0_vs_T0_immunome  %>% filter(name == "cellular response to virus") %>% 
               filter(Gene_significant=="Yes")     %>%  dplyr::select(Gene) # dim(F14R_aging) 172 

H2D_aging <-  GO_terms_H2D_C8_T0_vs_T0_immunome    %>% filter(name == "cellular response to virus") %>% 
              filter(Gene_significant=="Yes")      %>%  dplyr::select(Gene)  # dim(H2D_aging) 190

F14R_H2D_young <- GO_terms_F14R_C3_vs_H2D_C3_T0_immunome  %>% filter(name == "cellular response to virus") %>% 
                  filter(Gene_significant=="Yes") %>%  dplyr::select(Gene) # dim(F14R_H2D_young) 194

intersect_aging <- intersect(F14R_aging, H2D_aging)         # dim(intersect_aging) 104 common

exclusive_F14R_aging   <- setdiff(F14R_aging$Gene, H2D_aging$Gene) %>% 
                          as.data.frame()  %>% rename(Gene=".") # dim(exclusive_F14R_aging) 68

exclusive_H2D_aging   <- setdiff(H2D_aging$Gene, F14R_aging$Gene) %>% 
                         as.data.frame()  %>%  rename(Gene=".")  # dim(exclusive_H2D_aging) 86


# De los 68 cuanto tenias en joven  
#intersect_exclusive_A_C <- intersect(exclusive_A, C) # dim(intersect_exclusive_A_C) 31
# De los 104 cuantos tenias de joven 
# intersect_exclusive_A_C_intersect_A_B <- intersect(intersect_A_B, C) # dim(intersect_exclusive_A_C_intersect_A_B) 27

# De los 86 cuanto tenias en joven  
# intersect_exclusive_B_C <- intersect(exclusive_B, C) # dim(intersect_exclusive_B_C) 27
# De los 104 cuantos tenias de joven 
# intersect_exclusive_B_C_intersect_A_B <- intersect(intersect_A_B, C) # dim(intersect_exclusive_B_C_intersect_A_B) 27


# F14R_aging_FC heatmap 172 ---- 
F14R_aging_FC <- F14R_aging %>% 
                 left_join(F14R_C7_T0_vs_T0_FC) %>% rename(F14R_C7_FC = log2FoldChange)  %>%
                 left_join(F14R_C8_T0_vs_T0_FC) %>% rename(F14R_C8_FC = log2FoldChange)  %>%
                 mutate(F14R_C3_FC =0) %>%  dplyr::select(1,4,2,3)  %>%
                 mutate(across(everything(), ~ replace_na(., 0)))  %>% 
                 column_to_rownames(var = "Gene")

F14R_aging_FC_list <- F14R_aging %>% 
                      left_join(F14R_C7_T0_vs_T0_FC) %>% rename(F14R_C7_FC = log2FoldChange)  %>%
                      left_join(F14R_C8_T0_vs_T0_FC) %>% rename(F14R_C8_FC = log2FoldChange)  %>%
                      mutate(F14R_C3_FC =0) %>%  dplyr::select(1,4,2,3)  %>%
                      mutate(across(everything(), ~ replace_na(., 0)))   %>% left_join(Roseta)

write.table(F14R_aging_FC_list, "D:/Decicomp/R/MOFA_omics/F14R_aging_FC_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 

 
annotation_col <- data.frame(
  #Time = c("T0", "T0", "T0"),
  Age  = c(rep("4 months", 1), rep("16 months", 1), rep("28 months", 1)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
  Age    = c("4 months"   = "#dae3f3", "16 months"  = "#8faadc", "28 months"  ="#2f5597"))
  #Time   = c(T0 = "white"))

rownames(annotation_col) <- colnames(F14R_aging_FC) 

myBreaks <- c(seq(min(F14R_aging_FC), 0, length.out=ceiling(50/2) + 1),
              seq(max(F14R_aging_FC)/50, max(F14R_aging_FC), length.out=floor(50/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(50)

F14R_aging_FC_heatmap <- pheatmap(F14R_aging_FC, 
                                  cluster_cols = F, 
                                  #scale = "row",
                                  cluster_rows = T, 
                                  fontsize_row = 4, 
                                  color = mycolor, 
                                  breaks = myBreaks,
                                  border_color = "black",
                                  clustering_distance_rows = "euclidean",
                                  show_colnames = F, 
                                  show_rownames = F,
                                  cellwidth = 60,    # Adjust the cell width (you can change this value)
                                  #cellheight = 10, 
                                  annotation_col = annotation_col, 
                                  #annotation_row = annotation_row,
                                  annotation_names_col = FALSE,
                                  annotation_legend = FALSE,
                                  legend = F,
                                  legend_title = NULL,
                                  annotation_colors = annotation_colors)
#cutree_rows = 4,
# gaps_col =  c(5,10,15))
F14R_aging_FC_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/F14R_aging_FC_heatmap.tiff", F14R_aging_FC_heatmap, width = 8, height = 8, dpi = 300)



# H2D_aging_FC heatmap 190 ---- 
H2D_aging_FC <- H2D_aging %>% 
                left_join(H2D_C7_T0_vs_T0_FC) %>% rename(H2D_C7_FC = log2FoldChange)  %>%
                left_join(H2D_C8_T0_vs_T0_FC) %>% rename(H2D_C8_FC = log2FoldChange)  %>%
                mutate(H2D_C3_FC =0) %>%  dplyr::select(1,4,2,3)  %>%
                mutate(across(everything(), ~ replace_na(., 0)))  %>% 
                column_to_rownames(var = "Gene")

H2D_aging_FC_list <- H2D_aging %>% 
                     left_join(H2D_C7_T0_vs_T0_FC) %>% rename(H2D_C7_FC = log2FoldChange)  %>%
                     left_join(H2D_C8_T0_vs_T0_FC) %>% rename(H2D_C8_FC = log2FoldChange)  %>%
                     mutate(H2D_C3_FC =0) %>% dplyr::select(1,4,2,3)  %>%
                     mutate(across(everything(), ~ replace_na(., 0)))  %>% left_join(Roseta)

write.table(H2D_aging_FC_list, "D:/Decicomp/R/MOFA_omics/H2D_aging_FC_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 


annotation_col <- data.frame(
  #Time = c("T0", "T0", "T0"),
  Age  = c(rep("4 months", 1), rep("16 months", 1), rep("28 months", 1)))

annotation_colors  <-  list( #Family = c(H2D         = "#2f5597"),
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"))
  #Time   = c(T0 = "white"))

rownames(annotation_col) <- colnames(H2D_aging_FC) 

myBreaks <- c(seq(min(H2D_aging_FC), 0, length.out=ceiling(50/2) + 1),
              seq(max(H2D_aging_FC)/50, max(H2D_aging_FC), length.out=floor(50/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(50)

H2D_aging_FC_heatmap <- pheatmap(H2D_aging_FC, 
                                 cluster_cols = F, 
                                 #scale = "row",
                                 cluster_rows = T, 
                                 fontsize_row = 4, 
                                 color = mycolor, 
                                 breaks = myBreaks,
                                 border_color = T,
                                 clustering_distance_rows = "euclidean",
                                 show_colnames = F, 
                                 show_rownames = F,
                                 cellwidth = 60,    # Adjust the cell width (you can change this value)
                                 #cellheight = 10, 
                                 annotation_col = annotation_col, 
                                 #annotation_row = annotation_row,
                                 annotation_names_col = FALSE,
                                 annotation_legend = FALSE,
                                 legend = F,
                                 legend_title = F,
                                 annotation_colors = annotation_colors)
#cutree_rows = 4,
# gaps_col =  c(5,10,15))
H2D_aging_FC_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/H2D_aging_FC_heatmap.tiff", H2D_aging_FC_heatmap, width = 8, height = 8, dpi = 300)


# F14R and H2D aging_FC heatmap 104 ---- 
F14R_aging_FC <- F14R_aging %>% 
  left_join(F14R_C7_T0_vs_T0_FC) %>% rename(F14R_C7_FC = log2FoldChange)  %>%
  left_join(F14R_C8_T0_vs_T0_FC) %>% rename(F14R_C8_FC = log2FoldChange)  %>%
  mutate(F14R_C3_FC =0) %>%
  dplyr::select(1,4,2,3)  %>%
  mutate(across(everything(), ~ replace_na(., 0))) 

H2D_aging_FC <- H2D_aging %>% 
  left_join(H2D_C7_T0_vs_T0_FC) %>% rename(H2D_C7_FC = log2FoldChange)  %>%
  left_join(H2D_C8_T0_vs_T0_FC) %>% rename(H2D_C8_FC = log2FoldChange)  %>%
  mutate(H2D_C3_FC =0) %>%
  dplyr::select(1,4,2,3)  %>%
  mutate(across(everything(), ~ replace_na(., 0)))  

intersect_aging_FC_genes  <- intersect_aging %>% left_join(F14R_aging_FC) %>% left_join(H2D_aging_FC)  
# write.table(intersect_aging_FC_genes, "D:/Decicomp/R/MOFA_omics/intersect_aging_FC_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 
excluded_genes_F14R_genes <- setdiff( F14R_aging_FC$Gene, intersect_aging_FC_genes$Gene) %>% as.data.frame() %>%  dplyr::rename("Gene"=".") %>% left_join(F14R_aging_FC)
# write.table(excluded_genes_F14R_genes, "D:/Decicomp/R/MOFA_omics/excluded_genes_F14R_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 
excluded_genes_H2D_genes  <- setdiff( H2D_aging_FC$Gene,  intersect_aging_FC_genes$Gene) %>% as.data.frame() %>%  dplyr::rename("Gene"=".") %>% left_join(H2D_aging_FC)
# write.table(excluded_genes_H2D_genes, "D:/Decicomp/R/MOFA_omics/excluded_genes_H2D_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 


## Heatmap intersection 
intersect_aging_FC <- intersect_aging %>% left_join(F14R_aging_FC) %>% left_join(H2D_aging_FC) %>% column_to_rownames(var = "Gene")

# Define annotation for the columns
annotation_col <- data.frame( 
 # Time = rep(c("T0"), times = 6),
  Age  = c(rep("4 months_F14R", 1), rep("16 months_F14R", 1), rep("28 months_F14R", 1),
           rep("4 months_H2D", 1), rep("16 months_H2D", 1), rep("28 months_H2D", 1)))
# Define colors for annotations
annotation_colors <- list(
  Age = c(
    "4 months_F14R" = "#dae3f3", "16 months_F14R" = "#8faadc", "28 months_F14R" = "#2f5597",
    "4 months_H2D"  = "#E2F0D9", "16 months_H2D"  = "#A9D18E", "28 months_H2D"  = "#385700"))
 # Time = c(T0 ="white"))

# Assign rownames to annotation_col based on the colnames of the data
rownames(annotation_col) <- colnames(intersect_aging_FC)

myBreaks <- c(seq(min(intersect_aging_FC), 0, length.out=ceiling(50/2) + 1),
              seq(max(intersect_aging_FC)/50, max(intersect_aging_FC), length.out=floor(50/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(50)

intersect_aging_FC_heatmap <- pheatmap(intersect_aging_FC, 
                                       cluster_cols = F, 
                                       #scale = "row",
                                       cluster_rows = T, 
                                       fontsize_row = 6, 
                                       color = mycolor, 
                                       breaks = myBreaks,
                                       border_color = T,
                                       clustering_distance_rows = "euclidean",
                                       show_colnames = F, 
                                       show_rownames = F,
                                       cellwidth = 40,    # Adjust the cell width (you can change this value)
                                       #cellheight = 10, 
                                       annotation_col = annotation_col, 
                                       #annotation_row = annotation_row,
                                       annotation_colors = annotation_colors,
                                       annotation_names_col = FALSE,
                                       annotation_legend = FALSE,
                                       legend = F,
                                       legend_title = NULL,
                                       gaps_col =  c(3))
intersect_aging_FC_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/intersect_aging_FC_heatmap.tiff", intersect_aging_FC_heatmap, width = 8, height = 8, dpi = 300)







### EFFECTS OF AGEING in GENE EXPRESSION COMPARING FAMILY ----
GO_terms_F14R_C3_vs_H2D_C3_T0_immunome

Comparing_families_gene  <- GO_terms_F14R_C3_vs_H2D_C3_T0_immunome  %>% 
                            filter(name == "cellular response to virus") %>% filter(Gene_significant=="Yes") %>%  
                            select(Gene) # dim(F14R_aging)

Comparing_families_gene_FC <- Comparing_families_gene %>%   
  left_join(F14R_C3_vs_H2D_C3_T0_FC) %>% rename(Genetic_C3_FC = log2FoldChange)  %>%
  left_join(F14R_C7_vs_H2D_C7_T0_FC) %>% rename(Genetic_C7_FC = log2FoldChange)  %>% 
  left_join(F14R_C8_vs_H2D_C8_T0_FC) %>% rename(Genetic_C8_FC = log2FoldChange)  %>% 
  #filter(!Gene=="G9465") %>% 
  column_to_rownames(var = "Gene")

Comparing_families_gene_FC_list <- Comparing_families_gene %>% 
  left_join(F14R_C3_vs_H2D_C3_T0_FC) %>% rename(Genetic_C3_FC = log2FoldChange)  %>%
  left_join(F14R_C7_vs_H2D_C7_T0_FC) %>% rename(Genetic_C7_FC = log2FoldChange)  %>% 
  left_join(F14R_C8_vs_H2D_C8_T0_FC) %>% rename(Genetic_C8_FC = log2FoldChange)    %>% left_join(Roseta)

write.table(Comparing_families_gene_FC_list, "D:/Decicomp/R/MOFA_omics/Comparing_families_gene_FC_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 

annotation_col           <- data.frame( Age = rep(c("4 months", "16 months", "28 months"), each = 1))
annotation_colors        <- list(Age = c("4 months" = "#FDC299","16 months" = "#EE820D", "28 months"  ="#FD6A00"))
rownames(annotation_col) <- colnames(Comparing_families_gene_FC) 

# Gradient enrichment                        
myBreaks <- c(seq(min(Comparing_families_gene_FC), 0, length.out=ceiling(100/2) + 1),
              seq(max(Comparing_families_gene_FC)/100, max(Comparing_families_gene_FC), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)


Comparing_families_gene_FC_heatmap <- pheatmap(Comparing_families_gene_FC, 
                                          cluster_cols = F, 
                                          #scale = "row",
                                          cluster_rows = T, 
                                          fontsize_row = 4, 
                                          color = mycolor, 
                                          breaks = myBreaks,
                                          border_color = T,
                                          clustering_distance_rows = "euclidean",
                                          show_colnames = F, 
                                          show_rownames = T,
                                          cellwidth = 40,    # Adjust the cell width (you can change this value)
                                          #cellheight = 10, 
                                          annotation_col = annotation_col, 
                                          #annotation_row = annotation_row,
                                          annotation_colors = annotation_colors)

Comparing_families_gene_FC_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/Comparing_families_gene_FC_heatmap.png", Comparing_families_gene_FC_heatmap, width = 8, height = 8, dpi = 300)



# UNDERSTANDIN GENES FAMIY and AFFECTED BY TINE
Comparing_families_gene_FC <- Comparing_families_gene %>%   
                              left_join(F14R_C3_vs_H2D_C3_T0_FC) %>% rename(Genetic_C3_FC = log2FoldChange)  %>%
                              left_join(F14R_C7_vs_H2D_C7_T0_FC) %>% rename(Genetic_C7_FC = log2FoldChange)  %>% 
                              left_join(F14R_C8_vs_H2D_C8_T0_FC) %>% rename(Genetic_C8_FC = log2FoldChange) 



Comparing_families_gene_FC_aging <- intersect_aging %>% left_join(Comparing_families_gene_FC) %>%
                                    filter(if_any(starts_with("Genetic_C"), ~ !is.na(.)))  %>% 
                                    column_to_rownames(var = "Gene")


Comparing_families_gene_FC_aging_list  <- intersect_aging %>% left_join(Comparing_families_gene_FC) %>%
                                           filter(if_any(starts_with("Genetic_C"), ~ !is.na(.)))  %>%  left_join(Roseta)

write.table(Comparing_families_gene_FC_aging_list, "D:/Decicomp/R/MOFA_omics/Comparing_families_gene_FC_agingl_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 



annotation_col           <- data.frame( Age = rep(c("4 months", "16 months", "28 months"), each = 1))
annotation_colors        <-  list(      Age = c(    "4 months" = "#FDC299",    "16 months" = "#EE820D", "28 months"  ="#FD6A00"))
rownames(annotation_col) <- colnames(Comparing_families_gene_FC_aging) 

# Gradient enrichment                        
myBreaks <- c(seq(min(Comparing_families_gene_FC_aging), 0, length.out=ceiling(100/2) + 1),
              seq(max(Comparing_families_gene_FC_aging)/100, max(Comparing_families_gene_FC_aging), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)


Comparing_families_gene_FC_aging_heatmap  <- pheatmap(Comparing_families_gene_FC_aging, 
                                             cluster_cols = F, 
                                             #scale = "row",
                                             cluster_rows = T, 
                                             fontsize_row = 14, 
                                             color = mycolor, 
                                             breaks = myBreaks,
                                             border_color = T,
                                             clustering_distance_rows = "euclidean",
                                             show_colnames = F, 
                                             show_rownames = T,
                                             cellwidth = 40,    # Adjust the cell width (you can change this value)
                                             #cellheight = 10, 
                                             annotation_col = annotation_col, 
                                             #annotation_row = annotation_row,
                                             annotation_colors = annotation_colors)
Comparing_families_gene_FC_aging_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/Comparing_families_gene_FC_aging_heatmap.png", Comparing_families_gene_FC_aging_heatmap, width = 8, height = 8, dpi = 300)





# KEGGS pathway ----

# This data came from the 172 DEGS found in the "Cellular response to virus" using DAVID.

# F14R
KEEGS_imunomne_F14R <- read.delim("D:/Decicomp/R/MOFA_omics/Data_compacted/KEGG_immunome/F14R/KEEGS_imunomne_F14R.txt")
KEEGS_imunomne_F14R <- KEEGS_imunomne_F14R %>%
                         select(5, 2, 3, 8,6, 10) %>%
                         rename(FDR = PValue,
                         Pathway = Term,
                         nGenes = Count,
                         Pathway_Genes = Pop.Hits) %>%
                         mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
                         label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
                         filter(FDR <= 0.05)

KEEGS_imunomne_F14R_plot <- ggplot(KEEGS_imunomne_F14R, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
  geom_text(aes(label = label), hjust = -0.3, size = 5) +
  scale_fill_gradient2(low = "#dae3f3", mid = "#8faadc", high = "#2f5597") +
  scale_x_continuous(expand=c(0,0), limits = c(0, 12), breaks = c(0, 4,8,12)) +
  labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
#  title = "KEGG Pathways -log10(p-value)") +
#annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
#annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
#annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
Style_format_theme

KEEGS_imunomne_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_imunomne_F14R_plot.png", KEEGS_imunomne_F14R_plot, width = 10, height = 5, dpi = 300)

# H2D
KEEGS_imunomne_H2D <- read.delim("D:/Decicomp/R/MOFA_omics/Data_compacted/KEGG_immunome/H2D/KEEGS_imunomne_H2D.txt")
KEEGS_imunomne_H2D <- KEEGS_imunomne_H2D %>%
                       select(5, 2, 3, 8,6, 10) %>%
                       rename(FDR = PValue,
                       Pathway = Term,
                       nGenes = Count,
                       Pathway_Genes = Pop.Hits) %>%
                       mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
                       label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
                       filter(FDR <= 0.05)

KEEGS_imunomne_H2D_plot <- ggplot(KEEGS_imunomne_H2D, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
  geom_text(aes(label = label), hjust = -0.3, size = 5) +
  scale_fill_gradient2(low = "#E2F0D9", mid = "#A9D18E", high ="#385700") +
  scale_x_continuous(expand=c(0,0), limits = c(0, 12), breaks = c(0, 4,8,12)) +
  labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
  Style_format_theme
  
  KEEGS_imunomne_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_imunomne_H2D_plot.png", KEEGS_imunomne_H2D_plot, width = 10, height = 5, dpi = 300)


# Families
KEEGS_imunomne_families <- read.delim("D:/Decicomp/R/MOFA_omics/Data_compacted/KEGG_immunome/families/KEEGS_imunomne_families.txt")
KEEGS_imunomne_families <- KEEGS_imunomne_families %>%
                            select(5, 2, 3, 8,6, 10) %>%
                            rename(FDR = PValue,
                            Pathway = Term,
                            nGenes = Count,
                            Pathway_Genes = Pop.Hits) %>%
                            mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
                            label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
                            filter(FDR <= 0.05)
  
KEEGS_imunomne_families_plot <- ggplot(KEEGS_imunomne_families, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
    #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
    #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
    geom_col(color="black") +
    #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
    geom_text(aes(label = label), hjust = -0.3, size = 5) +
    scale_fill_gradient2(low = "#FDC299", mid = "#EE820D", high = "#FD6A00") +
    scale_x_continuous(expand=c(0,0), limits = c(0, 12), breaks = c(0, 4,8,12)) +
    labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
    #  title = "KEGG Pathways -log10(p-value)") +
    #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
    #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
    #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
    Style_format_theme
  
KEEGS_imunomne_families_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_imunomne_families_plot.png", KEEGS_imunomne_families_plot, width = 10, height = 5, dpi = 300)
  

# FOLD CHANGES PLOTS ----

# To fix 
F14R_aging_FC_list # 172 
F14R_C7_T0_vs_T0_FC
F14R_C8_T0_vs_T0_FC

H2D_aging_FC_list  # 190
H2D_C7_T0_vs_T0_FC
H2D_C8_T0_vs_T0_FC

Comparing_families_gene_FC_list # 194

Comparing_families_gene_FC_aging_list # 27

counts_norm 
counts_norm <- counts_norm %>% rownames_to_column("Gene")

## Genes list to use F14R_C8_T0 vs F14R_C3_T0 GO_MWU
F14R_C7_T0_vs_T0          <- results(dds, contrast=c("group", "F14R_C7_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 

## Genes list to use F14R_C8_T0 vs F14R_C3_T0 GO_MWU
F14R_C8_T0_vs_T0          <- results(dds, contrast=c("group", "F14R_C8_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 

# 4 months between families (F14R C3  H2D C3 at T0 ) 
F14R_C3_vs_H2D_C3_T0          <- results(dds, contrast=c("group", "H2D_C3_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 


# 4 months between families (F14R C3  H2D C8 at T0 ) 
F14R_C3_vs_H2D_C7_T0          <- results(dds, contrast=c("group", "H2D_C7_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 

# 4 months between families (F14R C3  H2D C8 at T0 ) 
F14R_C3_vs_H2D_C8_T0          <- results(dds, contrast=c("group", "H2D_C8_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% rownames_to_column(var = "Gene") %>% dplyr::filter(!is.na(padj)) 


gene <- "G26695" 

F14R_C7_T0_vs_T0_filter     <- F14R_C7_T0_vs_T0      %>% filter(Gene==gene) %>% select(Gene, log2FoldChange) %>% rename(FC =log2FoldChange) %>% mutate(Family="F14R", Age="C7")
F14R_C8_T0_vs_T0_filter     <- F14R_C8_T0_vs_T0      %>% filter(Gene==gene) %>% select(Gene, log2FoldChange) %>% rename(FC =log2FoldChange) %>% mutate(Family="F14R", Age="C8")
F14R_C3_vs_H2D_C3_T0_filter <- F14R_C3_vs_H2D_C3_T0  %>% filter(Gene==gene) %>% select(Gene, log2FoldChange) %>% rename(FC =log2FoldChange) %>% mutate(Family="H2D",  Age="C3")
F14R_C3_vs_H2D_C7_T0_filter <- F14R_C3_vs_H2D_C7_T0  %>% filter(Gene==gene) %>% select(Gene, log2FoldChange) %>% rename(FC =log2FoldChange) %>% mutate(Family="H2D",  Age="C7")
F14R_C3_vs_H2D_C8_T0_filter <- F14R_C3_vs_H2D_C8_T0  %>% filter(Gene==gene) %>% select(Gene, log2FoldChange) %>% rename(FC =log2FoldChange) %>% mutate(Family="H2D",  Age="C8")

Dataframe <- data.frame(
  Gene = gene,
  FC = 0,                # All fold changes are set to 0
  Family = "F14R",       # All entries have Family as "F14R"
  Age = "C3" )            # All entries have Age as "C3"


Fold_change_gene <- rbind(Dataframe, F14R_C7_T0_vs_T0_filter,F14R_C8_T0_vs_T0_filter, F14R_C3_vs_H2D_C3_T0_filter, F14R_C3_vs_H2D_C7_T0_filter, F14R_C3_vs_H2D_C8_T0_filter) %>% mutate(Family_Age = paste(Family, Age, sep = "_"))

Fold_change_gene_plot <-  ggplot(Fold_change_gene, aes(x = Age, y = FC, fill = Family_Age)) +
  geom_line(aes(group = Family, color = Family), size = 1, linetype ='solid') +
  geom_point(size = 8, color = "black", shape = 21) +  # Use shape = 21 to allow fill color
  scale_x_discrete(limits = c("C3", "C7", "C8"),  # Control the order of the x-axis labels
                   labels = c("C3" = "4", "C7" = "16", "C8" = "28"),
                   expand = expansion(mult = c(0.1, 0.1))) + # Add some space around the x-axis categories
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) + # Add some space around the y-axis
  scale_fill_manual(values = c("F14R_C3" = "#dae3f3", "F14R_C7" = "#8faadc", "F14R_C8" = "#2f5597",
                               "H2D_C3"  = "#E2F0D9", "H2D_C7"  = "#A9D18E", "H2D_C8"  = "#385700")) + 
  scale_color_manual(values = c("F14R" = "#2f5597", "H2D" = "#385700")) +
  labs(x = "Age (months)", y = "Log2 FC") + # Correct the legend title
  Style_format_theme  
Fold_change_gene_plot

# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Fold_change_gene_plot.png", Fold_change_gene_plot, width = 4, height = 4, dpi = 300)



#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --
#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --
#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --
#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --



# TRANSCRIPTOMICS ----

# TRANSCRIPTOMICS F14R ----
dds_norm      <- vst(dds, blind = F)
F14R_dds_norm <- assay(dds_norm) %>% as.data.frame() %>% dplyr::select(contains("F14R_C3_T0"), contains("F14R_C7_T0"), contains("F14R_C8_T0")) 

# For PCA
F14R_dds_norm_col_data <- colnames(F14R_dds_norm) %>% as.data.frame() %>% 
                            dplyr::rename("sample"=".") %>%
                            mutate(Family="F14R") %>%
                            mutate(Age = case_when(
                            grepl("C3", sample) ~ 4,
                            grepl("C7", sample) ~ 16,
                            grepl("C8", sample) ~ 28,
                            TRUE ~ NA_real_  )) %>%
                            mutate(Group = paste(Family, Age, sep = "_"))  %>%
                            tibble::column_to_rownames("sample")

F14R_dds_norm_t  <- assay(dds_norm) %>% as.data.frame() %>%
                    dplyr::select(contains("F14R_C3_T0"), 
                                  contains("F14R_C7_T0"), 
                                  contains("F14R_C8_T0"))  %>% t()

sample_pca_transcriptome               <- prcomp(F14R_dds_norm_t) 
summary_sample_pca_transcriptome       <- summary(sample_pca_transcriptome)
pca_scores_transcriptome_F14R          <- data.frame(sample = rownames(F14R_dds_norm_col_data), 
                                                      PC1 = sample_pca_transcriptome$x[, 1], 
                                                      PC2 = sample_pca_transcriptome$x[, 2], 
                                                      PC3 = sample_pca_transcriptome$x[, 3])

pca_scores_transcriptome_F14R  <- cbind(pca_scores_transcriptome_F14R, F14R_dds_norm_col_data)

PCA_F14R_transcriptomics  <- ggplot(data = pca_scores_transcriptome_F14R, aes(x = PC1 , y = PC2)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = Group, color=Group), linetype = 1)  +
  geom_point(aes(fill = Group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (29.5%)", y = "PC2 (12.6%)") +
  #scale_y_continuous(limits=c(-700,700)) +
  scale_x_continuous(limits=c(-60,60)) +
  # scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "F14R_16"="#8faadc",  "F14R_28"="#2f5597", "F14R_4"="#dae3f3" )) +
  scale_color_manual(values=c(  "F14R_16"="#8faadc", "F14R_28"="#2f5597",  "F14R_4"="#dae3f3" )) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_F14R_transcriptomics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_F14R_transcriptomics.tiff", PCA_F14R_transcriptomics, width = 10, height = 8, dpi = 600)

# For WGCNA
F14R_dds_norm <- assay(dds_norm) %>% as.data.frame() %>%
                 dplyr::select(contains("F14R_C3_T0"), contains("F14R_C7_T0"), contains("F14R_C8_T0")) %>% t()

Transcriptome_F14R_in_genes_list <- assay(dds_norm) %>% as.data.frame() %>% dplyr::select(contains("F14R_C3_T0"), contains("F14R_C7_T0"), contains("F14R_C8_T0")) %>% tibble::rownames_to_column("Gene")

power <- c(c(1:10), seq(from =1, to =30, by =1 ))
sft <- pickSoftThreshold(F14R_dds_norm, 
                         powerVector = power,
                         networkType = "signed", 
                         RsquaredCut = 0.95,
                         verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "scale free topology model filt, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "Mean Connectivity")+
  theme_classic()
grid.arrange(a1, a2, nrow=2)


F14R_dds_norm[] <- sapply(F14R_dds_norm, as.numeric)
#We want a power that is above 0.8 in R^2 AND a low mean connectivity
soft_power <- 13 
temp_cor   <- cor
cor        <- WGCNA::cor

bwnet_F14R_transcriptomic  <- blockwiseModules(F14R_dds_norm,
                               maxBlockSize = 25000,
                               networkType = "signed",
                               TOMType = "signed",
                               power = soft_power,
                               mergeCutHeight = 0.25,
                               minModuleSize= 100,
                               numericLabels = FALSE,
                               randomSeed=1234,
                               nThreads = 4,
                               verbose =3)

cor <- temp_cor
# Savethe network
# saveRDS(object = bwnet_F14R_transcriptomic, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_F14R_transcriptomic.RDS")

bwnet_F14R_transcriptomic <- readRDS("D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_F14R_transcriptomic.RDS")

module_eigengene_F14R_transcriptomic <- bwnet_F14R_transcriptomic$MEs
table(bwnet_F14R_transcriptomic$colors)

Dendrogram_plot <- plotDendroAndColors(bwnet_F14R_transcriptomic$dendrograms[[1]], 
                                       cbind(bwnet_F14R_transcriptomic$unmergedColors, 
                                             bwnet_F14R_transcriptomic$colors), 
                                       c("unmerged", "merged"),
                                       dendroLabels = FALSE,
                                       addGuide = F,
                                       hang = 0.03,
                                       guideHang = 0.05)
Dendrogram_plot

traits <- F14R_dds_norm_col_data %>% dplyr::select(Age) %>% 
          mutate(Age = case_when(Age == 4  ~ 1,
                                 Age == 16 ~ 2,
                                 Age == 28 ~ 3, TRUE ~ NA_real_))  

nSamples     <- nrow(F14R_dds_norm)
nCompounds   <- ncol(F14R_dds_norm)

module.trait.cor          <- cor(module_eigengene_F14R_transcriptomic, traits, use = 'p')
module.trait.cor.pvalues  <- corPvalueStudent(module.trait.cor, nSamples)

textMatrix      = paste(signif(module.trait.cor, 2), "\n(", signif(module.trait.cor.pvalues, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.cor)
textMatrix 

heatmap.data              <- merge(module_eigengene_F14R_transcriptomic, traits, by ='row.names')
heatmap.data              <- heatmap.data %>%  column_to_rownames(var='Row.names')
heatmap.data

Corrrelations_color_permissive <-  CorLevelPlot(heatmap.data,
                                                x = names(heatmap.data)[40],
                                                y = names(heatmap.data)[1:39],
                                                col=c("blue1", "skyblue", "white", "pink", "red"))
Corrrelations_color_permissive
# tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Corrrelations_F14R_permissive.tiff", width = 6, height = 8, units = "in", res = 600)
# print(Corrrelations_color_permissive)
# dev.off()

ME_positive_transcriptome_F14R <- module_eigengene_F14R_transcriptomic  %>%  dplyr::select(MEturquoise, MEbrown) %>%
                                  tibble::rownames_to_column(var = "samples") %>%
                                  gather(module, eigen_value, c(MEturquoise, MEbrown), factor_key = TRUE) %>%
                                  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                                                      ifelse(grepl("C7", samples), "16", 
                                                      ifelse(grepl("C8", samples), "28", NA)))) %>%
                                 mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
                                 arrange(Age)

ME_negative_transcriptome_F14R <- module_eigengene_F14R_transcriptomic  %>%  dplyr::select(MEpurple, MEblue, MEyellow, MElightcyan) %>%
                                  tibble::rownames_to_column(var = "samples") %>%
                                  gather(module, eigen_value, c(MEpurple, MEblue, MEyellow, MElightcyan), factor_key = TRUE) %>%
                                  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                                  ifelse(grepl("C7", samples), "16", 
                                  ifelse(grepl("C8", samples), "28", NA)))) %>%
                                  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
                                  arrange(Age)

# TO check the module
ggplot(ME_negative_transcriptome_F14R, aes(x = as.factor(Age), y = eigen_value, fill = module)) +
  geom_boxplot() +
  #ggtitle("MEblue") +
  #scale_fill_manual(values = c("red", "yellow","turquoise")) +  # Set custom colors
  labs(x = "Condition", y = "Eigen Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) 

# Important
module.gene.mapping_transcriptomics_F14R <- data.frame(cluster = as.character(bwnet_F14R_transcriptomic$colors), 
                                            Gene = names(bwnet_F14R_transcriptomic$colors), stringsAsFactors = FALSE)

positive_F14R_transcriptome    <- module.gene.mapping_transcriptomics_F14R  %>%  filter(cluster == "turquoise" | cluster == "brown") 
table(positive_F14R_transcriptome$cluster)

negative_F14R_transcriptome    <- module.gene.mapping_transcriptomics_F14R  %>%  filter(cluster == "purple" | cluster == "blue"   |
                                                                                        cluster == "yellow" | cluster == "lightcyan")
table(negative_F14R_transcriptome$cluster)
                                                                                           
geneModuleMembership_F14R <- as.data.frame(cor(F14R_dds_norm, module_eigengene_F14R_transcriptomic, use='p')) 


# Extract the GO terms and KEGGs positive
kME.positive_F14R_turquoise  <-  geneModuleMembership_F14R %>% dplyr::select(., MEturquoise) %>% rename(kME=MEturquoise)
kME.positive_F14R_brown      <-  geneModuleMembership_F14R %>% dplyr::select(., MEbrown)     %>% rename(kME=MEbrown)

kME.positive_F14R_transcriptome <- rbind(kME.positive_F14R_turquoise, kME.positive_F14R_brown) %>% dplyr::filter(., rownames(.) %in% positive_F14R_transcriptome$Gene) %>% rownames_to_column("Gene") %>% dplyr::filter(!is.na(kME))
# write.table(kME.positive_F14R_transcriptome, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_F14R_transcriptome.txt", sep = "\t", row.names = T, col.names = T)

List_kME.positive_F14R_transcriptome_clean <- kME.positive_F14R_transcriptome   %>% dplyr::select(Gene)                  
# write(write.table(List_kME.positive_F14R_clean, "D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/List_kME.positive_F14R_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.positive_F14R_transcriptome <- left_join(List_genes, kME.positive_F14R_transcriptome, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.positive_F14R_transcriptome, "D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/List_genes_kME.positive_F14R.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/")
input="List_genes_kME.positive_F14R.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

GO_terms_genes_kME.positive_F14R_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/MWU_BP_List_genes_kME.positive_F14R.txt", sep="") %>%
                                                   filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
                                                   mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
                                                   dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
                                                   mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(GO_terms_genes_kME.positive_F14R_transcriptome, "D:/Decicomp/Paper/Paper Valdi/GO_terms_genes_kME.positive_F14R_transcriptome.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

category_counts_positive_F14R_transcriptome  <-    GO_terms_genes_kME.positive_F14R_transcriptome %>% group_by(term_parent) %>% summarise(Count = n()) %>% arrange(desc(Count))

# GO_terms_genes_kME.positive_F14R_plot <- ggplot(GO_terms_genes_kME.positive_F14R_transcriptome,
       # aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
# geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
# geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
#  scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
#  scale_fill_manual(values = c("cellular process"      = "tomato",  
# "metabolic process"     = "red4",  
#  "localization"          = "grey80",  
#    "reproductive process"  = "grey80",
#   "biological regulation" = "grey80", 
#   "biological process involved in interspecies interaction between organisms"= "red",
                               #    "developmental process"            = "grey80",
#   "response to stimulus"             = "orange",
#  "multicellular organismal process" = "grey80")) +
#labs(x = "Category and GO Term", y = "Frequency") +
# theme(legend.position = "none", 
#  panel.background = element_blank(),  
#  panel.grid = element_blank(),        
#  axis.title = element_blank(),     
#  axis.text.y = element_blank(),       
#  axis.ticks = element_blank(),     
#   axis.text.x = element_text(size = 14,  hjust = 1)) 
#GO_terms_genes_kME.positive_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_genes_kME.positive_F14R_plot.png", GO_terms_genes_kME.positive_F14R_plot, width = 10, height = 10, dpi = 400)

GO_terms_F14R_transcriptome_positive          <- GO_terms_genes_kME.positive_F14R_transcriptome$last_term  %>% as.vector()
matrix_GO_terms_F14R_transcriptome_positive   <- GO_similarity(GO_terms_F14R_transcriptome_positive)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_F14R_transcriptome_positive.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_F14R_transcriptome_positive, fontsize_range = c(1, 20), order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1),show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_F14R_transcriptome_positive_simply <- simplifyGO(matrix_GO_terms_F14R_transcriptome_positive, fontsize_range = c(1, 20), order_by_size = F, 
                                                          bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1),
                                                          show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_F14R_transcriptome_positive_simply <- GO_terms_F14R_transcriptome_positive_simply %>% dplyr::rename(last_term=id) %>% left_join(GO_terms_genes_kME.positive_F14R_transcriptome)

# F14R KEGGS positive transcriptome
# KEEGS_positive_F14R <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/KEGGS_positive_F14R.txt")
# KEEGS_positive_F14R <- KEEGS_positive_F14R %>%
#  dplyr::select(5, 2, 3, 8,6, 10) %>%
# rename(FDR = PValue,
#  Pathway = Term,
#    nGenes = Count,
#   Pathway_Genes = Pop.Hits) %>%
# mutate(Pathway = gsub("crg\\d{5}:", "", Pathway), label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%  filter(FDR <= 0.05)

# KEEGS_positive_F14R_plot <- ggplot(KEEGS_positive_F14R, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
# geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
# geom_text(aes(label = label), hjust = -0.3, size = 5) +
# scale_fill_gradient2(low = "#dae3f3", mid = "#8faadc", high = "#2f5597") +
#  scale_x_continuous(expand=c(0,0), limits = c(0, 4), breaks = c(0, 2,4)) +
# labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
# Style_format_theme

# KEEGS_positive_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_positive_F14R_plot.png", KEEGS_positive_F14R_plot, width = 14, height = 8, dpi = 600)


# Extract the GO terms and KEGGs negative
kME.negative_F14R_blue       <- geneModuleMembership_F14R %>% dplyr::select(MEblue)      %>% rename(kME=MEblue)      %>% dplyr::filter(!is.na(kME))
kME.negative_F14R_purple     <- geneModuleMembership_F14R %>% dplyr::select(MEpurple)    %>% rename(kME=MEpurple)    %>% dplyr::filter(!is.na(kME))
kME.negative_F14R_yellow     <- geneModuleMembership_F14R %>% dplyr::select(MEyellow)    %>% rename(kME=MEyellow)    %>% dplyr::filter(!is.na(kME))
kME.negative_F14R_lightcyan  <- geneModuleMembership_F14R %>% dplyr::select(MElightcyan) %>% rename(kME=MElightcyan) %>% dplyr::filter(!is.na(kME))

kME.negative_F14R_transcriptome  <- rbind(kME.negative_F14R_blue, kME.negative_F14R_purple, kME.negative_F14R_yellow, kME.negative_F14R_lightcyan) %>%
                                    filter(., rownames(.) %in% negative_F14R_transcriptome$Gene) %>% rownames_to_column("Gene")
# write.table(kME.negative_F14R_transcriptome, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_F14R_transcriptome.txt", sep = "\t", row.names = T, col.names = T)

List_kME.negative_F14R_transcriptome_clean <- kME.negative_F14R_transcriptome %>% dplyr::select(Gene)                  
# write(write.table(List_kME.negative_F14R_clean, "D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/List_kME.negative_F14R_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.negative_F14R_transcriptome    <- left_join(List_genes, kME.negative_F14R_transcriptome, by = "Gene") %>% mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.negative_F14R_transcriptome, "D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/List_genes_kME.negative_F14R.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/")
input="List_genes_kME.negative_F14R.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

GO_terms_genes_kME.negative_F14R_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/MWU_BP_List_genes_kME.negative_F14R.txt", sep="") %>%
                                                    filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>% filter(name != "unknown")  %>%
                                                    mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
                                                    dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
                                                    mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(GO_terms_genes_kME.negative_F14R_transcriptome, "D:/Decicomp/Paper/Paper Valdi/GO_terms_genes_kME.negative_F14R_transcriptome.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

category_counts_negative_F14R_transcriptome  <- GO_terms_genes_kME.negative_F14R_transcriptome %>% group_by(term_parent) %>% summarise(Count = n()) %>% arrange(desc(Count))

# GO_terms_genes_kME.negative_F14R_plot <- ggplot(category_counts_negative_F14R_transcriptome,
#  aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
#  geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
#    scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
#     scale_fill_manual(values = c( "cellular process"      = "lightblue", 
#  "localization"          = "grey80",  
#    "developmental process" = "grey80",
# "metabolic process"     = "blue4", 
#   "biological regulation" = "grey80", 
#    "positive regulation of biological process"= "grey80",
#    "reproductive process"  = "grey80",
#     "multicellular organismal process" = "grey80",
#     "response to stimulus"             = "orange",
#      "immune system process"            = "darkorange",
#        "biological process involved in interspecies interaction between organisms"= "blue" )) +
#        labs(x = "Category and GO Term", y = "Frequency") +
# theme(legend.position = "none", 
#      panel.background = element_blank(),  
#       panel.grid = element_blank(),        
#      axis.title = element_blank(),     
#      axis.text.y = element_blank(),       
#     axis.ticks = element_blank(),     
#     axis.text.x = element_text(size = 14,  hjust = 1)) 
# GO_terms_genes_kME.negative_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_genes_kME.negative_F14R_plot.png", GO_terms_genes_kME.negative_F14R_plot, width = 10, height = 10, dpi = 400)

GO_terms_F14R_transcriptome_negative         <- GO_terms_genes_kME.negative_F14R_transcriptome$last_term  %>% as.vector()
matrix_GO_terms_F14R_transcriptome_negative  <- GO_similarity(GO_terms_F14R_transcriptome_negative)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_F14R_transcriptome_negative.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_F14R_transcriptome_negative, fontsize_range = c(1, 20),order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1), show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_F14R_transcriptome_negative_simply <- simplifyGO(matrix_GO_terms_F14R_transcriptome_negative, fontsize_range = c(1, 20),order_by_size = F, 
                                                          bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1), show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_F14R_transcriptome_negative_simply <- GO_terms_F14R_transcriptome_negative_simply %>% dplyr::rename(last_term=id) %>% left_join(GO_terms_genes_kME.negative_F14R_transcriptome)

# F14R KEGGS negative transcriptome
# KEEGS_negative_F14R <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/KEGGS_negative_F14R.txt")
# KEEGS_negative_F14R <- KEEGS_negative_F14R %>%
# dplyr::select(5, 2, 3, 8,6, 10) %>%
#  rename(FDR = PValue,
# Pathway = Term,
#  nGenes = Count,
#   Pathway_Genes = Pop.Hits) %>%
# mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
#       label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>% filter(FDR <= 0.05)

#KEEGS_negative_F14R_plot <- ggplot(KEEGS_negative_F14R, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
# geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
# geom_text(aes(label = label), hjust = -0.3, size = 5) +
# scale_fill_gradient2(low = "#dae3f3", mid = "#8faadc", high = "#2f5597") +
# scale_x_continuous(expand=c(0,0), limits = c(0, 4), breaks = c(0, 2,4)) +
#  labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
# Style_format_theme

#KEEGS_negative_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_negative_F14R_plot.png", KEEGS_negative_F14R_plot, width = 14, height = 8, dpi = 600)

# Heat Map F14R transcriptome
Transcriptome_F14R_in_genes_list
positive_F14R_transcriptome
negative_F14R_transcriptome

F14R_list_transcriptome_associated <- rbind(positive_F14R_transcriptome, negative_F14R_transcriptome) %>% dplyr::select(Gene) %>% 
                                      left_join(Transcriptome_F14R_in_genes_list) %>% column_to_rownames(var = "Gene")

# Colum Colors heatmap   
annotation_col <- data.frame(Age  = c(rep("4 months", 6), rep("16 months", 5), rep("28 months", 6)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
                              Age    = c("4 months"   = "#dae3f3", "16 months"  = "#8faadc", "28 months"  ="#2f5597"))

rownames(annotation_col) <- colnames(F14R_list_transcriptome_associated) 

# Gradient enrichment                        
myBreaks <- c(seq(min(F14R_list_transcriptome_associated), 0, length.out=ceiling(100/2) + 1),
              seq(max(F14R_list_transcriptome_associated)/100, max(F14R_list_transcriptome_associated), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "red"))(100)

F14R_list_transcriptome_genes_associated_heatmap <- pheatmap(F14R_list_transcriptome_associated, 
                                                             cluster_cols =F, 
                                                             scale = "row",
                                                             cluster_rows = F, 
                                                             fontsize_row = 1, 
                                                             color = mycolor, 
                                                             # breaks = myBreaks,
                                                             border_color = "black",
                                                             clustering_distance_rows = "euclidean",
                                                             show_colnames = F, 
                                                             show_rownames = F,
                                                             annotation_col = annotation_col, 
                                                             #annotation_row = annotation_row,
                                                             annotation_colors = annotation_colors,
                                                             legend = FALSE)
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/F14R_list_transcriptome_genes_associated_heatmap.tiff", 
       F14R_list_transcriptome_genes_associated_heatmap, 
       width =  5000 / 600,  # Convert to inches
       height = 5000 / 600,  # Convert to inches
       dpi = 600)


# sacar lista de genes en los GO terms.

BP_List_genes_kME.positive_F14R <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/BP_List_genes_kME.positive_F14R.txt") %>% 
  filter(name=="response to virus" | 
           name=="defense response to virus" | 
           name=="defense response to symbiont" |
           name=="cellular response to endogenous stimulus") %>% 
  filter(value > 0) %>%  dplyr::select(seq) %>% dplyr::rename(Gene=seq) %>% distinct()

BP_List_genes_kME.negative_F14R <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/BP_List_genes_kME.negative_F14R.txt") %>% 
    filter(name=="immune response" | 
           name=="biological process involved in interspecies interaction between organisms" | 
           name=="defense response" |
           name=="toll-like receptor signaling pathway" |
           name=="defense response to other organism") %>% 
  filter(value > 0) %>%  dplyr::select(seq) %>% dplyr::rename(Gene=seq) %>% distinct()


Transcriptome_F14R_in_genes_list_counts <- assay(dds) %>% as.data.frame() %>%
  dplyr::select(contains("F14R_C3_T0"), contains("F14R_C7_T0"), contains("F14R_C8_T0")) %>%
  tibble::rownames_to_column("Gene")

Transcriptome_F14R_in_genes_list_wide <- Transcriptome_F14R_in_genes_list_counts  %>%
  pivot_longer(cols = starts_with("F14R"), names_to = "Sample", values_to = "Expression") %>% 
  mutate(Age = case_when(grepl("C3", Sample) ~ 4,
                         grepl("C7", Sample) ~ 16,
                         grepl("C8", Sample) ~ 28, TRUE ~ NA_real_ )) %>% 
  group_by(Gene, Age) %>%
  summarize(Mean_Expression = mean(Expression, na.rm = TRUE),  # Media de la metilación
            SEM = sd(Expression, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

Gene_expression_I_F14R <- Transcriptome_F14R_in_genes_list_wide %>% filter(Gene=="G13495")
max_expression_value <- max(Gene_expression_I_F14R$Mean_Expression + Gene_expression_I_F14R$SEM, na.rm = TRUE)
buffered_max_value   <- max_expression_value * 1.1

Gene_expression_I_plot_F14R <- ggplot(Gene_expression_I, aes(x = factor(Age), y = Mean_Expression, fill = factor(Age))) +
  geom_line(aes(group = 1), color = "#2f5597", linewidth = 1) +
  geom_errorbar(aes(ymin = Mean_Expression - SEM, ymax = Mean_Expression + SEM), width = 0.1) +  
  geom_point(size = 6, shape=21) +  # Añadir los puntos
  scale_fill_manual(values = c("4" = "#dae3f3", "16" = "#8faadc", "28" = "#2f5597")) +  # Relleno personalizado
  scale_y_continuous(limits = c(0, buffered_max_value), expand = expansion(mult = c(0, 0.1))) +  
  labs(x = "Age", y = "Counts")+
  Style_format_theme  

Gene_expression_I_plot_F14R
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Gene_expression_I_plot_F14R.png", Gene_expression_I_plot_F14R, width = 5, height = 4, dpi = 600)



# TRANSCRIPTOMICS H2D ----
dds_norm    <- vst(dds, blind = F)
H2D_dds_norm <- assay(dds_norm) %>% as.data.frame() %>% dplyr::select(contains("H2D_C3_T0"), contains("H2D_C7_T0"), contains("H2D_C8_T0")) 

# For PCA
H2D_dds_norm_col_data <- colnames(H2D_dds_norm) %>% as.data.frame() %>% 
                         dplyr::rename("sample"=".") %>% mutate(Family="H2D") %>% 
                         mutate(Age = case_when(
                              grepl("C3", sample) ~ 4,
                              grepl("C7", sample) ~ 16,
                              grepl("C8", sample) ~ 28,
                              TRUE ~ NA_real_  )) %>%
                              mutate(Group = paste(Family, Age, sep = "_"))  %>%
                             tibble::column_to_rownames("sample")

H2D_dds_norm_t  <- assay(dds_norm) %>% as.data.frame() %>%
  dplyr::select(contains("H2D_C3_T0"), contains("H2D_C7_T0"), contains("H2D_C8_T0"))  %>% t()

sample_pca_transcriptome              <- prcomp(H2D_dds_norm_t) 
summary_sample_pca_transcriptome      <- summary(sample_pca_transcriptome)
pca_scores_transcriptome_H2D          <- data.frame(sample = rownames(H2D_dds_norm_col_data), 
                                                    PC1 = sample_pca_transcriptome$x[, 1], 
                                                    PC2 = sample_pca_transcriptome$x[, 2], 
                                                    PC3 = sample_pca_transcriptome$x[, 3])

pca_scores_transcriptome_H2D  <- cbind(pca_scores_transcriptome_H2D, 
                                       H2D_dds_norm_col_data)

PCA_H2D_transcriptomics  <- ggplot(data = pca_scores_transcriptome_H2D, aes(x = PC1 , y = PC2)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = Group, color=Group), linetype = 1)  +
  geom_point(aes(fill = Group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (25.7%)", y = "PC2 (13.7%)") +
  #scale_y_continuous(limits=c(-700,700)) +
  scale_x_continuous(limits=c(-60,60)) +
  # scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "H2D_16"="#A9D18E",  "H2D_28"="#385700", "H2D_4"="#E2F0D9" )) +
  scale_color_manual(values=c(  "H2D_16"="#A9D18E", "H2D_28"="#385700",  "H2D_4"="#E2F0D9" )) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_H2D_transcriptomics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_H2D_transcriptomics.tiff", PCA_H2D_transcriptomics, width = 10, height = 8, dpi = 600)

# For WGCNA
H2D_dds_norm_col_data <- colnames(H2D_dds_norm) %>% as.data.frame() %>% 
  dplyr::rename("sample"=".") %>%
  mutate(Family="H2D") %>%
  mutate(Age = case_when(
    grepl("C3", sample) ~ 4,
    grepl("C7", sample) ~ 16,
    grepl("C8", sample) ~ 28,
    TRUE ~ NA_real_  )) %>%
  mutate(Group = paste(Family, Age, sep = "_"))  %>%
  tibble::column_to_rownames("sample")

H2D_dds_norm <- assay(dds_norm) %>% as.data.frame() %>%
                dplyr::select(contains("H2D_C3_T0"), contains("H2D_C7_T0"), contains("H2D_C8_T0")) %>% t()

Transcriptome_H2D_in_genes_list <- assay(dds_norm) %>% as.data.frame() %>%
  dplyr::select(contains("H2D_C3_T0"), contains("H2D_C7_T0"), contains("H2D_C8_T0")) %>% tibble::rownames_to_column("Gene")

power <- c(c(1:10), seq(from =1, to =30, by =1 ))
sft <- pickSoftThreshold(H2D_dds_norm, 
                         powerVector = power,
                         networkType = "signed", 
                         RsquaredCut = 0.95,
                         verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "scale free topology model filt, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "Mean Connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow=2)


H2D_dds_norm[] <- sapply(H2D_dds_norm, as.numeric)
#We want a power that is above 0.8 in R^2 AND a low mean connectivity
soft_power <- 13 
temp_cor   <- cor
cor        <- WGCNA::cor

bwnet_H2D_transcriptomic      <- blockwiseModules(H2D_dds_norm,
                                                  maxBlockSize = 30000,
                                                  networkType = "signed",
                                                  TOMType = "signed",
                                                  power = soft_power,
                                                  mergeCutHeight = 0.25,
                                                  minModuleSize= 100,
                                                  numericLabels = FALSE,
                                                  randomSeed=1234,
                                                  nThreads = 4,
                                                  verbose =3)

cor <- temp_cor
# Savethe network
# saveRDS(object = bwnet_H2D_transcriptomic, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_H2D_transcriptomic.RDS")

bwnet_H2D_transcriptomic <- readRDS("D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_H2D_transcriptomic.RDS")
module_eigengene_H2D_transcriptomic <- bwnet_H2D_transcriptomic$MEs
table(bwnet_H2D_transcriptomic$colors)

Dendrogram_plot <- plotDendroAndColors(bwnet_H2D_transcriptomic$dendrograms[[1]], 
                                       cbind(bwnet_H2D_transcriptomic$unmergedColors, 
                                             bwnet_H2D_transcriptomic$colors), 
                                       c("unmerged", "merged"),
                                       dendroLabels = FALSE,
                                       addGuide = F,
                                       hang = 0.03,
                                       guideHang = 0.05)
Dendrogram_plot

traits <- H2D_dds_norm_col_data %>% dplyr::select(Age) %>% 
  mutate(Age = case_when(Age == 4  ~ 1,
                         Age == 16 ~ 2,
                         Age == 28 ~ 3, TRUE ~ NA_real_))  

nSamples     <- nrow(H2D_dds_norm)
nCompounds   <- ncol(H2D_dds_norm)

module.trait.cor          <- cor(module_eigengene_H2D_transcriptomic, traits, use = 'p')
module.trait.cor.pvalues  <- corPvalueStudent(module.trait.cor, nSamples)

textMatrix =  paste(signif(module.trait.cor, 2), "\n(",
                    signif(module.trait.cor.pvalues, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.cor)
textMatrix 

heatmap.data              <- merge(module_eigengene_H2D_transcriptomic, traits, by ='row.names')
heatmap.data              <- heatmap.data %>%  tibble::column_to_rownames(var='Row.names')
heatmap.data

Corrrelations_color_permissive <-  CorLevelPlot(heatmap.data,
                                                x = names(heatmap.data)[39],
                                                y = names(heatmap.data)[1:38],
                                                col=c("blue1", "skyblue", "white", "pink", "red"))
Corrrelations_color_permissive

# tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Corrrelations_H2D_permissive.tiff", width = 6, height = 8, units = "in", res = 600)
# print(Corrrelations_color_permissive)
# dev.off()

ME_positive_transcriptome_H2D <- module_eigengene_H2D_transcriptomic  %>%  dplyr::select(MEblue) %>%
  tibble::rownames_to_column(var = "samples") %>%
  gather(module, eigen_value, c(MEblue), factor_key = TRUE) %>%
  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                             ifelse(grepl("C7", samples), "16", 
                                    ifelse(grepl("C8", samples), "28", NA)))) %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%  arrange(Age)

ME_negative_transcriptome_H2D <- module_eigengene_H2D_transcriptomic  %>%  dplyr::select(MEyellow, MEturquoise) %>%
  tibble::rownames_to_column(var = "samples") %>%
  gather(module, eigen_value, c(MEyellow:MEturquoise), factor_key = TRUE) %>%
  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                             ifelse(grepl("C7", samples), "16", 
                                    ifelse(grepl("C8", samples), "28", NA)))) %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%  arrange(Age)

# TO check the module
ggplot(ME_positive_transcriptome_H2D, aes(x = as.factor(Age), y = eigen_value, fill = module)) +
  geom_boxplot() +
  #ggtitle("MEblue") +
  #scale_fill_manual(values = c("red", "yellow","turquoise")) +  # Set custom colors
  labs(x = "Condition", y = "Eigen Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) 

module.gene.mapping_transcriptomics_H2D <- data.frame(cluster = as.character(bwnet_H2D_transcriptomic$colors), 
                                                      Gene = names(bwnet_H2D_transcriptomic$colors), stringsAsFactors = FALSE)

positive_H2D_transcriptome    <- module.gene.mapping_transcriptomics_H2D  %>%  filter(cluster == "blue"  | 
                                                                                      cluster == "brown"  |
                                                                                      cluster == "salmon" |
                                                                                      cluster == "midnightblue" |
                                                                                      cluster == "green")          # dim(positive_H2D_transcriptome)    6,340     

table(positive_H2D_transcriptome$cluster)

negative_H2D_transcriptome    <- module.gene.mapping_transcriptomics_H2D  %>%  filter(cluster == "turquoise"      |
                                                                                      cluster == "paleturquoise"  |
                                                                                      cluster == "skyblue3"  |
                                                                                      cluster == "yellow")        # dim(negative_H2D_transcriptome)    5193             

table(negative_H2D_transcriptome$cluster)

geneModuleMembership_H2D <- as.data.frame(cor(H2D_dds_norm, module_eigengene_H2D_transcriptomic, use='p')) 


# the code can make conflicts with The PCA and WGCNa when ask for H2D_dds_norm

# Extract the GO terms and KEGGs positive
kME.positive_H2D_blue          <-  geneModuleMembership_H2D %>% dplyr::select(., MEblue)         %>% rename(kME=MEblue)         %>% dplyr::filter(!is.na(kME))
kME.positive_H2D_brown         <-  geneModuleMembership_H2D %>% dplyr::select(., MEbrown)        %>% rename(kME=MEbrown)        %>% dplyr::filter(!is.na(kME))
kME.positive_H2D_salmon        <-  geneModuleMembership_H2D %>% dplyr::select(., MEsalmon)       %>% rename(kME=MEsalmon)       %>% dplyr::filter(!is.na(kME))
kME.positive_H2D_midnightblue  <-  geneModuleMembership_H2D %>% dplyr::select(., MEmidnightblue) %>% rename(kME=MEmidnightblue) %>% dplyr::filter(!is.na(kME))
kME.positive_H2D_green         <-  geneModuleMembership_H2D %>% dplyr::select(., MEgreen)        %>% rename(kME=MEgreen)        %>% dplyr::filter(!is.na(kME))

kME.positive_H2D_transcriptome <- rbind(kME.positive_H2D_blue, kME.positive_H2D_brown, kME.positive_H2D_salmon, kME.positive_H2D_midnightblue, kME.positive_H2D_green) %>% 
                                  dplyr::filter(., rownames(.) %in% positive_H2D_transcriptome$Gene)  %>% rownames_to_column("Gene")
# write.table(kME.positive_H2D_transcriptome, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_H2D_transcriptome.txt", sep = "\t", row.names = T, col.names = T)


List_kME.positive_H2D_transcriptome_clean <- kME.positive_H2D_transcriptome   %>% dplyr::select(Gene)                  
# write(write.table(List_kME.positive_H2D_transcritpome_clean, "D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/List_kME.positive_H2D_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.positive_H2D_transcritpome    <- left_join(List_genes, kME.positive_H2D_transcriptome, by = "Gene") %>% mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.positive_H2D_transcritpome, "D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/List_genes_kME.positive_H2D.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/")
input="List_genes_kME.positive_H2D.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

GO_terms_genes_kME.positive_H2D_transcriptome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/MWU_BP_List_genes_kME.positive_H2D.txt", sep="") %>%
                                          filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
                                          mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
                                          dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
                                          mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(GO_terms_genes_kME.positive_H2D_transcriptome, "D:/Decicomp/Paper/Paper Valdi/GO_terms_genes_kME.positive_H2D_transcriptome.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

category_counts_positive_H2D  <- GO_terms_genes_kME.positive_H2D_transcriptome %>% group_by(term_parent) %>% summarise(Count = n()) %>% arrange(desc(Count))

# GO_terms_genes_kME.positive_H2D_plot <- ggplot(GO_terms_genes_kME.positive_H2D, aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
# geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
# geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
# scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
# scale_fill_manual(values = c("cellular process"      = "tomato",  
# "metabolic process"     = "red4",  
#  "localization"          = "grey80",  
# "reproductive process"  = "grey80",
# "biological regulation" = "grey80", 
#  "biological process involved in interspecies interaction between organisms"= "red",
#  "developmental process"            = "grey80",
#  "response to stimulus"             = "orange",
#  "multicellular organismal process" = "grey80")) +
#labs(x = "Category and GO Term", y = "Frequency") +
# theme(legend.position = "none", 
#    panel.background = element_blank(),  
#    panel.grid = element_blank(),        
#    axis.title = element_blank(),     
#    axis.text.y = element_blank(),       
#    axis.ticks = element_blank(),     
#    axis.text.x = element_text(size = 14,  hjust = 1)) 
#GO_terms_genes_kME.positive_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_genes_kME.positive_H2D_plot.png", GO_terms_genes_kME.positive_H2D_plot, width = 10, height = 10, dpi = 400)

# H2D KEGGS positive transcriptome
# KEEGS_positive_H2D <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/KEGGS_positive_H2D.txt")
# KEEGS_positive_H2D <- KEEGS_positive_H2D %>%
#  dplyr::select(5, 2, 3, 8,6, 10) %>%
#  dplyr::rename(FDR = PValue,
#  Pathway = Term,
#   nGenes = Count,
#   Pathway_Genes = Pop.Hits) %>%
# mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
#        label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
# filter(FDR <= 0.05)

# KEEGS_positive_H2D_plot <- ggplot(KEEGS_positive_H2D, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
#  geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
#   geom_text(aes(label = label), hjust = -0.3, size = 5) +
#   scale_fill_gradient2(low = "#E2F0D9", mid ="#A9D18E", high ="#385700") +
#   scale_x_continuous(expand=c(0,0), limits = c(0, 6), breaks = c(0, 2,4,6)) +
#   labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
#  Style_format_theme

# KEEGS_positive_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_positive_H2D_plot.png", KEEGS_positive_H2D_plot, width = 10, height = 5, dpi = 300)

GO_terms_H2D_transcriptome_positive          <- GO_terms_genes_kME.positive_H2D_transcriptome$last_term  %>% as.vector()
matrix_GO_terms_H2D_transcriptome_positive   <- GO_similarity(GO_terms_H2D_transcriptome_positive)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_H2D_transcriptome_positive.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_H2D_transcriptome_positive, fontsize_range = c(1, 20), order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1), show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_H2D_transcriptome_positive_simply <- simplifyGO(matrix_GO_terms_H2D_transcriptome_positive, fontsize_range = c(1, 20), order_by_size = F, 
                                                          bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1),
                                                          show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_H2D_transcriptome_positive_simply <- GO_terms_H2D_transcriptome_positive_simply %>% dplyr::rename(last_term=id) %>%
                                               left_join(GO_terms_genes_kME.positive_H2D_transcriptome)


# Extract the GO terms and KEGGs negative
kME.negative_H2D_yellow         <-  geneModuleMembership_H2D %>% dplyr::select(., MEyellow)        %>% rename(kME=MEyellow)        %>% dplyr::filter(!is.na(kME))
kME.negative_H2D_skyblue3       <-  geneModuleMembership_H2D %>% dplyr::select(., MEskyblue3)      %>% rename(kME=MEskyblue3)      %>% dplyr::filter(!is.na(kME))
kME.negative_H2D_paleturquoise  <-  geneModuleMembership_H2D %>% dplyr::select(., MEpaleturquoise) %>% rename(kME=MEpaleturquoise) %>% dplyr::filter(!is.na(kME))
kME.negative_H2D_turquoise      <-  geneModuleMembership_H2D %>% dplyr::select(., MEturquoise)     %>% rename(kME=MEturquoise)     %>% dplyr::filter(!is.na(kME))

kME.negative_H2D_transcriptome <- rbind(kME.negative_H2D_yellow, kME.negative_H2D_skyblue3, kME.negative_H2D_paleturquoise, kME.negative_H2D_turquoise) %>% 
                                  dplyr::filter(., rownames(.) %in% negative_H2D_transcriptome$Gene)  %>% rownames_to_column("Gene")
# write.table(kME.negative_H2D_transcriptome, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_H2D_transcriptome.txt", sep = "\t", row.names = T, col.names = T)


List_kME.negative_H2D_transcriptome_clean <- kME.negative_H2D_transcriptome   %>% dplyr::select(Gene)                  
# write(write.table(List_kME.negative_H2D_transcriptome_clean, "D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/List_kME.negative_H2D_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.negative_H2D_transcriptome    <- left_join(List_genes, kME.negative_H2D_transcriptome, by = "Gene") %>%  
                                                mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.negative_H2D_transcriptome, "D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/List_genes_kME.negative_H2D.txt", row.names = F, quote = F, col.names=F, sep = ","))


setwd("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/")
input="List_genes_kME.negative_H2D.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

GO_terms_genes_kME.negative_H2D_transcriptome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/MWU_BP_List_genes_kME.negative_H2D.txt", sep="") %>%
                                                       filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>% filter(name != "unknown")  %>%
                                                       mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
                                                       dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
                                                       mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(GO_terms_genes_kME.negative_H2D_transcriptome, "D:/Decicomp/Paper/Paper Valdi/GO_terms_genes_kME.negative_H2D_transcriptome.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

category_counts_negative_H2D_transcriptome  <- GO_terms_genes_kME.negative_H2D_transcriptome %>% group_by(term_parent) %>% summarise(Count = n()) %>% arrange(desc(Count))

# GO_terms_genes_kME.negative_H2D_plot <- ggplot(GO_terms_genes_kME.negative_H2D, aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
  # geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
# geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
# scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
# scale_fill_manual(values = c( "cellular process"      = "lightblue", 
# "localization"          = "grey80",  
# "developmental process" = "grey80",
# "homeostatic process"   = "orange4",
# "metabolic process"     = "blue4", 
# "biological regulation" = "grey80", 
#  "positive regulation of biological process"= "grey80",
#  "reproductive process"  = "grey80",
#  "multicellular organismal process" = "grey80",
#  "response to stimulus"             = "orange",
# "immune system process"            = "darkorange",
# "biological process involved in interspecies interaction between organisms"= "blue" )) +
# labs(x = "Category and GO Term", y = "Frequency") +
# theme(legend.position = "none", 
# panel.background = element_blank(),  
# panel.grid = element_blank(),        
#  axis.title = element_blank(),     
# axis.text.y = element_blank(),       
# axis.ticks = element_blank(),     
# axis.text.x = element_text(size = 14,  hjust = 1)) 
#GO_terms_genes_kME.negative_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_genes_kME.negative_H2D_plot.png", GO_terms_genes_kME.negative_H2D_plot, width = 10, height = 10, dpi = 400)


GO_terms_H2D_transcriptome_negative         <- GO_terms_genes_kME.negative_H2D_transcriptome$last_term  %>% as.vector()
matrix_GO_terms_H2D_transcriptome_negative  <- GO_similarity(GO_terms_H2D_transcriptome_negative)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_H2D_transcriptome_negative.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_H2D_transcriptome_negative, fontsize_range = c(1, 20),order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1), show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_H2D_transcriptome_negative_simply <- simplifyGO(matrix_GO_terms_H2D_transcriptome_negative, fontsize_range = c(1, 20),order_by_size = F, 
                                                          bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1), show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_H2D_transcriptome_negative_simply <- GO_terms_H2D_transcriptome_negative_simply %>% dplyr::rename(last_term=id) %>% left_join(GO_terms_genes_kME.negative_H2D_transcriptome)



# H2D KEGGS negative transcriptome
# KEEGS_negative_H2D <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/KEGGS_negative_H2D.txt")
# KEEGS_negative_H2D <- KEEGS_negative_H2D %>%
# dplyr::select(5, 2, 3, 8,6, 10) %>%
# rename(FDR = PValue,
# Pathway = Term,
# nGenes = Count,
# Pathway_Genes = Pop.Hits) %>%
# mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
# label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
#  filter(FDR <= 0.05)

# KEEGS_negative_H2D_plot <- ggplot(KEEGS_negative_H2D, aes(x = Fold.Enrichment, y = reorder(Pathway, desc(Fold.Enrichment)), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
#  geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
#  geom_text(aes(label = label), hjust = -0.3, size = 5) +
#  scale_fill_gradient2(low = "#E2F0D9", mid ="#A9D18E", high ="#385700") +
# scale_x_continuous(expand=c(0,0), limits = c(0, 6), breaks = c(0, 2,4,6)) +
#  labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
#Style_format_theme

# KEEGS_negative_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_negative_H2D_plot.png", KEEGS_negative_H2D_plot, width = 15, height = 10, dpi = 300)

# Heat Map H2D transcriptome
Transcriptome_H2D_in_genes_list
positive_H2D_transcriptome
negative_H2D_transcriptome

H2D_list_transcriptome_associated <- rbind(positive_H2D_transcriptome, negative_H2D_transcriptome) %>% dplyr::select(Gene) %>% 
                                     left_join(Transcriptome_H2D_in_genes_list) %>% column_to_rownames(var = "Gene")

# Colum Colors heatmap   
annotation_col <- data.frame(Age  = c(rep("4 months", 6), rep("16 months", 6), rep("28 months", 6)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
                         Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"))

rownames(annotation_col) <- colnames(H2D_list_transcriptome_associated) 

# Gradient enrichment                        
myBreaks <- c(seq(min(H2D_list_transcriptome_associated), 0, length.out=ceiling(100/2) + 1),
              seq(max(H2D_list_transcriptome_associated)/100, max(H2D_list_transcriptome_associated), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue1", "black", "red1"))(100)

H2D_list_transcriptome_genes_associated_heatmap <- pheatmap(H2D_list_transcriptome_associated, 
                                                            cluster_cols =F, 
                                                            scale = "row",
                                                            cluster_rows = F, 
                                                            fontsize_row = 1, 
                                                            color = mycolor, 
                                                            # breaks = myBreaks,
                                                            border_color = "black",
                                                            clustering_distance_rows = "euclidean",
                                                            show_colnames = F, 
                                                            show_rownames = F,
                                                            annotation_col = annotation_col, 
                                                            #annotation_row = annotation_row,
                                                            annotation_colors = annotation_colors)
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/H2D_list_transcriptome_genes_associated_heatmap.tiff", 
       H2D_list_transcriptome_genes_associated_heatmap, 
       width =  5000 / 600,  # Convert to inches
       height = 5000 / 600,  # Convert to inches
       dpi = 600)

# sacar lista de genes en los GO terms.
BP_List_genes_kME.positive_H2D <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/BP_List_genes_kME.positive_H2D.txt") %>% 
  filter(name=="response to stimulus" | 
           name=="biological process involved in interspecies interaction between organisms") %>% 
  filter(value > 0) %>%  dplyr::select(seq) %>% rename(Gene=seq) %>% distinct()


BP_List_genes_kME.negative_H2D <- read.delim("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/BP_List_genes_kME.negative_H2D.txt") %>% 
  filter(name=="MyD88-dependent toll-like receptor signaling pathway" | 
           name=="toll-like receptor signaling pathway" | 
           name=="response to abiotic stimulus" |
           name=="hemokine-mediated signaling pathway" |
           name=="response to nutrient levels" |
           name=="response to osmotic stress" |
           name=="response to oxidative stress" |
           name=="response to chemical" |
           name=="detoxification") %>% 
  filter(value > 0) %>%  dplyr::select(seq) %>% rename(Gene=seq) %>% distinct()


Transcriptome_H2D_in_genes_list_counts <- assay(dds) %>% as.data.frame() %>%
  dplyr::select(contains("H2D_C3_T0"), contains("H2D_C7_T0"), contains("H2D_C8_T0")) %>%
  tibble::rownames_to_column("Gene")

Transcriptome_H2D_in_genes_list_wide <- Transcriptome_H2D_in_genes_list_counts  %>%
  pivot_longer(cols = starts_with("H2D"), names_to = "Sample", values_to = "Expression") %>% 
  mutate(Age = case_when(grepl("C3", Sample) ~ 4,
                         grepl("C7", Sample) ~ 16,
                         grepl("C8", Sample) ~ 28, TRUE ~ NA_real_ )) %>% 
  group_by(Gene, Age) %>%
  summarize(Mean_Expression = mean(Expression, na.rm = TRUE), 
            SEM = sd(Expression, na.rm = TRUE) / sqrt(n()), .groups = 'drop')


Gene_expression_I_H2D <- Transcriptome_H2D_in_genes_list_wide %>% filter(Gene=="G13495")
max_expression_value <- max(Gene_expression_I_H2D$Mean_Expression + Gene_expression_I_H2D$SEM, na.rm = TRUE)
buffered_max_value   <- max_expression_value * 1.1


Gene_expression_I_plot_H2D <- ggplot(Gene_expression_I, aes(x = factor(Age), y = Mean_Expression, fill = factor(Age))) +
  geom_line(aes(group = 1), color = "#385700", linewidth = 1) +
  geom_errorbar(aes(ymin = Mean_Expression - SEM, ymax = Mean_Expression + SEM), width = 0.1) +  
  geom_point(size = 6, shape=21) +  # Añadir los puntos
  scale_fill_manual(values = c("4" = "#E2F0D9", "16" = "#A9D18E", "28" = "#385700")) +  # Relleno personalizado
  scale_y_continuous(limits = c(0, buffered_max_value), expand = expansion(mult = c(0, 0.1))) +  
  labs(x = "Age", y = "Counts")+
  Style_format_theme  

Gene_expression_I_plot_H2D
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Gene_expression_I_plot_H2D.png", Gene_expression_I_plot_H2D, width = 5, height = 4, dpi = 600)



# # TRANSCRIPTOMICS TRANSVERSAL ANALYSIS ----

Name_GO_interest <- "response to stimulus"

# F14R positive
List_kME.positive_F14R_transcriptome_clean
dim(List_kME.positive_F14R_transcriptome_clean)

GO_terms_genes_kME.positive_F14R_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/MWU_BP_List_genes_kME.positive_F14R.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.positive_F14R_transcriptome)

category_counts_positive_F14R_transcriptome  <- GO_terms_genes_kME.positive_F14R_transcriptome %>%
  group_by(term_parent) %>% summarise(Count = n()) %>%
  arrange(desc(Count))

Genes_extract_GO_positive_F14R_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/BP_List_genes_kME.positive_F14R.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big)  %>%
  filter(term_parent== Name_GO_interest) %>% distinct(seq)

# F14R negative
List_kME.negative_F14R_clean_transcriptome
dim(List_kME.negative_F14R_clean_transcriptome)

GO_terms_genes_kME.negative_F14R_transcriptome   <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/MWU_BP_List_genes_kME.negative_F14R.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>% filter(name != "unknown")  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.negative_F14R_transcriptome)

category_counts_negative_F14R_transcriptome  <- GO_terms_genes_kME.negative_F14R_transcriptome %>%
  group_by(term_parent) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

Genes_extract_GO_negative_F14R_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/BP_List_genes_kME.negative_F14R.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big)  %>%
  filter(term_parent== Name_GO_interest) %>% distinct(seq)

### ### #### ### ### ### #### ###

# H2D positive
List_kME.positive_H2D_transcriptome_clean
dim(List_kME.positive_H2D_transcriptome_clean)

GO_terms_genes_kME.positive_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/MWU_BP_List_genes_kME.positive_H2D.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj, pval)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.positive_H2D_transcriptome)

category_counts_positive_H2D <- GO_terms_genes_kME.positive_H2D_transcriptome %>%
  group_by(term_parent) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

Genes_extract_GO_positive_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/BP_List_genes_kME.positive_H2D.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big)  %>%
  filter(term_parent== Name_GO_interest) %>% distinct(seq)


# H2D negative
List_kME.negative_H2D_transcriptome_clean
dim(List_kME.negative_H2D_transcriptome_clean)

GO_terms_genes_kME.negative_H2D_transcriptome       <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/MWU_BP_List_genes_kME.negative_H2D.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>% filter(name != "unknown")  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.negative_H2D_transcriptome)

category_counts_negative_H2D_transcriptome  <- GO_terms_genes_kME.negative_H2D_transcriptome %>%
  group_by(term_parent) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

Genes_extract_GO_negative_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/BP_List_genes_kME.negative_H2D.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big)  %>%
  filter(term_parent== Name_GO_interest) %>% distinct(seq)


Genes_extract_GO_positive_F14R_transcriptome 
dim(Genes_extract_GO_positive_F14R_transcriptome)

Genes_extract_GO_negative_F14R_transcriptome 
dim(Genes_extract_GO_negative_F14R_transcriptome)

Genes_extract_GO_positive_H2D_transcriptome 
dim(Genes_extract_GO_positive_H2D_transcriptome)

Genes_extract_GO_negative_H2D_transcriptome 
dim(Genes_extract_GO_negative_H2D_transcriptome)

intersect_GO_positive <- intersect(Genes_extract_GO_positive_F14R_transcriptome$seq, 
                                   Genes_extract_GO_positive_H2D_transcriptome$seq) %>% as.data.frame()

dim(intersect_GO_positive)



intersect_GO_negative <- intersect(Genes_extract_GO_negative_F14R_transcriptome$seq, 
                                   Genes_extract_GO_negative_H2D_transcriptome$seq) %>% as.data.frame()

dim(intersect_GO_negative)




#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --
#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --
#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --
#### -- #### -- #### -- #### -- #### -- ### -- #### -- #### -- ### -- #### -- #### --


# METABOLITES ----
# METABOLITES F14R ----
## Input files  raw data (all approaches in one shoot) 
General_metabolites <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "General_metabolites") %>% filter(!(row_number() %in% 2:15))
Nucleotides         <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "Nucleotides")         %>% filter(row_number() != 1) 
#PUFA                <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "PUFA")                %>% filter(row_number() != 1) 
TCA                 <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "TCA")                 %>% filter(!(row_number() %in% 1:3))

compounds_F14R  <- rbind(General_metabolites, Nucleotides, #PUFA, 
                         TCA) %>% dplyr::select(1, contains("F14R")) %>%  dplyr::select(-"F14R_C3.4") 

compounds_F14R <- compounds_F14R %>%  as.data.frame() %>% mutate_all(~ifelse(. == "0", "0.000001", .)) %>%
                  tibble::column_to_rownames("Compound_Method")

compounds_F14R[1,]          <-  rep(1:3, times = c(5, 6, 6))
rownames(compounds_F14R)[1] <- "Level"

write.csv(compounds_F14R, "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/compounds_F14R.csv")

setwd("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/compounds_F14R.csv", "colu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)
mSet <- SanityCheckData(mSet)
mSet <- FilterVariable(mSet, "F", 25, "iqr", 10, "mean", 0)
mSet <- PreparePrenormData(mSet)

F14R_normalization            <- Normalization(mSet, "QuantileNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
F14R_normalization_metabolite <- PlotNormSummary(F14R_normalization,         "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/Normalization_metabolite_F14R", "png", width=NA)
F14R_normalization_sample     <- PlotSampleNormSummary(F14R_normalization,   "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/Normalization_sample_F14R",    "png", width=NA)
F14R_normalized_transform     <- SaveTransformedData(F14R_normalization)

file.remove("complete_norm.qs")
file.remove("data_orig.qs")
file.remove("data_proc.qs")
file.remove("prenorm.qs")
file.remove("row_norm.qs")
file.remove("preproc.qs")
file.remove("data_original.csv")
file.remove("data_processed.csv")
file.remove("raw_dataview.csv")
file.remove("data_prefilter_iqr.csv")
file.remove("Normalization_metabolite_F14Rdpi72.png")
file.remove("Normalization_sample_F14Rdpi72.png")
file.rename("data_normalized.csv", "F14R_normalized.csv")

F14R_normalized <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_normalized.csv") %>%
                   tibble::column_to_rownames("X") %>%  dplyr::slice(-1)  %>% as.data.frame()

# dim(F14R_normalized) 76 17
F14R_normalized_list <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_normalized.csv") %>%
                        dplyr::rename(compound=X) %>%  dplyr::slice(-1)  %>% as.data.frame()

gsg <- goodSamplesGenes(t(F14R_normalized))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

F14R_normalized          <- F14R_normalized[gsg$goodGenes==TRUE,]
F14R_normalized_col_data <- colnames(F14R_normalized) %>% as.data.frame() %>% 
  dplyr::rename("sample"=".") %>%
  mutate(Family="F14R") %>%
  mutate(Age = case_when(
    grepl("C3", sample) ~ 4,
    grepl("C7", sample) ~ 16,
    grepl("C8", sample) ~ 28,
    TRUE ~ NA_real_  )) %>%
  mutate(Group = paste(Family, Age, sep = "_"))  %>%
  tibble::column_to_rownames("sample")

annotation_col     <- data.frame(Age = rep(c("F14R_4", "F14R_16", "F14R_28"), times = c(5, 6, 6)))
annotation_colors  <-  list(Age = c("F14R_4"   = "#dae3f3", "F14R_16" = "#8faadc", "F14R_28"   = "#2f5597"))

rownames(annotation_col) <- colnames(F14R_normalized) 
F14R_normalized_heatmap     <- pheatmap(F14R_normalized, 
                                        cluster_rows = T, 
                                        cluster_cols = F,
                                        border_color = F,
                                        fontsize_row = 8, 
                                        #breaks = myBreaks,
                                        color = colorRampPalette(c("black","white","#EF8944"))(50), 
                                        clustering_distance_rows = "euclidean", 
                                        #clustering_distance_cols = "euclidean",
                                        annotation_col    = annotation_col,
                                        annotation_colors = annotation_colors,
                                        show_colnames = F, 
                                        show_rownames = T,
                                        cutree_rows = NA, 
                                        cutree_cols = NA,
                                        cellwidth = 12,          # Adjust cell width
                                        #cellheight = 10 ,        # Adjust cell heigh
                                        scale = "row")
F14R_normalized_heatmap
# ggsave("D:/Decicomp/R/MOFA_omics/F14R_normalized_heatmap.png", F14R_normalized_heatmap, width = 9, height = 7, dpi = 600)


# PCA Metabolomics 
sample_pca_metaboplites_F14R <- prcomp(t(F14R_normalized))
summary(sample_pca_metaboplites_F14R)

pca_scores_metabolomics_F14R             <- data.frame(Sample = rownames(F14R_normalized_col_data), 
                                                  PC1 = sample_pca_metaboplites_F14R$x[, 1], 
                                                  PC2 = sample_pca_metaboplites_F14R$x[, 2], 
                                                  PC3 = sample_pca_metaboplites_F14R$x[, 3])

pca_scores_metabolomics_coldata_F14R    <- cbind(pca_scores_metabolomics_F14R, F14R_normalized_col_data)

PCA_F14R_metabolomics  <- ggplot(data = pca_scores_metabolomics_coldata_F14R, aes(x = PC1 , y = PC2)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = Group, color=Group), linetype = 1)  +
  geom_point(aes(fill = Group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (55.7%)", y = "PC2 (16.0%)") +
  scale_y_continuous(limits=c(-6,6)) +
  #scale_x_continuous(limits=c(-60,60)) +
  # scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "F14R_16"="#8faadc", "F14R_28"="#2f5597", "F14R_4"="#dae3f3" )) +
  scale_color_manual(values=c(  "F14R_16"="#8faadc", "F14R_28"="#2f5597", "F14R_4"="#dae3f3" )) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_F14R_metabolomics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_F14R_metabolomics.tiff", PCA_F14R_metabolomics, width = 10, height = 8, dpi = 600)

# WGCNA 
F14R_normalized_trasposed <- t(F14R_normalized) 

power <- c(c(1:10), seq(from =1, to =50, by =1 ))
sft   <- pickSoftThreshold(F14R_normalized_trasposed, 
                           powerVector = power,
                           networkType = "signed", 
                           RsquaredCut = 0.95,
                           verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "scale free topology model filt, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "Mean Connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow=2)


F14R_normalized_trasposed[] <- base::sapply(F14R_normalized_trasposed, as.numeric)
#We want a power that is above 0.8 in R^2 AND a low mean connectivity
soft_power <- 8
temp_cor   <- cor
cor        <- WGCNA::cor

bwnet_F14R  <- blockwiseModules(F14R_normalized_trasposed,
                                maxBlockSize = 25000,
                                networkType = "signed",
                                TOMType = "signed",
                                power = soft_power,
                                mergeCutHeight = 0.25,
                                minModuleSize= 5,
                                numericLabels = FALSE,
                                randomSeed=1234,
                                nThreads = 4,
                                verbose =3)
cor <- temp_cor

# Savethe network
# saveRDS(object = bwnet_F14R, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/bwnet_F14R.RDS")

bwnet_F14R <- readRDS("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/bwnet_F14R.RDS")

# Modules of the network
module_eigengenes <- bwnet_F14R$MEs
table(bwnet_F14R$colors)

Dendrogram_plot <- plotDendroAndColors(bwnet_F14R$dendrograms[[1]], 
                                       cbind(bwnet_F14R$unmergedColors, 
                                             bwnet_F14R$colors), 
                                       c("unmerged", "merged"),
                                       dendroLabels = FALSE,
                                       addGuide = F,
                                       hang = 0.03,
                                       guideHang = 0.05)

traits <- F14R_normalized_col_data %>% dplyr::select(Age) %>% 
  mutate(Age = case_when(
    Age == 4  ~ 1,
    Age == 16 ~ 2,
    Age == 28 ~ 3, TRUE ~ NA_real_))  


nSamples     <- nrow(F14R_normalized)
nCompounds   <- ncol(F14R_normalized)

module.trait.cor          <- cor(module_eigengenes, traits, use = 'p')
module.trait.cor.pvalues  <- corPvalueStudent(module.trait.cor, nSamples)

textMatrix =  paste(signif(module.trait.cor, 2), "\n(",
                    signif(module.trait.cor.pvalues, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.cor)
textMatrix 

heatmap.data              <- merge(module_eigengenes, traits, by ='row.names')
heatmap.data              <- heatmap.data %>%  column_to_rownames(var='Row.names')
heatmap.data

Corrrelations_color_permissive <-  CorLevelPlot(heatmap.data,
                                                x = names(heatmap.data)[7],
                                                y = names(heatmap.data)[1:6],
                                                col=c("blue1", "skyblue", "white", "pink", "red"))

Corrrelations_color_permissive
# tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Corrrelations_F14R_permissive.tiff", width = 6, height = 8, units = "in", res = 600)
# print(Corrrelations_color_permissive)
# dev.off()


ME_positive <- module_eigengenes  %>%  dplyr::select(c(MEblue)) %>% tibble::rownames_to_column(var = "samples")

ME_positive <- gather(ME_positive, module, eigen_value, c(MEblue), factor_key=TRUE)  %>%
               dplyr::mutate( Age= ifelse(grepl("C3",   samples), "4",
                                   ifelse(grepl("C7",   samples), "16", 
                                   ifelse(grepl("C8",   samples), "28", NA))))  %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%  arrange(Age)

ggplot(ME_positive, aes(x = as.factor(Age), y = eigen_value, fill = module)) +
  geom_boxplot() +
  #ggtitle("MEblue") +
  #scale_fill_manual(values = c("red", "yellow","turquoise")) +  # Set custom colors
  labs(x = "Condition", y = "Eigen Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) +
  geom_vline(xintercept = 2.5, color = "black", linetype = "dashed")


module.gene.mapping <- data.frame(cluster = as.character(bwnet_F14R$colors), 
                                  Gene = names(bwnet_F14R$colors), 
                                  stringsAsFactors = FALSE)

positive_F14R_metabolites    <- module.gene.mapping  %>%  filter(cluster == "turquoise")                         # dim(positive_F14R_metabolites) 21
table(positive_F14R_metabolites$clust)

negative_F14R_metabolites    <- module.gene.mapping  %>%  filter(cluster == "brown" |   cluster=="blue")         # dim(negative_F14R_metabolites) 17  
table(negative_F14R_metabolites$clust)

# Calculate the module membership and the associated p values
# This quantifies the similarity of all genes on the array to every module
module.membership.measure         <- cor(module_eigengenes, F14R_normalized_trasposed, use='p')
module.membership.measure.pvalues <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvalues <- as.data.frame(t(module.membership.measure.pvalues))

# Calculate the genes significance and associated pvalues
gene.signf.corr       <- cor(F14R_normalized_trasposed, traits$Age, use='p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

# Scaterplots 
geneModuleMembership <- as.data.frame(cor(F14R_normalized_trasposed, module_eigengenes, use='p')) 


#HEATMAP filter 
F14R_normalized_list
significant_compounds_F14R <- rbind(positive_F14R_metabolites, negative_F14R_metabolites) %>% 
                              dplyr::rename(compound=Gene)  %>% 
                              dplyr::select(compound)  %>% 
                              as.data.frame()

filtered_F14R_normalized   <- significant_compounds_F14R %>% left_join(F14R_normalized_list, by = "compound") %>% tibble::column_to_rownames("compound")

annotation_col     <- data.frame(Age = rep(c("F14R_4", "F14R_16", "F14R_28"), times = c(5, 6, 6)))

annotation_colors  <-  list(Age = c("F14R_4"   = "#dae3f3",
                                    "F14R_16"   = "#8faadc", 
                                    "F14R_28"   = "#2f5597"))

rownames(annotation_col) <- colnames(filtered_F14R_normalized) 

F14R_metabolite_plot  <- pheatmap( filtered_F14R_normalized, 
                                   cluster_rows = F, 
                                   cluster_cols = F,
                                   border_color = F,
                                   fontsize_row = 12, 
                                   #breaks = myBreaks,
                                   color = colorRampPalette(c("black","white","#EF8944"))(50), 
                                   clustering_distance_rows = "euclidean", 
                                   #clustering_distance_cols = "euclidean",
                                   annotation_col    = annotation_col,
                                   annotation_colors = annotation_colors,
                                   show_colnames = F, 
                                   show_rownames = T,
                                   cutree_rows = NA, 
                                   cutree_cols = NA,
                                   cellwidth = 12,          # Adjust cell width
                                   #cellheight = 10 ,        # Adjust cell heigh
                                   scale = "row")
F14R_metabolite_plot
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/F14R_metabolite_plot_heatmap.tiff", 
       F14R_metabolite_plot, 
       width =  5000 / 600,  # Convert to inches
       height = 5000 / 600,  # Convert to inches
       dpi = 600)


Compounds <- read_excel("D:/Decicomp/R/MOFA_omics/Data/data_base_metabolites/Conversion_Compound_to_KEGGID_GM.xlsx", sheet = "compounds_KEGG") %>% 
                        dplyr::rename(compound=Compound) 

significant_compounds_F14R_code  <- significant_compounds_F14R %>% left_join(Compounds)

graph  <- buildGraphFromKEGGREST(organism = "crg")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
                   databaseDir = tmpdir,
                   internalDir = FALSE,
                   matrices = c("diffusion"),
                   normality = c("diffusion"),
                   niter = 1e3) 

fella.data <- loadKEGGdata(databaseDir = tmpdir,
                           internalDir = FALSE,
                           loadMatrix = c("diffusion"))

# Release 109.0+/02-14, Feb 24
cat(getInfo(fella.data))

# We need the compounds 
code_keeg_compounds <- significant_compounds_F14R_code %>% filter(grepl("^C00\\d+", code)) %>% 
                       dplyr::select(code)  %>% pull(code)

analysis <- defineCompounds(compounds = code_keeg_compounds, data = fella.data)

getInput(analysis)
getExcluded(analysis)

analysis <- runDiffusion(
  object = analysis,
  data = fella.data,
  approx = "normality" ,  #"simulation" for permutation, "normality" for the parametric alternatives
  niter = 1000)           #The number of iteration if "simulation" is used

# analysis <- runPagerank(
# object = analysis,
# data = fella.data,
# approx = "simulation",  #"simulation" for permutation, "normality" for the parametric alternatives
# niter = 1000, #The number of iteration if "simulation" is used
# dampingFactor = 0.85 #it is to add some randomness i.e 1-d is the probability of choosing a random link. It is usually 0.85
# )

nlimit           <- 1000
vertex.label.cex <- 0.6

plot(analysis,
     method = "diffusion",
     data = fella.data,
     nlimit = nlimit,
     vertex.label.cex = vertex.label.cex)

exportResults(format = "csv",
              file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact.csv",
              method = "diffusion",
              object = analysis,
              data = fella.data,
              nlimit = nlimit)

nodes_conections <- generateResultsGraph(
  object = analysis,
  method = "diffusion",
  nlimit = nlimit,
  data = fella.data,
  format = "igraph")

F14R_metabo_results_compact   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact.csv")

enzyme_results                        <- F14R_metabo_results_compact %>% filter(Entry.type=="enzyme")
Keeg_metabo_F14R_results              <- F14R_metabo_results_compact %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
compound_results                      <- F14R_metabo_results_compact %>% filter(Entry.type=="compound")
module_results                        <- F14R_metabo_results_compact %>% filter(Entry.type=="module")

#Keeg_metabo_F14R_results <- Keeg_metabo_F14R_results  %>% 
#mutate(log_pvalue = -log10(Keeg_metabo_F14R_results$p.score)) %>%
#dplyr::select(!2) %>% 
#dplyr::rename(code=KEGG.id)

# Create the plot
#Keeg_metabo_F14R_results_plot <- ggplot(Keeg_metabo_F14R_results, aes(x = log_pvalue, y = reorder(KEGG.name, -log_pvalue), fill = log_pvalue)) +
  #geom_col(color="black") +
  #scale_fill_gradient(low = "white", high = "#2f5597") +
  #scale_y_discrete(expand = c(0, 0)) +
#scale_x_continuous(expand=c(0,0), limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
#labs(x = "-Log10(FDR)",   y = " ")+
  #  title = "KEGG Pathways -log10(p-value)") +
# Style_format_theme

#Keeg_metabo_F14R_results_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Keeg_metabo_F14R_results_plot.png", Keeg_metabo_F14R_results_plot, width = 12, height = 10, dpi = 300)


Keeg_metabo_F14R_results_long <- Keeg_metabo_F14R_results    %>% dplyr::select(!2) %>% 
                                 dplyr::rename(code=KEGG.id) %>% mutate(code = str_remove(code, "crg")) %>% 
                                 mutate(code = as.numeric(code)) %>% 
                                 left_join(KEGGS_Magallana_gigas)

Keeg_metabo_F14R_plot <- ggplot(Keeg_metabo_F14R_results_long, aes(x = p.score, y = reorder(KEGG.name, p.score), fill = p.score))+
  geom_point(shape=21, color="black", size= 3) +
  scale_fill_gradientn(colours= c("#2f5597")) +
  #scale_x_continuous(limits=c(0-1e-4,max(KEGG_M_R$p.score)+1e-4), expand = c(0,0)) +
  labs(x='p.score', y=' ', color='-Log (p-value)', size='Gene counts') +
  facet_nested( Global + General~., scales = "free_y", space = "free") +
  theme(
    strip.text.y = element_text(angle = 0, size = 9),  # Rotate Y facet labels
    legend.position = "none",  
    axis.line          = element_line(color="black", size= 0.5),
    axis.text          = element_text(color="black", size=10),
    #panel.background = element_rect(fill = "#f0f0f0"),  
    plot.background  = element_rect(fill = "#ffffff"),  
    strip.background = element_rect(fill = "#dae3f3", color = "black"),
    panel.background = element_blank()) 

Keeg_metabo_F14R_plot
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Keeg_metabo_F14R_plot.tiff", Keeg_metabo_F14R_plot, width = 10, height = 7, dpi = 600)



 
###### METABOLITES F14R POSITIVE ----
Compounds <- read_excel("D:/Decicomp/R/MOFA_omics/Data/data_base_metabolites/Conversion_Compound_to_KEGGID_GM.xlsx", sheet = "compounds_KEGG") %>% 
  dplyr::rename(compound=Compound) 

F14R_normalized_list 
significant_compounds_F14R_positive <- rbind(positive_F14R_metabolites) %>% 
  dplyr::rename(compound=Gene)  %>% 
  dplyr::select(compound)  %>% 
  as.data.frame()

significant_compounds_F14R_code_positive  <- significant_compounds_F14R_positive %>% left_join(Compounds)

graph  <- buildGraphFromKEGGREST(organism = "crg")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
                   databaseDir = tmpdir,
                   internalDir = FALSE,
                   matrices = c("diffusion"),
                   normality = c("diffusion"),
                   niter = 1e3) 

fella.data <- loadKEGGdata(databaseDir = tmpdir,
                           internalDir = FALSE,
                           loadMatrix = c("diffusion"))

# Release 109.0+/02-14, Feb 24
cat(getInfo(fella.data))
code_keeg_compounds <- significant_compounds_F14R_code_positive %>% filter(grepl("^C00\\d+", code)) %>% dplyr::select(code)  %>% pull(code)
analysis            <- defineCompounds(compounds = code_keeg_compounds, data = fella.data)

getInput(analysis)
getExcluded(analysis)

analysis <- runDiffusion(
  object = analysis,
  data = fella.data,
  approx = "normality" , 
  niter = 1000)          

nlimit           <- 1000
vertex.label.cex <- 0.6

plot(analysis,
     method = "diffusion",
     data = fella.data,
     nlimit = nlimit,
     vertex.label.cex = vertex.label.cex)

exportResults(format = "csv",
              file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact_positive.csv",
              method = "diffusion",
              object = analysis,
              data = fella.data,
              nlimit = nlimit)

F14R_metabo_results_compact_positive   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact_positive.csv")

enzyme_results_positive                        <- F14R_metabo_results_compact_positive %>% filter(Entry.type=="enzyme")
Keeg_metabo_F14R_results_positive              <- F14R_metabo_results_compact_positive %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
compound_results_positive                      <- F14R_metabo_results_compact_positive %>% filter(Entry.type=="compound")
module_results_positive                        <- F14R_metabo_results_compact_positive %>% filter(Entry.type=="module")

Keeg_metabo_F14R_results_positive <- Keeg_metabo_F14R_results_positive  %>% 
  mutate(log_pvalue = -log10(Keeg_metabo_F14R_results_positive$p.score)) %>%
  dplyr::select(!2) %>% 
  dplyr::rename(code=KEGG.id)

# Create the plot
Keeg_metabo_F14R_results_positive_plot <- ggplot(Keeg_metabo_F14R_results_positive, aes(x = log_pvalue, y = reorder(KEGG.name, -log_pvalue), fill = log_pvalue)) +
  geom_col(color="black") +
  scale_fill_gradient(low = "white", high = "#2f5597") +
  #scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0), limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "-Log10(FDR)",   y = " ")+
  #  title = "KEGG Pathways -log10(p-value)") +
  Style_format_theme

Keeg_metabo_F14R_results_positive_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Keeg_metabo_F14R_results_positive_plot.png", Keeg_metabo_F14R_results_positive_plot, width = 12, height = 10, dpi = 300)


###### METABOLITES F14R NEGATIVE ----
Compounds <- read_excel("D:/Decicomp/R/MOFA_omics/Data/data_base_metabolites/Conversion_Compound_to_KEGGID_GM.xlsx", sheet = "compounds_KEGG") %>% 
  dplyr::rename(compound=Compound) 

F14R_normalized_list 
significant_compounds_F14R_negative <- rbind(negative_F14R_metabolites) %>% 
  dplyr::rename(compound=Gene)  %>% 
  dplyr::select(compound)  %>% 
  as.data.frame()

significant_compounds_F14R_code_negative  <- significant_compounds_F14R_negative %>% left_join(Compounds)

graph  <- buildGraphFromKEGGREST(organism = "crg")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
                   databaseDir = tmpdir,
                   internalDir = FALSE,
                   matrices = c("diffusion"),
                   normality = c("diffusion"),
                   niter = 1e3) 

fella.data <- loadKEGGdata(databaseDir = tmpdir,
                           internalDir = FALSE,
                           loadMatrix = c("diffusion"))

# Release 109.0+/02-14, Feb 24
cat(getInfo(fella.data))
code_keeg_compounds <- significant_compounds_F14R_code_negative %>% filter(grepl("^C00\\d+", code)) %>% dplyr::select(code)  %>% pull(code)
analysis            <- defineCompounds(compounds = code_keeg_compounds, data = fella.data)

getInput(analysis)
getExcluded(analysis)

analysis <- runDiffusion(
  object = analysis,
  data = fella.data,
  approx = "normality" , 
  niter = 1000)          

nlimit           <- 1000
vertex.label.cex <- 0.6

plot(analysis,
     method = "diffusion",
     data = fella.data,
     nlimit = nlimit,
     vertex.label.cex = vertex.label.cex)

exportResults(format = "csv",
              file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact_negative.csv",
              method = "diffusion",
              object = analysis,
              data = fella.data,
              nlimit = nlimit)

F14R_metabo_results_compact_negative   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact_negative.csv")

enzyme_results_negative                        <- F14R_metabo_results_compact_negative %>% filter(Entry.type=="enzyme")
Keeg_metabo_F14R_results_negative              <- F14R_metabo_results_compact_negative %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
compound_results_negative                      <- F14R_metabo_results_compact_negative %>% filter(Entry.type=="compound")
module_results_negative                        <- F14R_metabo_results_compact_negative %>% filter(Entry.type=="module")

Keeg_metabo_F14R_results_negative <- Keeg_metabo_F14R_results_negative  %>% 
  mutate(log_pvalue = -log10(Keeg_metabo_F14R_results_negative$p.score)) %>%
  dplyr::select(!2) %>% 
  dplyr::rename(code=KEGG.id)

# Create the plot
Keeg_metabo_F14R_results_negative_plot <- ggplot(Keeg_metabo_F14R_results_negative, aes(x = log_pvalue, y = reorder(KEGG.name, -log_pvalue), fill = log_pvalue)) +
  geom_col(color="black") +
  scale_fill_gradient(low = "white", high = "#2f5597") +
  #scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0), limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "-Log10(FDR)",   y = " ")+
  #  title = "KEGG Pathways -log10(p-value)") +
  Style_format_theme

Keeg_metabo_F14R_results_negative_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Keeg_metabo_F14R_results_negative_plot.png", Keeg_metabo_F14R_results_negative_plot, width = 12, height = 10, dpi = 300)

# METABOLITES H2D ----
## Input files  raw data (all approaches in one shoot) 
General_metabolites <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "General_metabolites") %>% filter(!(row_number() %in% 2:15))
Nucleotides         <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "Nucleotides")         %>% filter(row_number() != 1) 
#PUFA                <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "PUFA")                %>% filter(row_number() != 1) 
TCA                 <- read_excel("D:/Decicomp/Matebolomics/Data_raw/Metabolites.xlsx", sheet = "TCA")                 %>% filter(!(row_number() %in% 1:3))

compounds_H2D  <- rbind(General_metabolites, Nucleotides, #PUFA, 
                        TCA) %>% 
  dplyr::select(1, contains("H2D")) %>% 
  dplyr::select(-c("H2D_C3.3", "H2D_C7.6", "H2D_C8.2"))

compounds_H2D <- compounds_H2D %>%  as.data.frame() %>% 
  mutate_all(~ifelse(. == "0", "0.000001", .)) %>%
  tibble::column_to_rownames("Compound_Method")

compounds_H2D[1,]          <-  rep(1:3, times = c(5, 5, 5))
rownames(compounds_H2D)[1] <- "Level"

write.csv(compounds_H2D, "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/compounds_H2D.csv")

setwd("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/compounds_H2D.csv", "colu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)
mSet <- SanityCheckData(mSet)
mSet <- FilterVariable(mSet, "F", 25, "iqr", 10, "mean", 0)
mSet <- PreparePrenormData(mSet)

H2D_normalization            <- Normalization(mSet, "QuantileNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
H2D_normalization_metabolite <- PlotNormSummary(H2D_normalization,         "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/Normalization_metabolite_H2D", "png", width=NA)
H2D_normalization_sample     <- PlotSampleNormSummary(H2D_normalization,   "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/Normalization_sample_H2D",    "png", width=NA)
H2D_normalized_transform     <- SaveTransformedData(H2D_normalization)

file.remove("complete_norm.qs")
file.remove("data_orig.qs")
file.remove("data_proc.qs")
file.remove("prenorm.qs")
file.remove("row_norm.qs")
file.remove("preproc.qs")
file.remove("data_original.csv")
file.remove("data_processed.csv")
file.remove("raw_dataview.csv")
file.remove("data_prefilter_iqr.csv")
file.remove("Normalization_metabolite_H2Ddpi72.png")
file.remove("Normalization_sample_H2Ddpi72.png")
file.rename("data_normalized.csv", "H2D_normalized.csv")

H2D_normalized <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_normalized.csv") %>%
  tibble::column_to_rownames("X") %>%  
  dplyr::slice(-1)  %>% as.data.frame()

H2D_normalized_list <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_normalized.csv") %>%
  dplyr::rename(compound=X) %>%  
  dplyr::slice(-1)  %>% as.data.frame()

gsg <- goodSamplesGenes(t(H2D_normalized))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

H2D_normalized          <- H2D_normalized[gsg$goodGenes==TRUE,]
H2D_normalized_col_data <- colnames(H2D_normalized) %>% as.data.frame() %>% 
  dplyr::rename("sample"=".") %>%
  mutate(Family="H2D") %>%
  mutate(Age = case_when(
    grepl("C3", sample) ~ 4,
    grepl("C7", sample) ~ 16,
    grepl("C8", sample) ~ 28,
    TRUE ~ NA_real_  )) %>%
  mutate(Group = paste(Family, Age, sep = "_"))  %>%
  tibble::column_to_rownames("sample")


annotation_col     <- data.frame(Age = rep(c("H2D_4", "H2D_16", "H2D_28"), times = c(5, 5, 5)))
annotation_colors  <- list(Age = c("H2D_4"   = "#E2F0D9", "H2D_16" = "#A9D18E", "H2D_28"   = "#385700"))

rownames(annotation_col) <- colnames(H2D_normalized) 
H2D_normalized_heatmap     <- pheatmap(H2D_normalized, 
                                       cluster_rows = F, 
                                       cluster_cols = F,
                                       border_color = F,
                                       fontsize_row = 12, 
                                       #breaks = myBreaks,
                                       color = colorRampPalette(c("black","white","#EF8944"))(50), 
                                       clustering_distance_rows = "euclidean", 
                                       #clustering_distance_cols = "euclidean",
                                       annotation_col    = annotation_col,
                                       annotation_colors = annotation_colors,
                                       show_colnames = F, 
                                       show_rownames = T,
                                       cutree_rows = NA, 
                                       cutree_cols = NA,
                                       cellwidth = 12,          # Adjust cell width
                                       #cellheight = 10 ,        # Adjust cell heigh
                                       scale = "row")
H2D_normalized_heatmap
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/H2D_metabolite_plot_heatmap.tiff", 
       H2D_normalized_heatmap, 
       width =  5000 / 600,  # Convert to inches
       height = 5000 / 600,  # Convert to inches
       dpi = 600)

# PCA Metabolomics 
sample_pca_metaboplites_H2D <- prcomp(t(H2D_normalized))
summary(sample_pca_metaboplites_H2D)

pca_scores_metabolomics_H2D             <- data.frame(Sample = rownames(H2D_normalized_col_data), 
                                                       PC1 = sample_pca_metaboplites_H2D$x[, 1], 
                                                       PC2 = sample_pca_metaboplites_H2D$x[, 2], 
                                                       PC3 = sample_pca_metaboplites_H2D$x[, 3])

pca_scores_metabolomics_coldata_H2D    <- cbind(pca_scores_metabolomics_H2D, H2D_normalized_col_data)

PCA_H2D_metabolomics  <- ggplot(data = pca_scores_metabolomics_coldata_H2D, aes(x = PC1 , y = PC2)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = Group, color=Group), linetype = 1)  +
  geom_point(aes(fill = Group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (56.1%)", y = "PC2 (13.9%)") +
  #scale_y_continuous(limits=c(-6,6)) +
  #scale_x_continuous(limits=c(-60,60)) +
  # scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "H2D_16"="#A9D18E", "H2D_28"="#385700", "H2D_4"="#E2F0D9" )) +
  scale_color_manual(values=c(  "H2D_16"="#A9D18E", "H2D_28"="#385700", "H2D_4"="#E2F0D9" )) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_H2D_metabolomics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_H2D_metabolomics.tiff", PCA_H2D_metabolomics, width = 10, height = 8, dpi = 600)


# WGCNA
H2D_normalized_trasposed <- t(H2D_normalized) 

power <- c(c(1:10), seq(from =1, to =50, by =1 ))
sft   <- pickSoftThreshold(H2D_normalized_trasposed, 
                           powerVector = power,
                           networkType = "signed", 
                           RsquaredCut = 0.95,
                           verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "scale free topology model filt, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "Mean Connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow=2)


H2D_normalized_trasposed[] <- base::sapply(H2D_normalized_trasposed, as.numeric)
#We want a power that is above 0.8 in R^2 AND a low mean connectivity
soft_power <- 22
temp_cor   <- cor
cor        <- WGCNA::cor

bwnet_H2D  <- blockwiseModules(H2D_normalized_trasposed,
                               maxBlockSize = 25000,
                               networkType = "signed",
                               TOMType = "signed",
                               power = soft_power,
                               mergeCutHeight = 0.25,
                               minModuleSize= 5,
                               numericLabels = FALSE,
                               randomSeed=1234,
                               nThreads = 4,
                               verbose =3)
cor <- temp_cor

# Savethe network
#saveRDS(object = bwnet_H2D, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/bwnet_H2D.RDS")

bwnet_H2D <- readRDS("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/bwnet_H2D.RDS")

# Modules of the network
module_eigengenes <- bwnet_H2D$MEs
table(bwnet_H2D$colors)

Dendrogram_plot <- plotDendroAndColors(bwnet_H2D$dendrograms[[1]], 
                                       cbind(bwnet_H2D$unmergedColors, 
                                             bwnet_H2D$colors), 
                                       c("unmerged", "merged"),
                                       dendroLabels = FALSE,
                                       addGuide = F,
                                       hang = 0.03,
                                       guideHang = 0.05)

traits <- H2D_normalized_col_data %>% dplyr::select(Age) %>% 
  mutate(Age = case_when(
    Age == 4  ~ 1,
    Age == 16 ~ 2,
    Age == 28 ~ 3, TRUE ~ NA_real_))  

nSamples     <- nrow(H2D_normalized)
nCompounds   <- ncol(H2D_normalized)

module.trait.cor          <- cor(module_eigengenes, traits, use = 'p')
module.trait.cor.pvalues  <- corPvalueStudent(module.trait.cor, nSamples)

textMatrix =  paste(signif(module.trait.cor, 2), "\n(",
                    signif(module.trait.cor.pvalues, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.cor)
textMatrix 

heatmap.data              <- merge(module_eigengenes, traits, by ='row.names')
heatmap.data              <- heatmap.data %>%  column_to_rownames(var='Row.names')
heatmap.data

Corrrelations_color_permissive <-  CorLevelPlot(heatmap.data,
                                                x = names(heatmap.data)[4],
                                                y = names(heatmap.data)[1:3],
                                                col=c("blue1", "skyblue", "white", "pink", "red"))

Corrrelations_color_permissive
# tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Corrrelations_H2D_permissive.tiff", width = 6, height = 8, units = "in", res = 600)
# print(Corrrelations_color_permissive)
# dev.off()


ME_positive <- module_eigengenes  %>%  dplyr::select(c(MEblue)) %>%
  tibble::rownames_to_column(var = "samples")

ME_negative <- module_eigengenes  %>%  dplyr::select(c(MEturquoise)) %>%
  tibble::rownames_to_column(var = "samples")

ME_positive <- gather(ME_positive, module, eigen_value,
                      c(MEblue), factor_key=TRUE)  %>%
  dplyr::mutate( Age= ifelse(grepl("C3",   samples), "4",
                             ifelse(grepl("C7",   samples), "16", 
                                    ifelse(grepl("C8",   samples), "28", NA))))  %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
  arrange(Age)

ME_negative<- gather(ME_negative, module, eigen_value,
                     c(MEturquoise), factor_key=TRUE)  %>%
  dplyr::mutate( Age= ifelse(grepl("C3",   samples), "4",
                             ifelse(grepl("C7",   samples), "16", 
                                    ifelse(grepl("C8",   samples), "28", NA))))  %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
  arrange(Age)


ggplot(ME_negative, aes(x = as.factor(Age), y = eigen_value, fill = module)) +
  geom_boxplot() +
  #ggtitle("MEblue") +
  #scale_fill_manual(values = c("red", "yellow","turquoise")) +  # Set custom colors
  labs(x = "Condition", y = "Eigen Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) +
  geom_vline(xintercept = 2.5, color = "black", linetype = "dashed")


module.gene.mapping <- data.frame(cluster = as.character(bwnet_H2D$colors), 
                                  Gene = names(bwnet_H2D$colors), 
                                  stringsAsFactors = FALSE)

positive_H2D_metabolites    <- module.gene.mapping  %>%  filter(cluster == "blue")        # dim(positive_H2D_metabolites)   23    
negative_H2D_metabolites    <- module.gene.mapping  %>%  filter(cluster == "turquoise")   # dim(negative_H2D_metabolites)        26  


# Calculate the module membership and the associated p values
# This quantifies the similarity of all genes on the array to every module
module.membership.measure         <- cor(module_eigengenes, H2D_normalized_trasposed, use='p')
module.membership.measure.pvalues <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvalues <- as.data.frame(t(module.membership.measure.pvalues))

# Calculate the genes significance and associated pvalues
gene.signf.corr       <- cor(H2D_normalized_trasposed, traits$Age, use='p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

# Scaterplots 
geneModuleMembership <- as.data.frame(cor(H2D_normalized_trasposed, module_eigengenes, use='p')) 


### ### ### ### ### ### ### ###
##       Enrichment analysis -
### ### ### ### ### ### ### ###

kME.turquoise <-  geneModuleMembership %>% dplyr::select(.,MEturquoise)  %>% filter(., rownames(.) %in% positive_H2D_metabolites$Gene) 
colnames(kME.turquoise)[1] <- "kME"

kME.blue <-  geneModuleMembership %>% dplyr::select(.,MEblue)  %>%  filter(., rownames(.) %in% negative_H2D_metabolites$Gene) 
colnames(kME.blue )[1] <- "kME"

kME.permissiveness <- rbind(kME.turquoise, kME.blue) #kME.green, #kME.blue)

#HEATMAP filter 
H2D_normalized_list

significant_compounds_H2D <- rbind(positive_H2D_metabolites, negative_H2D_metabolites) %>% 
  dplyr::rename(compound=Gene)  %>% 
  dplyr::select(compound)  %>% 
  as.data.frame()

filtered_H2D_normalized   <- significant_compounds_H2D %>% left_join(H2D_normalized_list, by = "compound") %>% tibble::column_to_rownames("compound")

annotation_col     <- data.frame(Age = rep(c("H2D_4", "H2D_16", "H2D_28"), times = c(5, 5, 5)))

annotation_colors  <-  list(Age = c("H2D_4"    = "#E2F0D9",   
                                    "H2D_16"   = "#A9D18E", 
                                    "H2D_28"   = "#385700"))

rownames(annotation_col) <- colnames(filtered_H2D_normalized) 

H2D_metabolite_plot  <- pheatmap( filtered_H2D_normalized, 
                                  cluster_rows = F, 
                                  cluster_cols = F,
                                  border_color = F,
                                  fontsize_row = 12, 
                                  #breaks = myBreaks,
                                  color = colorRampPalette(c("black","white","#EF8944"))(50), 
                                  clustering_distance_rows = "euclidean", 
                                  #clustering_distance_cols = "euclidean",
                                  annotation_col    = annotation_col,
                                  annotation_colors = annotation_colors,
                                  show_colnames = F, 
                                  show_rownames = T,
                                  cutree_rows = NA, 
                                  cutree_cols = NA,
                                  cellwidth = 12,          # Adjust cell width
                                  #cellheight = 10 ,        # Adjust cell heigh
                                  scale = "row")
H2D_metabolite_plot
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/H2D_metabolite_plot_heatmap.tiff", 
       H2D_metabolite_plot, 
       width =  5000 / 600,  # Convert to inches
       height = 5000 / 600,  # Convert to inches
       dpi = 600)




Compounds <- read_excel("D:/Decicomp/R/MOFA_omics/Data/data_base_metabolites/Conversion_Compound_to_KEGGID_GM.xlsx", sheet = "compounds_KEGG") %>% 
  dplyr::rename(compound=Compound) 

significant_compounds_H2D_code  <- significant_compounds_H2D %>% left_join(Compounds)

graph  <- buildGraphFromKEGGREST(organism = "crg")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
                   databaseDir = tmpdir,
                   internalDir = FALSE,
                   matrices = c("diffusion"),
                   normality = c("diffusion"),
                   niter = 1e3) 

fella.data <- loadKEGGdata(databaseDir = tmpdir,
                           internalDir = FALSE,
                           loadMatrix = c("diffusion"))

# Release 109.0+/02-14, Feb 24
cat(getInfo(fella.data))

# We need the compounds 
code_keeg_compounds <- significant_compounds_H2D_code %>% filter(grepl("^C00\\d+", code)) %>% 
                       dplyr::select(code)  %>% pull(code)

analysis <- defineCompounds(compounds = code_keeg_compounds, data = fella.data)

getInput(analysis)
getExcluded(analysis)

analysis <- runDiffusion(
  object = analysis,
  data = fella.data,
  approx = "normality" ,  #"simulation" for permutation, "normality" for the parametric alternatives
  niter = 1000)           #The number of iteration if "simulation" is used

# analysis <- runPagerank(
# object = analysis,
# data = fella.data,
# approx = "simulation",  #"simulation" for permutation, "normality" for the parametric alternatives
# niter = 1000, #The number of iteration if "simulation" is used
# dampingFactor = 0.85 #it is to add some randomness i.e 1-d is the probability of choosing a random link. It is usually 0.85
# )

nlimit           <- 1000
vertex.label.cex <- 0.6

plot(analysis,
     method = "diffusion",
     data = fella.data,
     nlimit = nlimit,
     vertex.label.cex = vertex.label.cex)

exportResults(format = "csv",
              file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact.csv",
              method = "diffusion",
              object = analysis,
              data = fella.data,
              nlimit = nlimit)


nodes_conections <- generateResultsGraph(
  object = analysis,
  method = "diffusion",
  nlimit = nlimit,
  data = fella.data,
  format = "igraph")

H2D_metabo_results_compact   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact.csv")

enzyme_results                        <- H2D_metabo_results_compact %>% filter(Entry.type=="enzyme")
Keeg_metabo_H2D_results               <- H2D_metabo_results_compact %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
compound_results                      <- H2D_metabo_results_compact %>% filter(Entry.type=="compound")
module_results                        <- H2D_metabo_results_compact %>% filter(Entry.type=="module")

#Keeg_metabo_H2D_results <- Keeg_metabo_H2D_results  %>% 
#mutate(log_pvalue = -log10(Keeg_metabo_H2D_results$p.score)) %>%
#dplyr::select(!2) %>% 
#dplyr::rename(code=KEGG.id)

# Create the plot
#Keeg_metabo_H2D_results_plot <- ggplot(Keeg_metabo_H2D_results, aes(x = log_pvalue, y = reorder(KEGG.name, -log_pvalue), fill = log_pvalue)) +
#geom_col(color="black") +
# scale_fill_gradient(low = "white", high = "#5E3C99") +
  #scale_y_discrete(expand = c(0, 0)) +
#scale_x_continuous(expand=c(0,0), limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
#labs(x = "-Log10(FDR)",   y = " ")+
  #  title = "KEGG Pathways -log10(p-value)") +
# Style_format_theme

#Keeg_metabo_H2D_results_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Keeg_metabo_H2D_results_plot.png", Keeg_metabo_H2D_results_plot, width = 12, height = 10, dpi = 300)

Keeg_metabo_H2D_results_long <- Keeg_metabo_H2D_results    %>% dplyr::select(!2) %>% 
                                dplyr::rename(code=KEGG.id) %>% mutate(code = str_remove(code, "crg")) %>% 
                                mutate(code = as.numeric(code)) %>% 
                                left_join(KEGGS_Magallana_gigas)

Keeg_metabo_H2D_plot <- ggplot(Keeg_metabo_H2D_results_long, aes(x = p.score, y = reorder(KEGG.name, p.score), fill = p.score))+
  geom_point(shape=21, color="black", size= 3) +
  scale_fill_gradientn(colours= c("#385700")) +
  #scale_x_continuous(limits=c(0-1e-4,max(KEGG_M_R$p.score)+1e-4), expand = c(0,0)) +
  labs(x='p.score', y=' ', color='-Log (p-value)', size='Gene counts') +
  facet_nested( Global + General~., scales = "free_y", space = "free") +
  theme(
    strip.text.y = element_text(angle = 0, size = 9),  # Rotate Y facet labels
    legend.position = "none",  
    axis.line          = element_line(color="black", size= 0.5),
    axis.text          = element_text(color="black", size=10),
    #panel.background = element_rect(fill = "#f0f0f0"),  
    plot.background  = element_rect(fill = "#ffffff"),  
    strip.background = element_rect(fill = "#A9D18E", color = "black"),
    panel.background = element_blank()) 

Keeg_metabo_H2D_plot
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Keeg_metabo_H2D_plot.tiff", Keeg_metabo_H2D_plot, width = 10, height = 7, dpi = 600)



###### METABOLITES H2D POSITIVE ----
Compounds <- read_excel("D:/Decicomp/R/MOFA_omics/Data/data_base_metabolites/Conversion_Compound_to_KEGGID_GM.xlsx", sheet = "compounds_KEGG") %>% 
  dplyr::rename(compound=Compound) 

H2D_normalized_list 
significant_compounds_H2D_positive <- rbind(positive_H2D_metabolites) %>% 
  dplyr::rename(compound=Gene)  %>% 
  dplyr::select(compound)  %>% 
  as.data.frame()

significant_compounds_H2D_code_positive  <- significant_compounds_H2D_positive %>% left_join(Compounds)

graph  <- buildGraphFromKEGGREST(organism = "crg")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
                   databaseDir = tmpdir,
                   internalDir = FALSE,
                   matrices = c("diffusion"),
                   normality = c("diffusion"),
                   niter = 1e3) 

fella.data <- loadKEGGdata(databaseDir = tmpdir,
                           internalDir = FALSE,
                           loadMatrix = c("diffusion"))

# Release 109.0+/02-14, Feb 24
cat(getInfo(fella.data))
code_keeg_compounds <- significant_compounds_H2D_code_positive %>% filter(grepl("^C00\\d+", code)) %>% dplyr::select(code)  %>% pull(code)
analysis            <- defineCompounds(compounds = code_keeg_compounds, data = fella.data)

getInput(analysis)
getExcluded(analysis)

analysis <- runDiffusion(
  object = analysis,
  data = fella.data,
  approx = "normality" , 
  niter = 1000)          

nlimit           <- 1000
vertex.label.cex <- 0.6

plot(analysis,
     method = "diffusion",
     data = fella.data,
     nlimit = nlimit,
     vertex.label.cex = vertex.label.cex)

exportResults(format = "csv",
              file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact_positive.csv",
              method = "diffusion",
              object = analysis,
              data = fella.data,
              nlimit = nlimit)

H2D_metabo_results_compact_positive   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact_positive.csv")

enzyme_results_positive                        <- H2D_metabo_results_compact_positive %>% filter(Entry.type=="enzyme")
Keeg_metabo_H2D_results_positive              <- H2D_metabo_results_compact_positive %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
compound_results_positive                      <- H2D_metabo_results_compact_positive %>% filter(Entry.type=="compound")
module_results_positive                        <- H2D_metabo_results_compact_positive %>% filter(Entry.type=="module")

Keeg_metabo_H2D_results_positive <- Keeg_metabo_H2D_results_positive  %>% 
  mutate(log_pvalue = -log10(Keeg_metabo_H2D_results_positive$p.score)) %>%
  dplyr::select(!2) %>% 
  dplyr::rename(code=KEGG.id)

# Create the plot
Keeg_metabo_H2D_results_positive_plot <- ggplot(Keeg_metabo_H2D_results_positive, aes(x = log_pvalue, y = reorder(KEGG.name, -log_pvalue), fill = log_pvalue)) +
  geom_col(color="black") +
  scale_fill_gradient(low = "white", high = "#385700") +
  #scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0), limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "-Log10(FDR)",   y = " ")+
  #  title = "KEGG Pathways -log10(p-value)") +
  Style_format_theme

Keeg_metabo_H2D_results_positive_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Keeg_metabo_H2D_results_positive_plot.png", Keeg_metabo_H2D_results_positive_plot, width = 12, height = 10, dpi = 300)



###### METABOLITES H2D NEGATIVE ----
Compounds <- read_excel("D:/Decicomp/R/MOFA_omics/Data/data_base_metabolites/Conversion_Compound_to_KEGGID_GM.xlsx", sheet = "compounds_KEGG") %>% 
  dplyr::rename(compound=Compound) 

H2D_normalized_list 
significant_compounds_H2D_negative <- rbind(negative_H2D_metabolites) %>% 
  dplyr::rename(compound=Gene)  %>% 
  dplyr::select(compound)  %>% 
  as.data.frame()

significant_compounds_H2D_code_negative  <- significant_compounds_H2D_negative %>% left_join(Compounds)

graph  <- buildGraphFromKEGGREST(organism = "crg")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
                   databaseDir = tmpdir,
                   internalDir = FALSE,
                   matrices = c("diffusion"),
                   normality = c("diffusion"),
                   niter = 1e3) 

fella.data <- loadKEGGdata(databaseDir = tmpdir,
                           internalDir = FALSE,
                           loadMatrix = c("diffusion"))

# Release 109.0+/02-14, Feb 24
cat(getInfo(fella.data))
code_keeg_compounds <- significant_compounds_H2D_code_negative %>% filter(grepl("^C00\\d+", code)) %>% dplyr::select(code)  %>% pull(code)
analysis            <- defineCompounds(compounds = code_keeg_compounds, data = fella.data)

getInput(analysis)
getExcluded(analysis)

analysis <- runDiffusion(
  object = analysis,
  data = fella.data,
  approx = "normality" , 
  niter = 1000)          

nlimit           <- 1000
vertex.label.cex <- 0.6

plot(analysis,
     method = "diffusion",
     data = fella.data,
     nlimit = nlimit,
     vertex.label.cex = vertex.label.cex)

exportResults(format = "csv",
              file = "D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact_negative.csv",
              method = "diffusion",
              object = analysis,
              data = fella.data,
              nlimit = nlimit)

H2D_metabo_results_compact_negative   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact_negative.csv")

enzyme_results_negative                        <- H2D_metabo_results_compact_negative %>% filter(Entry.type=="enzyme")
Keeg_metabo_H2D_results_negative               <- H2D_metabo_results_compact_negative %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
compound_results_negative                      <- H2D_metabo_results_compact_negative %>% filter(Entry.type=="compound")
module_results_negative                        <- H2D_metabo_results_compact_negative %>% filter(Entry.type=="module")

Keeg_metabo_H2D_results_negative <- Keeg_metabo_H2D_results_negative  %>% 
  mutate(log_pvalue = -log10(Keeg_metabo_H2D_results_negative$p.score)) %>%
  dplyr::select(!2) %>% 
  dplyr::rename(code=KEGG.id)

# Create the plot
Keeg_metabo_H2D_results_negative_plot <- ggplot(Keeg_metabo_H2D_results_negative, aes(x = log_pvalue, y = reorder(KEGG.name, -log_pvalue), fill = log_pvalue)) +
  geom_col(color="black") +
  scale_fill_gradient(low = "white", high = "#385700") +
  #scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0), limits = c(0, 8), breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "-Log10(FDR)",   y = " ")+
  #  title = "KEGG Pathways -log10(p-value)") +
  Style_format_theme

Keeg_metabo_H2D_results_negative_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Keeg_metabo_H2D_results_negative_plot.png", Keeg_metabo_H2D_results_negative_plot, width = 12, height = 10, dpi = 300)

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
# save.image("D:/Decicomp/R/Melthilkit projects/coverage8X/meth_F14R.RData")

###### PCA Epigenetics 
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
                                              PC3 = sample_pca_epigenetics$x[, 3])

pca_scores_epigenetic_coldata_F14R <- cbind(pca_scores_epigenetics, coldata_epigenetics) %>% 
  tibble::rownames_to_column(var = "sample")  %>%   
  dplyr::mutate(sample = gsub("_DNA", "", sample))

PCA_F14R_epigenetics  <- ggplot(data = pca_scores_epigenetic_coldata_F14R, aes(x = PC1 , y = PC2)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = group, color=group), linetype = 1)  +
  geom_point(aes(fill = group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (14.0%)", y = "PC2 (9.2%)") +
  #scale_y_continuous(limits=c(-700,700)) +
  #scale_x_continuous(limits=c(-500,500)) +
  #scale_shape_manual(values  = c(22,21, 23, 24, 25)) +
  scale_fill_manual (values=c(  "F14R_C7"="#8faadc", "F14R_C8"="#2f5597", "F14R_C3"="#dae3f3" )) +
  scale_color_manual(values=c(  "F14R_C7"="#8faadc", "F14R_C8"="#2f5597", "F14R_C3"="#dae3f3" )) +
  #geom_text_repel(aes(label = sample), size = 1.8, max.overlaps = Inf) +
  Style_format_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
PCA_F14R_epigenetics 
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/PCA_F14R_epigenetics.tiff", PCA_F14R_epigenetics, width = 10, height = 8, dpi = 600)


Coordenates_all_cpg_meth_F14R        <- data.frame(chr   = meth_F14R$chr, 
                                                   start = meth_F14R$start, 
                                                   end   = meth_F14R$end)

Methylation_meth_F14R                <- percMethylation(meth_F14R, rowids = F)
Betas_all_cpg_meth_F14R              <- cbind(Coordenates_all_cpg_meth_F14R, Methylation_meth_F14R)
setDT(Betas_all_cpg_meth_F14R,  key = c("chr", "start", "end"))

Genes_body_region                              <- read.delim("D:/Gestinnov/1) Oyster proyect/R in datarmor/Immuno_genes_list/Genes_body_region.txt")
Genes_body_region                              <- Genes_body_region[,c(4,1,2,3)]
setDT(Genes_body_region,  key = c("chr", "start", "end"))

Meth_F14R_in_genes  <- data.table::foverlaps(Betas_all_cpg_meth_F14R, Genes_body_region,  nomatch = NULL) %>% 
  dplyr::select(1,2, 7:24) %>% 
  gather(sample, methylation, 3:20, factor_key=TRUE) %>%
  mutate(sample = str_replace(sample, "_T0", ""))    %>%
  mutate(sample = str_replace(sample, "_DNA", ""))   %>%
  group_by(Gene, sample) %>%
  summarize(methylation = mean(methylation, na.rm = TRUE)) %>%
  ungroup()   %>%
  tidyr::pivot_wider(names_from = sample, values_from = methylation) %>% # dim(Meth_F14R_in_genes) 26767    
  filter(rowSums(dplyr::select(., -Gene) != 0) > 0)   %>%                # dim(Meth_F14R_in_genes) 26104     
  tibble::column_to_rownames("Gene") 

write.table(Meth_F14R_in_genes, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/Meth_F14R_in_genes_table.tsv", 
            sep = "\t",  row.names = TRUE, col.names = TRUE, quote = FALSE)

Meth_F14R_in_genes_t    <- Meth_F14R_in_genes  %>% t()
Meth_F14R_in_genes_list <- Meth_F14R_in_genes  %>% tibble::rownames_to_column("Gene") 


power <- c(c(1:10), seq(from =1, to =30, by =1 ))
sft <- pickSoftThreshold(Meth_F14R_in_genes_t, 
                         powerVector = power,
                         networkType = "signed", 
                         RsquaredCut = 0.95,
                         verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "scale free topology model filt, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "Mean Connectivity")+
  theme_classic()
grid.arrange(a1, a2, nrow=2)



Meth_F14R_in_genes_t[] <- sapply(Meth_F14R_in_genes_t, as.numeric)
#We want a power that is above 0.8 in R^2 AND a low mean connectivity
soft_power <- 25 #12 es el bueno
temp_cor   <- cor
cor        <- WGCNA::cor

bwnet_F14R_methylation      <- blockwiseModules(Meth_F14R_in_genes_t,
                                                maxBlockSize = 30000,
                                                networkType = "signed",
                                                TOMType = "signed",
                                                power = soft_power,
                                                mergeCutHeight = 0.25,
                                                minModuleSize= 100,
                                                numericLabels = FALSE,
                                                randomSeed=1234,
                                                nThreads = 4,
                                                verbose =3)

cor <- temp_cor

# Savethe network
# saveRDS(object = bwnet_F14R_methylation, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_F14R_methylation_25.RDS")

bwnet_F14R_methylation <- readRDS("D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_F14R_methylation.RDS")

module_eigengene_F14R_methylation <- bwnet_F14R_methylation$MEs
table(bwnet_F14R_methylation$colors)

Dendrogram_plot <- plotDendroAndColors(bwnet_F14R_methylation$dendrograms[[1]], 
                                       cbind(bwnet_F14R_methylation$unmergedColors, 
                                             bwnet_F14R_methylation$colors), 
                                       c("unmerged", "merged"),
                                       dendroLabels = FALSE,
                                       addGuide = F,
                                       hang = 0.03,
                                       guideHang = 0.05)
Dendrogram_plot

traits <- coldata_epigenetics %>%
  tibble::rownames_to_column("sample") %>%
  mutate(sample = str_replace(sample, "_T0", ""))    %>%
  mutate(sample = str_replace(sample, "_DNA", ""))   %>%
  mutate(age = case_when(
    grepl("C3", group) ~ 1,
    grepl("C7", group) ~ 2,
    grepl("C8", group) ~ 3,
    TRUE ~ NA_real_ )) %>%
  tibble:: column_to_rownames("sample") %>% 
  dplyr::select(age)

nSamples            <- nrow(Meth_F14R_in_genes)
nmethylates_genes   <- ncol(Meth_F14R_in_genes)

module.trait.cor          <- cor(module_eigengene_F14R_methylation, traits, use = 'p')
module.trait.cor.pvalues  <- corPvalueStudent(module.trait.cor, nSamples)

textMatrix =  paste(signif(module.trait.cor, 2), "\n(",
                    signif(module.trait.cor.pvalues, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.cor)
textMatrix 

heatmap.data              <- merge(module_eigengene_F14R_methylation, traits, by ='row.names')
heatmap.data              <- heatmap.data %>%  column_to_rownames(var='Row.names')
heatmap.data

Corrrelations_color_permissive <-  CorLevelPlot(heatmap.data,
                                                x = names(heatmap.data)[25],
                                                y = names(heatmap.data)[1:24],
                                                col=c("blue1", "skyblue", "white", "pink", "red"))
Corrrelations_color_permissive
# tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Corrrelations_F14R_permissive_methylation.tiff", width = 6, height = 8, units = "in", res = 600)
# print(Corrrelations_color_permissive)
# dev.off()


ME_positive_methyaltion_F14R <- module_eigengene_F14R_methylation  %>%  dplyr::select(MEtan) %>%
  tibble::rownames_to_column(var = "samples") %>%
  gather(module, eigen_value, c(MEtan), factor_key = TRUE) %>%
  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                             ifelse(grepl("C7", samples), "16", 
                                    ifelse(grepl("C8", samples), "28", NA)))) %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
  arrange(Age)

ME_negative_methyaltion_F14R <- module_eigengene_F14R_methylation  %>%  dplyr::select(MElightcyan) %>%
  tibble::rownames_to_column(var = "samples") %>%
  gather(module, eigen_value, c(MElightcyan), factor_key = TRUE) %>%
  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                             ifelse(grepl("C7", samples), "16", 
                                    ifelse(grepl("C8", samples), "28", NA)))) %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
  arrange(Age)


# TO check the module
ggplot(ME_negative_methyaltion_F14R, aes(x = as.factor(Age), y = eigen_value, fill = module)) +
  geom_boxplot() +
  #ggtitle("MEblue") +
  #scale_fill_manual(values = c("red", "yellow","turquoise")) +  # Set custom colors
  labs(x = "Condition", y = "Eigen Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) 

# Important
module.gene.mapping_methyaltion_F14R <- data.frame(cluster = as.character(bwnet_F14R_methylation$colors), 
                                                   Gene = names(bwnet_F14R_methylation$colors), stringsAsFactors = FALSE)

positive_F14R_methylation    <- module.gene.mapping_methyaltion_F14R  %>%  filter(cluster == "tan")                                    # dim(positive_F14R_methylation)    693       
negative_F14R_methylation    <- module.gene.mapping_methyaltion_F14R  %>%  filter(cluster == "lightcyan" | cluster == "darkred")       # dim(negative_F14R_methylation)    999    

geneModuleMembership_F14R_methylation <- as.data.frame(cor(Meth_F14R_in_genes_t, module_eigengene_F14R_methylation, use='p')) 

table(positive_F14R_methylation$cluster)
table(negative_F14R_methylation$cluster)

# write.table(positive_F14R_methylation, "D:/Decicomp/Paper/Paper Valdi/positive_F14R_methylation.tsv", row.names = F, quote = F, col.names=T, sep = "\t")
# write.table(negative_F14R_methylation, "D:/Decicomp/Paper/Paper Valdi/negative_F14R_methylation.tsv", row.names = F, quote = F, col.names=T, sep = "\t")



###### EPIGENETICS F14R POSITIVE ----
kME.positive_F14R_tan         <- geneModuleMembership_F14R_methylation %>% dplyr::select(., MEtan) %>% rename(kME=MEtan)

kME.positive_F14R_methylation <- rbind(kME.positive_F14R_tan) %>% dplyr::filter(., rownames(.) %in% positive_F14R_methylation$Gene) %>% rownames_to_column("Gene") %>% dplyr::filter(!is.na(kME))
# write.table(kME.positive_F14R_methylation, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_F14R_methylation.txt", sep = "\t", row.names = T, col.names = T)

List_kME.positive_F14R_methylation_clean  <- kME.positive_F14R_methylation %>% dplyr::select(Gene)                  
# write(write.table(List_kME.positive_F14R_methylation_clean, "D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/List_kME.positive_F14R_methylation_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.positive_F14R_methylation  <- left_join(List_genes, kME.positive_F14R_methylation, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.positive_F14R_methylation, "D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/List_genes_kME.positive_F14R_methylation.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/")
input="List_genes_kME.positive_F14R_methylation.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

David_GO_BP_kME.positive_F14R_methylation  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/David_GO_BP_List_genes_kME.positive_F14R.txt") %>%
                                               dplyr::select(2,5) %>% dplyr::rename(p.adj=PValue) %>% separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big) %>% 
                                               mutate(Count=1) %>% arrange(p.adj) %>% mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(David_GO_BP_kME.positive_F14R_methylation, "D:/Decicomp/Paper/Paper Valdi/David_GO_BP_kME.positive_F14R_methylation.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

category_counts_positive_F14R_methylation  <- David_GO_BP_kME.positive_F14R_methylation %>% group_by(term_parent) %>%
                                              summarise(Count = n()) %>%  arrange(desc(Count))

GO_terms_F14R_methylation_positive          <- David_GO_BP_kME.positive_F14R_methylation$last_term  %>% as.vector()
matrix_GO_terms_F14R_methylation_positive   <- GO_similarity(GO_terms_F14R_methylation_positive)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_F14R_methylation_positive.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_F14R_methylation_positive, fontsize_range = c(10, 30), order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1),show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_F14R_methylation_positive_simply <- simplifyGO(matrix_GO_terms_F14R_methylation_positive, fontsize_range = c(10, 30), order_by_size = F, 
                                                          bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1),
                                                          show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_F14R_methylation_positive_simply <- GO_terms_F14R_methylation_positive_simply %>% dplyr::rename(last_term=id) %>% left_join(David_GO_BP_kME.positive_F14R_methylation)





#GO_terms_methylation_kME.positive_F14R_plot <- ggplot(David_GO_BP_kME.positive_F14R, aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
#geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
#geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
#scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
#scale_fill_manual(values = c("cellular process"      = "tomato",  
# "metabolic process"     = "red4",  
# "localization"          = "grey80",  
# "reproductive process"  = "grey80",
# "biological regulation" = "grey80", 
# "biological process involved in interspecies interaction between organisms"= "red",
# "developmental process"            = "grey80",
# "response to stimulus"             = "orange",
# "immune system process"            = "orange4",
#  "multicellular organismal process" = "grey80")) +
# labs(x = "Category and GO Term", y = "Frequency") +
#theme(legend.position = "none", 
# panel.background = element_blank(),  
# panel.grid = element_blank(),        
# axis.title = element_blank(),     
#axis.text.y = element_blank(),       
# axis.ticks = element_blank(),     
# axis.text.x = element_text(size = 14,  hjust = 1)) 
#GO_terms_methylation_kME.positive_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_methylation_kME.positive_F14R_plot.png", GO_terms_methylation_kME.positive_F14R_plot, width = 10, height = 10, dpi = 600)


# F14R KEGGS positive methyaltion
#KEEGS_positive_F14R_methylation  <- read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/KEGGS_positive_methylation_F14R.txt")
#KEEGS_positive_F14R_methylation  <- KEEGS_positive_F14R_methylation %>%
# dplyr::select(5, 2, 3, 8,6, 10) %>%
#dplyr::rename(FDR = PValue,
#Pathway = Term,
#nGenes = Count,
# Pathway_Genes = Pop.Hits) %>%
#mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
#label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
# filter(FDR <= 0.05)

#KEEGS_positive_F14R_methylation_plot <- ggplot(KEEGS_positive_F14R_methylation, aes(x = Fold.Enrichment, y = Pathway, fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
# geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
# geom_text(aes(label = label), hjust = -0.3, size = 5) +
# scale_fill_gradient2(low = "#dae3f3", mid = "#8faadc", high = "#2f5597") +
# scale_x_continuous(expand=c(0,0), limits = c(0, 4), breaks = c(0, 2,4)) +
#labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
# Style_format_theme

#KEEGS_positive_F14R_methylation_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_positive_F14R_methylation_plot.png", KEEGS_positive_F14R_methylation_plot, width = 14, height = 8, dpi = 600)


###### EPIGENETICS F14R NEGATIVE ----
kME.negative_F14R_lightcyan   <-  geneModuleMembership_F14R_methylation %>% dplyr::select(., MElightcyan) %>% dplyr::rename(kME=MElightcyan)
kME.negative_F14R_darkred     <-  geneModuleMembership_F14R_methylation %>% dplyr::select(., MEdarkred)   %>% dplyr::rename(kME=MEdarkred)

kME.negative_F14R_methylation <-  rbind(kME.negative_F14R_lightcyan, kME.negative_F14R_darkred) %>% dplyr::filter(., rownames(.) %in% negative_F14R_methylation$Gene) %>% rownames_to_column("Gene") %>% dplyr::filter(!is.na(kME))
# write.table(kME.negative_F14R_methylation, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_F14R_methylation.txt", sep = "\t", row.names = T, col.names = T)

List_kME.negative_F14R_methylation_clean  <- kME.negative_F14R_methylation   %>% dplyr::select(Gene)                  
# write(write.table(List_kME.negative_F14R_methylation_clean, "D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/List_kME.negative_F14R_methylation_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.negative_F14R_methylation  <- left_join(List_genes, kME.negative_F14R_methylation, by = "Gene") %>% mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.negative_F14R_methylation, "D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/List_genes_kME.negative_F14R_methylation.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/")
input="List_genes_kME.negative_F14R_methylation.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

David_GO_BP_kME.negative_F14R_methylation  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/David_GO_BP_List_genes_kME.negative_F14R.txt") %>%
                                                dplyr::select(2,5) %>% dplyr::rename(p.adj=PValue) %>%
                                                separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big) %>% 
                                                mutate(Count=1) %>% arrange(p.adj)  %>%
                                                mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(David_GO_BP_kME.negative_F14R_methylation, "D:/Decicomp/Paper/Paper Valdi/David_GO_BP_kME.negative_F14R_methylation.tsv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

category_counts_negative_F14R_methylation  <- David_GO_BP_kME.negative_F14R_methylation %>% group_by(term_parent) %>% summarise(Count = n()) %>% arrange(desc(Count))


GO_terms_F14R_methylation_negative          <- David_GO_BP_kME.negative_F14R_methylation$last_term  %>% as.vector()
matrix_GO_terms_F14R_methylation_negative   <- GO_similarity(GO_terms_F14R_methylation_negative)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_F14R_methylation_negative.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_F14R_methylation_negative, fontsize_range = c(10, 30), order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1),show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_F14R_methylation_negative_simply <- simplifyGO(matrix_GO_terms_F14R_methylation_negative, fontsize_range = c(10, 30), order_by_size = F, 
                                                        bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1),
                                                        show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_F14R_methylation_negative_simply <- GO_terms_F14R_methylation_negative_simply %>% dplyr::rename(last_term=id) %>% left_join(David_GO_BP_kME.negative_F14R_methylation)


# GO_terms_methylation_kME.negative_F14R_plot <- ggplot(David_GO_BP_kME.negative_F14R,
                                                      # aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
#geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
#geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
#scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
# scale_fill_manual(values = c( "cellular process"      = "lightblue", 
#"localization"          = "grey80",  
# "developmental process" = "grey80",
# "metabolic process"     = "blue4", 
# "biological regulation" = "grey80", 
# "positive regulation of biological process"= "grey80",
# "reproductive process"  = "grey80",
# "multicellular organismal process" = "grey80",
# "response to stimulus"             = "orange",
# "immune system process"            = "darkorange",
# "biological process involved in interspecies interaction between organisms"= "blue" )) +
#labs(x = "Category and GO Term", y = "Frequency") +
# theme(legend.position = "none", 
# panel.background = element_blank(),  
# panel.grid = element_blank(),        
# axis.title = element_blank(),     
# axis.text.y = element_blank(),       
# axis.ticks = element_blank(),     
#  axis.text.x = element_text(size = 14,  hjust = 1)) 
#GO_terms_methylation_kME.negative_F14R_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_methylation_kME.negative_F14R_plot.png", GO_terms_methylation_kME.negative_F14R_plot, width = 10, height = 10, dpi = 400)

# F14R KEGGS negative methyaltion
#KEEGS_negative_F14R_methylation  <- read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/KEGGS_negative_methylation_F14R.txt")
#KEEGS_negative_F14R_methylation  <- KEEGS_negative_F14R_methylation %>%
#dplyr::select(5, 2, 3, 8,6, 10) %>%
# dplyr::rename(FDR = PValue,
#Pathway = Term,
# nGenes = Count,
# Pathway_Genes = Pop.Hits) %>%
#mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
#label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
# filter(FDR <= 0.05)

#KEEGS_negative_F14R_methylation_plot <- ggplot(KEEGS_negative_F14R_methylation, aes(x = Fold.Enrichment, y = Pathway, fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
#geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
#geom_text(aes(label = label), hjust = -0.3, size = 5) +
#scale_fill_gradient2(low = "#dae3f3", mid = "#8faadc", high = "#2f5597") +
#scale_x_continuous(expand=c(0,0), limits = c(0, 6), breaks = c(0, 2,4,6)) +
#labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
# Style_format_theme

#KEEGS_negative_F14R_methylation_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_negative_F14R_methylation_plot.png", KEEGS_negative_F14R_methylation_plot, width = 14, height = 8, dpi = 600)


# Heat Map F14R methylation
Meth_F14R_in_genes_list
positive_F14R_methylation
negative_F14R_methylation

F14R_list_methylated_genes_associated <- rbind(positive_F14R_methylation, negative_F14R_methylation) %>% dplyr::select(Gene) %>% 
                                         left_join(Meth_F14R_in_genes_list) %>% column_to_rownames(var = "Gene")

# Colum Colors heatmap   
annotation_col <- data.frame(
  Age  = c(rep("4 months", 6), rep("16 months", 6), rep("28 months", 6)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
  Age    = c("4 months"   = "#dae3f3", "16 months"  = "#8faadc", "28 months"  ="#2f5597"))

rownames(annotation_col) <- colnames(F14R_list_methylated_genes_associated) 


# Gradient enrichment                        
myBreaks <- c(seq(min(F14R_list_methylated_genes_associated), 0, length.out=ceiling(100/2) + 1),
              seq(max(F14R_list_methylated_genes_associated)/100, max(F14R_list_methylated_genes_associated), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

F14R_list_methylated_genes_associated_heatmap <- pheatmap(F14R_list_methylated_genes_associated, 
                                                          cluster_cols = F, 
                                                          scale = "row",
                                                          cluster_rows = F, 
                                                          fontsize_row = 1, 
                                                          color = mycolor, 
                                                          #breaks = myBreaks,
                                                          border_color = "black",
                                                          clustering_distance_rows = "euclidean",
                                                          show_colnames = F, 
                                                          show_rownames = F,
                                                          annotation_col = annotation_col, 
                                                          #annotation_row = annotation_row,
                                                          annotation_colors = annotation_colors)

# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/F14R_list_methylated_genes_associated_heatmap.tiff", 
  # F14R_list_methylated_genes_associated_heatmap, width =  5000 / 600,  # Convert to inches height = 5000 / 600, dpi = 600)


F14R_list_methylated_genes_associated_wide <- rbind(positive_F14R_methylation, negative_F14R_methylation) %>% dplyr::select(Gene) %>% 
  left_join(Meth_F14R_in_genes_list) %>%
  pivot_longer(cols = starts_with("F14R"), 
               names_to = "Sample",
               values_to = "Methylation") %>%
  mutate(Age = case_when(grepl("C3", Sample) ~ 4,
                         grepl("C7", Sample) ~ 16,
                         grepl("C8", Sample) ~ 28, TRUE ~ NA_real_ )) %>%
  group_by(Gene, Age) %>%
  summarize(Mean_Methylation = mean(Methylation, na.rm = TRUE),  # Media de la metilación
            SEM = sd(Methylation, na.rm = TRUE) / sqrt(n()),     # Error estándar de la media (SEM)
            .groups = 'drop')


# sacar lista de genes en los GO terms.
David_GO_BP_kME.positive_F14R_extract_genes  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/David_GO_BP_List_genes_kME.positive_F14R.txt") %>%   
  dplyr::rename(Gene=Genes) %>% separate(Term, into = c("last_term", "name"), sep = "~") %>% dplyr::select(3,7) %>% 
  separate_rows(Gene, sep = ", ") %>% mutate(Gene = trimws(Gene)) 


GO_temrs_interest_extract_genes <- David_GO_BP_kME.positive_F14R_extract_genes %>% 
                                   filter(name=="immune response") %>% 
                                          #name=="response to biotic stimulus" | 
                                          #name=="response to other organism"  | 
                                          #name=="response to external stimulus") %>%
                                   dplyr::select(Gene) %>% distinct()  %>% 
                                   left_join(Roseta) %>% 
                                   #left_join(Gene_conversion_Manu)

Gene_methylation_I <- F14R_list_methylated_genes_associated_wide %>% filter(Gene=="G5681")

#max_methylation_value <- max(Gene_methylation_I$Mean_Methylation + Gene_methylation_I$SEM, na.rm = TRUE)
#buffered_max_value    <- max_methylation_value * 1.1

Gene_methylation_I_plot_F14R <- ggplot(Gene_methylation_I, aes(x = factor(Age), y = Mean_Methylation, fill = factor(Age))) +
  geom_line(aes(group = 1), color = "#2f5597", linewidth = 1) +
  geom_errorbar(aes(ymin = Mean_Methylation - SEM, ymax = Mean_Methylation + SEM), width = 0.1) +  
  geom_point(size = 6, shape=21) +  # Añadir los puntos
  scale_fill_manual(values = c("4" = "#dae3f3", "16" = "#8faadc", "28" = "#2f5597")) +  # Relleno personalizado
 # scale_y_continuous(limits = c(0, buffered_max_value), expand = expansion(mult = c(0, 0.1))) +  
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) + 
  labs(x = " ", y = " ")+
  #labs(x = "Age", y = "Methylation")+
  Style_format_theme  

Gene_methylation_I_plot_F14R
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Gene_methylation_I_plot_F14R.tiff", Gene_methylation_I_plot_F14R, width = 5, height = 4, dpi = 600)

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
# save.image("D:/Decicomp/R/Melthilkit projects/coverage8X/meth_H2D.RData")


###### PCA Epigenetics 
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
                                              PC3 = sample_pca_epigenetics$x[, 3])

pca_scores_epigenetic_coldata_H2D <- cbind(pca_scores_epigenetics, coldata_epigenetics) %>% 
  tibble::rownames_to_column(var = "sample")  %>%   
  dplyr::mutate(sample = gsub("_DNA", "", sample))

PCA_H2D_epigenetics  <- ggplot(data = pca_scores_epigenetic_coldata_H2D, aes(x = PC1 , y = PC2)) +
  geom_hline(yintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, color= "grey", size = 0.25) +
  stat_ellipse(geom="polygon", alpha = 0.2, level = 0.95, size = 1, aes(fill = group, color=group), linetype = 1)  +
  geom_point(aes(fill = group ), size = 6, shape = 21, stroke = 1) +
  labs(x = "PC1 (10.8%)", y = "PC2 (9.9%)") +
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

Coordenates_all_cpg_meth_H2D        <- data.frame(chr   = meth_H2D$chr, 
                                                  start = meth_H2D$start, 
                                                  end   = meth_H2D$end)

Methylation_meth_H2D                <- percMethylation(meth_H2D, rowids = F)
Betas_all_cpg_meth_H2D              <- cbind(Coordenates_all_cpg_meth_H2D, Methylation_meth_H2D)
setDT(Betas_all_cpg_meth_H2D,  key = c("chr", "start", "end"))

Genes_body_region                              <- read.delim("D:/Gestinnov/1) Oyster proyect/R in datarmor/Immuno_genes_list/Genes_body_region.txt")
Genes_body_region                              <- Genes_body_region[,c(4,1,2,3)]
setDT(Genes_body_region,  key = c("chr", "start", "end"))

Meth_H2D_in_genes  <- data.table::foverlaps(Betas_all_cpg_meth_H2D, Genes_body_region,  nomatch = NULL) %>% 
  dplyr::select(1,2, 7:21) %>% 
  gather(sample, methylation, 3:17, factor_key=TRUE) %>%
  mutate(sample = str_replace(sample, "_T0", ""))    %>%
  mutate(sample = str_replace(sample, "_DNA", ""))   %>%
  group_by(Gene, sample) %>%
  summarize(methylation = mean(methylation, na.rm = TRUE)) %>%
  ungroup()   %>%
  tidyr::pivot_wider(names_from = sample, values_from = methylation) %>% # dim(Meth_H2D_in_genes) 26952        
  filter(rowSums(dplyr::select(., -Gene) != 0) > 0)   %>%                # dim(Meth_H2D_in_genes) 26308         
  tibble::column_to_rownames("Gene") 


write.table(Meth_H2D_in_genes, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/Meth_H2D_in_genes_table.tsv", 
            sep = "\t",  row.names = TRUE, col.names = TRUE, quote = FALSE)


Meth_H2D_in_genes_t    <- Meth_H2D_in_genes  %>% t()
Meth_H2D_in_genes_list <- Meth_H2D_in_genes  %>% tibble::rownames_to_column("Gene") 

power <- c(c(1:10), seq(from =1, to =30, by =1 ))
sft <- pickSoftThreshold(Meth_H2D_in_genes_t, 
                         powerVector = power,
                         networkType = "signed", 
                         RsquaredCut = 0.95,
                         verbose = 5)


sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "scale free topology model filt, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point()+
  geom_text(nudge_y =0.1) +
  geom_hline(yintercept = 0.8, color="red")+
  labs(x="Power", y= "Mean Connectivity")+
  theme_classic()
grid.arrange(a1, a2, nrow=2)


Meth_H2D_in_genes_t[] <- sapply(Meth_H2D_in_genes_t, as.numeric)
#We want a power that is above 0.8 in R^2 AND a low mean connectivity
soft_power <- 30 
temp_cor   <- cor
cor        <- WGCNA::cor

bwnet_H2D_methylation      <- blockwiseModules(Meth_H2D_in_genes_t,
                                               maxBlockSize = 30000,
                                               networkType = "signed",
                                               TOMType = "signed",
                                               power = soft_power,
                                               mergeCutHeight = 0.25,
                                               minModuleSize= 100,
                                               numericLabels = FALSE,
                                               randomSeed=1234,
                                               nThreads = 4,
                                               verbose =3)

cor <- temp_cor

# Savethe network
# saveRDS(object = bwnet_H2D_methylation, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_H2D_methylation_test_30.RDS")

bwnet_H2D_methylation <- readRDS("D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/bwnet_H2D_methylation_test_30.RDS")

module_eigengene_H2D_methylation <- bwnet_H2D_methylation$MEs
table(bwnet_H2D_methylation$colors)

Dendrogram_plot <- plotDendroAndColors(bwnet_H2D_methylation$dendrograms[[1]], 
                                       cbind(bwnet_H2D_methylation$unmergedColors, 
                                             bwnet_H2D_methylation$colors), 
                                       c("unmerged", "merged"),
                                       dendroLabels = FALSE,
                                       addGuide = F,
                                       hang = 0.03,
                                       guideHang = 0.05)
Dendrogram_plot


traits <- coldata_epigenetics %>%
  tibble::rownames_to_column("sample") %>%
  mutate(sample = str_replace(sample, "_T0", ""))    %>%
  mutate(sample = str_replace(sample, "_DNA", ""))   %>%
  mutate(age = case_when(
    grepl("C3", group) ~ 1,
    grepl("C7", group) ~ 2,
    grepl("C8", group) ~ 3,
    TRUE ~ NA_real_ )) %>%
  tibble:: column_to_rownames("sample") %>% 
  dplyr::select(age)

nSamples            <- nrow(Meth_H2D_in_genes)
nmethylates_genes   <- ncol(Meth_H2D_in_genes)

module.trait.cor          <- cor(module_eigengene_H2D_methylation, traits, use = 'p')
module.trait.cor.pvalues  <- corPvalueStudent(module.trait.cor, nSamples)

textMatrix =  paste(signif(module.trait.cor, 2), "\n(",
                    signif(module.trait.cor.pvalues, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.cor)
textMatrix 

heatmap.data              <- merge(module_eigengene_H2D_methylation, traits, by ='row.names')
heatmap.data              <- heatmap.data %>% column_to_rownames(var='Row.names')
heatmap.data

Corrrelations_color_permissive <-  CorLevelPlot(heatmap.data,
                                                x = names(heatmap.data)[17],
                                                y = names(heatmap.data)[1:16],
                                                col=c("blue1", "skyblue", "white", "pink", "red"))
Corrrelations_color_permissive

# tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/Corrrelations_H2D_permissive_methylation.tiff", width = 6, height = 8, units = "in", res = 600)
# print(Corrrelations_color_permissive)
# dev.off()


ME_negative_methyaltion_H2D <- module_eigengene_H2D_methylation  %>%  dplyr::select(MEgreenyellow) %>%
  tibble::rownames_to_column(var = "samples") %>%
  gather(module, eigen_value, c(MEgreenyellow), factor_key = TRUE) %>%
  dplyr::mutate(Age = ifelse(grepl("C3", samples), "4", 
                             ifelse(grepl("C7", samples), "16", 
                                    ifelse(grepl("C8", samples), "28", NA)))) %>%
  mutate(Age = factor(Age, levels = c("4", "16", "28"))) %>%
  arrange(Age)


# TO check the module
ggplot(ME_negative_methyaltion_H2D, aes(x = as.factor(Age), y = eigen_value, fill = module)) +
  geom_boxplot() +
  #ggtitle("MEblue") +
  #scale_fill_manual(values = c("red", "yellow","turquoise")) +  # Set custom colors
  labs(x = "Condition", y = "Eigen Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) 

# Important
module.gene.mapping_methyaltion_H2D <- data.frame(cluster = as.character(bwnet_H2D_methylation$colors), 
                                                  Gene = names(bwnet_H2D_methylation$colors), stringsAsFactors = FALSE)

positive_H2D_methylation    <- module.gene.mapping_methyaltion_H2D  %>%  filter(cluster == "yellow" |  cluster == "pink")         # dim(positive_H2D_methylation)   823  
table(positive_H2D_methylation$cluster)
                                                                                           
negative_H2D_methylation    <- module.gene.mapping_methyaltion_H2D  %>%  filter(cluster == "blue")                   # dim(negative_H2D_methylation)   2085          
table(negative_H2D_methylation$cluster)

geneModuleMembership_H2D <- as.data.frame(cor(Meth_H2D_in_genes_t, module_eigengene_H2D_methylation, use='p')) 

# write.table(positive_H2D_methylation, "D:/Decicomp/Paper/Paper Valdi/positive_H2D_methylation.tsv", row.names = F, quote = F, col.names=T, sep = "\t")
# write.table(negative_H2D_methylation, "D:/Decicomp/Paper/Paper Valdi/negative_H2D_methylation.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


###### EPIGENETICS H2D POSITIVE ----
kME.positive_H2D_yellow  <- geneModuleMembership_H2D %>% dplyr::select(., MEyellow) %>% dplyr::rename(kME=MEyellow)
kME.positive_H2D_pink    <- geneModuleMembership_H2D %>% dplyr::select(., MEpink)   %>% dplyr::rename(kME=MEpink)

kME.positive_H2D_methylation         <- rbind(kME.positive_H2D_yellow, kME.positive_H2D_pink) %>% dplyr::filter(., rownames(.) %in% positive_H2D_methylation$Gene) %>% rownames_to_column("Gene") %>% dplyr::filter(!is.na(kME))
# write.table(kME.positive_H2D_methylation, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_H2D_methylation.txt", sep = "\t", row.names = T, col.names = T)


kME.positive_H2D_methylation_clean  <- kME.positive_H2D_methylation   %>% dplyr::select(Gene)                  
# write(write.table(kME.positive_H2D_methylation_clean, "D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/List_kME.positive_H2D_methylation_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.positive_H2D_methylation    <- left_join(List_genes, kME.positive_H2D_methylation, by = "Gene") %>% mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.positive_H2D_methylation, "D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/List_genes_kME.positive_H2D_methylation.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/")
input="List_genes_kME.positive_H2D_methylation.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

David_GO_BP_kME.positive_H2D_methylation  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/David_GO_BP_List_genes_kME.positive_H2D.txt") %>%
                                               dplyr::select(2,5) %>% dplyr::rename(p.adj=PValue) %>%
                                               separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big) %>% 
                                               mutate(Count=1) %>% arrange(p.adj)  %>%
                                               mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(David_GO_BP_kME.positive_H2D_methylation, "D:/Decicomp/Paper/Paper Valdi/David_GO_BP_kME.positive_H2D_methylation.tsv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

category_counts_positive_H2D_methylation  <- David_GO_BP_kME.positive_H2D_methylation %>% group_by(term_parent) %>% summarise(Count = n()) %>% arrange(desc(Count))


GO_terms_H2D_methylation_positive          <- David_GO_BP_kME.positive_H2D_methylation$last_term %>% as.vector()
matrix_GO_terms_H2D_methylation_positive   <- GO_similarity(GO_terms_H2D_methylation_positive)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_H2D_methylation_positive.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_H2D_methylation_positive, fontsize_range = c(10, 30), order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1),show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_H2D_methylation_positive_simply <- simplifyGO(matrix_GO_terms_H2D_methylation_positive, fontsize_range = c(10, 30), order_by_size = F, 
                                                        bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1),
                                                        show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_H2D_methylation_positive_simply <- GO_terms_H2D_methylation_positive_simply %>% dplyr::rename(last_term=id) %>% left_join(David_GO_BP_kME.positive_H2D_methylation)


# H2D KEGGS positive transcriptome
#KEEGS_positive_H2D <- read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/KEGGS_positive_methylation_H2D.txt")
#KEEGS_positive_H2D <- KEEGS_positive_H2D %>%
#dplyr::select(5, 2, 3, 8,6, 10) %>%
#dplyr::rename(FDR = PValue,
#Pathway = Term,
#nGenes = Count,
#Pathway_Genes = Pop.Hits) %>%
#mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
#label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
#filter(FDR <= 0.05)

#KEEGS_positive_methyaltion_H2D_plot <- ggplot(KEEGS_positive_H2D, aes(x = Fold.Enrichment, y = reorder(Pathway, Fold.Enrichment), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
#geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
# geom_text(aes(label = label), hjust = -0.3, size = 5) +
#scale_fill_gradient2(low = "#E2F0D9", mid = "#A9D18E", high = "#385700") +
#scale_x_continuous(expand=c(0,0), limits = c(0, 6), breaks = c(0, 2,4,6)) +
#labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
# Style_format_theme

#KEEGS_positive_methyaltion_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_positive_methyaltion_H2D_plot.png", KEEGS_positive_methyaltion_H2D_plot, width = 14, height = 8, dpi = 600)



###### EPIGENETICS H2D NEGATIVE ----
kME.negative_H2D_blue         <-  geneModuleMembership_H2D %>% dplyr::select(., MEblue) %>% dplyr::rename(kME=MEblue)
kME.negative_H2D_methylation  <-  rbind(kME.negative_H2D_blue) %>% dplyr::filter(., rownames(.) %in% negative_H2D_methylation$Gene) %>% rownames_to_column("Gene") %>% dplyr::filter(!is.na(kME))
# write.table(kME.negative_H2D_methylation, file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_H2D_methylation.txt", sep = "\t", row.names = T, col.names = T)

List_kME.negative_H2D_methylation_clean <- kME.negative_H2D_methylation %>% dplyr::select(Gene)                  
# write(write.table(List_kME.negative_H2D_methylation_clean, "D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/List_kME.negative_H2D_methylation_clean.txt", row.names = F, quote = F, col.names=F, sep = ","))

List_genes_kME.negative_H2D_methylation    <- left_join(List_genes, kME.negative_H2D_methylation, by = "Gene") %>% mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_genes_kME.negative_H2D_methylation, "D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/List_genes_kME.negative_H2D_methylation.txt", row.names = F, quote = F, col.names=F, sep = ","))

setwd("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/")
input="List_genes_kME.negative_H2D_methylation.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

David_GO_BP_kME.negative_H2D_methylation  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/David_GO_BP_List_genes_kME.negative_H2D.txt") %>%
                                              dplyr::select(2,5) %>% dplyr::rename(p.adj=PValue) %>%
                                              separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big) %>% 
                                              mutate(Count=1) %>% arrange(p.adj)  %>%
                                              mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# write.table(David_GO_BP_kME.negative_H2D_methylation, "D:/Decicomp/Paper/Paper Valdi/David_GO_BP_kME.negative_H2D_methylation.tsv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

category_counts_negative_H2D_methylation  <- David_GO_BP_kME.negative_H2D_methylation %>% group_by(term_parent) %>% summarise(Count = n()) %>%  arrange(desc(Count))


GO_terms_H2D_methylation_negative          <- David_GO_BP_kME.positive_H2D_methylation$last_term %>% as.vector()
matrix_GO_terms_H2D_methylation_negative   <- GO_similarity(GO_terms_H2D_methylation_negative)

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/matrix_GO_terms_H2D_methylation_negative.tiff", width = 4800, height = 3600, res = 600)
simplifyGO(matrix_GO_terms_H2D_methylation_negative, fontsize_range = c(10, 30), order_by_size = T, 
           bg_gp = gpar(fill = "white", col = "black", lwd = 0.75, lty = 1),show_heatmap_legend = FALSE,  column_title = NULL)
dev.off()

GO_terms_H2D_methylation_negative_simply <- simplifyGO(matrix_GO_terms_H2D_methylation_negative, fontsize_range = c(10, 30), order_by_size = F, 
                                                       bg_gp = gpar(fill = "white", col = "white", lwd = 0.75, lty = 1),
                                                       show_heatmap_legend = FALSE,  column_title = NULL)

GO_terms_H2D_methylation_negative_simply <- GO_terms_H2D_methylation_negative_simply %>% dplyr::rename(last_term=id) %>% left_join(David_GO_BP_kME.negative_H2D_methylation)

#GO_terms_methylation_kME.negative_H2D_plot <- ggplot(GO_terms_genes_kME.negative_H2D,
#aes(axis2 = term_parent, axis1 = last_term, y = Count)) +
#geom_alluvium(aes(fill = term_parent), width = 1/12) +  geom_stratum(width = 1/12, fill = "grey95") +
#geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2,   nudge_x = -0.2, hjust = -0.3) +
#scale_x_discrete(limits = c("GO Term", "Category"), expand = c(0.2, 0.2)) +
# scale_fill_manual(values = c( "cellular process"      = "lightblue", 
# "localization"          = "grey80",  
#"developmental process" = "grey80",
#"metabolic process"     = "blue4", 
# "biological regulation" = "grey80", 
#"positive regulation of biological process"= "grey80",
# "reproductive process"  = "grey80",
# "multicellular organismal process" = "grey80",
#"response to stimulus"             = "orange",
# "immune system process"            = "darkorange",
#"biological process involved in interspecies interaction between organisms"= "blue" )) +
# labs(x = "Category and GO Term", y = "Frequency") +
# labs(x = "Category and GO Term", y = "Frequency") +
#theme(legend.position = "none", 
# panel.background = element_blank(),  
#panel.grid = element_blank(),        
# axis.title = element_blank(),     
# axis.text.y = element_blank(),       
# axis.ticks = element_blank(),     
#axis.text.x = element_text(size = 14,  hjust = 1)) 
#GO_terms_methylation_kME.negative_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_methylation_kME.negative_H2D_plot.png", GO_terms_methylation_kME.negative_H2D_plot, width = 10, height = 10, dpi = 400)

# H2D KEGGS negative transcriptome
#KEEGS_negative_H2D <- read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/KEGGS_negative_methylation_H2D.txt")
#KEEGS_negative_H2D <- KEEGS_negative_H2D %>%
#dplyr::select(5, 2, 3, 8,6, 10) %>%
#dplyr::rename(FDR = PValue,
#Pathway = Term,
#nGenes = Count,
#Pathway_Genes = Pop.Hits) %>%
# mutate(Pathway = gsub("crg\\d{5}:", "", Pathway),
#label = paste(nGenes, "/", Pathway_Genes, sep = "")) %>%
# filter(FDR <= 0.05)


#KEEGS_negative_methyaltion_H2D_plot <- ggplot(KEEGS_negative_H2D, aes(x = Fold.Enrichment, y = reorder(Pathway, Fold.Enrichment), fill = -log10(FDR))) +
  #geom_hline(yintercept = 4.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
  #geom_hline(yintercept = 1.5, linetype = "dashed", color = "#CA0020", size = 0.5) +
# geom_col(color="black") +
  #scale_fill_gradient2(high = "darkred", low = "white", midpoint = 0) + 
    #geom_text(aes(label = label), hjust = -0.3, size = 5) +
#scale_fill_gradient2(low = "#E2F0D9", mid = "#A9D18E", high = "#385700") +
# scale_x_continuous(expand=c(0,0), limits = c(0, 6), breaks = c(0, 2,4,6)) +
#labs(fill = expression(-Log[10](FDR)),  y = "", x = "Fold Enrichment" ) +
  #  title = "KEGG Pathways -log10(p-value)") +
  #annotate("text", x = 10, y = 1, label = bquote(italic(P) < 0.001), size = 9, alpha = 0.9,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 3, label = bquote(italic(P) < 0.01),  size = 9, alpha = 0.6,  color = "darkgrey", angle = 0, fontface = "bold") +
  #annotate("text", x = 10, y = 7, label = bquote(italic(P) < 0.05),  size = 9, alpha = 0.3,  color = "darkgrey", angle = 0, fontface = "bold") +
 # Style_format_theme

# KEEGS_negative_methyaltion_H2D_plot
# ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/KEEGS_negative_methyaltion_H2D_plot.png", KEEGS_negative_methyaltion_H2D_plot, width = 14, height = 8, dpi = 600)


# Heat Map H2D methylation
Meth_H2D_in_genes_list
positive_H2D_methylation 
negative_H2D_methylation 

H2D_list_methylated_genes_associated <- rbind( positive_H2D_methylation, negative_H2D_methylation) %>% dplyr::select(Gene) %>% 
  left_join(Meth_H2D_in_genes_list) %>%
  column_to_rownames(var = "Gene")

# Colum Colors heatmap   
annotation_col <- data.frame(
  Age  = c(rep("4 months", 5), rep("16 months", 6), rep("28 months", 4)))

annotation_colors  <-  list( #Family = c(F14R         = "#2f5597"),
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"))

rownames(annotation_col) <- colnames(H2D_list_methylated_genes_associated) 


# Gradient enrichment                        
myBreaks <- c(seq(min(H2D_list_methylated_genes_associated), 0, length.out=ceiling(100/2) + 1),
              seq(max(H2D_list_methylated_genes_associated)/100, max(H2D_list_methylated_genes_associated), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

H2D_list_methylated_genes_associated_heatmap <- pheatmap(H2D_list_methylated_genes_associated, 
                                                         cluster_cols = F, 
                                                         scale = "row",
                                                         cluster_rows = F, 
                                                         fontsize_row = 1, 
                                                         color = mycolor, 
                                                         #breaks = myBreaks,
                                                         border_color = "black",
                                                         clustering_distance_rows = "euclidean",
                                                         show_colnames = F, 
                                                         show_rownames = F,
                                                         annotation_col = annotation_col, 
                                                         #annotation_row = annotation_row,
                                                         annotation_colors = annotation_colors)
ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/H2D_list_methylated_genes_associated_heatmap.tiff", 
       H2D_list_methylated_genes_associated_heatmap, 
       width =  5000 / 600,  # Convert to inches
       height = 5000 / 600,  # Convert to inches
       dpi = 600)


H2D_list_methylated_genes_associated_wide <- rbind(positive_H2D_methylation, negative_H2D_methylation) %>% dplyr::select(Gene) %>% 
  left_join(Meth_H2D_in_genes_list) %>%
  pivot_longer(cols = starts_with("H2D"), 
               names_to = "Sample",
               values_to = "Methylation") %>%
  mutate(Age = case_when(grepl("C3", Sample) ~ 4,
                         grepl("C7", Sample) ~ 16,
                         grepl("C8", Sample) ~ 28, TRUE ~ NA_real_ )) %>%
  group_by(Gene, Age) %>%
  summarize(Mean_Methylation = mean(Methylation, na.rm = TRUE),  # Media de la metilación
            SEM = sd(Methylation, na.rm = TRUE) / sqrt(n()),     # Error estándar de la media (SEM)
            .groups = 'drop')

BP_List_genes_kME.negative_H2D <- read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/BP_List_genes_kME.negative_H2D_methylation.txt") %>% 
  filter(name=="response to virus") %>% 
  filter(value > 0) %>%  dplyr::select(seq) %>% dplyr::rename(Gene=seq) %>% distinct()

Gene_methylation_I <- H2D_list_methylated_genes_associated_wide %>% filter(Gene=="G26336")

max_methylation_value <- max(Gene_methylation_I$Mean_Methylation + Gene_methylation_I$SEM, na.rm = TRUE)
buffered_max_value    <- max_methylation_value * 1.1

Gene_methylation_I_plot_H2D <- ggplot(Gene_methylation_I, aes(x = factor(Age), y = Mean_Methylation, fill = factor(Age))) +
  geom_line(aes(group = 1), color = "#385700", linewidth = 1) +
  geom_errorbar(aes(ymin = Mean_Methylation - SEM, ymax = Mean_Methylation + SEM), width = 0.1) +  
  geom_point(size = 6, shape=21) +  # Añadir los puntos
  scale_fill_manual(values = c("4" = "#E2F0D9", "16" =  "#A9D18E",  "28" ="#385700")) +  # Relleno personalizado
  scale_y_continuous(limits = c(0, buffered_max_value), expand = expansion(mult = c(0, 0.1))) +  
  labs(x = "Age", y = "Methylation")+
  Style_format_theme  

Gene_methylation_I_plot_H2D
#ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Gene_methylation_I_plot_H2D.png", Gene_methylation_I_plot_H2D, width = 5, height = 4, dpi = 600)


## INTEGRATION ----

# Download the pathways of Magalla gigas
crg_pathways            <- downloadPathways("crg") #cached version from 2024-04-20
crg_pathways_individual <- getGenesets(org = "crg", db = "kegg", cache = TRUE, return.type="list") #cached version from 2024-04-25

pathway_dfs <- list()
for (i in seq_along(crg_pathways_individual)) {
  kegg_pathway     <- names(crg_pathways_individual)[i]
  pathway_genes    <- crg_pathways_individual[[i]]
  pathway_df       <- data.frame(Gene_LOC = pathway_genes, kegg_pathway = kegg_pathway)
  pathway_df$code  <- substring(kegg_pathway, 1, 8)  # Corrected code assignment
  pathway_dfs[[i]] <- pathway_df}

All_data_kegg_pathways_genes <- bind_rows(pathway_dfs) %>%  mutate(kegg_pathway = substr(kegg_pathway, 9, nchar(kegg_pathway)),  # Remove first 8 characters
                                                                   kegg_pathway = gsub("_", " ", kegg_pathway)) %>%  
                                                                   # %>% left_join(Genes_conversion, by="Gene_LOC")
                                                                   mutate(Gene_LOC = str_replace(as.character(Gene_LOC), "^(\\d+)", "LOC\\1")) %>%  
                                                                   left_join(Gene_conversion_Manu) %>%
                                                                   filter(!is.na(roslin))  %>%
                                                                   separate_rows(roslin, sep = ",\\s*") %>%
                                                                   dplyr::rename("Gene"="roslin")


######  INTEGRATION F14R (Recuperamos la lista de genes que hemos conseguido antes para cada -omic) ---- 
## Epigenetics and Transcriptome 

positive_F14R_methylation     <- read.table("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/List_kME.positive_F14R_methylation_clean.txt", quote="\"", comment.char="") %>%  dplyr::rename ("Gene"="V1") %>% mutate(Tendency="positive", Family="F14R")
negative_F14R_methylation     <- read.table("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/List_kME.negative_F14R_methylation_clean.txt", quote="\"", comment.char="") %>%  dplyr::rename ("Gene"="V1") %>% mutate(Tendency="negative", Family="F14R")
F14R_methylation_integrate    <- rbind(positive_F14R_methylation, negative_F14R_methylation)  %>% mutate(Omic="Methylation")

positive_F14R_transcriptome   <- read.table("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/List_kME.positive_F14R_clean.txt", quote="\"", comment.char="") %>% dplyr::rename ("Gene"="V1") %>% mutate(Tendency="positive", Family="F14R")
negative_F14R_transcriptome   <- read.table("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/List_kME.negative_F14R_clean.txt", quote="\"", comment.char="") %>% dplyr::rename ("Gene"="V1") %>% mutate(Tendency="negative", Family="F14R")
F14R_transcriptome_integrate  <- rbind(positive_F14R_transcriptome, negative_F14R_transcriptome)  %>% mutate(Omic="Transcriptome")

F14R_metabo_results_compact   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/F14R_metabo_results_compact.csv") %>%
                                 filter(p.score < 0.05) %>% filter(Entry.type=="pathway") %>% mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) %>% dplyr::select(!2) %>% 
                                 dplyr::rename(code=KEGG.id)

F14R_metabo_results_compact_KEGG_data <- F14R_metabo_results_compact %>% left_join(All_data_kegg_pathways_genes) %>% dplyr::select(1,2,3,4,7)

F14R_methylation_integrate_to_metabo         <-  F14R_methylation_integrate %>%  merge(F14R_metabo_results_compact_KEGG_data)
F14R_methylation_integrate_to_metabo_table   <-  table(F14R_methylation_integrate_to_metabo$KEGG.name, F14R_methylation_integrate_to_metabo$Tendency) %>% as.data.frame() %>% 
                                                 dplyr::rename("KEGG.name"="Var1", "Methylation"="Var2", "freq_methylation"="Freq") %>%
                                                 merge(F14R_metabo_results_compact)

F14R_transcriptome_integrate_to_metabo       <- F14R_transcriptome_integrate %>%  merge(F14R_metabo_results_compact_KEGG_data)
F14R_transcriptome_integrate_to_metabo_table <- table(F14R_transcriptome_integrate_to_metabo$KEGG.name, F14R_transcriptome_integrate_to_metabo$Tendency) %>% as.data.frame() %>% 
                                                dplyr::rename("KEGG.name"="Var1", "Transcriptome"="Var2", "freq_transcriptome"="Freq") %>%
                                                merge(F14R_metabo_results_compact)


all_KEGG <- union(F14R_methylation_integrate_to_metabo_table$KEGG.name, 
                  F14R_transcriptome_integrate_to_metabo_table$KEGG.name)

meth_complete <- full_join(data.frame(KEGG.name = all_KEGG), 
                           F14R_methylation_integrate_to_metabo_table, 
                           by = "KEGG.name")  %>%
  bind_rows(data.frame(  KEGG.name = "Pyrimidine metabolism",
      freq_methylation = 0,
      code = "crg00240",
      p.score = 3.960853e-05,
      Methylation = c("negative", "positive"))) %>% 
  filter(!(KEGG.name == "Pyrimidine metabolism" & is.na(freq_methylation)))

  
trans_complete <- full_join(data.frame(KEGG.name = all_KEGG), 
                            F14R_transcriptome_integrate_to_metabo_table, 
                            by = "KEGG.name")

Integration_F14R_metabo_with_trans_meth <- bind_cols(meth_complete, trans_complete %>% dplyr::select(-KEGG.name)) %>% dplyr::select(1,2,3,6,7,5) %>% dplyr::rename("p.score"="p.score...5")



# Associate_age_trasncriptome_F14R             <- rbind(positive_F14R_transcriptome_tendency, negative_F14R_transcriptome_tendency)
# Associate_age_methylation_F14R               <- rbind(positive_F14R_methylation_tendency,   negative_F14R_methylation_tendency)  
#Intersect_Age_F14R_transcriptomic_methylation <- intersect(Associate_age_trasncriptome_F14R$Gene, Associate_age_methylation_F14R$Gene) %>% as.data.frame() %>% dplyr::rename(Gene=".")

# List_tendency_transcriptome_Genes_F14R   <- Associate_age_trasncriptome_F14R %>% rename(Transcriptome = Tendency)
# List_tendency_methylation_Genes_F14R     <- Associate_age_methylation_F14R  %>% rename(Methylation = Tendency)
# Omics_F14R_tendency <- merge(List_tendency_transcriptome_Genes_F14R, List_tendency_methylation_Genes_F14R, by = "Gene", all = TRUE) #%>% left_join(Genes_conversion)




#Omics_F14R_tendency_with_metabolites <- Omics_F14R_tendency %>% left_join(Keeg_metabo_F14R_results_KEGG_database) %>% filter %>%
# filter(!is.na(Global))

#Omics_F14R_tendency_with_metabolites_sumarize <- Omics_F14R_tendency_with_metabolites %>%
#group_by(Global, Genereal, kegg_pathway) %>%
# summarize( positive_transcriptome = sum(Transcriptome == "positive", na.rm = TRUE),  
# negative_transcriptome = sum(Transcriptome == "negative", na.rm = TRUE),
                                                 #
#positive_methylation = sum(Methylation == "positive", na.rm = TRUE),
# negative_methylation = sum(Methylation == "negative", na.rm = TRUE),
# total_genes = n()) %>% arrange(desc(total_genes))
# write.table(Omics_F14R_tendency_with_metabolites_sumarize, "D:/Decicomp/R/MOFA_omics/Data_compacted/Omics_F14R_tendency_with_metabolites_sumarize.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 

# Omics_F14R_transcriptome_positive_Genes <- Omics_F14R_tendency_with_metabolites %>% filter(Transcriptome=="positive") %>% distinct(Gene) # 535
# Omics_F14R_transcriptome_negative_Genes <- Omics_F14R_tendency_with_metabolites %>% filter(Transcriptome=="negative") %>% distinct(Gene) # 325

# Omics_F14R_methylation_positive_Genes <- Omics_F14R_tendency_with_metabolites %>% filter(Methylation=="positive") %>% distinct(Gene) # 34
# Omics_F14R_methylation_negative_Genes <- Omics_F14R_tendency_with_metabolites %>% filter(Methylation=="negative") %>% distinct(Gene) # 75

# Needs to load  the modules.membership of Kmer for trascriptome and methylation
# kME.positive_F14R_transcriptome <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_F14R_transcriptome.txt", sep = "\t", header = TRUE, row.names = 1)
#kME.negative_F14R_transcriptome <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_F14R_transcriptome.txt", sep = "\t", header = TRUE, row.names = 1)

#kME.positive_F14R_methylation   <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_F14R_methylation.txt", sep = "\t", header = TRUE, row.names = 1)
#kME.negative_F14R_methylation   <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_F14R_methylation.txt", sep = "\t", header = TRUE, row.names = 1)

# --

#Omics_F14R_transcriptome_positive_Genes_kmer        <- Omics_F14R_transcriptome_positive_Genes %>% left_join(kME.positive_F14R_transcriptome)
#List_Omics_F14R_transcriptome_positive_Genes_kmer   <- left_join(List_genes, Omics_F14R_transcriptome_positive_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_F14R_transcriptome_positive_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/positive/List_Omics_F14R_transcriptome_positive_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

#setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/positive/")
#input="List_Omics_F14R_transcriptome_positive_Genes_kmer.txt" 
#goAnnotations="all_go.tab" 
#goDatabase="go.obo" 
#goDivision="BP" 
#source("gomwu.functions.R")
#gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

#GO_terms_Omics_F14R_transcriptome_positive_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/positive/MWU_BP_List_Omics_F14R_transcriptome_positive_Genes_kmer.txt", sep="") %>%
#filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
#mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
#dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
# mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

#category_counts_GO_terms_Omics_F14R_transcriptome_positive_Genes_kmer  <- GO_terms_Omics_F14R_transcriptome_positive_Genes_kmer %>% group_by(term_parent) %>% 
# summarise(Count = n()) %>% arrange(desc(Count))

# --


#Omics_F14R_transcriptome_negative_Genes_kmer        <- Omics_F14R_transcriptome_negative_Genes %>% left_join(kME.negative_F14R_transcriptome)
#List_Omics_F14R_transcriptome_negative_Genes_kmer   <- left_join(List_genes, Omics_F14R_transcriptome_negative_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_F14R_transcriptome_negative_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/negative/List_Omics_F14R_transcriptome_negative_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

#setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/negative/")
#input="List_Omics_F14R_transcriptome_negative_Genes_kmer.txt" 
#goAnnotations="all_go.tab" 
#goDatabase="go.obo" 
#goDivision="BP" 
#source("gomwu.functions.R")
#gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

#GO_terms_Omics_F14R_transcriptome_negative_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/negative/MWU_BP_List_Omics_F14R_transcriptome_negative_Genes_kmer.txt", sep="") %>%
# filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
#mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
# dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
# mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

#category_counts_GO_terms_Omics_F14R_transcriptome_negative_Genes_kmer  <- GO_terms_Omics_F14R_transcriptome_negative_Genes_kmer %>% group_by(term_parent) %>% 
#summarise(Count = n()) %>% arrange(desc(Count))

# --

#Omics_F14R_methylation_positive_Genes_kmer      <- Omics_F14R_methylation_positive_Genes %>% left_join(kME.positive_F14R_methylation)
#List_Omics_F14R_methylation_positive_Genes_kmer <- left_join(List_genes, Omics_F14R_methylation_positive_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_F14R_methylation_positive_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/positive/List_Omics_F14R_methylation_positive_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

#setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/positive/")
#input="List_Omics_F14R_methylation_positive_Genes_kmer.txt" 
#goAnnotations="all_go.tab" 
#goDatabase="go.obo" 
#goDivision="BP" 
#source("gomwu.functions.R")
#gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

#GO_terms_Omics_F14R_methylation_positive_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/positive/MWU_BP_List_Omics_F14R_methylation_positive_Genes_kmer.txt", sep="") %>%
# filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
# mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
# dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
# mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

#category_counts_GO_terms_Omics_F14R_methylation_positive_Genes_kmer  <- GO_terms_Omics_F14R_methylation_positive_Genes_kmer %>% group_by(term_parent) %>% 
#summarise(Count = n()) %>% arrange(desc(Count))

# --

#Omics_F14R_methylation_negative_Genes_kmer      <- Omics_F14R_methylation_negative_Genes %>% left_join(kME.negative_F14R_methylation)
#List_Omics_F14R_methylation_negative_Genes_kmer <- left_join(List_genes, Omics_F14R_methylation_negative_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_F14R_methylation_negative_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/negative/List_Omics_F14R_methylation_negative_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

#setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/negative/")
#input="List_Omics_F14R_methylation_negative_Genes_kmer.txt" 
#goAnnotations="all_go.tab" 
#goDatabase="go.obo" 
#goDivision="BP" 
#source("gomwu.functions.R")
#gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

#GO_terms_Omics_F14R_methylation_negative_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/F14R/negative/MWU_BP_List_Omics_F14R_methylation_negative_Genes_kmer.txt", sep="") %>%
# filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
# mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
# dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
#  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

#category_counts_GO_terms_Omics_F14R_methylation_negative_Genes_kmer  <- GO_terms_Omics_F14R_methylation_negative_Genes_kmer %>% group_by(term_parent) %>% 
#summarise(Count = n()) %>% arrange(desc(Count))


######  INTEGRATION H2D ----
positive_H2D_methylation     <- read.table("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/List_kME.positive_H2D_methylation_clean.txt", quote="\"", comment.char="") %>%  dplyr::rename ("Gene"="V1") %>% mutate(Tendency="positive", Family="H2D")
negative_H2D_methylation     <- read.table("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/List_kME.negative_H2D_methylation_clean.txt", quote="\"", comment.char="") %>%  dplyr::rename ("Gene"="V1") %>% mutate(Tendency="negative", Family="H2D")
H2D_methylation_integrate    <- rbind(positive_H2D_methylation, negative_H2D_methylation)  %>% mutate(Omic="Methylation")

positive_H2D_transcriptome   <- read.table("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/List_kME.positive_H2D_clean.txt", quote="\"", comment.char="") %>% dplyr::rename ("Gene"="V1") %>% mutate(Tendency="positive", Family="H2D")
negative_H2D_transcriptome   <- read.table("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/List_kME.negative_H2D_clean.txt", quote="\"", comment.char="") %>% dplyr::rename ("Gene"="V1") %>% mutate(Tendency="negative", Family="H2D")
H2D_transcriptome_integrate  <- rbind(positive_H2D_transcriptome, negative_H2D_transcriptome)  %>% mutate(Omic="Transcriptome")

H2D_metabo_results_compact   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact.csv") %>%
                                filter(p.score < 0.05) %>% filter(Entry.type=="pathway") %>% mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) %>% dplyr::select(!2) %>% 
                                dplyr::rename(code=KEGG.id)

H2D_metabo_results_compact_KEGG_data <- H2D_metabo_results_compact %>% left_join(All_data_kegg_pathways_genes) %>% dplyr::select(1,2,3,4,7)

H2D_methylation_integrate_to_metabo         <-  H2D_methylation_integrate %>%  merge(H2D_metabo_results_compact_KEGG_data)
H2D_methylation_integrate_to_metabo_table   <-  table(H2D_methylation_integrate_to_metabo$KEGG.name, H2D_methylation_integrate_to_metabo$Tendency) %>% as.data.frame() %>% 
                                                dplyr::rename("KEGG.name"="Var1", "Methylation"="Var2", "freq_methylation"="Freq") %>%
                                                merge(H2D_metabo_results_compact)

H2D_transcriptome_integrate_to_metabo       <- H2D_transcriptome_integrate %>%  merge(H2D_metabo_results_compact_KEGG_data)
H2D_transcriptome_integrate_to_metabo_table <- table(H2D_transcriptome_integrate_to_metabo$KEGG.name, H2D_transcriptome_integrate_to_metabo$Tendency) %>% as.data.frame() %>% 
                                               dplyr::rename("KEGG.name"="Var1", "Transcriptome"="Var2", "freq_transcriptome"="Freq") %>%
                                               merge(H2D_metabo_results_compact)

all_KEGG <- union(H2D_methylation_integrate_to_metabo_table$KEGG.name, 
                  H2D_transcriptome_integrate_to_metabo_table$KEGG.name)

meth_complete <- full_join(data.frame(KEGG.name = all_KEGG), 
                           H2D_methylation_integrate_to_metabo_table, 
                           by = "KEGG.name")  %>%
  bind_rows(data.frame(  KEGG.name = "Pyrimidine metabolism",
                         freq_methylation = 0,
                         code = "crg00240",
                         p.score = 3.960853e-05,
                         Methylation = c("negative", "positive"))) %>% 
  filter(!(KEGG.name == "Pyrimidine metabolism" & is.na(freq_methylation)))


trans_complete <- full_join(data.frame(KEGG.name = all_KEGG), 
                            H2D_transcriptome_integrate_to_metabo_table, 
                            by = "KEGG.name")

Integration_H2D_metabo_with_trans_meth <- bind_cols(meth_complete, trans_complete %>% dplyr::select(-KEGG.name)) %>% dplyr::select(1,2,3,6,7,5) %>% dplyr::rename("p.score"="p.score...5")



#Associate_age_trasncriptome_H2D <- rbind(positive_H2D_transcriptome_tendency, negative_H2D_transcriptome_tendency)
#Associate_age_methylation_H2D   <- rbind(positive_H2D_methylation_tendency,   negative_H2D_methylation_tendency)  
# Intersect_Age_H2D_transcriptomic_methylation <- intersect(Associate_age_trasncriptome_H2D$Gene, Associate_age_methylation_H2D$Gene) %>% as.data.frame() %>% dplyr::rename(Gene=".")


# Transcriptomic and Epigenetic and with Metabolites
# List_tendency_transcriptome_Genes_H2D   <- Associate_age_trasncriptome_H2D %>% rename(Transcriptome = Tendency)
# List_tendency_methylation_Genes_H2D     <- Associate_age_methylation_H2D   %>% rename(Methylation = Tendency)
#Omics_H2D_tendency <- merge(List_tendency_transcriptome_Genes_H2D, List_tendency_methylation_Genes_H2D, by = "Gene", all = TRUE) #%>% left_join(Genes_conversion)


H2D_metabo_results_compact   <- read.csv("D:/Decicomp/R/MOFA_omics/Data_compacted/metabolomic/landscape/H2D_metabo_results_compact.csv")
Keeg_metabo_H2D_results      <- H2D_metabo_results_compact %>% filter(Entry.type=="pathway") %>%  mutate(KEGG.name = str_remove(KEGG.name, "- Cr.*| -.*")) 
Keeg_metabo_H2D_results      <- Keeg_metabo_H2D_results  %>% mutate(log_pvalue = -log10(Keeg_metabo_H2D_results$p.score)) %>% 
  dplyr::select(!2) %>% dplyr::rename(code=KEGG.id)


Keeg_metabo_H2D_results_KEGG_database <- left_join(All_data_kegg_pathways_genes, Keeg_metabo_H2D_results, by="code") %>%
                                         filter(!is.na(KEGG.name)) %>% filter(!is.na(Gene)) %>% dplyr::select(1:5) %>%
                                         left_join(KEGGS_Magallana_gigas) %>% dplyr::select(6,7,2,3,1,4,5)


#Omics_H2D_tendency_with_metabolites <- Omics_H2D_tendency %>% left_join(Keeg_metabo_H2D_results_KEGG_database) %>% filter %>%
                                       filter(!is.na(Global))

                                       #Omics_H2D_tendency_with_metabolites_sumarize <- Omics_H2D_tendency_with_metabolites %>%
                                                group_by(Global, Genereal, kegg_pathway) %>%
                                                summarize(positive_transcriptome = sum(Transcriptome == "positive", na.rm = TRUE),  
                                                          negative_transcriptome = sum(Transcriptome == "negative", na.rm = TRUE),
             
                                                          positive_methylation = sum(Methylation == "positive", na.rm = TRUE),
                                                          negative_methylation = sum(Methylation == "negative", na.rm = TRUE),
                                                          total_genes = n()) %>% arrange(desc(total_genes))
# write.table(Omics_H2D_tendency_with_metabolites_sumarize, "D:/Decicomp/R/MOFA_omics/Data_compacted/Omics_H2D_tendency_with_metabolites_sumarize.tsv", sep = "\t", row.names = FALSE, quote = FALSE) 


# Omics_H2D_transcriptome_positive_Genes <- Omics_H2D_tendency_with_metabolites %>% filter(Transcriptome=="positive") %>% distinct(Gene) # 290
# Omics_H2D_transcriptome_negative_Genes <- Omics_H2D_tendency_with_metabolites %>% filter(Transcriptome=="negative") %>% distinct(Gene) # 265

# Omics_H2D_methylation_positive_Genes <- Omics_H2D_tendency_with_metabolites %>% filter(Methylation=="positive") %>% distinct(Gene) # 34
# Omics_H2D_methylation_negative_Genes <- Omics_H2D_tendency_with_metabolites %>% filter(Methylation=="negative") %>% distinct(Gene) # 75

# Needs to load  the modules.membership of Kmer for trascriptome and methylation
# kME.positive_H2D_transcriptome <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_H2D_transcriptome.txt", sep = "\t", header = TRUE, row.names = 1)
# kME.negative_H2D_transcriptome <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_H2D_transcriptome.txt", sep = "\t", header = TRUE, row.names = 1)

# kME.positive_H2D_methylation   <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.positive_H2D_methylation.txt", sep = "\t", header = TRUE, row.names = 1)
# kME.negative_H2D_methylation   <- read.table(file = "D:/Decicomp/R/MOFA_omics/Data_compacted/transcriptome/landscape/kME.negative_H2D_methylation.txt", sep = "\t", header = TRUE, row.names = 1)

# --

# Omics_H2D_transcriptome_positive_Genes_kmer        <- Omics_H2D_transcriptome_positive_Genes %>% left_join(kME.positive_H2D_transcriptome)
# List_Omics_H2D_transcriptome_positive_Genes_kmer   <- left_join(List_genes, Omics_H2D_transcriptome_positive_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_H2D_transcriptome_positive_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/positive/List_Omics_H2D_transcriptome_positive_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

# setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/positive/")
# input="List_Omics_H2D_transcriptome_positive_Genes_kmer.txt" 
# goAnnotations="all_go.tab" 
# goDatabase="go.obo" 
# goDivision="BP" 
# source("gomwu.functions.R")
# gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

# GO_terms_Omics_H2D_transcriptome_positive_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/positive/MWU_BP_List_Omics_H2D_transcriptome_positive_Genes_kmer.txt", sep="") %>%
#  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
# mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
# dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
#  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# category_counts_GO_terms_Omics_H2D_transcriptome_positive_Genes_kmer  <- GO_terms_Omics_H2D_transcriptome_positive_Genes_kmer %>% group_by(term_parent) %>% 
#  summarise(Count = n()) %>% arrange(desc(Count))

# --

# Omics_H2D_transcriptome_negative_Genes_kmer        <- Omics_H2D_transcriptome_negative_Genes %>% left_join(kME.negative_H2D_transcriptome)
# List_Omics_H2D_transcriptome_negative_Genes_kmer   <- left_join(List_genes, Omics_H2D_transcriptome_negative_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_H2D_transcriptome_negative_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/negative/List_Omics_H2D_transcriptome_negative_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

# setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/negative/")
# input="List_Omics_H2D_transcriptome_negative_Genes_kmer.txt" 
# goAnnotations="all_go.tab" 
# goDatabase="go.obo" 
# goDivision="BP" 
# source("gomwu.functions.R")
# gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

# GO_terms_Omics_H2D_transcriptome_negative_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/negative/MWU_BP_List_Omics_H2D_transcriptome_negative_Genes_kmer.txt", sep="") %>%
#  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
    #  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
#  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
#   mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# category_counts_GO_terms_Omics_H2D_transcriptome_negative_Genes_kmer  <- GO_terms_Omics_H2D_transcriptome_negative_Genes_kmer %>% group_by(term_parent) %>% 
#  summarise(Count = n()) %>% arrange(desc(Count))

# --

# Omics_H2D_methylation_positive_Genes_kmer      <- Omics_H2D_methylation_positive_Genes %>% left_join(kME.positive_H2D_methylation)
# List_Omics_H2D_methylation_positive_Genes_kmer <- left_join(List_genes, Omics_H2D_methylation_positive_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_H2D_methylation_positive_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/positive/List_Omics_H2D_methylation_positive_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

# setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/positive/")
# input="List_Omics_H2D_methylation_positive_Genes_kmer.txt" 
# goAnnotations="all_go.tab" 
# goDatabase="go.obo" 
# goDivision="BP" 
# source("gomwu.functions.R")
# gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

# GO_terms_Omics_H2D_methylation_positive_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/positive/MWU_BP_List_Omics_H2D_methylation_positive_Genes_kmer.txt", sep="") %>%
#  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
#  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
#  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
# mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# category_counts_GO_terms_Omics_H2D_methylation_positive_Genes_kmer  <- GO_terms_Omics_H2D_methylation_positive_Genes_kmer %>% group_by(term_parent) %>% 
# summarise(Count = n()) %>% arrange(desc(Count))

# --



# Omics_H2D_methylation_negative_Genes_kmer      <- Omics_H2D_methylation_negative_Genes %>% left_join(kME.negative_H2D_methylation)
# List_Omics_H2D_methylation_negative_Genes_kmer <- left_join(List_genes, Omics_H2D_methylation_negative_Genes_kmer, by = "Gene") %>%  mutate(across(everything(), ~ replace_na(., 0)))
# write(write.table(List_Omics_H2D_methylation_negative_Genes_kmer, "D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/negative/List_Omics_H2D_methylation_negative_Genes_kmer.txt", row.names = F, quote = F, col.names=F, sep = ","))

# setwd("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/negative/")
# input="List_Omics_H2D_methylation_negative_Genes_kmer.txt" 
# goAnnotations="all_go.tab" 
# goDatabase="go.obo" 
# goDivision="BP" 
# source("gomwu.functions.R")
# gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25,  Module=TRUE, Alternative="g")

# GO_terms_Omics_H2D_methylation_negative_Genes_kmer  <- read.csv("D:/Decicomp/R/GO_terms/Metabolites/WGCNA/H2D/negative/MWU_BP_List_Omics_H2D_methylation_negative_Genes_kmer.txt", sep="") %>%
# filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
# mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
# dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
# mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

# category_counts_GO_terms_Omics_H2D_methylation_negative_Genes_kmer  <- GO_terms_Omics_H2D_methylation_negative_Genes_kmer %>% group_by(term_parent) %>% 
# summarise(Count = n()) %>% arrange(desc(Count))


library(UpSetR)
library(gridExtra)

#### VEN diagram  ----

# F14R 
positive_F14R_methylation         # dim(positive_F14R_methylation)               693
positive_F14R_methylation_list    <- positive_F14R_methylation$Gene

positive_F14R_transcriptome       # dim(positive_F14R_transcriptome_tendency)    7888
positive_F14R_transcriptome_list   <- positive_F14R_transcriptome$Gene

Omics_F14R_transcriptome_positive_Genes # dim(Omics_F14R_transcriptome_positive_Genes) 535
Omics_F14R_methylation_positive_Genes   # dim(Omics_F14R_methylation_positive_Genes)   34
positive_F14R_metabolites      <- rbind(Omics_F14R_transcriptome_positive_Genes, Omics_F14R_methylation_positive_Genes) # dim(positive_F14R_metabolites) # 569
positive_F14R_metabolites_list <- unique(positive_F14R_metabolites$Gene)  # length(positive_F14R_metabolites_list) # 558



negative_F14R_methylation         # dim(negative_F14R_methylation)               999
negative_F14R_methylation_list    <- negative_F14R_methylation$Gene

negative_F14R_transcriptome       # dim(negative_F14R_transcriptome_tendency)    6658
negative_F14R_transcriptome_list   <- negative_F14R_transcriptome$Gene

Omics_F14R_transcriptome_negative_Genes # dim(Omics_F14R_transcriptome_negative_Genes) 325
Omics_F14R_methylation_negative_Genes   # dim(Omics_F14R_methylation_negative_Genes)   75
negative_F14R_metabolites      <- rbind(Omics_F14R_transcriptome_negative_Genes, Omics_F14R_methylation_negative_Genes) # dim(negative_F14R_metabolites) # 400
negative_F14R_metabolites_list <- unique(negative_F14R_metabolites$Gene)  # length(negative_F14R_metabolites_list) # 387







# H2D 
positive_H2D_methylation         # dim(positive_H2D_methylation)               823
positive_H2D_methylation_list    <- positive_H2D_methylation$Gene

positive_H2D_transcriptome       # dim(positive_H2D_transcriptome_tendency)    6340
positive_H2D_transcriptome_list   <- positive_H2D_transcriptome$Gene

Omics_H2D_transcriptome_positive_Genes # dim(Omics_H2D_transcriptome_positive_Genes) 290
Omics_H2D_methylation_positive_Genes   # dim(Omics_H2D_methylation_positive_Genes)   33
positive_H2D_metabolites      <- rbind(Omics_H2D_transcriptome_positive_Genes, Omics_H2D_methylation_positive_Genes) # dim(positive_H2D_metabolites) # 323
positive_H2D_metabolites_list <- unique(positive_H2D_metabolites$Gene)  # length(positive_H2D_metabolites_list) # 323



negative_H2D_methylation         # dim(negative_H2D_methylation)               2085
negative_H2D_methylation_list    <- negative_H2D_methylation$Gene

negative_H2D_transcriptome       # dim(negative_H2D_transcriptome_tendency)    5193
negative_H2D_transcriptome_list   <- negative_H2D_transcriptome$Gene

Omics_H2D_transcriptome_negative_Genes # dim(Omics_H2D_transcriptome_negative_Genes) 265
Omics_H2D_methylation_negative_Genes   # dim(Omics_H2D_methylation_negative_Genes)   121
negative_H2D_metabolites      <- rbind(Omics_H2D_transcriptome_negative_Genes, Omics_H2D_methylation_negative_Genes) # dim(negative_H2D_metabolites) # 386
negative_H2D_metabolites_list <- unique(negative_H2D_metabolites$Gene)  # length(negative_H2D_metabolites_list) # 386



# Plots 
Gene_sets_F14R_positive <- list(
  Epigenetics     = positive_F14R_methylation_list,
  Transcriptomics = positive_F14R_transcriptome_list,
  Metabolism      = positive_F14R_metabolites_list) 

gene_matrix_F14R <- fromList(Gene_sets_F14R_positive)

F14R_positive  <- upset(
  gene_matrix_F14R, 
  nsets = 6,
  set_size.show = F,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#2f5597","#2f5597","#2f5597"),
  empty.intersections = T,
  show.numbers = "yes",
  keep.order = F)  

F14R_positive

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/F14R_positive.tiff", 
     width = 7, height = 6, units = "in", res = 600)
print(F14R_positive)
dev.off()


Gene_sets_F14R_negative <- list(
  Epigenetics     = negative_F14R_methylation_list,
  Transcriptomics = negative_F14R_transcriptome_list,
  Metabolism      = negative_F14R_metabolites_list) 

gene_matrix_F14R <- fromList(Gene_sets_F14R_negative)

F14R_negative  <- upset(
  gene_matrix_F14R, 
  nsets = 6,
  set_size.show = F,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#2f5597","#2f5597","#2f5597"),
  empty.intersections = T,
  show.numbers = "yes",
  keep.order = F)  

F14R_negative

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/F14R_negative.tiff", 
     width = 7, height = 6, units = "in", res = 600)
print(F14R_negative)
dev.off()



Gene_sets_H2D_positive <- list(
  Epigenetics     = positive_H2D_methylation_list,
  Transcriptomics = positive_H2D_transcriptome_list,
  Metabolism      = positive_H2D_metabolites_list) 

gene_matrix_H2D <- fromList(Gene_sets_H2D_positive)

H2D_positive  <- upset(
  gene_matrix_H2D, 
  nsets = 6,
  set_size.show = F,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#385700","#385700","#385700"),
  empty.intersections = T,
  show.numbers = "yes",
  keep.order = F)  

H2D_positive

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/H2D_positive.tiff", 
     width = 7, height = 6, units = "in", res = 600)
print(H2D_positive)
dev.off()



Gene_sets_H2D_negative <- list(
  Epigenetics     = negative_H2D_methylation_list,
  Transcriptomics = negative_H2D_transcriptome_list,
  Metabolism      = negative_H2D_metabolites_list) 

gene_matrix_H2D <- fromList(Gene_sets_H2D_negative)

H2D_negative  <- upset(
  gene_matrix_H2D, 
  nsets = 6,
  set_size.show = F,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#385700","#385700","#385700"),
  empty.intersections = T,
  show.numbers = "yes",
  keep.order = F)  

H2D_negative

tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/H2D_negative.tiff", 
     width = 7, height = 6, units = "in", res = 600)
print(H2D_negative)
dev.off()


###  


Gene_sets_families_positive <- list(
  Epigenetics_F14R     = positive_F14R_methylation_list,
  Transcriptomics_F14R = positive_F14R_transcriptome_list,
  Metabolism_F14R      = positive_F14R_metabolites_list, 
  Epigenetics_H2D      = positive_H2D_methylation_list,
  Transcriptomics_H2D = positive_H2D_transcriptome_list,
  Metabolism_H2D      = positive_H2D_metabolites_list) 

gene_matrix_families_positive <- fromList(Gene_sets_families_positive)


families_positive  <- upset(
  gene_matrix_families_positive, 
  nsets = 6,
  set_size.show = F,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=2.5,
  sets.bar.color=c("#2f5597","#385700", "#385700", "#2f5597", "#2f5597","#385700"),
  empty.intersections = NULL,
  show.numbers = "yes",
  keep.order =F) 

families_positive


tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/families_positive.tiff", 
     width = 17, height = 8, units = "in", res = 600)
print(families_positive)
dev.off()



Gene_sets_families_negative <- list(
  Epigenetics_F14R     = negative_F14R_methylation_list,
  Transcriptomics_F14R = negative_F14R_transcriptome_list,
  Metabolism_F14R      = negative_F14R_metabolites_list, 
  Epigenetics_H2D      = negative_H2D_methylation_list,
  Transcriptomics_H2D = negative_H2D_transcriptome_list,
  Metabolism_H2D      = negative_H2D_metabolites_list) 

gene_matrix_families_negative <- fromList(Gene_sets_families_negative)


families_negative  <- upset(
  gene_matrix_families_negative, 
  nsets = 6,
  set_size.show = F,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=2.5,
  sets.bar.color=c("#2f5597","#385700", "#385700", "#2f5597", "#2f5597","#385700"),
  empty.intersections = NULL,
  show.numbers = "yes",
  keep.order = F) 

families_negative


tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Final_figures/families_negative.tiff", 
     width = 17, height = 8, units = "in", res = 600)
print(families_negative)
dev.off()

















# Transcriptome comparative Families GO terms ----
Name_GO_interest <- "immune system process"

# F14R positive
List_kME.positive_F14R_transcriptome_clean
dim(List_kME.positive_F14R_transcriptome_clean)

GO_terms_genes_kME.positive_F14R_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/MWU_BP_List_genes_kME.positive_F14R.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.positive_F14R_transcriptome)

category_counts_positive_F14R_transcriptome  <- GO_terms_genes_kME.positive_F14R_transcriptome %>%
  group_by(term_parent) %>% summarise(Count = n()) %>%
  arrange(desc(Count))

GO_extract_GO_positive_F14R_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/BP_List_genes_kME.positive_F14R.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(name)  %>%  right_join(GO_terms_genes_kME.positive_F14R_transcriptome)  %>%
  filter(term_parent=="cellular process")
dim(GO_extract_GO_positive_F14R_transcriptome)

GO_terms_kME.positive_F14R_transcriptome_2  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/MWU_BP_List_genes_kME.positive_F14R.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, pval)  %>% left_join(GO_terms_big) %>% filter(term_parent==Name_GO_interest) %>%
  dplyr::select(last_term, pval)
# write(write.table(GO_terms_kME.positive_F14R_transcriptome_2, "D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_kME.positive_F14R_transcriptome_2.txt", row.names = F, quote = F, col.names=F, sep = "\t"))


Genes_extract_GO_positive_F14R_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/positive/BP_List_genes_kME.positive_F14R.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(seq) %>% rename(Gene=seq)

dim(Genes_extract_GO_positive_F14R_transcriptome)
# write(write.table(Genes_extract_GO_positive_F14R_transcriptome, "D:/Decicomp/R/MOFA_omics/Data_compacted/Genes_extract_GO_positive_F14R_transcriptome.txt", row.names = F, quote = F, col.names=F, sep = ","))


# F14R negative
List_kME.negative_F14R_transcriptome_clean
dim(List_kME.negative_F14R_transcriptome_clean)

GO_terms_genes_kME.negative_F14R_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/MWU_BP_List_genes_kME.negative_F14R.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.negative_F14R_transcriptome)

category_counts_negative_F14R_transcriptome  <- GO_terms_genes_kME.negative_F14R_transcriptome %>%
  group_by(term_parent) %>% summarise(Count = n()) %>%
  arrange(desc(Count))

GO_extract_GO_negative_F14R_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/BP_List_genes_kME.negative_F14R.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(name)  %>%  right_join(GO_terms_genes_kME.negative_F14R_transcriptome)  %>%
  filter(term_parent=="immune system process")
dim(GO_extract_GO_negative_F14R_transcriptome)

GO_terms_kME.negative_F14R_transcriptome_2  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/MWU_BP_List_genes_kME.negative_F14R.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, pval)  %>% left_join(GO_terms_big) %>% filter(term_parent==Name_GO_interest) %>%
  dplyr::select(last_term, pval)
write(write.table(GO_terms_kME.negative_F14R_transcriptome_2, "D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_kME.negative_F14R_transcriptome_2.txt", row.names = F, quote = F, col.names=F, sep = "\t"))

Genes_extract_GO_negative_F14R_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/F14R/negative/BP_List_genes_kME.negative_F14R.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(seq) %>% rename(Gene=seq) %>% left_join(Annot_Cg_Ros_IHPEV1)

dim(Genes_extract_GO_negative_F14R_transcriptome)
# write(write.table(Genes_extract_GO_negative_F14R_transcriptome, "D:/Decicomp/R/MOFA_omics/Data_compacted/Genes_extract_GO_negative_F14R_transcriptome.txt", row.names = F, quote = F, col.names=F, sep = ","))


# H2D positive
List_kME.positive_H2D_transcriptome_clean
dim(List_kME.positive_H2D_transcriptome_clean)

GO_terms_genes_kME.positive_H2D_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/MWU_BP_List_genes_kME.positive_H2D.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.positive_H2D_transcriptome)

category_counts_positive_H2D_transcriptome  <- GO_terms_genes_kME.positive_H2D_transcriptome %>%
  group_by(term_parent) %>% summarise(Count = n()) %>%
  arrange(desc(Count))

GO_extract_GO_positive_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/BP_List_genes_kME.positive_H2D.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(name)  %>%  right_join(GO_terms_genes_kME.positive_H2D_transcriptome)  %>%
  filter(term_parent=="immune system process")
dim(GO_extract_GO_positive_H2D_transcriptome)

GO_terms_kME.positive_H2D_transcriptome_2  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/MWU_BP_List_genes_kME.positive_H2D.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, pval)  %>% left_join(GO_terms_big) %>% filter(term_parent==Name_GO_interest) %>%
  dplyr::select(last_term, pval)
write(write.table(GO_terms_kME.positive_H2D_transcriptome_2, "D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_kME.positive_H2D_transcriptome_2.txt", row.names = F, quote = F, col.names=F, sep = "\t"))


Genes_extract_GO_positive_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/positive/BP_List_genes_kME.positive_H2D.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(seq) %>% rename(Gene=seq)

dim(Genes_extract_GO_positive_H2D_transcriptome)
write(write.table(Genes_extract_GO_positive_H2D_transcriptome, "D:/Decicomp/R/MOFA_omics/Data_compacted/Genes_extract_GO_positive_H2D_transcriptome.txt", row.names = F, quote = F, col.names=F, sep = ","))

# H2D negative
List_kME.negative_H2D_transcriptome_clean
dim(List_kME.negative_H2D_transcriptome_clean)

GO_terms_genes_kME.negative_H2D_transcriptome  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/MWU_BP_List_genes_kME.negative_H2D.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, p.adj)  %>% left_join(GO_terms_big) %>% mutate(Count=1) %>% arrange(p.adj)  %>%
  mutate(term_parent = factor(term_parent, levels = unique(term_parent)),last_term = factor(last_term, levels = unique(last_term)))

dim(GO_terms_genes_kME.negative_H2D_transcriptome)

category_counts_negative_H2D_transcriptome  <- GO_terms_genes_kME.negative_H2D_transcriptome %>%
  group_by(term_parent) %>% summarise(Count = n()) %>%
  arrange(desc(Count))

GO_extract_GO_negative_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/BP_List_genes_kME.negative_H2D.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(name)  %>%  right_join(GO_terms_genes_kME.negative_H2D_transcriptome)  %>%
  filter(term_parent=="immune system process")
dim(GO_extract_GO_negative_H2D_transcriptome)

GO_terms_kME.negative_H2D_transcriptome_2  <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/MWU_BP_List_genes_kME.negative_H2D.txt", sep="") %>%
  filter(level != -1) %>% filter(p.adj <= 0.05) %>% filter(level >= 2)  %>%
  mutate(last_term = sapply(strsplit(as.character(term), ";"), tail, 1)) %>%
  dplyr::select( name, last_term, pval)  %>% left_join(GO_terms_big) %>% filter(term_parent==Name_GO_interest) %>%
  dplyr::select(last_term, pval)
write(write.table(GO_terms_kME.negative_H2D_transcriptome_2, "D:/Decicomp/R/MOFA_omics/Data_compacted/GO_terms_kME.negative_H2D_transcriptome_2.txt", row.names = F, quote = F, col.names=F, sep = "\t"))

Genes_extract_GO_negative_H2D_transcriptome <- read.csv("D:/Decicomp/R/GO_terms/RNA/WGCN/H2D/negative/BP_List_genes_kME.negative_H2D.txt", sep="") %>%
  filter(lev != -1) %>%  filter(lev >= 2) %>%  filter(value > 0) %>%
  left_join(GO_terms_big) %>% filter(term_parent== Name_GO_interest) %>% distinct(seq) %>% rename(Gene=seq) %>% left_join(Annot_Cg_Ros_IHPEV1)

dim(Genes_extract_GO_negative_H2D_transcriptome)
# write(write.table(Genes_extract_GO_negative_H2D_transcriptome, "D:/Decicomp/R/MOFA_omics/Data_compacted/Genes_extract_GO_negative_H2D_transcriptome.txt", row.names = F, quote = F, col.names=F, sep = ","))

### Intersect Transcriptome comparative Families ----
Genes_extract_GO_positive_F14R_transcriptome 
dim(Genes_extract_GO_positive_F14R_transcriptome)

Genes_extract_GO_negative_F14R_transcriptome 
dim(Genes_extract_GO_negative_F14R_transcriptome)

Genes_extract_GO_positive_H2D_transcriptome 
dim(Genes_extract_GO_positive_H2D_transcriptome)

Genes_extract_GO_negative_H2D_transcriptome 
dim(Genes_extract_GO_negative_H2D_transcriptome)

intersect_GO_positive_transcriptome <- intersect(Genes_extract_GO_positive_F14R_transcriptome$Gene, 
                                                 Genes_extract_GO_positive_H2D_transcriptome$Gene) %>% as.data.frame()
dim(intersect_GO_positive_transcriptome)
# write(write.table(intersect_GO_positive_transcriptome, "D:/Decicomp/R/MOFA_omics/Data_compacted/intersect_GO_positive_transcriptome.txt", row.names = F, quote = F, col.names=F, sep = ","))


intersect_GO_negative_transcriptome <- intersect(Genes_extract_GO_negative_F14R_transcriptome$Gene, 
                                                 Genes_extract_GO_negative_H2D_transcriptome$Gene) %>% as.data.frame()
dim(intersect_GO_negative_transcriptome)
# write(write.table(intersect_GO_negative_transcriptome, "D:/Decicomp/R/MOFA_omics/Data_compacted/intersect_GO_negative_transcriptome.txt", row.names = F, quote = F, col.names=F, sep = ","))



# Mehtylation comparative Families ----
Name_GO_interest <- "immune system process"

GO_positive_F14R_methylation_genes  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/positive/David_GO_BP_List_genes_kME.positive_F14R.txt") %>%
  separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big)  %>%
  dplyr::select(last_term, name, Genes, term_parent) %>%
  rename(Gene=Genes) %>% 
  separate_rows(Gene, sep = ",") %>%  mutate(Gene = str_trim(Gene)) %>%
  filter(term_parent== Name_GO_interest) %>% 
  distinct(Gene) 

GO_positive_F14R_methylation_genes_list_name <- GO_positive_F14R_methylation_genes %>% left_join(Annot_Cg_Ros_IHPEV1)
dim(GO_positive_F14R_methylation_genes)




GO_negative_F14R_methylation_genes  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/F14R/negative/David_GO_BP_List_genes_kME.negative_F14R.txt") %>%
  separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big)  %>%
  dplyr::select(last_term, name, Genes, term_parent) %>%
  rename(Gene=Genes) %>%
  separate_rows(Gene, sep = ",") %>%  mutate(Gene = str_trim(Gene)) %>%
  filter(term_parent== Name_GO_interest) %>% 
  distinct(Gene)  

dim(GO_negative_F14R_methylation_genes)



GO_positive_H2D_methylation_genes  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/positive/David_GO_BP_List_genes_kME.positive_H2D.txt") %>%
  separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big)  %>%
  dplyr::select(last_term, name, Genes, term_parent) %>%
  rename(Gene=Genes) %>%
  separate_rows(Gene, sep = ",") %>%  mutate(Gene = str_trim(Gene)) %>%
  filter(term_parent== Name_GO_interest) %>% 
  distinct(Gene) 

GO_positive_H2D_methylation_genes_list_name <- GO_positive_H2D_methylation_genes %>% left_join(Annot_Cg_Ros_IHPEV1)
dim(GO_positive_H2D_methylation_genes)

dim(GO_positive_H2D_methylation_genes)

GO_negative_H2D_methylation_genes  <-  read.delim("D:/Decicomp/R/GO_terms/DNA/WGCN/H2D/negative/David_GO_BP_List_genes_kME.negative_H2D.txt") %>%
  separate(Term, into = c("last_term", "name"), sep = "~") %>% left_join(GO_terms_big)  %>%
  dplyr::select(last_term, name, Genes, term_parent) %>%
  rename(Gene=Genes) %>%
  separate_rows(Gene, sep = ",") %>%  mutate(Gene = str_trim(Gene)) %>%
  filter(term_parent== Name_GO_interest) %>% 
  distinct(Gene)  

dim(GO_negative_H2D_methylation_genes)



intersect_GO_positive_methylation <- intersect(GO_positive_F14R_methylation_genes$Gene, 
                                               GO_positive_H2D_methylation_genes$Gene) %>% as.data.frame()
dim(intersect_GO_positive_methylation)


intersect_GO_negative_methylation <- intersect(GO_negative_F14R_methylation_genes$Gene, 
                                               GO_negative_H2D_methylation_genes$Gene) %>% as.data.frame()
dim(intersect_GO_negative_methylation)


library(treemap) 	# se pueso 0.5 uni prot standarized y luego sacamos de tree map de R y abajo cambiammos.					

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006139","nucleobase-containing compound metabolic process",18.39314412781983,23.053189440138908,0.9309068335558746,0.0781618,"nucleobase-containing compound metabolic process"),
                     c("GO:0006396","RNA processing",4.402491558667012,16.205043383586016,0.8241131716197491,0,"RNA processing"),
                     c("GO:0000209","protein polyubiquitination",0.4445080348651629,3.5477742906833183,0.8729881733531571,0.4361715,"RNA processing"),
                     c("GO:0000966","RNA 5'-end processing",0.16293664216619538,3.4152958643229816,0.821593789626112,0.33907537,"RNA processing"),
                     c("GO:0001510","RNA methylation",0.712861962256531,6.792280395298492,0.8236606994510961,0.39383712,"RNA processing"),
                     c("GO:0006259","DNA metabolic process",5.572970721566783,5.413548660986488,0.8644035893583336,0.21831618,"RNA processing"),
                     c("GO:0006270","DNA replication initiation",0.19247287749043998,4.017815129827215,0.8513490216573609,0.21831618,"RNA processing"),
                     c("GO:0006399","tRNA metabolic process",2.3934466773236682,5.788233955140595,0.8529087476910917,0.27280945,"RNA processing"),
                     c("GO:0006400","tRNA modification",1.0258525252816637,4.952629433755633,0.7994200576478064,0.43849162,"RNA processing"),
                     c("GO:0006468","protein phosphorylation",3.756826993126081,7.651908030358493,0.8440340052964356,0.24781368,"RNA processing"),
                     c("GO:0006513","protein monoubiquitination",0.044016374691955525,3.5104253132710617,0.8729881733531571,0.47593805,"RNA processing"),
                     c("GO:0006793","phosphorus metabolic process",12.714975181439863,9.130215034208002,0.9116170965399192,0.10046538,"RNA processing"),
                     c("GO:0008033","tRNA processing",1.3800535206429487,7.294860879540894,0.819644491853204,0.29637339,"RNA processing"),
                     c("GO:0009057","macromolecule catabolic process",3.5074771718129205,7.03252681023289,0.9008922302072931,0.11926201,"RNA processing"),
                     c("GO:0009059","macromolecule biosynthetic process",16.264263355707197,2.8812692475128694,0.8744044016374207,0.15753374,"RNA processing"),
                     c("GO:0009187","cyclic nucleotide metabolic process",0.160918332752372,5.171402250706756,0.8476715251394624,0.18383601,"RNA processing"),
                     c("GO:0009451","RNA modification",1.762112108620884,8.115538590648507,0.835623361313993,0.27280945,"RNA processing"),
                     c("GO:0016070","RNA metabolic process",8.105171248384922,21.217057881356958,0.8644035893583336,0.21831618,"RNA processing"),
                     c("GO:0016071","mRNA metabolic process",1.6631731043434954,6.8555320718191854,0.8529087476910917,0.27280945,"RNA processing"),
                     c("GO:0016072","rRNA metabolic process",1.615169337882767,3.079483351985235,0.8529087476910917,0.27280945,"RNA processing"),
                     c("GO:0016073","snRNA metabolic process",0.16282341993078575,5.138145366062886,0.8529087476910917,0.27280945,"RNA processing"),
                     c("GO:0016310","phosphorylation",6.253736641437635,5.608853998860012,0.9001210413683305,0.22444672,"RNA processing"),
                     c("GO:0018193","peptidyl-amino acid modification",2.338036699855386,4.519866903060484,0.8786330158678246,0.27057677,"RNA processing"),
                     c("GO:0018212","peptidyl-tyrosine modification",0.03475922627074853,5.781112833871167,0.876542963464948,0.27057677,"RNA processing"),
                     c("GO:0019637","organophosphate metabolic process",6.773690066732693,3.4535693431683026,0.9001534746605937,0.22392112,"RNA processing"),
                     c("GO:0031123","RNA 3'-end processing",0.3622643876048038,4.801858268077015,0.821593789626112,0.33907537,"RNA processing"),
                     c("GO:0032774","RNA biosynthetic process",6.5238873564778945,3.3882965176181927,0.8251737981680586,0.28440176,"RNA processing"),
                     c("GO:0034654","nucleobase-containing compound biosynthetic process",10.209893841355267,3.7046309009569254,0.8822491399151907,0.18383601,"RNA processing"),
                     c("GO:0043412","macromolecule modification",10.203408176348653,17.82640604327904,0.908390855894966,0.11926201,"RNA processing"),
                     c("GO:0043414","macromolecule methylation",0.926246494356536,7.994345512752895,0.877954446694526,0.11926201,"RNA processing"),
                     c("GO:0043687","post-translational protein modification",1.8860683042179072,4.502720793989289,0.8786330158678246,0.27057677,"RNA processing"),
                     c("GO:0044249","cellular biosynthetic process",23.43500411119782,4.724822598822131,0.9025717282161574,0.13625948,"RNA processing"),
                     c("GO:0051603","proteolysis involved in protein catabolic process",1.7812048234050641,7.1436557966200755,0.8890755514357329,0.2607896,"RNA processing"),
                     c("GO:0052652","cyclic purine nucleotide metabolic process",0.122159407948132,2.98273086775907,0.8476715251394624,0.36823636,"RNA processing"),
                     c("GO:0070646","protein modification by small protein removal",0.39107698516360445,4.203648046051087,0.8735924309346343,0.27057677,"RNA processing"),
                     c("GO:0090407","organophosphate biosynthetic process",4.68253445118795,3.863813706007878,0.888568106209005,0.29229367,"RNA processing"),
                     c("GO:1901293","nucleoside phosphate biosynthetic process",2.839852415308414,3.033536658759147,0.8376240097166288,0.33652359,"RNA processing"),
                     c("GO:1901657","glycosyl compound metabolic process",0.5856321672444037,3.260265350703324,0.9358540056350532,0.29472003,"RNA processing"),
                     c("GO:1901658","glycosyl compound catabolic process",0.16646129958155523,4.174088881363266,0.9263548955576935,0.29472003,"RNA processing"),
                     c("GO:1903510","mucopolysaccharide metabolic process",0.06785211659862092,4.74237199569436,0.9010988339312026,0.11926201,"RNA processing"),
                     c("GO:0006520","amino acid metabolic process",5.5080894579711925,2.8023395405175293,0.9309068335558746,0.0781618,"amino acid metabolic process"),
                     c("GO:0006887","exocytosis",0.35007084512177816,4.047465037932256,0.9527803805091507,0.03189197,"exocytosis"),
                     c("GO:0016192","vesicle-mediated transport",2.591260690701909,3.6879742802298576,0.9527803805091507,0.1831236,"exocytosis"),
                     c("GO:0051169","nuclear transport",0.5619367224555284,4.005467423797187,0.9474626618210207,0.1831236,"exocytosis"),
                     c("GO:0051170","import into nucleus",0.2587448054992138,4.167310278626721,0.9474626618210207,0.38606433,"exocytosis"),
                     c("GO:0051641","cellular localization",5.613085851842994,4.274443584038131,0.9532372471697568,0.1767677,"exocytosis"),
                     c("GO:0051649","establishment of localization in cell",3.464570867298568,3.434502168291362,0.9487666151767916,0.31269958,"exocytosis"),
                     c("GO:0055085","transmembrane transport",12.953593503918547,5.729620597239131,0.9527803805091507,0.1831236,"exocytosis"),
                     c("GO:0007049","cell cycle",2.051845347681118,5.043990879569662,0.9670869746379425,0.03189197,"cell cycle"),
                     c("GO:0007155","cell adhesion",1.0284049482842672,6.360006194414478,0.9670869746379425,0.03189197,"cell adhesion"),
                     c("GO:0008152","metabolic process",58.08643542977757,51.39392012426401,1,-0,"metabolic process"),
                     c("GO:0009056","catabolic process",9.05402280817628,5.723137394973818,0.9432850422786514,0.05898132,"catabolic process"),
                     c("GO:0009058","biosynthetic process",28.507692540195496,7.777630091302633,0.9432850422786514,0.05898132,"biosynthetic process"),
                     c("GO:0009311","oligosaccharide metabolic process",0.35118829935821205,3.171547448521602,0.9232884827700948,0.0781618,"oligosaccharide metabolic process"),
                     c("GO:0005996","monosaccharide metabolic process",1.0398527008253557,3.163031745112628,0.9172364861359323,0.321463,"oligosaccharide metabolic process"),
                     c("GO:0016051","carbohydrate biosynthetic process",1.1871720585638565,2.9176838833428516,0.9097314912034135,0.321463,"oligosaccharide metabolic process"),
                     c("GO:0016477","cell migration",0.48278207313949656,3.0842650488074,0.9670869746379425,0.03189197,"cell migration"),
                     c("GO:0022402","cell cycle process",1.8491578554743764,8.570483892540278,0.9484653462541209,0.03189197,"cell cycle process"),
                     c("GO:1903047","mitotic cell cycle process",0.8772532640124455,3.1833132928267798,0.9481077269983269,0.42196384,"cell cycle process"),
                     c("GO:0022412","cellular process involved in reproduction in multicellular organism",0.14250987408653654,4.143346697409058,0.9596869787243146,0.03189197,"cellular process involved in reproduction in multicellular organism"),
                     c("GO:0070192","chromosome organization involved in meiotic cell cycle",0.08269653620199752,4.0130219385826305,0.9081795929675467,0.49555985,"cellular process involved in reproduction in multicellular organism"),
                     c("GO:0030029","actin filament-based process",0.8043086081733359,3.6810216556838604,0.9670869746379425,0.03189197,"actin filament-based process"),
                     c("GO:0032259","methylation",2.6585910004764686,6.409325004842066,0.9432850422786514,0.05898132,"methylation"),
                     c("GO:0033554","cellular response to stress",4.425081855984172,14.296402560153844,0.9108350068111117,0.03189197,"cellular response to stress"),
                     c("GO:0000725","recombinational repair",0.5276648440676299,4.214770773541057,0.7956140654198518,0.38585642,"cellular response to stress"),
                     c("GO:0006298","mismatch repair",0.21791588286933303,2.8535517511054858,0.7956140654198518,0.38585642,"cellular response to stress"),
                     c("GO:0006301","postreplication repair",0.11389418476323088,2.829959039070865,0.7956140654198518,0.38585642,"cellular response to stress"),
                     c("GO:0006302","double-strand break repair",0.7311276624516326,5.1590481543384055,0.7956140654198518,0.33851967,"cellular response to stress"),
                     c("GO:0006974","DNA damage response",3.1870828591331772,14.655923206067566,0.9061316982275448,0.29516635,"cellular response to stress"),
                     c("GO:0007165","signal transduction",8.851839893322994,7.727827470817362,0.8636959415565606,0.25923296,"cellular response to stress"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7427092633657988,9.854963034705403,0.863332521094354,0.21493909,"cellular response to stress"),
                     c("GO:0007167","enzyme-linked receptor protein signaling pathway",0.5096157429315727,6.6493580541708965,0.8553533975181411,0.26324161,"cellular response to stress"),
                     c("GO:0007169","cell surface receptor protein tyrosine kinase signaling pathway",0.36239976201670654,5.062216230456837,0.8533470346251024,0.43969377,"cellular response to stress"),
                     c("GO:0034976","response to endoplasmic reticulum stress",0.2589072547934972,2.7963372857678013,0.9061316982275448,0.33851967,"cellular response to stress"),
                     c("GO:0035556","intracellular signal transduction",4.385343312708344,5.15730638328922,0.863332521094354,0.26324161,"cellular response to stress"),
                     c("GO:0036297","interstrand cross-link repair",0.07709695825510948,2.8065508562358517,0.7956140654198518,0.38585642,"cellular response to stress"),
                     c("GO:0042276","error-prone translesion synthesis",0.07315879354521021,3.1106630686724728,0.7825872320716518,0.38585642,"cellular response to stress"),
                     c("GO:0043401","steroid hormone receptor signaling pathway",0.053017542407019065,3.0353540196364315,0.863332521094354,0.26324161,"cellular response to stress"),
                     c("GO:0071526","semaphorin-plexin signaling pathway",0.0656516670669647,4.8906718665272635,0.8553533975181411,0.43969377,"cellular response to stress"),
                     c("GO:0043170","macromolecule metabolic process",33.33895178165279,35.41144656575706,0.9432850422786514,0.05898132,"macromolecule metabolic process"),
                     c("GO:0044237","cellular metabolic process",39.640441131550894,14.008136106311058,0.9232466595090747,0.05898132,"cellular metabolic process"),
                     c("GO:0044281","small molecule metabolic process",14.792283225323425,3.3462902957525524,0.9432850422786514,0.05898132,"small molecule metabolic process"),
                     c("GO:0048869","cellular developmental process",2.113307792037927,3.4093954137594467,0.9619403224547716,0.03189197,"cellular developmental process"),
                     c("GO:0048468","cell development",1.0529397144269397,3.6617660246689243,0.9619403224547716,0.35418145,"cellular developmental process"),
                     c("GO:0051056","regulation of small GTPase mediated signal transduction",0.24989870301960257,7.705770720393395,0.9131464819466458,-0,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0000184","nuclear-transcribed mRNA catabolic process, nonsense-mediated decay",0.09922452121885607,3.4289771274641887,0.7536925399088289,0.49980419,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0006357","regulation of transcription by RNA polymerase II",4.1876991326832185,4.358643541832591,0.9148563084099328,0.15132909,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0009891","positive regulation of biosynthetic process",2.016283720350727,4.633747925431611,0.907748236873826,0.3386863,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0009892","negative regulation of metabolic process",2.9251604469836634,4.398733148404224,0.9084573205841351,0.33279065,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0009893","positive regulation of metabolic process",2.881710183468756,6.514971846082802,0.9101408799094697,0.3386863,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0010608","post-transcriptional regulation of gene expression",1.6327261684300864,2.8916215475071922,0.9163313768677227,0.21476949,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0010629","negative regulation of gene expression",1.0018051070218412,5.097057509599465,0.9028141429710156,0.38346258,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0019219","regulation of nucleobase-containing compound metabolic process",11.820748427527752,19.091140355587363,0.9189925339582946,0.21066885,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0031047","regulatory ncRNA-mediated gene silencing",0.2015528084997015,5.091927486856027,0.9013663468301626,0.40605715,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0031324","negative regulation of cellular metabolic process",2.62213590202752,5.366114331514383,0.9040490008508437,0.33279065,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0033044","regulation of chromosome organization",0.1861028960721779,2.8277879472845475,0.9286525435868986,0.15132909,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0045935","positive regulation of nucleobase-containing compound metabolic process",1.8080311091383094,4.353101821117515,0.9075790300750759,0.38508742,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0048519","negative regulation of biological process",4.664847168934615,6.333916010549782,0.9312434880267227,0.14346073,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0048522","positive regulation of cellular process",3.6764170538087693,6.943046343618236,0.9164566044562067,0.15132909,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0048523","negative regulation of cellular process",4.3068778422165455,6.533678134013103,0.9143676714778236,0.15132909,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0051052","regulation of DNA metabolic process",0.2657522773298909,3.8802131605845487,0.9177265591605476,0.23183876,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0051054","positive regulation of DNA metabolic process",0.07547738801816341,3.556914233585299,0.9073043295597335,0.38508742,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0080090","regulation of primary metabolic process",14.365627383358822,17.969555274339392,0.9196186982523421,0.20062793,"regulation of small GTPase mediated signal transduction"),
                     c("GO:1901987","regulation of cell cycle phase transition",0.4407249353907159,3.5915471225394153,0.9286525435868986,0.15132909,"regulation of small GTPase mediated signal transduction"),
                     c("GO:1902531","regulation of intracellular signal transduction",0.9735955409342439,4.2666476917657015,0.9150426406687414,0.44192018,"regulation of small GTPase mediated signal transduction"),
                     c("GO:0051276","chromosome organization",1.530548023678806,5.75931319606703,0.9298854038655432,0.03189197,"chromosome organization"),
                     c("GO:0000226","microtubule cytoskeleton organization",0.8833918782540008,3.2822535347691986,0.922417783492871,0.295867,"chromosome organization"),
                     c("GO:0000280","nuclear division",0.43374453844241945,3.70214957855307,0.9298854038655434,0.295867,"chromosome organization"),
                     c("GO:0006996","organelle organization",6.55438598080312,3.5103035032648844,0.9363557930755378,0.22720981,"chromosome organization"),
                     c("GO:0007015","actin filament organization",0.47463991760177976,2.9127877234094046,0.9220788017564445,0.43043603,"chromosome organization"),
                     c("GO:0043933","protein-containing complex organization",2.844366536607136,3.029532402148241,0.9363557930755378,0.22720981,"chromosome organization"),
                     c("GO:0044085","cellular component biogenesis",6.870847051478853,2.913190419477842,0.9385950596988646,0.21140611,"chromosome organization"),
                     c("GO:0051225","spindle assembly",0.14606652909016435,3.284053705637745,0.907252869975767,0.46719168,"chromosome organization"),
                     c("GO:0097435","supramolecular fiber organization",0.789545413217101,4.8772742089186,0.9363557930755378,0.22720981,"chromosome organization"),
                     c("GO:0051301","cell division",1.377833380287743,5.896508149298385,0.9670869746379425,0.03189197,"cell division"),
                     c("GO:0071840","cellular component organization or biogenesis",14.26840886478907,7.260886717037079,0.9670869746379425,0.03189197,"cellular component organization or biogenesis"),
                     c("GO:0097502","mannosylation",0.08825919385473027,2.875197620336557,0.9432850422786514,0.05898132,"mannosylation"));

stuff        <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value          <- as.numeric( 1 );
stuff$frequency      <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness     <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Revigo_TreeMap_F14R_positive.tiff", width = 12, height = 8, units = 'in', res = 600)

treemap(
  stuff,
  index = c("representative", "description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  bg.labels = "#CCCCCCAA",
  position.legend = "none")

dev.off()


#################



revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000281","mitotic cytokinesis",0.10825030246335648,5.108691076121528,0.9429425096959133,0.03189197,"mitotic cytokinesis"),
                     c("GO:0000910","cytokinesis",0.33092398057283656,3.109885706717696,0.9457294559969499,0.46520082,"mitotic cytokinesis"),
                     c("GO:0006357","regulation of transcription by RNA polymerase II",4.1876991326832185,10.328691843307576,0.8933650129860267,-0,"regulation of transcription by RNA polymerase II"),
                     c("GO:0007165","signal transduction",8.851839893322994,17.3128113597647,0.8447721319310091,0.25923296,"regulation of transcription by RNA polymerase II"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7427092633657988,4.880000644310527,0.8442263811677575,0.26324161,"regulation of transcription by RNA polymerase II"),
                     c("GO:0007186","G protein-coupled receptor signaling pathway",1.464221945905271,3.514609595558031,0.8442263811677575,0.26324161,"regulation of transcription by RNA polymerase II"),
                     c("GO:0007218","neuropeptide signaling pathway",0.15345304927416917,6.176546519333313,0.8442263811677575,0.26324161,"regulation of transcription by RNA polymerase II"),
                     c("GO:0008063","Toll signaling pathway",0.00457319376937054,3.6492464088542125,0.8382323400382036,0.43969377,"regulation of transcription by RNA polymerase II"),
                     c("GO:0016055","Wnt signaling pathway",0.19157448366599422,9.42210173926811,0.8382323400382036,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0019722","calcium-mediated signaling",0.1524365105084264,6.494678731097938,0.8429293473778766,0.26324161,"regulation of transcription by RNA polymerase II"),
                     c("GO:0023056","positive regulation of signaling",0.6428118574791976,7.153470138934448,0.8454641096664763,0.14346073,"regulation of transcription by RNA polymerase II"),
                     c("GO:0030111","regulation of Wnt signaling pathway",0.18630718836650392,5.2472213994056185,0.8576380811540493,0.44192018,"regulation of transcription by RNA polymerase II"),
                     c("GO:0030178","negative regulation of Wnt signaling pathway",0.1054320533428348,3.8976024442342334,0.8543941900103375,0.44192018,"regulation of transcription by RNA polymerase II"),
                     c("GO:0031323","regulation of cellular metabolic process",14.72286076554673,9.861897799546599,0.8937165466298352,0.20062793,"regulation of transcription by RNA polymerase II"),
                     c("GO:0032411","positive regulation of transporter activity",0.051941931170627835,6.4854087050962,0.8708847089043874,0.3386863,"regulation of transcription by RNA polymerase II"),
                     c("GO:0032886","regulation of microtubule-based process",0.17552892382609833,6.044904544965419,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0035556","intracellular signal transduction",4.385343312708344,5.084467871858955,0.8442263811677575,0.26324161,"regulation of transcription by RNA polymerase II"),
                     c("GO:0040017","positive regulation of locomotion",0.16887096411342487,3.1047128946640856,0.8777245887424803,0.3386863,"regulation of transcription by RNA polymerase II"),
                     c("GO:0045787","positive regulation of cell cycle",0.12797558495406453,3.8002472465724964,0.8732786845066448,0.35864381,"regulation of transcription by RNA polymerase II"),
                     c("GO:0048519","negative regulation of biological process",4.664847168934615,4.216906681520611,0.9015940760367622,0.14346073,"regulation of transcription by RNA polymerase II"),
                     c("GO:0048522","positive regulation of cellular process",3.6764170538087693,6.262151359307686,0.8743151121957952,0.3386863,"regulation of transcription by RNA polymerase II"),
                     c("GO:0050678","regulation of epithelial cell proliferation",0.059951173649385486,3.1232998933171747,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051128","regulation of cellular component organization",1.4586420587819326,2.7805341793634004,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051246","regulation of protein metabolic process",2.477140061467859,5.002136016581367,0.8972254088379591,0.21066885,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051248","negative regulation of protein metabolic process",0.5470258463226723,3.0579846021777293,0.8939074254684327,0.21066885,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051302","regulation of cell division",0.16835161864230688,6.971879492244718,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051493","regulation of cytoskeleton organization",0.47119894618650526,4.819432382818121,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051726","regulation of cell cycle",0.9785034787039558,3.6712487085790935,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0051781","positive regulation of cell division",0.05264095540663495,4.428505991371453,0.8732786845066448,0.3386863,"regulation of transcription by RNA polymerase II"),
                     c("GO:0060070","canonical Wnt signaling pathway",0.10527206540149517,7.025973668578032,0.8382323400382036,0.43969377,"regulation of transcription by RNA polymerase II"),
                     c("GO:0060632","regulation of microtubule-based movement",0.032425863680133206,4.915427986389654,0.8975862302209535,0.15132909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0062197","cellular response to chemical stress",0.37857577356261785,3.173099634013462,0.9284510630574763,0.21493909,"regulation of transcription by RNA polymerase II"),
                     c("GO:0080090","regulation of primary metabolic process",14.365627383358822,14.31409986949664,0.8975777011094312,0.20062793,"regulation of transcription by RNA polymerase II"),
                     c("GO:0090092","regulation of transmembrane receptor protein serine/threonine kinase signaling pathway",0.09221212668229166,2.789922693877269,0.8576380811540493,0.44192018,"regulation of transcription by RNA polymerase II"),
                     c("GO:0099177","regulation of trans-synaptic signaling",0.24673094178107735,6.812755209344133,0.8591506625169564,0.42687478,"regulation of transcription by RNA polymerase II"),
                     c("GO:0099536","synaptic signaling",0.4961324515060551,3.5565342328990064,0.8505761412745563,0.25923296,"regulation of transcription by RNA polymerase II"),
                     c("GO:0141124","intracellular signaling cassette",0.7898358528644561,3.4741270192363674,0.8429293473778766,0.3394991,"regulation of transcription by RNA polymerase II"),
                     c("GO:1902531","regulation of intracellular signal transduction",0.9735955409342439,3.7953708250478058,0.8576380811540493,0.44192018,"regulation of transcription by RNA polymerase II"),
                     c("GO:1902533","positive regulation of intracellular signal transduction",0.412717200244501,5.463944396418106,0.8376585836470826,0.42722693,"regulation of transcription by RNA polymerase II"),
                     c("GO:0006412","translation",4.432288697403288,9.416629642092984,0.9228784258553101,0.03189197,"translation"),
                     c("GO:0001522","pseudouridine synthesis",0.27709911440028817,3.4028215192017877,0.9551904316645192,0.24781368,"translation"),
                     c("GO:0006022","aminoglycan metabolic process",1.0723474823879122,4.568992705435275,0.964071399758336,0.11926201,"translation"),
                     c("GO:0006040","amino sugar metabolic process",0.3589095635425583,4.412820525606365,0.9739612825912263,0.29472003,"translation"),
                     c("GO:0006468","protein phosphorylation",3.756826993126081,6.023394931013309,0.9141041258867636,0.18083678,"translation"),
                     c("GO:0006508","proteolysis",5.329242630377004,4.344498100793902,0.9617149585994127,0.18083678,"translation"),
                     c("GO:0006644","phospholipid metabolic process",1.7199195964602003,4.532208304635846,0.9291221054374792,0.22444672,"translation"),
                     c("GO:0009059","macromolecule biosynthetic process",16.264263355707197,5.3953873667529955,0.9287059262123749,0.15753374,"translation"),
                     c("GO:0016071","mRNA metabolic process",1.6631731043434954,3.2685937241553216,0.9645196551794198,0.27280945,"translation"),
                     c("GO:0016310","phosphorylation",6.253736641437635,4.437401249230823,0.9351943037281288,0.22444672,"translation"),
                     c("GO:0018193","peptidyl-amino acid modification",2.338036699855386,5.044606329041383,0.954504764692118,0.27057677,"translation"),
                     c("GO:0018200","peptidyl-glutamic acid modification",0.04871263610851042,3.7803393419070104,0.9520239184151046,0.27057677,"translation"),
                     c("GO:0018208","peptidyl-proline modification",0.22083750881348954,3.4632495862768566,0.9520239184151046,0.40778717,"translation"),
                     c("GO:0018958","phenol-containing compound metabolic process",0.1696069086435873,2.894673428288157,0.9711947540280369,0.44620405,"translation"),
                     c("GO:0044249","cellular biosynthetic process",23.43500411119782,4.338427304901265,0.9369698033947731,0.13625948,"translation"),
                     c("GO:0044255","cellular lipid metabolic process",4.451927832540966,3.5384332460958055,0.933553555536321,0.30582377,"translation"),
                     c("GO:1901617","organic hydroxy compound biosynthetic process",0.8333107299088005,4.154810173367419,0.9669684162083689,0.13625948,"translation"),
                     c("GO:0006629","lipid metabolic process",5.980051423570241,3.7258376037878267,0.9747643773547617,0.0781618,"lipid metabolic process"),
                     c("GO:0007017","microtubule-based process",1.352798959497502,15.831177313993082,0.9614007662692289,0.03189197,"microtubule-based process"),
                     c("GO:0007018","microtubule-based movement",0.5954751176662081,13.09331022186228,0.8764224700228849,0.03189197,"microtubule-based movement"),
                     c("GO:0000226","microtubule cytoskeleton organization",0.8833918782540008,3.9200300937722097,0.8237442985573106,0.46719168,"microtubule-based movement"),
                     c("GO:0003341","cilium movement",0.13264477148823883,NaN,0.8702551335683347,0.46719168,"microtubule-based movement"),
                     c("GO:0098534","centriole assembly",0.03715658403789972,3.6981064494176685,0.8124327826989505,0.46719168,"microtubule-based movement"),
                     c("GO:0007155","cell adhesion",1.0284049482842672,3.990701227732296,0.9614007662692289,0.03189197,"cell adhesion"),
                     c("GO:0008152","metabolic process",58.08643542977757,5.01957980696974,1,-0,"metabolic process"),
                     c("GO:0009058","biosynthetic process",28.507692540195496,4.389775824914547,0.9782743675028893,0.05898132,"biosynthetic process"),
                     c("GO:0030031","cell projection assembly",0.5724417768191847,24.48026877402233,0.8586277224721697,0,"cell projection assembly"),
                     c("GO:0006996","organelle organization",6.55438598080312,8.986939498178593,0.8917789664015157,0.22720981,"cell projection assembly"),
                     c("GO:0007010","cytoskeleton organization",1.8978237258769566,4.628940646043917,0.8762386545373954,0.22720981,"cell projection assembly"),
                     c("GO:0016050","vesicle organization",0.4876752427915464,3.694163996413555,0.8762386545373954,0.295867,"cell projection assembly"),
                     c("GO:0022607","cellular component assembly",4.343086805371125,13.172758538920963,0.8829120579384916,0.29074743,"cell projection assembly"),
                     c("GO:0030030","cell projection organization",1.214633373356573,20.57679028522921,0.8917789664015157,0.22720981,"cell projection assembly"),
                     c("GO:0045010","actin nucleation",0.07539616337102174,3.4447029706979158,0.8622982269319094,0.295867,"cell projection assembly"),
                     c("GO:0048284","organelle fusion",0.2981486047746998,4.101799639086141,0.8762386545373954,0.295867,"cell projection assembly"),
                     c("GO:0061024","membrane organization",1.185564795091629,4.291057901792766,0.8917789664015157,0.22720981,"cell projection assembly"),
                     c("GO:0061025","membrane fusion",0.29388554147623386,3.333775755501191,0.8917789664015157,0.22720981,"cell projection assembly"),
                     c("GO:0065003","protein-containing complex assembly",2.452194249383591,3.250703976539501,0.875263247327097,0.34055037,"cell projection assembly"),
                     c("GO:0070925","organelle assembly",1.5672344893044614,17.45688132323049,0.8673575304407358,0.34055037,"cell projection assembly"),
                     c("GO:1905515","non-motile cilium assembly",0.05417437829055198,7.852745962667691,0.8400409356204637,0.4788887,"cell projection assembly"),
                     c("GO:0044237","cellular metabolic process",39.640441131550894,5.038423227920164,0.9441135709925537,0.05898132,"cellular metabolic process"),
                     c("GO:0045165","cell fate commitment",0.1569284796306553,7.817085249184954,0.9485672953969492,0.03189197,"cell fate commitment"),
                     c("GO:0030154","cell differentiation",2.06266791657451,3.482188484569003,0.9485672953969492,0.41875931,"cell fate commitment"),
                     c("GO:0048869","cellular developmental process",2.113307792037927,4.772321153234118,0.9497323794174287,0.35418145,"cell fate commitment"),
                     c("GO:0048870","cell motility",0.8376747086779577,3.480040096279695,0.9614007662692289,0.03189197,"cell motility"),
                     c("GO:0051649","establishment of localization in cell",3.464570867298568,12.241285863177106,0.9181921924214835,0.03189197,"establishment of localization in cell"),
                     c("GO:0006886","intracellular protein transport",1.529423185383541,6.024584042645504,0.9112768323951745,0.31269958,"establishment of localization in cell"),
                     c("GO:0006888","endoplasmic reticulum to Golgi vesicle-mediated transport",0.416106483247958,3.1129120556493373,0.9119593999833391,0.38606433,"establishment of localization in cell"),
                     c("GO:0035725","sodium ion transmembrane transport",0.5221267999443341,3.5409110649950373,0.9344716733688867,0.18012579,"establishment of localization in cell"),
                     c("GO:0051641","cellular localization",5.613085851842994,16.682480501879215,0.9352747558612293,0.1767677,"establishment of localization in cell"),
                     c("GO:0051905","establishment of pigment granule localization",0.005405131064336763,3.546443862289474,0.9135499760717493,0.36508764,"establishment of localization in cell"),
                     c("GO:0055085","transmembrane transport",12.953593503918547,16.385761989535613,0.9344716733688867,0.1831236,"establishment of localization in cell"),
                     c("GO:0061512","protein localization to cilium",0.06874804907012301,5.900358205370607,0.9176264329319752,0.35130364,"establishment of localization in cell"),
                     c("GO:0060285","cilium-dependent cell motility",0.09471286127307771,5.942327795020085,0.9419101596908703,0.03189197,"cilium-dependent cell motility"),
                     c("GO:0071840","cellular component organization or biogenesis",14.26840886478907,12.905763149108259,0.9614007662692289,0.03189197,"cellular component organization or biogenesis"),
                     c("GO:0098609","cell-cell adhesion",0.542898157436109,4.296117447994224,0.9529748085167895,0.03189197,"cell-cell adhesion"),
                     c("GO:0098742","cell-cell adhesion via plasma-membrane adhesion molecules",0.2759619693403047,4.432309017687842,0.9529748085167895,0.49695896,"cell-cell adhesion"),
                     c("GO:1901615","organic hydroxy compound metabolic process",1.6412842926152864,3.112563656903549,0.9782743675028893,0.05898132,"organic hydroxy compound metabolic process"));

stuff        <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value          <- as.numeric( 1 );
stuff$frequency      <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness     <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Revigo_TreeMap_F14R_negative.tiff", width = 12, height = 8, units = 'in', res = 600)

treemap(
  stuff,
  index = c("representative", "description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  bg.labels = "#CCCCCCAA",
  position.legend = "none")

dev.off()


##########

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006139","nucleobase-containing compound metabolic process",18.39314412781983,4.755596088238773,0.9523491139099532,0.0781618,"nucleobase-containing compound metabolic process"),
                     c("GO:0006887","exocytosis",0.35007084512177816,4.174947218745009,0.9597063959175196,0.03189197,"exocytosis"),
                     c("GO:0016192","vesicle-mediated transport",2.591260690701909,3.1679093732414096,0.9597063959175196,0.1831236,"exocytosis"),
                     c("GO:0055085","transmembrane transport",12.953593503918547,3.3274587843765118,0.9597063959175196,0.1831236,"exocytosis"),
                     c("GO:0007015","actin filament organization",0.47463991760177976,3.4391452152347273,0.962938792770451,0.03189197,"actin filament organization"),
                     c("GO:0034330","cell junction organization",0.326665839980258,3.0889280975624045,0.962938792770451,0.22720981,"actin filament organization"),
                     c("GO:0007155","cell adhesion",1.0284049482842672,3.775221437307347,0.9688451853946736,0.03189197,"cell adhesion"),
                     c("GO:0007186","G protein-coupled receptor signaling pathway",1.464221945905271,14.635247152692484,0.8135252214348936,0,"G protein-coupled receptor signaling pathway"),
                     c("GO:0001817","regulation of cytokine production",0.14516567391277488,5.362625960804394,0.8616982117691914,0.15132909,"G protein-coupled receptor signaling pathway"),
                     c("GO:0006357","regulation of transcription by RNA polymerase II",4.1876991326832185,4.260458928561313,0.8606469217178796,0.21476949,"G protein-coupled receptor signaling pathway"),
                     c("GO:0006974","DNA damage response",3.1870828591331772,3.7575007846441904,0.9011168705611772,0.21493909,"G protein-coupled receptor signaling pathway"),
                     c("GO:0007165","signal transduction",8.851839893322994,24.665466525433146,0.814526780949329,0.25923296,"G protein-coupled receptor signaling pathway"),
                     c("GO:0007166","cell surface receptor signaling pathway",1.7427092633657988,14.285775708689929,0.8135252214348936,0.26324161,"G protein-coupled receptor signaling pathway"),
                     c("GO:0007167","enzyme-linked receptor protein signaling pathway",0.5096157429315727,6.404022354868767,0.7940551752992367,0.26324161,"G protein-coupled receptor signaling pathway"),
                     c("GO:0007169","cell surface receptor protein tyrosine kinase signaling pathway",0.36239976201670654,5.572009070303073,0.7940551752992367,0.43969377,"G protein-coupled receptor signaling pathway"),
                     c("GO:0007264","small GTPase-mediated signal transduction",0.3940010724607046,4.957803341597308,0.8114102616930879,0.3394991,"G protein-coupled receptor signaling pathway"),
                     c("GO:0008063","Toll signaling pathway",0.00457319376937054,3.4682011955900247,0.7940551752992367,0.43969377,"G protein-coupled receptor signaling pathway"),
                     c("GO:0009891","positive regulation of biosynthetic process",2.016283720350727,3.6901841019759165,0.8512939973542935,0.3386863,"G protein-coupled receptor signaling pathway"),
                     c("GO:0009892","negative regulation of metabolic process",2.9251604469836634,3.23255510441052,0.8578826166968443,0.33279065,"G protein-coupled receptor signaling pathway"),
                     c("GO:0009893","positive regulation of metabolic process",2.881710183468756,6.057505044098711,0.8519929729666075,0.3386863,"G protein-coupled receptor signaling pathway"),
                     c("GO:0019219","regulation of nucleobase-containing compound metabolic process",11.820748427527752,22.82028343749162,0.8664940443641042,0.20062793,"G protein-coupled receptor signaling pathway"),
                     c("GO:0031324","negative regulation of cellular metabolic process",2.62213590202752,4.273421873973661,0.8534317592845057,0.33279065,"G protein-coupled receptor signaling pathway"),
                     c("GO:0033554","cellular response to stress",4.425081855984172,3.249157384730851,0.9011168705611772,0.29516635,"G protein-coupled receptor signaling pathway"),
                     c("GO:0035556","intracellular signal transduction",4.385343312708344,9.11371177781204,0.8135252214348936,0.26324161,"G protein-coupled receptor signaling pathway"),
                     c("GO:0043068","positive regulation of programmed cell death",0.14822759697472157,3.4708480616810036,0.8586784328628402,0.3386863,"G protein-coupled receptor signaling pathway"),
                     c("GO:0048519","negative regulation of biological process",4.664847168934615,5.204309254531074,0.880222118020537,0.14346073,"G protein-coupled receptor signaling pathway"),
                     c("GO:0048522","positive regulation of cellular process",3.6764170538087693,7.698307859930513,0.8586784328628402,0.15132909,"G protein-coupled receptor signaling pathway"),
                     c("GO:0048523","negative regulation of cellular process",4.3068778422165455,5.998565687016824,0.8645910992918725,0.15132909,"G protein-coupled receptor signaling pathway"),
                     c("GO:0051056","regulation of small GTPase mediated signal transduction",0.24989870301960257,12.623163807511558,0.8322397126606976,0.15132909,"G protein-coupled receptor signaling pathway"),
                     c("GO:0071526","semaphorin-plexin signaling pathway",0.0656516670669647,4.525873763651901,0.7940551752992367,0.43969377,"G protein-coupled receptor signaling pathway"),
                     c("GO:0080090","regulation of primary metabolic process",14.365627383358822,20.083344110418352,0.8670693186527986,0.20062793,"G protein-coupled receptor signaling pathway"),
                     c("GO:0080135","regulation of cellular response to stress",0.17710418971005804,3.8615274822252226,0.8333993235918875,0.4106238,"G protein-coupled receptor signaling pathway"),
                     c("GO:0097696","cell surface receptor signaling pathway via STAT",0.025120568143270038,3.1361271988288553,0.7940551752992367,0.43969377,"G protein-coupled receptor signaling pathway"),
                     c("GO:0141124","intracellular signaling cassette",0.7898358528644561,5.948321528997969,0.8114102616930882,0.26324161,"G protein-coupled receptor signaling pathway"),
                     c("GO:1902531","regulation of intracellular signal transduction",0.9735955409342439,6.36744593658207,0.8356669294828272,0.44192018,"G protein-coupled receptor signaling pathway"),
                     c("GO:0008152","metabolic process",58.08643542977757,25.116555748956394,1,-0,"metabolic process"),
                     c("GO:0016311","dephosphorylation",0.8318929906132368,5.410270374892913,0.9049615572929474,0.05898132,"dephosphorylation"),
                     c("GO:0006793","phosphorus metabolic process",12.714975181439863,12.56870640129184,0.926866584877115,0.10046538,"dephosphorylation"),
                     c("GO:0009187","cyclic nucleotide metabolic process",0.160918332752372,4.782880522970748,0.8888301330325752,0.22444672,"dephosphorylation"),
                     c("GO:0016310","phosphorylation",6.253736641437635,11.326112765063167,0.9049615572929474,0.22444672,"dephosphorylation"),
                     c("GO:0052652","cyclic purine nucleotide metabolic process",0.122159407948132,3.473265926609984,0.8888301330325752,0.36823636,"dephosphorylation"),
                     c("GO:0030029","actin filament-based process",0.8043086081733359,3.4157681263508337,0.9688451853946736,0.03189197,"actin filament-based process"),
                     c("GO:0032259","methylation",2.6585910004764686,3.348860070252668,0.9586949445685589,0.05898132,"methylation"),
                     c("GO:0043170","macromolecule metabolic process",33.33895178165279,23.650136567770026,0.9586949445685589,0.05898132,"macromolecule metabolic process"),
                     c("GO:0043687","post-translational protein modification",1.8860683042179072,12.5730968976993,0.9036275513457297,-0,"post-translational protein modification"),
                     c("GO:0001510","RNA methylation",0.712861962256531,4.114552255511823,0.8891736788974962,0.24781368,"post-translational protein modification"),
                     c("GO:0006396","RNA processing",4.402491558667012,3.0509892878530747,0.8946343075852659,0.27280945,"post-translational protein modification"),
                     c("GO:0006468","protein phosphorylation",3.756826993126081,10.300514136453048,0.8631857027684412,0.27057677,"post-translational protein modification"),
                     c("GO:0006470","protein dephosphorylation",0.49314683038536267,4.391022340694499,0.8560706428343213,0.27057677,"post-translational protein modification"),
                     c("GO:0009057","macromolecule catabolic process",3.5074771718129205,6.057003313000018,0.9250428080498713,0.11926201,"post-translational protein modification"),
                     c("GO:0016070","RNA metabolic process",8.105171248384922,4.277010382068444,0.9213868165569002,0.11926201,"post-translational protein modification"),
                     c("GO:0018212","peptidyl-tyrosine modification",0.03475922627074853,6.765000985480657,0.9036275513457297,0.27057677,"post-translational protein modification"),
                     c("GO:0043412","macromolecule modification",10.203408176348653,24.548340774993598,0.9352312473397119,0.11926201,"post-translational protein modification"),
                     c("GO:0043414","macromolecule methylation",0.926246494356536,3.6029318006708384,0.9043686809181112,0.39383712,"post-translational protein modification"),
                     c("GO:0043632","modification-dependent macromolecule catabolic process",1.2894658869034905,3.7078162046083953,0.9250428080498713,0.2607896,"post-translational protein modification"),
                     c("GO:0070646","protein modification by small protein removal",0.39107698516360445,3.510113170713573,0.9036275513457297,0.27057677,"post-translational protein modification"),
                     c("GO:1901658","glycosyl compound catabolic process",0.16646129958155523,3.400990469046989,0.9395297963724982,0.29472003,"post-translational protein modification"),
                     c("GO:1903510","mucopolysaccharide metabolic process",0.06785211659862092,5.1911769548451305,0.9282297337231078,0.11926201,"post-translational protein modification"),
                     c("GO:0044237","cellular metabolic process",39.640441131550894,9.044098291153489,0.9367236265475163,0.05898132,"cellular metabolic process"),
                     c("GO:0048869","cellular developmental process",2.113307792037927,5.386282952032796,0.9591088678659405,0.03189197,"cellular developmental process"),
                     c("GO:0048468","cell development",1.0529397144269397,4.5342315282675525,0.9591088678659405,0.35418145,"cellular developmental process"),
                     c("GO:0061919","process utilizing autophagic mechanism",0.4385663288591023,3.180554108933776,0.9688451853946736,0.03189197,"process utilizing autophagic mechanism"));

stuff        <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value          <- as.numeric( 1 );
stuff$frequency      <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness     <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Revigo_TreeMap_H2D_positive.tiff", width = 12, height = 8, units = 'in', res = 600)

treemap(
  stuff,
  index = c("representative", "description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  bg.labels = "#CCCCCCAA",
  position.legend = "none")

dev.off()


###############

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0003341","cilium movement",0.13264477148823883,8.900902116125637,0.8982089726440059,0.03189197,"cilium movement"),
                     c("GO:0007018","microtubule-based movement",0.5954751176662081,8.042430571895673,0.9019656186292613,0.46719168,"cilium movement"),
                     c("GO:0005975","carbohydrate metabolic process",5.177832504046095,4.062661295219734,0.921765237180423,0.0781618,"carbohydrate metabolic process"),
                     c("GO:0006412","translation",4.432288697403288,14.746412498584776,0.8597585747518881,0,"translation"),
                     c("GO:0006022","aminoglycan metabolic process",1.0723474823879122,7.492731307323866,0.8994429071425116,0.29472003,"translation"),
                     c("GO:0006040","amino sugar metabolic process",0.3589095635425583,10.572275115311026,0.9053813749211701,0.29472003,"translation"),
                     c("GO:0006081","cellular aldehyde metabolic process",0.5676347545201639,3.984368920511328,0.8857077208933207,0.10046538,"translation"),
                     c("GO:0006082","organic acid metabolic process",8.445551746967404,18.065533697854153,0.8436167888283749,0.2074912,"translation"),
                     c("GO:0006091","generation of precursor metabolites and energy",2.4286809447125486,10.749957883212996,0.8857077208933207,0.10046538,"translation"),
                     c("GO:0006099","tricarboxylic acid cycle",0.4936760212682554,5.598471938694151,0.8587885093511024,0.46826111,"translation"),
                     c("GO:0006457","protein folding",1.1925033490398829,3.203841957040198,0.8597585747518881,0.2212988,"translation"),
                     c("GO:0006508","proteolysis",5.329242630377004,4.076623081102472,0.9127704967248759,0.18083678,"translation"),
                     c("GO:0006525","arginine metabolic process",0.40362742382346467,2.9320080016375294,0.8229619954968694,0.35100309,"translation"),
                     c("GO:0006575","cellular modified amino acid metabolic process",1.0867980855202988,4.143526125705428,0.8857077208933207,0.10046538,"translation"),
                     c("GO:0006631","fatty acid metabolic process",1.7695330877461002,8.133287118539895,0.8115126116419102,0.27289201,"translation"),
                     c("GO:0006749","glutathione metabolic process",0.32127301568063965,3.0269201660499387,0.8857077208933207,0.10046538,"translation"),
                     c("GO:0006790","sulfur compound metabolic process",2.3822573668416673,4.663269051180058,0.8857077208933207,0.10046538,"translation"),
                     c("GO:0008610","lipid biosynthetic process",3.4846382778484486,4.0061541157885365,0.8829451868326795,0.30582377,"translation"),
                     c("GO:0008652","amino acid biosynthetic process",2.752513767454623,3.0760408210226546,0.8830084131255491,0.31474975,"translation"),
                     c("GO:0009059","macromolecule biosynthetic process",16.264263355707197,6.504885338255428,0.8656086366866137,0.15753374,"translation"),
                     c("GO:0009060","aerobic respiration",1.0799752151603983,9.808753345138143,0.8622151981926233,0.46826111,"translation"),
                     c("GO:0009063","amino acid catabolic process",1.1404334275161834,3.7811017318143305,0.8799157774360129,0.31474975,"translation"),
                     c("GO:0009066","aspartate family amino acid metabolic process",0.9159186574048253,3.9816797773761934,0.8227560488496893,0.31474975,"translation"),
                     c("GO:0009067","aspartate family amino acid biosynthetic process",0.7209376612647682,3.078029201011742,0.8134399017690428,0.36632792,"translation"),
                     c("GO:0009083","branched-chain amino acid catabolic process",0.18455224371765505,3.220257520166444,0.8056962152820935,0.42728573,"translation"),
                     c("GO:0009123","nucleoside monophosphate metabolic process",1.1288059962102055,3.927763950408506,0.8088332151625492,0.33652359,"translation"),
                     c("GO:0009132","nucleoside diphosphate metabolic process",0.5593006134528395,11.432509957724589,0.8088332151625492,0.33652359,"translation"),
                     c("GO:0009141","nucleoside triphosphate metabolic process",1.1643627008347086,11.081642343643683,0.8088332151625492,0.33652359,"translation"),
                     c("GO:0009142","nucleoside triphosphate biosynthetic process",0.5874806433051126,12.696133078357022,0.7944742250188941,0.13625948,"translation"),
                     c("GO:0009145","purine nucleoside triphosphate biosynthetic process",0.43565947103260794,13.124927129213289,0.7934937488891489,0.48347793,"translation"),
                     c("GO:0009152","purine ribonucleotide biosynthetic process",1.0720988857405997,5.647973182591526,0.7849447536144242,0.38667606,"translation"),
                     c("GO:0009161","ribonucleoside monophosphate metabolic process",0.8874924922581835,2.8063584430878863,0.8068316113345374,0.48684517,"translation"),
                     c("GO:0009208","pyrimidine ribonucleoside triphosphate metabolic process",0.13881538331806229,2.8637481139446934,0.8027192555601158,0.48347793,"translation"),
                     c("GO:0009218","pyrimidine ribonucleotide metabolic process",0.35888987271900885,4.128324103284932,0.7875247731819812,0.42641897,"translation"),
                     c("GO:0009226","nucleotide-sugar biosynthetic process",0.2857384632326297,3.258199995926268,0.7875233005107019,0.38667606,"translation"),
                     c("GO:0015980","energy derivation by oxidation of organic compounds",1.5167964447824265,5.3652120997725605,0.8646975078769747,0.40365738,"translation"),
                     c("GO:0016042","lipid catabolic process",1.1775407844952093,3.8775383405168813,0.8806601938293377,0.30582377,"translation"),
                     c("GO:0016052","carbohydrate catabolic process",1.673968598354507,5.007510940546821,0.8971887529161809,0.2607896,"translation"),
                     c("GO:0016053","organic acid biosynthetic process",4.243739216505085,5.803187400884035,0.8257302052769059,0.26834299,"translation"),
                     c("GO:0018193","peptidyl-amino acid modification",2.338036699855386,4.39729828580048,0.9102319172691189,0.18083678,"translation"),
                     c("GO:0018200","peptidyl-glutamic acid modification",0.04871263610851042,3.625938567338854,0.9082935864630958,0.40778717,"translation"),
                     c("GO:0018208","peptidyl-proline modification",0.22083750881348954,3.8723776056095978,0.9082935864630958,0.27057677,"translation"),
                     c("GO:0019637","organophosphate metabolic process",6.773690066732693,10.113220614261117,0.8583812535729545,0.22392112,"translation"),
                     c("GO:0019693","ribose phosphate metabolic process",2.078218744472724,10.752662478280435,0.8314430233824718,0.29229367,"translation"),
                     c("GO:0022900","electron transport chain",0.8110773287684753,12.758993002935902,0.8646975078769747,0.10046538,"translation"),
                     c("GO:0022904","respiratory electron transport chain",0.6574322939646986,10.408857746984792,0.8622151981926233,0.40365738,"translation"),
                     c("GO:0030258","lipid modification",0.6585497482011325,4.483345879131648,0.863234254316889,0.33786297,"translation"),
                     c("GO:0032787","monocarboxylic acid metabolic process",3.2972997825985404,12.528904409886332,0.8356881068685523,0.2074912,"translation"),
                     c("GO:0034440","lipid oxidation",0.4393564231540259,6.897167814330846,0.863234254316889,0.33786297,"translation"),
                     c("GO:0034654","nucleobase-containing compound biosynthetic process",10.209893841355267,2.7016347120677167,0.8870317052018232,0.18383601,"translation"),
                     c("GO:0044242","cellular lipid catabolic process",0.8441899099198973,5.342597143708266,0.840988252859924,0.33786297,"translation"),
                     c("GO:0044248","cellular catabolic process",3.339716277877099,5.886847690482252,0.8628472810042849,0.2607896,"translation"),
                     c("GO:0044249","cellular biosynthetic process",23.43500411119782,10.405384605414426,0.8712382917345722,0.13625948,"translation"),
                     c("GO:0044255","cellular lipid metabolic process",4.451927832540966,7.358162692801412,0.8645577723359904,0.30582377,"translation"),
                     c("GO:0044282","small molecule catabolic process",2.800286166738645,10.436156864376496,0.8547544103685616,0.2074912,"translation"),
                     c("GO:0044283","small molecule biosynthetic process",6.191735160786158,5.863633821672078,0.8639447655882346,0.2074912,"translation"),
                     c("GO:0045333","cellular respiration",1.3395396511898598,8.829456725910692,0.862586465868802,0.45476817,"translation"),
                     c("GO:0046434","organophosphate catabolic process",0.8689954248863754,4.17732421154897,0.8220450320297776,0.29229367,"translation"),
                     c("GO:0046940","nucleoside monophosphate phosphorylation",0.19151048848945834,4.480011584039445,0.7957970675805766,0.41586979,"translation"),
                     c("GO:0055086","nucleobase-containing small molecule metabolic process",5.031190018367107,11.482958200711197,0.8691847521083499,0.2074912,"translation"),
                     c("GO:0072329","monocarboxylic acid catabolic process",0.7868330022731579,7.103195479650423,0.8114973961004885,0.37046039,"translation"),
                     c("GO:0072330","monocarboxylic acid biosynthetic process",1.1203143285544854,3.6824793872604986,0.8221189612004854,0.37046039,"translation"),
                     c("GO:0090407","organophosphate biosynthetic process",4.68253445118795,3.684620135848789,0.8336214298834116,0.29229367,"translation"),
                     c("GO:1901136","carbohydrate derivative catabolic process",1.245173840681842,4.171369693661025,0.8799592256646781,0.29472003,"translation"),
                     c("GO:1901137","carbohydrate derivative biosynthetic process",4.348346716611784,4.188749418460256,0.8858791805134951,0.29472003,"translation"),
                     c("GO:1901293","nucleoside phosphate biosynthetic process",2.839852415308414,5.50793792855431,0.8023017560389937,0.33652359,"translation"),
                     c("GO:1901605","alpha-amino acid metabolic process",3.9444608507292327,6.628403604414698,0.8239367440609553,0.27289201,"translation"),
                     c("GO:0006629","lipid metabolic process",5.980051423570241,10.005908597374901,0.921765237180423,0.0781618,"lipid metabolic process"),
                     c("GO:0007017","microtubule-based process",1.352798959497502,9.519122516600552,0.9579294153106881,0.03189197,"microtubule-based process"),
                     c("GO:0008152","metabolic process",58.08643542977757,13.926320160901598,1,-0,"metabolic process"),
                     c("GO:0009056","catabolic process",9.05402280817628,5.737503234793042,0.9321612972451532,0.05898132,"catabolic process"),
                     c("GO:0009058","biosynthetic process",28.507692540195496,20.531511671620787,0.9321612972451532,0.05898132,"biosynthetic process"),
                     c("GO:0018958","phenol-containing compound metabolic process",0.1696069086435873,8.721433919969968,0.9211169398876139,0.05898132,"phenol-containing compound metabolic process"),
                     c("GO:0006066","alcohol metabolic process",1.085552640930793,2.9414369038281047,0.8649863016688212,0.44620405,"phenol-containing compound metabolic process"),
                     c("GO:1901617","organic hydroxy compound biosynthetic process",0.8333107299088005,5.098489460718436,0.898187083177686,0.44620405,"phenol-containing compound metabolic process"),
                     c("GO:0030031","cell projection assembly",0.5724417768191847,12.814513523070655,0.9105772350313577,0.03189197,"cell projection assembly"),
                     c("GO:0017004","cytochrome complex assembly",0.2808428322276361,2.7181174265506725,0.9163008184477747,0.40261128,"cell projection assembly"),
                     c("GO:0022607","cellular component assembly",4.343086805371125,5.475089488822062,0.9216000096300085,0.29074743,"cell projection assembly"),
                     c("GO:0030030","cell projection organization",1.214633373356573,10.846127382988803,0.9270298613388435,0.22720981,"cell projection assembly"),
                     c("GO:0033108","mitochondrial respiratory chain complex assembly",0.23663693335901673,4.24682343639792,0.9124103708398202,0.34055037,"cell projection assembly"),
                     c("GO:0070286","axonemal dynein complex assembly",0.04858464575543869,5.7485372971780695,0.8676301384092949,0.4788887,"cell projection assembly"),
                     c("GO:0070925","organelle assembly",1.5672344893044614,8.421786963384166,0.9141676810472618,0.34055037,"cell projection assembly"),
                     c("GO:0032886","regulation of microtubule-based process",0.17552892382609833,7.124100000867983,0.9584277881002151,-0,"regulation of microtubule-based process"),
                     c("GO:0006357","regulation of transcription by RNA polymerase II",4.1876991326832185,4.613381793941384,0.9549373105179325,0.21066885,"regulation of microtubule-based process"),
                     c("GO:0007165","signal transduction",8.851839893322994,6.824281768656675,0.900468170644192,0.25923296,"regulation of microtubule-based process"),
                     c("GO:0007186","G protein-coupled receptor signaling pathway",1.464221945905271,3.349853555957328,0.9001864199900058,0.26324161,"regulation of microtubule-based process"),
                     c("GO:0007218","neuropeptide signaling pathway",0.15345304927416917,6.479995154880226,0.9001864199900058,0.26324161,"regulation of microtubule-based process"),
                     c("GO:0008063","Toll signaling pathway",0.00457319376937054,3.0416079183385336,0.8952328225490553,0.43969377,"regulation of microtubule-based process"),
                     c("GO:0016055","Wnt signaling pathway",0.19157448366599422,6.910349487753496,0.8952328225490553,0.15132909,"regulation of microtubule-based process"),
                     c("GO:0019722","calcium-mediated signaling",0.1524365105084264,3.972203573729965,0.9001864199900058,0.26324161,"regulation of microtubule-based process"),
                     c("GO:0023056","positive regulation of signaling",0.6428118574791976,5.142097201594627,0.9436408419283313,0.42687478,"regulation of microtubule-based process"),
                     c("GO:0031323","regulation of cellular metabolic process",14.72286076554673,2.7664368534102697,0.9552280608241247,0.20062793,"regulation of microtubule-based process"),
                     c("GO:0043266","regulation of potassium ion transport",0.06280634306406248,4.67109797441867,0.960026129702478,0.14346073,"regulation of microtubule-based process"),
                     c("GO:0051246","regulation of protein metabolic process",2.477140061467859,5.5367271976848285,0.9564187367716392,0.14346073,"regulation of microtubule-based process"),
                     c("GO:0051248","negative regulation of protein metabolic process",0.5470258463226723,3.1829496068045917,0.9564187367716392,0.21066885,"regulation of microtubule-based process"),
                     c("GO:0060070","canonical Wnt signaling pathway",0.10527206540149517,5.616319664799205,0.8952328225490553,0.43969377,"regulation of microtubule-based process"),
                     c("GO:0060632","regulation of microtubule-based movement",0.032425863680133206,7.365752134749796,0.9584277881002151,0.15132909,"regulation of microtubule-based process"),
                     c("GO:0062197","cellular response to chemical stress",0.37857577356261785,6.064935582596512,0.9268975166601516,0.21493909,"regulation of microtubule-based process"),
                     c("GO:0070887","cellular response to chemical stimulus",2.1752846591599666,5.323041009770405,0.9279541417576689,0.34156074,"regulation of microtubule-based process"),
                     c("GO:0071470","cellular response to osmotic stress",0.03679722650812141,6.1549669611925655,0.9268975166601516,0.41562098,"regulation of microtubule-based process"),
                     c("GO:0080090","regulation of primary metabolic process",14.365627383358822,5.330822520464328,0.9567097124986228,0.20062793,"regulation of microtubule-based process"),
                     c("GO:0099177","regulation of trans-synaptic signaling",0.24673094178107735,5.29715484917911,0.9423804643032545,0.42722693,"regulation of microtubule-based process"),
                     c("GO:0104004","cellular response to environmental stimulus",0.15989933263368555,5.005019438414177,0.9315717379996121,0.21493909,"regulation of microtubule-based process"),
                     c("GO:1902533","positive regulation of intracellular signal transduction",0.412717200244501,5.859549171850739,0.9402049566939823,0.15132909,"regulation of microtubule-based process"),
                     c("GO:0043603","amide metabolic process",1.8259719107448442,4.409725434110874,0.9321612972451532,0.05898132,"amide metabolic process"),
                     c("GO:0044237","cellular metabolic process",39.640441131550894,17.15803775200746,0.9131007609732111,0.05898132,"cellular metabolic process"),
                     c("GO:0044281","small molecule metabolic process",14.792283225323425,28.720622864428247,0.9321612972451532,0.05898132,"small molecule metabolic process"),
                     c("GO:0045165","cell fate commitment",0.1569284796306553,6.6131450162573495,0.9476743544034729,0.03189197,"cell fate commitment"),
                     c("GO:0030154","cell differentiation",2.06266791657451,3.3990475572919383,0.9476743544034729,0.41875931,"cell fate commitment"),
                     c("GO:0048869","cellular developmental process",2.113307792037927,3.8276824506552036,0.9486059257412051,0.35418145,"cell fate commitment"),
                     c("GO:0048870","cell motility",0.8376747086779577,3.005729504803782,0.9579294153106881,0.03189197,"cell motility"),
                     c("GO:0060285","cilium-dependent cell motility",0.09471286127307771,7.9473325518807805,0.9423492821382274,0.03189197,"cilium-dependent cell motility"),
                     c("GO:0061640","cytoskeleton-dependent cytokinesis",0.16478019552101697,3.774351637989755,0.9458650177181279,0.03189197,"cytoskeleton-dependent cytokinesis"),
                     c("GO:0000910","cytokinesis",0.33092398057283656,2.932771107077871,0.9458650177181279,0.46520082,"cytoskeleton-dependent cytokinesis"),
                     c("GO:0071840","cellular component organization or biogenesis",14.26840886478907,3.723001520735392,0.9579294153106881,0.03189197,"cellular component organization or biogenesis"),
                     c("GO:0072524","pyridine-containing compound metabolic process",1.1765094776118044,4.374057615146443,0.9321612972451532,0.05898132,"pyridine-containing compound metabolic process"),
                     c("GO:0072527","pyrimidine-containing compound metabolic process",0.8949011646186815,3.18764650423773,0.9321612972451532,0.05898132,"pyrimidine-containing compound metabolic process"),
                     c("GO:1901135","carbohydrate derivative metabolic process",6.623993042050591,14.172857204930027,0.9321612972451532,0.05898132,"carbohydrate derivative metabolic process"),
                     c("GO:1901615","organic hydroxy compound metabolic process",1.6412842926152864,7.569378445797743,0.9321612972451532,0.05898132,"organic hydroxy compound metabolic process"),
                     c("GO:1902600","proton transmembrane transport",1.4650292696708005,7.37141978491058,0.9449954665764069,0.03189197,"proton transmembrane transport"),
                     c("GO:0051641","cellular localization",5.613085851842994,6.7911843428216265,0.9453668322864652,0.1767677,"proton transmembrane transport"),
                     c("GO:0051649","establishment of localization in cell",3.464570867298568,3.7637430504521276,0.9394054491373154,0.18012579,"proton transmembrane transport"),
                     c("GO:0051905","establishment of pigment granule localization",0.005405131064336763,3.1826730631700704,0.9379015761798366,0.31269958,"proton transmembrane transport"),
                     c("GO:0055085","transmembrane transport",12.953593503918547,24.342720064279316,0.9449954665764069,0.1831236,"proton transmembrane transport"));

stuff        <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value          <- as.numeric( 1 );
stuff$frequency      <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness     <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );


tiff("D:/Decicomp/R/MOFA_omics/Data_compacted/Revigo_TreeMap_H2D_negative.tiff", width = 12, height = 8, units = 'in', res = 600)

treemap(
  stuff,
  index = c("representative", "description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = TRUE,
  lowerbound.cex.labels = 0,
  bg.labels = "#CCCCCCAA",
  position.legend = "none")

dev.off()



# Sacar los bloques Ruta metabolica  ----


library(ggplot2)
library(reshape2)
library(patchwork)
library(svglite)

# F14R Family ----
### F14R_C3_T0_vs_C7_T0_FC GO_MWU  ----
F14R_C7_T0_vs_C3_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame() %>% rownames_to_column(var = "gene")  %>% dplyr::select("gene", "log2FoldChange", "padj" ) %>% 
  dplyr::rename(Gene=gene) %>% dplyr::rename(F14R_C7=log2FoldChange, F14R_C7_padj=padj)

F14R_C8_T0_vs_C3_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T0", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame() %>% rownames_to_column(var = "gene")  %>% dplyr::select("gene", "log2FoldChange", "padj" ) %>% 
  dplyr::rename(Gene=gene) %>% dplyr::rename(F14R_C8=log2FoldChange, F14R_C8_padj=padj)

F14R_C3_C7_C8 <- merge(F14R_C7_T0_vs_C3_T0_FC, F14R_C8_T0_vs_C3_T0_FC, by="Gene" )  %>%
  mutate(F14R_C3=0) %>% dplyr::select(Gene, F14R_C3, F14R_C7, F14R_C8)

positive_F14R_transcriptome_clean <- positive_F14R_transcriptome %>% dplyr::select(Gene)
negative_F14R_transcriptome_clean <- negative_F14R_transcriptome %>% dplyr::select(Gene)

List_FC_genes_F14R_age_transcriptome <- rbind(positive_F14R_transcriptome_clean, negative_F14R_transcriptome_clean)  %>% 
                                        left_join(F14R_C3_C7_C8) %>%   pivot_longer(cols = starts_with("F14R"),    names_to = "Condition", values_to = "FoldChange")    

# H2D Family ----
### H2D_C3_T0_vs_C7_T0_FC GO_MWU  ----
H2D_C7_T0_vs_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T0", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame() %>% rownames_to_column(var = "gene")  %>% dplyr::select("gene", "log2FoldChange", "padj" ) %>% 
                            dplyr::rename(Gene=gene) %>% dplyr::rename(H2D_C7=log2FoldChange, H2D_C7_padj=padj)

H2D_C8_T0_vs_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T0", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame() %>% rownames_to_column(var = "gene")  %>% dplyr::select("gene", "log2FoldChange", "padj" ) %>% 
                            dplyr::rename(Gene=gene) %>% dplyr::rename(H2D_C8=log2FoldChange, H2D_C8_padj=padj)

H2D_C3_C7_C8 <- merge(H2D_C7_T0_vs_C3_T0_FC, H2D_C8_T0_vs_C3_T0_FC, by="Gene" )  %>%
                mutate(H2D_C3=0) %>% dplyr::select(Gene, H2D_C3, H2D_C7, H2D_C8)

positive_H2D_transcriptome_clean <- positive_H2D_transcriptome %>% dplyr::select(Gene)
negative_H2D_transcriptome_clean <- negative_H2D_transcriptome %>% dplyr::select(Gene)

List_FC_genes_H2D_age_transcriptome <- rbind(positive_H2D_transcriptome_clean, negative_H2D_transcriptome_clean)  %>% 
                                       left_join(H2D_C3_C7_C8) %>%  pivot_longer(cols = starts_with("H2D"),  names_to = "Condition",     
                                       values_to = "FoldChange")


Gene_to_plot <- c("G22104") 


Gene_List_FC_genes_F14R_age_transcriptome_filter <- filter(List_FC_genes_F14R_age_transcriptome, Gene==Gene_to_plot)

Gene_to_plot_heatmap_F14R <- ggplot(Gene_List_FC_genes_F14R_age_transcriptome_filter, aes(Condition, Gene, fill = FoldChange)) +
  geom_tile(color = "black", size = 0.2, width = 1, height = 1) +
  scale_fill_gradientn(
    colors = c("darkblue", "blue", "white", "red", "darkred"),  # Colores del gradiente
    values = scales::rescale(c(-5, -2.5, 0, 2.5, 5)),  # Asignación de los valores del rango
    limits = c(-5, 5)  ) +
  theme_minimal() +  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 20, color = "#2f5597", face = "bold"), # Activate gene names here
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none") +
  coord_fixed()

Gene_to_plot_heatmap_F14R

#ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Gene_to_plot_heatmap.tiff", plot = Gene_to_plot_heatmap_F14R, width = 2, height = 1, dpi = 600, bg = "transparent")



Gene_List_FC_genes_H2D_age_transcriptome_filter <- filter(List_FC_genes_H2D_age_transcriptome, Gene==Gene_to_plot)

Gene_to_plot_heatmap_H2D <- ggplot(Gene_List_FC_genes_H2D_age_transcriptome_filter, aes(Condition, Gene, fill = FoldChange)) +
  geom_tile(color = "black", size = 0.2, width = 1, height = 1) +
  scale_fill_gradientn(
    colors = c("darkblue", "blue", "white", "red", "darkred"),
    values = scales::rescale(c(-5, -2.5, 0, 2.5, 5)),
    limits = c(-5, 5)) +
  theme_minimal() +  
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 20, color = "#385700", face = "bold"),# Activate gene names here
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none") +
  coord_fixed()

Gene_to_plot_heatmap_H2D
#ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Gene_to_plot_heatmap.tiff", plot = Gene_to_plot_heatmap_H2D, width = 2, height = 1, dpi = 600, bg = "transparent")


combined_heatmaps <- Gene_to_plot_heatmap_F14R / Gene_to_plot_heatmap_H2D 
combined_heatmaps

ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/Combined_Gene_to_plot_heatmap.tiff",
       plot = combined_heatmaps,    width = 2, height = 1, dpi = 600, bg = "transparent")


























To integrate



script_path    <- "D:/Decicomp/R/RNA_decicomp_GO_terms_infection_kinetics_reference_fix.R"
script_lines   <- readLines(script_path)

specific_lines <- script_lines[78]
eval(parse(text = specific_lines))

specific_lines <- script_lines[3798:3802]
eval(parse(text = specific_lines))

specific_lines <- script_lines[3892:3896]
eval(parse(text = specific_lines))

specific_lines <- script_lines[6120:6124]
eval(parse(text = specific_lines))

specific_lines <- script_lines[6210:6214]
eval(parse(text = specific_lines))

specific_lines <- script_lines[4288:4292]
eval(parse(text = specific_lines))

specific_lines <- script_lines[4382:4386]
eval(parse(text = specific_lines))

specific_lines <- script_lines[6634:6638]
eval(parse(text = specific_lines))

specific_lines <- script_lines[6697:6701]
eval(parse(text = specific_lines))

specific_lines <- script_lines[6962:6966]
eval(parse(text = specific_lines))

specific_lines <- script_lines[6986:6990]
eval(parse(text = specific_lines))

specific_lines <- script_lines[7009:7013]
eval(parse(text = specific_lines))

specific_lines <- script_lines[7032:7036]
eval(parse(text = specific_lines))

specific_lines <- script_lines[7136:7140]
eval(parse(text = specific_lines))

specific_lines <- script_lines[7159:7163]
eval(parse(text = specific_lines))

specific_lines <- script_lines[7182:7186]
eval(parse(text = specific_lines))

specific_lines <- script_lines[7207:7211]
eval(parse(text = specific_lines))


# GO TERMS VEN DIAGRAM

# F14R # EPI
David_GO_BP_kME.positive_F14R_methylation
dim(David_GO_BP_kME.positive_F14R_methylation)  # 10 
#write.table(David_GO_BP_kME.positive_F14R_methylation, file = "D:/Decicomp/temp/aa.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_F14R_positive_methylation        <- David_GO_BP_kME.positive_F14R_methylation$last_term # %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
Revigo_GO_terms_list_F14R_positive_methylation <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_F14R_positive_methylation.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_F14R_positive_methylation) # 8
Revigo_GO_terms_list_F14R_positive_methylation_terms <- Revigo_GO_terms_list_F14R_positive_methylation$TermID

David_GO_BP_kME.negative_F14R_methylation
dim(David_GO_BP_kME.negative_F14R_methylation) # 9 
#write.table(David_GO_BP_kME.negative_F14R_methylation, file = "D:/Decicomp/temp/bb.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_F14R_negative_methylation  <- David_GO_BP_kME.negative_F14R_methylation$last_term # %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
Revigo_GO_terms_list_F14R_negative_methylation <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_F14R_negative_methylation.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_F14R_negative_methylation) # 8
Revigo_GO_terms_list_F14R_negative_methylation_terms <- Revigo_GO_terms_list_F14R_negative_methylation$TermID


# F14R # TRANS
GO_terms_genes_kME.positive_F14R_transcriptome
dim(GO_terms_genes_kME.positive_F14R_transcriptome) # 167
#write.table(GO_terms_genes_kME.positive_F14R_transcriptome, file = "D:/Decicomp/temp/a.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_F14R_postive_transcriptome        <- GO_terms_genes_kME.positive_F14R_transcriptome$last_term
Revigo_GO_terms_list_F14R_positive_transcriptome <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_F14R_positive_transcriptome.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_F14R_positive_transcriptome) # 74 
Revigo_GO_terms_list_F14R_positive_transcriptome_terms <- Revigo_GO_terms_list_F14R_positive_transcriptome$TermID


GO_terms_genes_kME.negative_F14R_transcriptome
dim(GO_terms_genes_kME.negative_F14R_transcriptome) # 175
#write.table(GO_terms_genes_kME.negative_F14R_transcriptome, file = "D:/Decicomp/temp/b.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_F14R_negative_transcriptome         <- GO_terms_genes_kME.negative_F14R_transcriptome$last_term
Revigo_GO_terms_list_F14R_negative_transcriptome  <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_F14R_negative_transcriptome.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_F14R_negative_transcriptome) # 91
Revigo_GO_terms_list_F14R_negative_transcriptome_terms <- Revigo_GO_terms_list_F14R_negative_transcriptome$TermID


# F14R # METABO
GO_terms_Omics_F14R_transcriptome_positive_Genes_kmer
dim(GO_terms_Omics_F14R_transcriptome_positive_Genes_kmer) # 165

GO_terms_Omics_F14R_methylation_positive_Genes_kmer
dim(GO_terms_Omics_F14R_methylation_positive_Genes_kmer)   #  9



GO_terms_Omics_F14R_transcriptome_negative_Genes_kmer
dim(GO_terms_Omics_F14R_transcriptome_negative_Genes_kmer) # 96

GO_terms_Omics_F14R_methylation_negative_Genes_kmer      
dim(GO_terms_Omics_F14R_methylation_negative_Genes_kmer)   # 41


# Para los GO+terms de metabolics, como viene de la mezlca de genes tenemos que juntarlos, porque entendemos que la informacion es redundant.
GO_terms_Omics_F14R_join_positive <- rbind(GO_terms_Omics_F14R_transcriptome_positive_Genes_kmer,
                                           GO_terms_Omics_F14R_methylation_positive_Genes_kmer)
#write.table(GO_terms_Omics_F14R_join_positive, file = "D:/Decicomp/temp/c.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
Revigo_GO_terms_list_F14R_positive_metabolics  <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_F14R_positive_metabolics.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_F14R_positive_metabolics) # 68
Revigo_GO_terms_list_F14R_positive_metabolics_terms <- Revigo_GO_terms_list_F14R_positive_metabolics$TermID



GO_terms_Omics_F14R_join_negative <- rbind(GO_terms_Omics_F14R_transcriptome_negative_Genes_kmer,
                                           GO_terms_Omics_F14R_methylation_negative_Genes_kmer)
#write.table(GO_terms_Omics_F14R_join_negative, file = "D:/Decicomp/temp/d.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
Revigo_GO_terms_list_F14R_negative_metabolics  <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_F14R_negative_metabolics.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_F14R_negative_metabolics) # 49
Revigo_GO_terms_list_F14R_negative_metabolics_terms <- Revigo_GO_terms_list_F14R_negative_metabolics$TermID



# H2D # EPI
David_GO_BP_kME.positive_H2D_methylation
dim(David_GO_BP_kME.positive_H2D_methylation) # 10 
#write.table(David_GO_BP_kME.positive_H2D_methylation, file = "D:/Decicomp/temp/e.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_H2D_postive_methylation        <- David_GO_BP_kME.positive_H2D_methylation$last_term # %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
Revigo_GO_terms_list_H2D_positive_methylation <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_H2D_positive_methylation.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_H2D_positive_methylation) # 8
Revigo_GO_terms_list_H2D_positive_methylation_terms <- Revigo_GO_terms_list_H2D_positive_methylation$TermID


David_GO_BP_kME.negative_H2D_methylation
dim(David_GO_BP_kME.negative_H2D_methylation) # 23
#write.table(David_GO_BP_kME.negative_H2D_methylation, file = "D:/Decicomp/temp/f.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_H2D_negative_methylation   <- David_GO_BP_kME.negative_H2D_methylation$last_term # %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
Revigo_GO_terms_list_H2D_negative_methylation <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_H2D_negative_methylation.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_H2D_negative_methylation) # 18
Revigo_GO_terms_list_H2D_negative_methylation_terms <- Revigo_GO_terms_list_H2D_negative_methylation$TermID

# H2D # TRANS
GO_terms_genes_kME.positive_H2D_transcriptome
dim(GO_terms_genes_kME.positive_H2D_transcriptome) # 95
#write.table(GO_terms_genes_kME.positive_H2D_transcriptome, file = "D:/Decicomp/temp/g.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_list_H2D_postive_transcriptome   <- GO_terms_genes_kME.positive_H2D_transcriptome$last_term # %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
Revigo_GO_terms_list_H2D_positive_transcriptome <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_H2D_positive_transcriptome.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_list_H2D_positive_transcriptome) # 51
Revigo_GO_terms_list_H2D_positive_transcriptome_terms <- Revigo_GO_terms_list_H2D_positive_transcriptome$TermID


GO_terms_genes_kME.negative_H2D_transcriptome
dim(GO_terms_genes_kME.negative_H2D_transcriptome) # 195
#write.table(GO_terms_genes_kME.negative_H2D_transcriptome, file = "D:/Decicomp/temp/h.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
GO_terms_genes_negative_H2D_transcriptome        <- GO_terms_genes_kME.negative_H2D_transcriptome$last_term 
Revigo_GO_terms_genes_negative_H2D_transcriptome <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_H2D_negative_transcriptome.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_genes_negative_H2D_transcriptome) # 89
Revigo_GO_terms_list_H2D_negative_transcriptome_terms <- Revigo_GO_terms_genes_negative_H2D_transcriptome$TermID



# H2D # METABO
GO_terms_Omics_H2D_transcriptome_positive_Genes_kmer
dim(GO_terms_Omics_H2D_transcriptome_positive_Genes_kmer) # 84

GO_terms_Omics_H2D_methylation_positive_Genes_kmer
dim(GO_terms_Omics_H2D_methylation_positive_Genes_kmer)   #  4


GO_terms_Omics_H2D_transcriptome_negative_Genes_kmer
dim(GO_terms_Omics_H2D_transcriptome_negative_Genes_kmer) # 119

GO_terms_Omics_H2D_methylation_negative_Genes_kmer
dim(GO_terms_Omics_H2D_methylation_negative_Genes_kmer)   #  89


GO_terms_Omics_H2D_join_positive <- rbind(GO_terms_Omics_H2D_transcriptome_positive_Genes_kmer,
                                          GO_terms_Omics_H2D_methylation_positive_Genes_kmer)
#write.table(GO_terms_Omics_H2D_join_positive, file = "D:/Decicomp/temp/i.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
Revigo_GO_terms_Omics_H2D_join_positive_metabolics  <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_H2D_positive_metabolics.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_Omics_H2D_join_positive_metabolics) # 44
Revigo_GO_terms_list_H2D_positive_metabolics_terms <- Revigo_GO_terms_Omics_H2D_join_positive_metabolics$TermID


GO_terms_Omics_H2D_join_negative <- rbind(GO_terms_Omics_H2D_transcriptome_negative_Genes_kmer,
                                          GO_terms_Omics_H2D_methylation_negative_Genes_kmer)
#write.table(GO_terms_Omics_H2D_join_negative, file = "D:/Decicomp/temp/j.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
Revigo_GO_terms_Omics_H2D_join_negative_metabolics  <- read_delim("D:/Decicomp/R/Revigo_compacting_GO/Revigo_GO_terms_list_H2D_negative_metabolics.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
dim(Revigo_GO_terms_Omics_H2D_join_negative_metabolics) # 64
Revigo_GO_terms_list_H2D_negative_metabolics_terms <- Revigo_GO_terms_Omics_H2D_join_negative_metabolics$TermID 


GO_terms_F14R_positive <- list(
  Epigenetics     = Revigo_GO_terms_list_F14R_positive_methylation_terms,
  Transcriptomics = Revigo_GO_terms_list_F14R_positive_transcriptome_terms,
  Metabolism      = Revigo_GO_terms_list_F14R_positive_metabolics_terms) 

GO_terms_matrix_F14R_positive <- fromList(GO_terms_F14R_positive)

F14R_positive  <- upset(
  GO_terms_matrix_F14R_positive, 
  nsets = 6,
  set_size.show = T,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#2f5597","#2f5597","#2f5597"),
  empty.intersections = F,
  show.numbers = "yes",
  keep.order = F)

F14R_positive


GO_terms_F14R_negative <- list(
  Epigenetics     = Revigo_GO_terms_list_F14R_negative_methylation_terms,
  Transcriptomics = Revigo_GO_terms_list_F14R_negative_transcriptome_terms,
  Metabolism      = Revigo_GO_terms_list_F14R_negative_metabolics_terms) 

GO_terms_matrix_F14R_negative <- fromList(GO_terms_F14R_negative)

F14R_negative  <- upset(
  GO_terms_matrix_F14R_negative, 
  nsets = 6,
  set_size.show = T,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#2f5597","#2f5597","#2f5597"),
  empty.intersections = F,
  show.numbers = "yes",
  keep.order = F)

F14R_negative




GO_terms_H2D_positive <- list(
  Epigenetics     = Revigo_GO_terms_list_H2D_positive_methylation_terms,
  Transcriptomics = Revigo_GO_terms_list_H2D_positive_transcriptome_terms,
  Metabolism      = Revigo_GO_terms_list_H2D_positive_metabolics_terms) 

GO_terms_matrix_H2D_positive <- fromList(GO_terms_H2D_positive)




H2D_positive  <- upset(
  GO_terms_matrix_H2D_positive, 
  nsets = 6,
  set_size.show = T,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#385700","#385700","#385700"),
  empty.intersections = F,
  show.numbers = "yes",
  keep.order = F)

H2D_positive


GO_terms_H2D_negative <- list(
  Epigenetics     = Revigo_GO_terms_list_H2D_negative_methylation_terms,
  Transcriptomics = Revigo_GO_terms_list_H2D_negative_transcriptome_terms,
  Metabolism      = Revigo_GO_terms_list_H2D_negative_metabolics_terms) 

GO_terms_matrix_H2D_negative <- fromList(GO_terms_H2D_negative)

H2D_negative  <- upset(
  GO_terms_matrix_H2D_negative, 
  nsets = 6,
  set_size.show = T,
  set_size.numbers_size=8,
  point.size=6,
  line.size =2,
  text.scale=3,
  sets.bar.color=c("#385700","#385700","#385700"),
  empty.intersections = F,
  show.numbers = "yes",
  keep.order = F)

H2D_negative




