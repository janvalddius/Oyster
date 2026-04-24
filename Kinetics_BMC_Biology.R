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
library(tibble)
library(pheatmap)
library(readxl)
library(tidyr)


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
countData            <- countData %>% #dplyr::select(matches(grep("T0|T0|T6|T12", colnames(countData), value = TRUE))) %>% 
  dplyr::select(-c("F14R_C3_T24_S5", 
                   "F14R_C3_T24_S6",
                   
                   "F14R_C7_T0_S3",  
                   "F14R_C7_T24_S4",
                   
                   "F14R_C8_T6_S3",  
                   "F14R_C8_T24_S1",
                  
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


##### ANALYSIS KINETICS ## ----

# F14R ----

# F14R_C3_T3_vs_C3_T0 ----
F14R_C3_T3_vs_C3_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_F14R_C3_T3_vs_C3_T0_FC               <- left_join(List_genes, F14R_C3_T3_vs_C3_T0_FC)
write.table(List_genes_F14R_C3_T3_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T3_vs_C3_T0/List_genes_F14R_C3_T3_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T3_vs_C3_T0/")
input="List_genes_F14R_C3_T3_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C3_T3_vs_C3_T0  <- results(dds, contrast=c("group", "F14R_C3_T3", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C3_T3_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T3_vs_C3_T0/BP_List_genes_F14R_C3_T3_vs_C3_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)    %>% 
                                      left_join(F14R_C3_T3_vs_C3_T0, by = "Gene") 
  
gene_counts_F14R_C3_T3_vs_C3_T0   <- BP_List_genes_F14R_C3_T3_vs_C3_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                      nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T3_vs_C3_T0  <- BP_List_genes_F14R_C3_T3_vs_C3_T0 %>% 
                                      left_join(gene_counts_F14R_C3_T3_vs_C3_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C3_T3_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T3_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T3_vs_C3_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_F14R_C3_T3_vs_C3_T0_sig)

GO_terms_F14R_C3_T3_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T3_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T3_vs_C3_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_F14R_C3_T3_vs_C3_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "F14R_C3_T3_vs_C3_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C3_T3_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T3_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C3_T6_vs_C3_T0 ----
F14R_C3_T6_vs_C3_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_F14R_C3_T6_vs_C3_T0_FC               <- left_join(List_genes, F14R_C3_T6_vs_C3_T0_FC)
write.table(List_genes_F14R_C3_T6_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T6_vs_C3_T0/List_genes_F14R_C3_T6_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T6_vs_C3_T0/")
input="List_genes_F14R_C3_T6_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C3_T6_vs_C3_T0  <- results(dds, contrast=c("group", "F14R_C3_T6", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C3_T6_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T6_vs_C3_T0/BP_List_genes_F14R_C3_T6_vs_C3_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)    %>% 
                                      left_join(F14R_C3_T6_vs_C3_T0, by = "Gene") 

gene_counts_F14R_C3_T6_vs_C3_T0   <- BP_List_genes_F14R_C3_T6_vs_C3_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T6_vs_C3_T0  <- BP_List_genes_F14R_C3_T6_vs_C3_T0 %>% 
                                      left_join(gene_counts_F14R_C3_T6_vs_C3_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C3_T6_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T6_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T6_vs_C3_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_F14R_C3_T6_vs_C3_T0_sig)

GO_terms_F14R_C3_T6_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T6_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T6_vs_C3_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_F14R_C3_T6_vs_C3_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "F14R_C3_T6_vs_C3_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C3_T6_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T6_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C3_T12_vs_C3_T0 ----
F14R_C3_T12_vs_C3_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()                          %>% 
                              rownames_to_column(var = "gene")         %>% 
                              filter(!is.na(padj))                     %>% 
                              dplyr::select("gene", "log2FoldChange")  %>% 
                              dplyr::rename(Gene=gene)

List_genes_F14R_C3_T12_vs_C3_T0_FC               <- left_join(List_genes, F14R_C3_T12_vs_C3_T0_FC)
write.table(List_genes_F14R_C3_T12_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T12_vs_C3_T0/List_genes_F14R_C3_T12_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T12_vs_C3_T0/")
input="List_genes_F14R_C3_T12_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C3_T12_vs_C3_T0  <- results(dds, contrast=c("group", "F14R_C3_T12", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()  %>% 
                         rownames_to_column(var = "Gene") %>% 
                         dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C3_T12_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T12_vs_C3_T0/BP_List_genes_F14R_C3_T12_vs_C3_T0_FC.txt") %>% 
                                       dplyr::rename("Gene" = "seq") %>% 
                                       filter(lev != -1)             %>% 
                                       dplyr::select(-value)         %>% 
                                       left_join(F14R_C3_T12_vs_C3_T0, by = "Gene") 

gene_counts_F14R_C3_T12_vs_C3_T0   <- BP_List_genes_F14R_C3_T12_vs_C3_T0 %>%
                                      dplyr::group_by(term)             %>% 
                                      dplyr::summarize(nseqs_relative = dplyr::n(),
                                                       nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T12_vs_C3_T0  <- BP_List_genes_F14R_C3_T12_vs_C3_T0 %>% 
                                       left_join(gene_counts_F14R_C3_T12_vs_C3_T0, by = "term") %>% 
                                       mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C3_T12_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T12_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T12_vs_C3_T0_FC.txt", sep="") %>% 
                                      filter(p.adj <= 0.05)
dim(GO_terms_F14R_C3_T12_vs_C3_T0_sig)

GO_terms_F14R_C3_T12_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T12_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T12_vs_C3_T0_FC.txt", sep="") %>% 
                                  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                  left_join(BP_List_genes_F14R_C3_T12_vs_C3_T0, by = c("term", "name")) %>% 
                                  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                  mutate(Comparison = "F14R_C3_T12_vs_C3_T0") %>% 
                                  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C3_T12_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T12_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C3_T24_vs_C3_T0 ----
F14R_C3_T24_vs_C3_T0_FC    <- results(dds, contrast=c("group", "F14R_C3_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()                          %>% 
                              rownames_to_column(var = "gene")         %>% 
                              filter(!is.na(padj))                     %>% 
                              dplyr::select("gene", "log2FoldChange")  %>% 
                              dplyr::rename(Gene=gene)

List_genes_F14R_C3_T24_vs_C3_T0_FC               <- left_join(List_genes, F14R_C3_T24_vs_C3_T0_FC)
write.table(List_genes_F14R_C3_T24_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T24_vs_C3_T0/List_genes_F14R_C3_T24_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T24_vs_C3_T0/")
input="List_genes_F14R_C3_T24_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C3_T24_vs_C3_T0  <- results(dds, contrast=c("group", "F14R_C3_T24", "F14R_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()  %>% 
                         rownames_to_column(var = "Gene") %>% 
                         dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C3_T24_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T24_vs_C3_T0/BP_List_genes_F14R_C3_T24_vs_C3_T0_FC.txt") %>% 
                                       dplyr::rename("Gene" = "seq") %>% 
                                       filter(lev != -1)             %>% 
                                       dplyr::select(-value)    %>% 
                                       left_join(F14R_C3_T24_vs_C3_T0, by = "Gene") 

gene_counts_F14R_C3_T24_vs_C3_T0   <- BP_List_genes_F14R_C3_T24_vs_C3_T0 %>%
                                      dplyr::group_by(term)             %>% 
                                      dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C3_T24_vs_C3_T0  <- BP_List_genes_F14R_C3_T24_vs_C3_T0 %>% 
                                       left_join(gene_counts_F14R_C3_T24_vs_C3_T0, by = "term") %>% 
                                       mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C3_T24_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T24_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T24_vs_C3_T0_FC.txt", sep="") %>% 
                                      filter(p.adj <= 0.05)
dim(GO_terms_F14R_C3_T24_vs_C3_T0_sig)

GO_terms_F14R_C3_T24_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C3_T24_vs_C3_T0/MWU_BP_List_genes_F14R_C3_T24_vs_C3_T0_FC.txt", sep="") %>% 
                                  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                  left_join(BP_List_genes_F14R_C3_T24_vs_C3_T0, by = c("term", "name")) %>% 
                                  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                  mutate(Comparison = "F14R_C3_T24_vs_C3_T0") %>% 
                                  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C3_T24_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T24_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C7_T3_vs_C7_T0 ----
F14R_C7_T3_vs_C7_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T3", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_F14R_C7_T3_vs_C7_T0_FC               <- left_join(List_genes, F14R_C7_T3_vs_C7_T0_FC)
write.table(List_genes_F14R_C7_T3_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T3_vs_C7_T0/List_genes_F14R_C7_T3_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T3_vs_C7_T0/")
input="List_genes_F14R_C7_T3_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C7_T3_vs_C7_T0  <- results(dds, contrast=c("group", "F14R_C7_T3", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C7_T3_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T3_vs_C7_T0/BP_List_genes_F14R_C7_T3_vs_C7_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)    %>% 
                                      left_join(F14R_C7_T3_vs_C7_T0, by = "Gene") 

gene_counts_F14R_C7_T3_vs_C7_T0   <- BP_List_genes_F14R_C7_T3_vs_C7_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T3_vs_C7_T0  <- BP_List_genes_F14R_C7_T3_vs_C7_T0 %>% 
                                      left_join(gene_counts_F14R_C7_T3_vs_C7_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C7_T3_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T3_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T3_vs_C7_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_F14R_C7_T3_vs_C7_T0_sig)

GO_terms_F14R_C7_T3_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T3_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T3_vs_C7_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_F14R_C7_T3_vs_C7_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "F14R_C7_T3_vs_C7_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C7_T3_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T3_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C7_T6_vs_C7_T0 ----
F14R_C7_T6_vs_C7_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T6", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_F14R_C7_T6_vs_C7_T0_FC               <- left_join(List_genes, F14R_C7_T6_vs_C7_T0_FC)
write.table(List_genes_F14R_C7_T6_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T6_vs_C7_T0/List_genes_F14R_C7_T6_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T6_vs_C7_T0/")
input="List_genes_F14R_C7_T6_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C7_T6_vs_C7_T0  <- results(dds, contrast=c("group", "F14R_C7_T6", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C7_T6_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T6_vs_C7_T0/BP_List_genes_F14R_C7_T6_vs_C7_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)    %>% 
                                      left_join(F14R_C7_T6_vs_C7_T0, by = "Gene") 

gene_counts_F14R_C7_T6_vs_C7_T0   <- BP_List_genes_F14R_C7_T6_vs_C7_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T6_vs_C7_T0  <- BP_List_genes_F14R_C7_T6_vs_C7_T0 %>% 
                                      left_join(gene_counts_F14R_C7_T6_vs_C7_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C7_T6_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T6_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T6_vs_C7_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_F14R_C7_T6_vs_C7_T0_sig)

GO_terms_F14R_C7_T6_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T6_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T6_vs_C7_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_F14R_C7_T6_vs_C7_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "F14R_C7_T6_vs_C7_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C7_T6_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T6_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C7_T12_vs_C7_T0 ----
F14R_C7_T12_vs_C7_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T12", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()                          %>% 
                              rownames_to_column(var = "gene")         %>% 
                              filter(!is.na(padj))                     %>% 
                              dplyr::select("gene", "log2FoldChange")  %>% 
                              dplyr::rename(Gene=gene)

List_genes_F14R_C7_T12_vs_C7_T0_FC               <- left_join(List_genes, F14R_C7_T12_vs_C7_T0_FC)
write.table(List_genes_F14R_C7_T12_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T12_vs_C7_T0/List_genes_F14R_C7_T12_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T12_vs_C7_T0/")
input="List_genes_F14R_C7_T12_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C7_T12_vs_C7_T0  <- results(dds, contrast=c("group", "F14R_C7_T12", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()  %>% 
                         rownames_to_column(var = "Gene") %>% 
                         dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C7_T12_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T12_vs_C7_T0/BP_List_genes_F14R_C7_T12_vs_C7_T0_FC.txt") %>% 
                                       dplyr::rename("Gene" = "seq") %>% 
                                       filter(lev != -1)             %>% 
                                       dplyr::select(-value)    %>% 
                                       left_join(F14R_C7_T12_vs_C7_T0, by = "Gene") 

gene_counts_F14R_C7_T12_vs_C7_T0   <- BP_List_genes_F14R_C7_T12_vs_C7_T0 %>%
                                      dplyr::group_by(term)             %>% 
                                      dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T12_vs_C7_T0  <- BP_List_genes_F14R_C7_T12_vs_C7_T0 %>% 
                                       left_join(gene_counts_F14R_C7_T12_vs_C7_T0, by = "term") %>% 
                                       mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C7_T12_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T12_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T12_vs_C7_T0_FC.txt", sep="") %>% 
                                      filter(p.adj <= 0.05)
dim(GO_terms_F14R_C7_T12_vs_C7_T0_sig)

GO_terms_F14R_C7_T12_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T12_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T12_vs_C7_T0_FC.txt", sep="") %>% 
                                  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                  left_join(BP_List_genes_F14R_C7_T12_vs_C7_T0, by = c("term", "name")) %>% 
                                  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                  mutate(Comparison = "F14R_C7_T12_vs_C7_T0") %>% 
                                  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C7_T12_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T12_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C7_T24_vs_C7_T0 ----
F14R_C7_T24_vs_C7_T0_FC    <- results(dds, contrast=c("group", "F14R_C7_T24", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()                          %>% 
                              rownames_to_column(var = "gene")         %>% 
                              filter(!is.na(padj))                     %>% 
                              dplyr::select("gene", "log2FoldChange")  %>% 
                              dplyr::rename(Gene=gene)

List_genes_F14R_C7_T24_vs_C7_T0_FC               <- left_join(List_genes, F14R_C7_T24_vs_C7_T0_FC)
write.table(List_genes_F14R_C7_T24_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T24_vs_C7_T0/List_genes_F14R_C7_T24_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T24_vs_C7_T0/")
input="List_genes_F14R_C7_T24_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C7_T24_vs_C7_T0  <- results(dds, contrast=c("group", "F14R_C7_T24", "F14R_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()  %>% 
                         rownames_to_column(var = "Gene") %>% 
                         dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C7_T24_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T24_vs_C7_T0/BP_List_genes_F14R_C7_T24_vs_C7_T0_FC.txt") %>% 
                                       dplyr::rename("Gene" = "seq") %>% 
                                       filter(lev != -1)             %>% 
                                       dplyr::select(-value)    %>% 
                                       left_join(F14R_C7_T24_vs_C7_T0, by = "Gene") 

gene_counts_F14R_C7_T24_vs_C7_T0   <- BP_List_genes_F14R_C7_T24_vs_C7_T0 %>%
                                      dplyr::group_by(term)             %>% 
                                      dplyr::summarize(nseqs_relative = dplyr::n(),
                                                       nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C7_T24_vs_C7_T0  <- BP_List_genes_F14R_C7_T24_vs_C7_T0 %>% 
                                       left_join(gene_counts_F14R_C7_T24_vs_C7_T0, by = "term") %>% 
                                       mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C7_T24_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T24_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T24_vs_C7_T0_FC.txt", sep="") %>% 
                                      filter(p.adj <= 0.05)
dim(GO_terms_F14R_C7_T24_vs_C7_T0_sig)

GO_terms_F14R_C7_T24_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C7_T24_vs_C7_T0/MWU_BP_List_genes_F14R_C7_T24_vs_C7_T0_FC.txt", sep="") %>% 
                                  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                  left_join(BP_List_genes_F14R_C7_T24_vs_C7_T0, by = c("term", "name")) %>% 
                                  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                  mutate(Comparison = "F14R_C7_T24_vs_C7_T0") %>% 
                                  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C7_T24_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T24_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C8_T3_vs_C8_T0 ----
F14R_C8_T3_vs_C8_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T3", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_F14R_C8_T3_vs_C8_T0_FC               <- left_join(List_genes, F14R_C8_T3_vs_C8_T0_FC)
write.table(List_genes_F14R_C8_T3_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T3_vs_C8_T0/List_genes_F14R_C8_T3_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T3_vs_C8_T0/")
input="List_genes_F14R_C8_T3_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C8_T3_vs_C8_T0  <- results(dds, contrast=c("group", "F14R_C8_T3", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C8_T3_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T3_vs_C8_T0/BP_List_genes_F14R_C8_T3_vs_C8_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)    %>% 
                                      left_join(F14R_C8_T3_vs_C8_T0, by = "Gene") 

gene_counts_F14R_C8_T3_vs_C8_T0   <- BP_List_genes_F14R_C8_T3_vs_C8_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T3_vs_C8_T0  <- BP_List_genes_F14R_C8_T3_vs_C8_T0 %>% 
                                      left_join(gene_counts_F14R_C8_T3_vs_C8_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C8_T3_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T3_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T3_vs_C8_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_F14R_C8_T3_vs_C8_T0_sig)

GO_terms_F14R_C8_T3_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T3_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T3_vs_C8_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_F14R_C8_T3_vs_C8_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "F14R_C8_T3_vs_C8_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C8_T3_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T3_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C8_T6_vs_C8_T0 ----
F14R_C8_T6_vs_C8_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T6", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_F14R_C8_T6_vs_C8_T0_FC               <- left_join(List_genes, F14R_C8_T6_vs_C8_T0_FC)
write.table(List_genes_F14R_C8_T6_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T6_vs_C8_T0/List_genes_F14R_C8_T6_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T6_vs_C8_T0/")
input="List_genes_F14R_C8_T6_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C8_T6_vs_C8_T0  <- results(dds, contrast=c("group", "F14R_C8_T6", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C8_T6_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T6_vs_C8_T0/BP_List_genes_F14R_C8_T6_vs_C8_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)    %>% 
                                      left_join(F14R_C8_T6_vs_C8_T0, by = "Gene") 

gene_counts_F14R_C8_T6_vs_C8_T0   <- BP_List_genes_F14R_C8_T6_vs_C8_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T6_vs_C8_T0  <- BP_List_genes_F14R_C8_T6_vs_C8_T0 %>% 
                                      left_join(gene_counts_F14R_C8_T6_vs_C8_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C8_T6_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T6_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T6_vs_C8_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_F14R_C8_T6_vs_C8_T0_sig)

GO_terms_F14R_C8_T6_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T6_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T6_vs_C8_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_F14R_C8_T6_vs_C8_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "F14R_C8_T6_vs_C8_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C8_T6_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T6_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C8_T12_vs_C8_T0 ----
F14R_C8_T12_vs_C8_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T12", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()                          %>% 
                              rownames_to_column(var = "gene")         %>% 
                              filter(!is.na(padj))                     %>% 
                              dplyr::select("gene", "log2FoldChange")  %>% 
                              dplyr::rename(Gene=gene)

List_genes_F14R_C8_T12_vs_C8_T0_FC               <- left_join(List_genes, F14R_C8_T12_vs_C8_T0_FC)
write.table(List_genes_F14R_C8_T12_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T12_vs_C8_T0/List_genes_F14R_C8_T12_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T12_vs_C8_T0/")
input="List_genes_F14R_C8_T12_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C8_T12_vs_C8_T0  <- results(dds, contrast=c("group", "F14R_C8_T12", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()  %>% 
                         rownames_to_column(var = "Gene") %>% 
                         dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C8_T12_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T12_vs_C8_T0/BP_List_genes_F14R_C8_T12_vs_C8_T0_FC.txt") %>% 
                                       dplyr::rename("Gene" = "seq") %>% 
                                       filter(lev != -1)             %>% 
                                       dplyr::select(-value)    %>% 
                                       left_join(F14R_C8_T12_vs_C8_T0, by = "Gene") 

gene_counts_F14R_C8_T12_vs_C8_T0   <- BP_List_genes_F14R_C8_T12_vs_C8_T0 %>%
                                      dplyr::group_by(term)             %>% 
                                      dplyr::summarize(nseqs_relative = dplyr::n(),
                                                       nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T12_vs_C8_T0  <- BP_List_genes_F14R_C8_T12_vs_C8_T0 %>% 
                                       left_join(gene_counts_F14R_C8_T12_vs_C8_T0, by = "term") %>% 
                                       mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C8_T12_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T12_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T12_vs_C8_T0_FC.txt", sep="") %>% 
                                      filter(p.adj <= 0.05)
dim(GO_terms_F14R_C8_T12_vs_C8_T0_sig)

GO_terms_F14R_C8_T12_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T12_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T12_vs_C8_T0_FC.txt", sep="") %>% 
                                  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                  left_join(BP_List_genes_F14R_C8_T12_vs_C8_T0, by = c("term", "name")) %>% 
                                  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                  mutate(Comparison = "F14R_C8_T12_vs_C8_T0") %>% 
                                  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C8_T12_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T12_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# F14R_C8_T24_vs_C8_T0 ----
F14R_C8_T24_vs_C8_T0_FC    <- results(dds, contrast=c("group", "F14R_C8_T24", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                              as.data.frame()                          %>% 
                              rownames_to_column(var = "gene")         %>% 
                              filter(!is.na(padj))                     %>% 
                              dplyr::select("gene", "log2FoldChange")  %>% 
                              dplyr::rename(Gene=gene)

List_genes_F14R_C8_T24_vs_C8_T0_FC               <- left_join(List_genes, F14R_C8_T24_vs_C8_T0_FC)
write.table(List_genes_F14R_C8_T24_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T24_vs_C8_T0/List_genes_F14R_C8_T24_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T24_vs_C8_T0/")
input="List_genes_F14R_C8_T24_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

F14R_C8_T24_vs_C8_T0  <- results(dds, contrast=c("group", "F14R_C8_T24", "F14R_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                         as.data.frame()  %>% 
                         rownames_to_column(var = "Gene") %>% 
                         dplyr::filter(!is.na(padj)) 

BP_List_genes_F14R_C8_T24_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T24_vs_C8_T0/BP_List_genes_F14R_C8_T24_vs_C8_T0_FC.txt") %>% 
                                       dplyr::rename("Gene" = "seq") %>% 
                                       filter(lev != -1)             %>% 
                                       dplyr::select(-value)    %>% 
                                       left_join(F14R_C8_T24_vs_C8_T0, by = "Gene") 

gene_counts_F14R_C8_T24_vs_C8_T0   <- BP_List_genes_F14R_C8_T24_vs_C8_T0 %>%
                                      dplyr::group_by(term)             %>% 
                                      dplyr::summarize(nseqs_relative = dplyr::n(),
                                                       nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_F14R_C8_T24_vs_C8_T0  <- BP_List_genes_F14R_C8_T24_vs_C8_T0 %>% 
                                       left_join(gene_counts_F14R_C8_T24_vs_C8_T0, by = "term") %>% 
                                       mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_F14R_C8_T24_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T24_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T24_vs_C8_T0_FC.txt", sep="") %>% 
                                      filter(p.adj <= 0.05)
dim(GO_terms_F14R_C8_T24_vs_C8_T0_sig)

GO_terms_F14R_C8_T24_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/F14R_C8_T24_vs_C8_T0/MWU_BP_List_genes_F14R_C8_T24_vs_C8_T0_FC.txt", sep="") %>% 
                                  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                  left_join(BP_List_genes_F14R_C8_T24_vs_C8_T0, by = c("term", "name")) %>% 
                                  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                  mutate(Comparison = "F14R_C8_T24_vs_C8_T0") %>% 
                                  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_F14R_C8_T24_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T24_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


dim(GO_terms_F14R_C3_T3_vs_C3_T0_sig)  # 82
dim(GO_terms_F14R_C3_T6_vs_C3_T0_sig)  # 91
dim(GO_terms_F14R_C3_T12_vs_C3_T0_sig) # 126
dim(GO_terms_F14R_C3_T24_vs_C3_T0_sig) # 156

dim(GO_terms_F14R_C7_T3_vs_C7_T0_sig)  # 13
dim(GO_terms_F14R_C7_T6_vs_C7_T0_sig)  # 45
dim(GO_terms_F14R_C7_T12_vs_C7_T0_sig) # 93
dim(GO_terms_F14R_C7_T24_vs_C7_T0_sig) # 110

dim(GO_terms_F14R_C8_T3_vs_C8_T0_sig)  # 21
dim(GO_terms_F14R_C8_T6_vs_C8_T0_sig)  # 19
dim(GO_terms_F14R_C8_T12_vs_C8_T0_sig) # 69
dim(GO_terms_F14R_C8_T24_vs_C8_T0_sig) # 98


#### Data frames and heatmpas F14R ----
GO_terms_F14R_C3_T3_vs_C3_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T3_vs_C3_T0.tsv")
GO_terms_F14R_C3_T6_vs_C3_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T6_vs_C3_T0.tsv")
GO_terms_F14R_C3_T12_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T12_vs_C3_T0.tsv")
GO_terms_F14R_C3_T24_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C3_T24_vs_C3_T0.tsv")

GO_terms_F14R_C7_T3_vs_C7_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T3_vs_C7_T0.tsv")
GO_terms_F14R_C7_T6_vs_C7_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T6_vs_C7_T0.tsv")
GO_terms_F14R_C7_T12_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T12_vs_C7_T0.tsv")
GO_terms_F14R_C7_T24_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C7_T24_vs_C7_T0.tsv")

GO_terms_F14R_C8_T3_vs_C8_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T3_vs_C8_T0.tsv")
GO_terms_F14R_C8_T6_vs_C8_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T6_vs_C8_T0.tsv")
GO_terms_F14R_C8_T12_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T12_vs_C8_T0.tsv")
GO_terms_F14R_C8_T24_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R/GO_data_frames/GO_terms_F14R_C8_T24_vs_C8_T0.tsv")


GO_terms_F14R_C3_T3_vs_C3_T0_sig  <- GO_terms_F14R_C3_T3_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C3_T3_vs_C3_T0_sig)  # 82

GO_terms_F14R_C3_T6_vs_C3_T0_sig  <- GO_terms_F14R_C3_T6_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C3_T6_vs_C3_T0_sig)  # 91

GO_terms_F14R_C3_T12_vs_C3_T0_sig  <- GO_terms_F14R_C3_T12_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C3_T12_vs_C3_T0_sig)  # 126

GO_terms_F14R_C3_T24_vs_C3_T0_sig  <- GO_terms_F14R_C3_T24_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C3_T24_vs_C3_T0_sig)  # 156


GO_terms_F14R_C7_T3_vs_C7_T0_sig  <- GO_terms_F14R_C7_T3_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C7_T3_vs_C7_T0_sig)   # 13

GO_terms_F14R_C7_T6_vs_C7_T0_sig  <- GO_terms_F14R_C7_T6_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C7_T6_vs_C7_T0_sig)   # 45

GO_terms_F14R_C7_T12_vs_C7_T0_sig  <- GO_terms_F14R_C7_T12_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C7_T12_vs_C7_T0_sig)  # 93

GO_terms_F14R_C7_T24_vs_C7_T0_sig  <- GO_terms_F14R_C7_T24_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C7_T24_vs_C7_T0_sig)  # 110


GO_terms_F14R_C8_T3_vs_C8_T0_sig  <- GO_terms_F14R_C8_T3_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C8_T3_vs_C8_T0_sig)   # 21

GO_terms_F14R_C8_T6_vs_C8_T0_sig  <- GO_terms_F14R_C8_T6_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C8_T6_vs_C8_T0_sig)   # 19

GO_terms_F14R_C8_T12_vs_C8_T0_sig  <- GO_terms_F14R_C8_T12_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C8_T12_vs_C8_T0_sig)  # 69

GO_terms_F14R_C8_T24_vs_C8_T0_sig  <- GO_terms_F14R_C8_T24_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_F14R_C8_T24_vs_C8_T0_sig)  # 98


F14R_transcriptome_GO <- rbind(GO_terms_F14R_C3_T3_vs_C3_T0,
                               GO_terms_F14R_C3_T6_vs_C3_T0,
                               GO_terms_F14R_C3_T12_vs_C3_T0,
                               GO_terms_F14R_C3_T24_vs_C3_T0,
                            
                               GO_terms_F14R_C7_T3_vs_C7_T0,
                               GO_terms_F14R_C7_T6_vs_C7_T0,
                               GO_terms_F14R_C7_T12_vs_C7_T0,
                               GO_terms_F14R_C7_T24_vs_C7_T0,
                            
                               GO_terms_F14R_C8_T3_vs_C8_T0,
                               GO_terms_F14R_C8_T6_vs_C8_T0,
                               GO_terms_F14R_C8_T12_vs_C8_T0,
                               GO_terms_F14R_C8_T24_vs_C8_T0)  %>% 
  
                           filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% distinct(name, .keep_all = TRUE) %>%
                           dplyr::select(Comparison, name, Enrichment) %>% 
                           pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
                           as.data.frame() %>% 
                           mutate_all(~ replace(., is.na(.), 0)) %>%
                           mutate(F14R_C3_T0_vs_C3_T0  = 0,
                                  F14R_C7_T0_vs_C7_T0  = 0,
                                  F14R_C8_T0_vs_C8_T0  = 0) %>%
  
                           filter(name != "unknown") %>%
                           filter(!grepl("obsolete", name, ignore.case = TRUE)) %>% 
  
                           column_to_rownames(var = "name")  %>%
                           dplyr::select(13, 1:4, 14, 5:8, 15, 9:12)

GO_list_F14R_kinetics   <- row.names(F14R_transcriptome_GO) %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
# write.table(GO_list_F14R_kinetics, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_list_F14R_kinetics.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

GO_classified_full <- read_excel("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_classified_full.xlsx")


GO_list_F14R_kinetics_groups <- GO_list_F14R_kinetics %>% left_join(GO_classified_full) %>%
                                group_by(Category) %>% summarise(n = n())

GO_list_F14R_kinetics_table <- F14R_transcriptome_GO %>% tibble::rownames_to_column("GO_term" ) %>% left_join(GO_classified_full)
write.table(GO_list_F14R_kinetics_table, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/GO_list_F14R_kinetics_table.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# Heat maps global F14R  ----
# Colum Colors heatmap   
annotation_col <- data.frame(
  Time = c("T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24"),
  Age  = c(rep("4 months", 5), rep("16 months", 5), rep("28 months", 5)))

annotation_colors  <-  list(
  Age  = c("4 months"   = "#dae3f3", "16 months"  = "#8faadc", "28 months" = "#2f5597"),
  Time = c(T0 = "white",  T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred",  T24 = "red"))

rownames(annotation_col) <- colnames(F14R_transcriptome_GO) 

# Gradient enrichment                        
myBreaks <- c(seq(min(F14R_transcriptome_GO), 0, length.out=ceiling(100/2) + 1),
              seq(max(F14R_transcriptome_GO)/100, max(F14R_transcriptome_GO), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

F14R_transcriptome_GO_heatmap <- pheatmap(F14R_transcriptome_GO, 
                                                cluster_cols = F, 
                                                scale = "none",
                                                cluster_rows = T, 
                                                fontsize_row = 5, 
                                                color = mycolor, 
                                                breaks = myBreaks,
                                                border_color = "black",
                                                clustering_distance_rows = "euclidean",
                                                show_colnames = F, 
                                                show_rownames = F,
                                                annotation_col = annotation_col, 
                                                #annotation_row = annotation_row,
                                                annotation_colors = annotation_colors,
                                                #cutree_rows = 4,
                                                gaps_col =  c(5, 10, 15))
F14R_transcriptome_GO_heatmap
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R_transcriptome_GO_heatmap.tiff", plot = F14R_transcriptome_GO_heatmap$gtable,  width = 10,    height = 12,  dpi = 600)




# =========================================================
# Category
# =========================================================

mat <- F14R_transcriptome_GO
mat[is.na(mat)] <- 0


annotation_col <- data.frame(
  Time = c("T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24"),
  Age = c(rep("4 months",5),
          rep("16 months",5),
          rep("28 months",5)))

rownames(annotation_col) <- colnames(mat)

annotation_colors <- list(
  Age = c("4 months" = "#dae3f3", "16 months" = "#8faadc",  "28 months" = "#2f5597"),
  Time = c(T0 = "white",  T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred", T24 = "red" ))

GO_classified_full <- read_excel("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_classified_full.xlsx")

GO_terms <- data.frame(GO_term = rownames(mat))
GO_map <- GO_classified_full %>%  dplyr::select(GO_term, Category)

annotation_row <- GO_terms %>%  left_join(GO_map, by = "GO_term") %>%
  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>%
  column_to_rownames("GO_term")

annotation_row <- annotation_row[rownames(mat), , drop = FALSE]

category_order <- c(
  "Immune system procces",
  "Response to stimulus",
  "Metabolic procces",
  "Biological regulation",
  "Cellular procces",
  "Communication",
  "Transport",
  "Localisation",
  "Developmental procces",
  "Reproduction")
#"Unknown")

annotation_row$Category <- factor(annotation_row$Category, levels = category_order)
ord                     <- order(annotation_row$Category)

mat_ord <- mat[ord, ]
annotation_row_ord <- annotation_row[ord, , drop = FALSE]
gaps_row           <- cumsum(table(annotation_row_ord$Category))


category_colors <- c(
  # 🔴 TOP 3
  "Immune system procces" = "purple",
  "Response to stimulus"  = "darkorange",
  "Metabolic procces"     = "#4DAF4A",
  # 🌿 OTHERS
  "Biological regulation" = "#8C8C8C",
  "Cellular procces"      = "#A6A6A6",
  "Communication"         = "#BFBFBF",
  "Transport"             = "#D9D9D9",
  "Localisation"          = "#E0E0E0",
  "Developmental procces" = "#E6E6E6",
  "Reproduction"          = "#F2F2F2")
#"Unknown"               = "grey85")

annotation_colors$Category <- category_colors

mycolor  <- colorRampPalette(c("blue","black","yellow"))(100)
myBreaks <- c(seq(min(mat_ord), 0, length.out = 50), seq(max(mat_ord)/100, max(mat_ord), length.out = 50))

F14R_transcriptome_category_heatmap <- pheatmap(
  mat_ord,
  
  cluster_rows = T,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  
  scale = "none",
  color = mycolor,
  breaks = myBreaks,
  
  border_color = "black",
  
  annotation_col    = annotation_col,
  annotation_row    = annotation_row_ord,
  annotation_colors = annotation_colors,
  gaps_row = gaps_row,
  gaps_col = c(5,10,15),
  fontsize_row = 5,
  
  legend = T,              
  annotation_legend = FALSE)

F14R_transcriptome_category_heatmap
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R_transcriptome_category_heatmap.tiff", plot = F14R_transcriptome_category_heatmap$gtable,  width = 8,    height = 12,  dpi = 600)


# Only immunity and stress
F14R_transcriptome_GO_immune <- GO_list_F14R_kinetics_table %>% 
                                   filter(Category=="Immune system procces" | Category=="Response to stimulus") %>% 
                                   dplyr::select(-Category)  %>%
                                   tibble::column_to_rownames("GO_term")


mat_immune <- F14R_transcriptome_GO_immune
mat_immune[is.na(mat_immune)] <- 0

annotation_col <- data.frame(
  Time = c("T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24"),
  Age = c(rep("4 months",5),
          rep("16 months",5),
          rep("28 months",5)))

rownames(annotation_col) <- colnames(mat_immune)

annotation_colors <- list(
  Age = c("4 months" = "#dae3f3", "16 months" = "#8faadc",  "28 months" = "#2f5597"),
  Time = c(T0 = "white",  T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred", T24 = "red" ))

GO_classified_full <- read_excel("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_classified_full.xlsx")

GO_terms <- data.frame(GO_term = rownames(mat_immune))
GO_map <- GO_classified_full %>%  dplyr::select(GO_term, Category)

annotation_row <- GO_terms %>%  left_join(GO_map, by = "GO_term") %>%
  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>%
  column_to_rownames("GO_term")

annotation_row <- annotation_row[rownames(mat_immune), , drop = FALSE]

category_order <- c("Immune system procces",  "Response to stimulus")


annotation_row$Category <- factor(annotation_row$Category, levels = category_order)
ord                     <- order(annotation_row$Category)

mat_ord <- mat_immune[ord, ]
annotation_row_ord <- annotation_row[ord, , drop = FALSE]
gaps_row           <- cumsum(table(annotation_row_ord$Category))


category_colors <- c( "Immune system procces" = "purple", "Response to stimulus"  = "darkorange")

annotation_colors$Category <- category_colors

mycolor  <- colorRampPalette(c("blue","black","yellow"))(100)
myBreaks <- c(seq(min(mat_ord), 0, length.out = 50), seq(max(mat_ord)/100, max(mat_ord), length.out = 50))

F14R_transcriptome_category_Infection_heatmap <- pheatmap(
  mat_ord,
  
  cluster_rows = F,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = T,
  
  scale = "none",
  color = mycolor,
  breaks = myBreaks,
  
  border_color = "black",
  
  annotation_col    = annotation_col,
  annotation_row    = annotation_row_ord,
  annotation_colors = annotation_colors,
  gaps_row = gaps_row,
  gaps_col = c(5,10,15),
  fontsize_row = 14,
  
  legend = T,              
  annotation_legend = FALSE)

F14R_transcriptome_category_Infection_heatmap
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/F14R_transcriptome_category_Infection_heatmap.tiff", plot = F14R_transcriptome_category_Infection_heatmap$gtable,  width = 8,    height =6,  dpi = 600)






# H2D ----
# H2D_C3_T3_vs_C3_T0 ----
H2D_C3_T3_vs_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C3_T3_vs_C3_T0_FC               <- left_join(List_genes, H2D_C3_T3_vs_C3_T0_FC)
write.table(List_genes_H2D_C3_T3_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T3_vs_C3_T0/List_genes_H2D_C3_T3_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T3_vs_C3_T0/")
input="List_genes_H2D_C3_T3_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C3_T3_vs_C3_T0  <- results(dds, contrast=c("group", "H2D_C3_T3", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C3_T3_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T3_vs_C3_T0/BP_List_genes_H2D_C3_T3_vs_C3_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C3_T3_vs_C3_T0, by = "Gene") 

gene_counts_H2D_C3_T3_vs_C3_T0   <- BP_List_genes_H2D_C3_T3_vs_C3_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T3_vs_C3_T0  <- BP_List_genes_H2D_C3_T3_vs_C3_T0 %>% 
                                     left_join(gene_counts_H2D_C3_T3_vs_C3_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C3_T3_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T3_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T3_vs_C3_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C3_T3_vs_C3_T0_sig)

GO_terms_H2D_C3_T3_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T3_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T3_vs_C3_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C3_T3_vs_C3_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C3_T3_vs_C3_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C3_T3_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C3_T3_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C3_T6_vs_C3_T0 ----
H2D_C3_T6_vs_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C3_T6_vs_C3_T0_FC               <- left_join(List_genes, H2D_C3_T6_vs_C3_T0_FC)
write.table(List_genes_H2D_C3_T6_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T6_vs_C3_T0/List_genes_H2D_C3_T6_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T6_vs_C3_T0/")
input="List_genes_H2D_C3_T6_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C3_T6_vs_C3_T0  <- results(dds, contrast=c("group", "H2D_C3_T6", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C3_T6_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T6_vs_C3_T0/BP_List_genes_H2D_C3_T6_vs_C3_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C3_T6_vs_C3_T0, by = "Gene") 

gene_counts_H2D_C3_T6_vs_C3_T0   <- BP_List_genes_H2D_C3_T6_vs_C3_T0 %>%
                                    dplyr::group_by(term)            %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T6_vs_C3_T0  <- BP_List_genes_H2D_C3_T6_vs_C3_T0 %>% 
                                     left_join(gene_counts_H2D_C3_T6_vs_C3_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C3_T6_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T6_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T6_vs_C3_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C3_T6_vs_C3_T0_sig)

GO_terms_H2D_C3_T6_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T6_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T6_vs_C3_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C3_T6_vs_C3_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C3_T6_vs_C3_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C3_T6_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C3_T6_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C3_T12_vs_C3_T0 ----
H2D_C3_T12_vs_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_H2D_C3_T12_vs_C3_T0_FC               <- left_join(List_genes, H2D_C3_T12_vs_C3_T0_FC)
write.table(List_genes_H2D_C3_T12_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T12_vs_C3_T0/List_genes_H2D_C3_T12_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T12_vs_C3_T0/")
input="List_genes_H2D_C3_T12_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C3_T12_vs_C3_T0  <- results(dds, contrast=c("group", "H2D_C3_T12", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C3_T12_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T12_vs_C3_T0/BP_List_genes_H2D_C3_T12_vs_C3_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C3_T12_vs_C3_T0, by = "Gene") 

gene_counts_H2D_C3_T12_vs_C3_T0   <- BP_List_genes_H2D_C3_T12_vs_C3_T0 %>%
                                     dplyr::group_by(term)            %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T12_vs_C3_T0  <- BP_List_genes_H2D_C3_T12_vs_C3_T0 %>% 
                                      left_join(gene_counts_H2D_C3_T12_vs_C3_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C3_T12_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T12_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T12_vs_C3_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C3_T12_vs_C3_T0_sig)

GO_terms_H2D_C3_T12_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T12_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T12_vs_C3_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_H2D_C3_T12_vs_C3_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "H2D_C3_T12_vs_C3_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C3_T12_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C3_T12_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C3_T24_vs_C3_T0 ----
H2D_C3_T24_vs_C3_T0_FC    <- results(dds, contrast=c("group", "H2D_C3_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_H2D_C3_T24_vs_C3_T0_FC               <- left_join(List_genes, H2D_C3_T24_vs_C3_T0_FC)
write.table(List_genes_H2D_C3_T24_vs_C3_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T24_vs_C3_T0/List_genes_H2D_C3_T24_vs_C3_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T24_vs_C3_T0/")
input="List_genes_H2D_C3_T24_vs_C3_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C3_T24_vs_C3_T0  <- results(dds, contrast=c("group", "H2D_C3_T24", "H2D_C3_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C3_T24_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T24_vs_C3_T0/BP_List_genes_H2D_C3_T24_vs_C3_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C3_T24_vs_C3_T0, by = "Gene") 

gene_counts_H2D_C3_T24_vs_C3_T0   <- BP_List_genes_H2D_C3_T24_vs_C3_T0 %>%
                                     dplyr::group_by(term)            %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C3_T24_vs_C3_T0  <- BP_List_genes_H2D_C3_T24_vs_C3_T0 %>% 
                                      left_join(gene_counts_H2D_C3_T24_vs_C3_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C3_T24_vs_C3_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T24_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T24_vs_C3_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C3_T24_vs_C3_T0_sig)

GO_terms_H2D_C3_T24_vs_C3_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C3_T24_vs_C3_T0/MWU_BP_List_genes_H2D_C3_T24_vs_C3_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_H2D_C3_T24_vs_C3_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "H2D_C3_T24_vs_C3_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C3_T24_vs_C3_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C3_T24_vs_C3_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")



# H2D_C7_T3_vs_C7_T0 ----
H2D_C7_T3_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T3", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C7_T3_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T3_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T3_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/List_genes_H2D_C7_T3_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/")
input="List_genes_H2D_C7_T3_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T3_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T3", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T3_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/BP_List_genes_H2D_C7_T3_vs_C7_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C7_T3_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T3_vs_C7_T0   <- BP_List_genes_H2D_C7_T3_vs_C7_T0 %>%
                                    dplyr::group_by(term)             %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T3_vs_C7_T0  <- BP_List_genes_H2D_C7_T3_vs_C7_T0 %>% 
                                     left_join(gene_counts_H2D_C7_T3_vs_C7_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T3_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T3_vs_C7_T0_FC.txt", sep="") %>% 
                                    filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T3_vs_C7_T0_sig)

GO_terms_H2D_C7_T3_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T3_vs_C7_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C7_T3_vs_C7_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C7_T3_vs_C7_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T3_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T3_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C7_T6_vs_C7_T0 ----
H2D_C7_T6_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T6", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C7_T6_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T6_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T6_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/List_genes_H2D_C7_T6_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/")
input="List_genes_H2D_C7_T6_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T6_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T6", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T6_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/BP_List_genes_H2D_C7_T6_vs_C7_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C7_T6_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T6_vs_C7_T0   <- BP_List_genes_H2D_C7_T6_vs_C7_T0 %>%
                                    dplyr::group_by(term)             %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T6_vs_C7_T0  <- BP_List_genes_H2D_C7_T6_vs_C7_T0 %>% 
                                     left_join(gene_counts_H2D_C7_T6_vs_C7_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T6_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T6_vs_C7_T0_FC.txt", sep="") %>% 
                                    filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T6_vs_C7_T0_sig)

GO_terms_H2D_C7_T6_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T6_vs_C7_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C7_T6_vs_C7_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C7_T6_vs_C7_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T6_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T6_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C7_T12_vs_C7_T0 ----
H2D_C7_T12_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T12", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_H2D_C7_T12_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T12_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T12_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/List_genes_H2D_C7_T12_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/")
input="List_genes_H2D_C7_T12_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T12_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T12", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T12_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/BP_List_genes_H2D_C7_T12_vs_C7_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C7_T12_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T12_vs_C7_T0   <- BP_List_genes_H2D_C7_T12_vs_C7_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T12_vs_C7_T0  <- BP_List_genes_H2D_C7_T12_vs_C7_T0 %>% 
                                      left_join(gene_counts_H2D_C7_T12_vs_C7_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T12_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T12_vs_C7_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T12_vs_C7_T0_sig)

GO_terms_H2D_C7_T12_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T12_vs_C7_T0_FC.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_H2D_C7_T12_vs_C7_T0, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C7_T12_vs_C7_T0") %>% 
  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T12_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T12_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")



# H2D_C7_T24_vs_C7_T0 ----
H2D_C7_T24_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T24", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()                          %>% 
  rownames_to_column(var = "gene")         %>% 
  filter(!is.na(padj))                     %>% 
  dplyr::select("gene", "log2FoldChange")  %>% 
  dplyr::rename(Gene=gene)

List_genes_H2D_C7_T24_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T24_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T24_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/List_genes_H2D_C7_T24_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/")
input="List_genes_H2D_C7_T24_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T24_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T24", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()  %>% 
  rownames_to_column(var = "Gene") %>% 
  dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T24_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/BP_List_genes_H2D_C7_T24_vs_C7_T0_FC.txt") %>% 
  dplyr::rename("Gene" = "seq") %>% 
  filter(lev != -1)             %>% 
  dplyr::select(-value)         %>% 
  left_join(H2D_C7_T24_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T24_vs_C7_T0   <- BP_List_genes_H2D_C7_T24_vs_C7_T0 %>%
  dplyr::group_by(term)             %>% 
  dplyr::summarize(nseqs_relative = dplyr::n(),
                   nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T24_vs_C7_T0  <- BP_List_genes_H2D_C7_T24_vs_C7_T0 %>% 
  left_join(gene_counts_H2D_C7_T24_vs_C7_T0, by = "term") %>% 
  mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T24_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T24_vs_C7_T0_FC.txt", sep="") %>% 
  filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T24_vs_C7_T0_sig)

GO_terms_H2D_C7_T24_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T24_vs_C7_T0_FC.txt", sep="") %>% 
  mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
  mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
  left_join(BP_List_genes_H2D_C7_T24_vs_C7_T0, by = c("term", "name")) %>% 
  mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
  mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
  mutate(Comparison = "H2D_C7_T24_vs_C7_T0") %>% 
  left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T24_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T24_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C7_T3_vs_C7_T0 ----
H2D_C7_T3_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T3", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C7_T3_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T3_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T3_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/List_genes_H2D_C7_T3_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/")
input="List_genes_H2D_C7_T3_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T3_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T3", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T3_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/BP_List_genes_H2D_C7_T3_vs_C7_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C7_T3_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T3_vs_C7_T0   <- BP_List_genes_H2D_C7_T3_vs_C7_T0 %>%
                                    dplyr::group_by(term)             %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T3_vs_C7_T0  <- BP_List_genes_H2D_C7_T3_vs_C7_T0 %>% 
                                     left_join(gene_counts_H2D_C7_T3_vs_C7_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T3_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T3_vs_C7_T0_FC.txt", sep="") %>% 
                                    filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T3_vs_C7_T0_sig)

GO_terms_H2D_C7_T3_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T3_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T3_vs_C7_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C7_T3_vs_C7_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C7_T3_vs_C7_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T3_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T3_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C7_T6_vs_C7_T0 ----
H2D_C7_T6_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T6", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C7_T6_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T6_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T6_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/List_genes_H2D_C7_T6_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/")
input="List_genes_H2D_C7_T6_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T6_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T6", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T6_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/BP_List_genes_H2D_C7_T6_vs_C7_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C7_T6_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T6_vs_C7_T0   <- BP_List_genes_H2D_C7_T6_vs_C7_T0 %>%
                                    dplyr::group_by(term)             %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T6_vs_C7_T0  <- BP_List_genes_H2D_C7_T6_vs_C7_T0 %>% 
                                     left_join(gene_counts_H2D_C7_T6_vs_C7_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T6_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T6_vs_C7_T0_FC.txt", sep="") %>% 
                                    filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T6_vs_C7_T0_sig)

GO_terms_H2D_C7_T6_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T6_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T6_vs_C7_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C7_T6_vs_C7_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C7_T6_vs_C7_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T6_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T6_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C7_T12_vs_C7_T0 ----
H2D_C7_T12_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T12", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_H2D_C7_T12_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T12_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T12_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/List_genes_H2D_C7_T12_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/")
input="List_genes_H2D_C7_T12_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T12_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T12", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T12_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/BP_List_genes_H2D_C7_T12_vs_C7_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C7_T12_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T12_vs_C7_T0   <- BP_List_genes_H2D_C7_T12_vs_C7_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T12_vs_C7_T0  <- BP_List_genes_H2D_C7_T12_vs_C7_T0 %>% 
                                      left_join(gene_counts_H2D_C7_T12_vs_C7_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T12_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T12_vs_C7_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T12_vs_C7_T0_sig)

GO_terms_H2D_C7_T12_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T12_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T12_vs_C7_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_H2D_C7_T12_vs_C7_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "H2D_C7_T12_vs_C7_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T12_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T12_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")



# H2D_C7_T24_vs_C7_T0 ----
H2D_C7_T24_vs_C7_T0_FC    <- results(dds, contrast=c("group", "H2D_C7_T24", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
  as.data.frame()                          %>% 
  rownames_to_column(var = "gene")         %>% 
  filter(!is.na(padj))                     %>% 
  dplyr::select("gene", "log2FoldChange")  %>% 
  dplyr::rename(Gene=gene)

List_genes_H2D_C7_T24_vs_C7_T0_FC               <- left_join(List_genes, H2D_C7_T24_vs_C7_T0_FC)
write.table(List_genes_H2D_C7_T24_vs_C7_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/List_genes_H2D_C7_T24_vs_C7_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/")
input="List_genes_H2D_C7_T24_vs_C7_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C7_T24_vs_C7_T0  <- results(dds, contrast=c("group", "H2D_C7_T24", "H2D_C7_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C7_T24_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/BP_List_genes_H2D_C7_T24_vs_C7_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C7_T24_vs_C7_T0, by = "Gene") 

gene_counts_H2D_C7_T24_vs_C7_T0   <- BP_List_genes_H2D_C7_T24_vs_C7_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C7_T24_vs_C7_T0  <- BP_List_genes_H2D_C7_T24_vs_C7_T0 %>% 
                                      left_join(gene_counts_H2D_C7_T24_vs_C7_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C7_T24_vs_C7_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T24_vs_C7_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C7_T24_vs_C7_T0_sig)

GO_terms_H2D_C7_T24_vs_C7_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C7_T24_vs_C7_T0/MWU_BP_List_genes_H2D_C7_T24_vs_C7_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_H2D_C7_T24_vs_C7_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "H2D_C7_T24_vs_C7_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C7_T24_vs_C7_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C7_T24_vs_C7_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

# H2D_C8_T3_vs_C8_T0 ----
H2D_C8_T3_vs_C8_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T3", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C8_T3_vs_C8_T0_FC               <- left_join(List_genes, H2D_C8_T3_vs_C8_T0_FC)
write.table(List_genes_H2D_C8_T3_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T3_vs_C8_T0/List_genes_H2D_C8_T3_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T3_vs_C8_T0/")
input="List_genes_H2D_C8_T3_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C8_T3_vs_C8_T0  <- results(dds, contrast=c("group", "H2D_C8_T3", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C8_T3_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T3_vs_C8_T0/BP_List_genes_H2D_C8_T3_vs_C8_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C8_T3_vs_C8_T0, by = "Gene") 

gene_counts_H2D_C8_T3_vs_C8_T0   <- BP_List_genes_H2D_C8_T3_vs_C8_T0 %>%
                                    dplyr::group_by(term)             %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T3_vs_C8_T0  <- BP_List_genes_H2D_C8_T3_vs_C8_T0 %>% 
                                     left_join(gene_counts_H2D_C8_T3_vs_C8_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C8_T3_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T3_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T3_vs_C8_T0_FC.txt", sep="") %>% 
                                    filter(p.adj <= 0.05)
dim(GO_terms_H2D_C8_T3_vs_C8_T0_sig)

GO_terms_H2D_C8_T3_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T3_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T3_vs_C8_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C8_T3_vs_C8_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C8_T3_vs_C8_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C8_T3_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C8_T3_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C8_T6_vs_C8_T0 ----
H2D_C8_T6_vs_C8_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T6", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                            as.data.frame()                          %>% 
                            rownames_to_column(var = "gene")         %>% 
                            filter(!is.na(padj))                     %>% 
                            dplyr::select("gene", "log2FoldChange")  %>% 
                            dplyr::rename(Gene=gene)

List_genes_H2D_C8_T6_vs_C8_T0_FC               <- left_join(List_genes, H2D_C8_T6_vs_C8_T0_FC)
write.table(List_genes_H2D_C8_T6_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T6_vs_C8_T0/List_genes_H2D_C8_T6_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T6_vs_C8_T0/")
input="List_genes_H2D_C8_T6_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C8_T6_vs_C8_T0  <- results(dds, contrast=c("group", "H2D_C8_T6", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                       as.data.frame()  %>% 
                       rownames_to_column(var = "Gene") %>% 
                       dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C8_T6_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T6_vs_C8_T0/BP_List_genes_H2D_C8_T6_vs_C8_T0_FC.txt") %>% 
                                     dplyr::rename("Gene" = "seq") %>% 
                                     filter(lev != -1)             %>% 
                                     dplyr::select(-value)         %>% 
                                     left_join(H2D_C8_T6_vs_C8_T0, by = "Gene") 

gene_counts_H2D_C8_T6_vs_C8_T0   <- BP_List_genes_H2D_C8_T6_vs_C8_T0 %>%
                                    dplyr::group_by(term)             %>% 
                                    dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T6_vs_C8_T0  <- BP_List_genes_H2D_C8_T6_vs_C8_T0 %>% 
                                     left_join(gene_counts_H2D_C8_T6_vs_C8_T0, by = "term") %>% 
                                     mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C8_T6_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T6_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T6_vs_C8_T0_FC.txt", sep="") %>% 
                                    filter(p.adj <= 0.05)
dim(GO_terms_H2D_C8_T6_vs_C8_T0_sig)

GO_terms_H2D_C8_T6_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T6_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T6_vs_C8_T0_FC.txt", sep="") %>% 
                                mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                left_join(BP_List_genes_H2D_C8_T6_vs_C8_T0, by = c("term", "name")) %>% 
                                mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                mutate(Comparison = "H2D_C8_T6_vs_C8_T0") %>% 
                                left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C8_T6_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C8_T6_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C8_T12_vs_C8_T0 ----
H2D_C8_T12_vs_C8_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T12", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_H2D_C8_T12_vs_C8_T0_FC               <- left_join(List_genes, H2D_C8_T12_vs_C8_T0_FC)
write.table(List_genes_H2D_C8_T12_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T12_vs_C8_T0/List_genes_H2D_C8_T12_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T12_vs_C8_T0/")
input="List_genes_H2D_C8_T12_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C8_T12_vs_C8_T0  <- results(dds, contrast=c("group", "H2D_C8_T12", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C8_T12_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T12_vs_C8_T0/BP_List_genes_H2D_C8_T12_vs_C8_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C8_T12_vs_C8_T0, by = "Gene") 

gene_counts_H2D_C8_T12_vs_C8_T0   <- BP_List_genes_H2D_C8_T12_vs_C8_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T12_vs_C8_T0  <- BP_List_genes_H2D_C8_T12_vs_C8_T0 %>% 
                                      left_join(gene_counts_H2D_C8_T12_vs_C8_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C8_T12_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T12_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T12_vs_C8_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C8_T12_vs_C8_T0_sig)

GO_terms_H2D_C8_T12_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T12_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T12_vs_C8_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_H2D_C8_T12_vs_C8_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "H2D_C8_T12_vs_C8_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C8_T12_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C8_T12_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


# H2D_C8_T24_vs_C8_T0 ----
H2D_C8_T24_vs_C8_T0_FC    <- results(dds, contrast=c("group", "H2D_C8_T24", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                             as.data.frame()                          %>% 
                             rownames_to_column(var = "gene")         %>% 
                             filter(!is.na(padj))                     %>% 
                             dplyr::select("gene", "log2FoldChange")  %>% 
                             dplyr::rename(Gene=gene)

List_genes_H2D_C8_T24_vs_C8_T0_FC               <- left_join(List_genes, H2D_C8_T24_vs_C8_T0_FC)
write.table(List_genes_H2D_C8_T24_vs_C8_T0_FC, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T24_vs_C8_T0/List_genes_H2D_C8_T24_vs_C8_T0_FC.txt", row.names = F, quote = F, col.names=F, sep = ",")

### Find GO Immunome 
setwd("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T24_vs_C8_T0/")
input="List_genes_H2D_C8_T24_vs_C8_T0_FC.txt" 
goAnnotations="all_go.tab" 
goDatabase="go.obo" 
goDivision="BP" 
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision, perlPath="perl", largest=0.5, smallest=10, clusterCutHeight=0.25) 

H2D_C8_T24_vs_C8_T0  <- results(dds, contrast=c("group", "H2D_C8_T24", "H2D_C8_T0"), alpha = 0.05, lfcThreshold = 0) %>%
                        as.data.frame()  %>% 
                        rownames_to_column(var = "Gene") %>% 
                        dplyr::filter(!is.na(padj)) 

BP_List_genes_H2D_C8_T24_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T24_vs_C8_T0/BP_List_genes_H2D_C8_T24_vs_C8_T0_FC.txt") %>% 
                                      dplyr::rename("Gene" = "seq") %>% 
                                      filter(lev != -1)             %>% 
                                      dplyr::select(-value)         %>% 
                                      left_join(H2D_C8_T24_vs_C8_T0, by = "Gene") 

gene_counts_H2D_C8_T24_vs_C8_T0   <- BP_List_genes_H2D_C8_T24_vs_C8_T0 %>%
                                     dplyr::group_by(term)             %>% 
                                     dplyr::summarize(nseqs_relative = dplyr::n(),
                                                     nseqs_relative_significant = sum(padj <= 0.05, na.rm = TRUE))

BP_List_genes_H2D_C8_T24_vs_C8_T0  <- BP_List_genes_H2D_C8_T24_vs_C8_T0 %>% 
                                      left_join(gene_counts_H2D_C8_T24_vs_C8_T0, by = "term") %>% 
                                      mutate(Gene_significant = ifelse(padj <= 0.05, "Yes", "No"))

GO_terms_H2D_C8_T24_vs_C8_T0_sig  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T24_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T24_vs_C8_T0_FC.txt", sep="") %>% 
                                     filter(p.adj <= 0.05)
dim(GO_terms_H2D_C8_T24_vs_C8_T0_sig)

GO_terms_H2D_C8_T24_vs_C8_T0  <- read.csv("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/H2D_C8_T24_vs_C8_T0/MWU_BP_List_genes_H2D_C8_T24_vs_C8_T0_FC.txt", sep="") %>% 
                                 mutate(GO_significant = ifelse(p.adj <= 0.05, "Yes", "No")) %>% 
                                 mutate(GO_regulation  = ifelse(delta.rank > 0, 1, ifelse(delta.rank < 0, -1, 0))) %>% 
                                 left_join(BP_List_genes_H2D_C8_T24_vs_C8_T0, by = c("term", "name")) %>% 
                                 mutate(Enrichment = nseqs_relative_significant / nseqs) %>% 
                                 mutate(Enrichment = ifelse(Enrichment != 0, Enrichment * GO_regulation, 0)) %>% 
                                 mutate(Comparison = "H2D_C8_T24_vs_C8_T0") %>% 
                                 left_join(Roseta, by = c("Gene"))

write.table(GO_terms_H2D_C8_T24_vs_C8_T0, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D/GO_data_frames/GO_terms_H2D_C8_T24_vs_C8_T0.tsv", row.names = F, quote = F, col.names=T, sep = "\t")


#### Data frames and heatmpas H2D ----
GO_terms_H2D_C3_T3_vs_C3_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C3_T3_vs_C3_T0.tsv")
GO_terms_H2D_C3_T6_vs_C3_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C3_T6_vs_C3_T0.tsv")
GO_terms_H2D_C3_T12_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C3_T12_vs_C3_T0.tsv")
GO_terms_H2D_C3_T24_vs_C3_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C3_T24_vs_C3_T0.tsv")

GO_terms_H2D_C7_T3_vs_C7_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C7_T3_vs_C7_T0.tsv")
GO_terms_H2D_C7_T6_vs_C7_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C7_T6_vs_C7_T0.tsv")
GO_terms_H2D_C7_T12_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C7_T12_vs_C7_T0.tsv")
GO_terms_H2D_C7_T24_vs_C7_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C7_T24_vs_C7_T0.tsv")

GO_terms_H2D_C8_T3_vs_C8_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C8_T3_vs_C8_T0.tsv")
GO_terms_H2D_C8_T6_vs_C8_T0   <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C8_T6_vs_C8_T0.tsv")
GO_terms_H2D_C8_T12_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C8_T12_vs_C8_T0.tsv")
GO_terms_H2D_C8_T24_vs_C8_T0  <- read.delim("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/kinetics/H2D/GO_data_frames/GO_terms_H2D_C8_T24_vs_C8_T0.tsv")


GO_terms_H2D_C3_T3_vs_C3_T0_sig  <- GO_terms_H2D_C3_T3_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C3_T3_vs_C3_T0_sig)  # 101

GO_terms_H2D_C3_T6_vs_C3_T0_sig  <- GO_terms_H2D_C3_T6_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C3_T6_vs_C3_T0_sig)  # 154

GO_terms_H2D_C3_T12_vs_C3_T0_sig  <- GO_terms_H2D_C3_T12_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C3_T12_vs_C3_T0_sig)  # 187

GO_terms_H2D_C3_T24_vs_C3_T0_sig  <- GO_terms_H2D_C3_T24_vs_C3_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C3_T24_vs_C3_T0_sig)  # 161


GO_terms_H2D_C7_T3_vs_C7_T0_sig  <- GO_terms_H2D_C7_T3_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C7_T3_vs_C7_T0_sig)   # 45

GO_terms_H2D_C7_T6_vs_C7_T0_sig  <- GO_terms_H2D_C7_T6_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C7_T6_vs_C7_T0_sig)   # 43

GO_terms_H2D_C7_T12_vs_C7_T0_sig  <- GO_terms_H2D_C7_T12_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C7_T12_vs_C7_T0_sig)  # 48

GO_terms_H2D_C7_T24_vs_C7_T0_sig  <- GO_terms_H2D_C7_T24_vs_C7_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C7_T24_vs_C7_T0_sig)  # 89


GO_terms_H2D_C8_T3_vs_C8_T0_sig  <- GO_terms_H2D_C8_T3_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C8_T3_vs_C8_T0_sig)   # 14

GO_terms_H2D_C8_T6_vs_C8_T0_sig  <- GO_terms_H2D_C8_T6_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C8_T6_vs_C8_T0_sig)   # 35

GO_terms_H2D_C8_T12_vs_C8_T0_sig  <- GO_terms_H2D_C8_T12_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C8_T12_vs_C8_T0_sig)  # 12

GO_terms_H2D_C8_T24_vs_C8_T0_sig  <- GO_terms_H2D_C8_T24_vs_C8_T0 %>%  filter(GO_significant == "Yes") %>% distinct(term)  
dim(GO_terms_H2D_C8_T24_vs_C8_T0_sig)  # 56


H2D_transcriptome_GO <- rbind(GO_terms_H2D_C3_T3_vs_C3_T0,
                              GO_terms_H2D_C3_T6_vs_C3_T0,
                              GO_terms_H2D_C3_T12_vs_C3_T0,
                              GO_terms_H2D_C3_T24_vs_C3_T0,
                             
                              GO_terms_H2D_C7_T3_vs_C7_T0,
                              GO_terms_H2D_C7_T6_vs_C7_T0,
                              GO_terms_H2D_C7_T12_vs_C7_T0,
                              GO_terms_H2D_C7_T24_vs_C7_T0,
                               
                              GO_terms_H2D_C8_T3_vs_C8_T0,
                              GO_terms_H2D_C8_T6_vs_C8_T0,
                              GO_terms_H2D_C8_T12_vs_C8_T0,
                              GO_terms_H2D_C8_T24_vs_C8_T0)  %>% 
  
                      filter(GO_significant == "Yes") %>%  group_by(Comparison)  %>% distinct(name, .keep_all = TRUE) %>%
                      dplyr::select(Comparison, name, Enrichment) %>% 
                      pivot_wider( names_from = Comparison,  values_from = Enrichment) %>%
                      as.data.frame() %>% 
                      mutate_all(~ replace(., is.na(.), 0)) %>%
                      mutate(H2D_C3_T0_vs_C3_T0  = 0,
                             H2D_C7_T0_vs_C7_T0  = 0,
                             H2D_C8_T0_vs_C8_T0  = 0) %>%
                      
                      filter(name != "unknown") %>%
                      filter(!grepl("obsolete", name, ignore.case = TRUE)) %>% 
                      
                      column_to_rownames(var = "name")  %>%
                      dplyr::select(13, 1:4, 14, 5:8, 15, 9:12)



GO_list_H2D_kinetics       <- row.names(H2D_transcriptome_GO) %>% as.data.frame() %>% dplyr::rename("GO_term"=".")
#write.table(GO_list_H2D_kinetics, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_list_H2D.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

GO_list_H2D_kinetics_groups <- GO_list_H2D_kinetics %>% left_join(GO_classified_full) %>%
                               group_by(Category) %>% summarise(n = n())

GO_list_H2D_kinetics_table <- H2D_transcriptome_GO %>% tibble::rownames_to_column("GO_term" ) %>% left_join(GO_classified_full)
write.table(GO_list_H2D_kinetics_table, "C:/Users/avaldivi/Desktop/BMC_Biology_revision/GO_list_H2D_kinetics_table.tsv", row.names = F, quote = F, col.names=T, sep = "\t")

# Heat maps  kinetics ----
# Colum Colors heatmap   
annotation_col <- data.frame(
  Time = c("T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24", "T0", "T3", "T6", "T12", "T24"),
  Age  = c(rep("4 months", 5), rep("16 months", 5), rep("28 months", 5)))

annotation_colors  <-  list(
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"),
  Time = c(T0 = "white",  T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred",  T24 = "red"))

rownames(annotation_col) <- colnames(H2D_transcriptome_GO) 

# Gradient enrichment                        
myBreaks <- c(seq(min(H2D_transcriptome_GO), 0, length.out=ceiling(100/2) + 1),
              seq(max(H2D_transcriptome_GO)/100, max(H2D_transcriptome_GO), length.out=floor(100/2)))
mycolor  <- colorRampPalette(c("blue", "black", "yellow"))(100)

H2D_transcriptome_GO_heatmap <- pheatmap(H2D_transcriptome_GO, 
                                          cluster_cols = F, 
                                          scale = "none",
                                          cluster_rows = T, 
                                          fontsize_row = 5, 
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
                                          gaps_col =  c(5, 10, 15))
H2D_transcriptome_GO_heatmap


# =========================================================
# Category
# =========================================================

mat <- H2D_transcriptome_GO
mat[is.na(mat)] <- 0


annotation_col <- data.frame(
  Time = c("T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24"),
  Age = c(rep("4 months",5),
          rep("16 months",5),
          rep("28 months",5)))

rownames(annotation_col) <- colnames(mat)

annotation_colors <- list(
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"),
  Time = c(T0 = "white",  T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred", T24 = "red" ))

GO_classified_full <- read_excel("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_classified_full.xlsx")

GO_terms <- data.frame(GO_term = rownames(mat))
GO_map <- GO_classified_full %>%  dplyr::select(GO_term, Category)

annotation_row <- GO_terms %>%  left_join(GO_map, by = "GO_term") %>%
  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>%
  column_to_rownames("GO_term")

annotation_row <- annotation_row[rownames(mat), , drop = FALSE]

category_order <- c(
  "Immune system procces",
  "Response to stimulus",
  "Metabolic procces",
  "Biological regulation",
  "Cellular procces",
  "Communication",
  "Transport",
  "Localisation",
  "Developmental procces",
  "Reproduction")
#"Unknown")

annotation_row$Category <- factor(annotation_row$Category, levels = category_order)
ord                     <- order(annotation_row$Category)

mat_ord <- mat[ord, ]
annotation_row_ord <- annotation_row[ord, , drop = FALSE]
gaps_row           <- cumsum(table(annotation_row_ord$Category))


category_colors <- c(
  # 🔴 TOP 3
  "Immune system procces" = "purple",
  "Response to stimulus"  = "darkorange",
  "Metabolic procces"     = "#4DAF4A",
  # 🌿 OTHERS
  "Biological regulation" = "#8C8C8C",
  "Cellular procces"      = "#A6A6A6",
  "Communication"         = "#BFBFBF",
  "Transport"             = "#D9D9D9",
  "Localisation"          = "#E0E0E0",
  "Developmental procces" = "#E6E6E6",
  "Reproduction"          = "#F2F2F2")
#"Unknown"               = "grey85")

annotation_colors$Category <- category_colors

mycolor  <- colorRampPalette(c("blue","black","yellow"))(100)
myBreaks <- c(seq(min(mat_ord), 0, length.out = 50), seq(max(mat_ord)/100, max(mat_ord), length.out = 50))

H2D_transcriptome_category_heatmap <- pheatmap(
  mat_ord,
  
  cluster_rows = T,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  
  scale = "none",
  color = mycolor,
  breaks = myBreaks,
  
  border_color = "black",
  
  annotation_col    = annotation_col,
  annotation_row    = annotation_row_ord,
  annotation_colors = annotation_colors,
  gaps_row = gaps_row,
  gaps_col = c(5,10,15),
  fontsize_row = 5,
  
  legend = T,              
  annotation_legend = FALSE)

H2D_transcriptome_category_heatmap
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D_transcriptome_category_heatmap.tiff", plot = H2D_transcriptome_category_heatmap$gtable,  width = 8,    height = 12,  dpi = 600)



# Only immunity and stress
H2D_transcriptome_GO_immune <- GO_list_H2D_kinetics_table %>% 
                               filter(Category=="Immune system procces" | Category=="Response to stimulus") %>% 
                               dplyr::select(-Category)  %>%
                               tibble::column_to_rownames("GO_term")


mat_immune <- H2D_transcriptome_GO_immune
mat_immune[is.na(mat_immune)] <- 0

annotation_col <- data.frame(
  Time = c("T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24",
           "T0","T3","T6","T12","T24"),
  Age = c(rep("4 months",5),
          rep("16 months",5),
          rep("28 months",5)))

rownames(annotation_col) <- colnames(mat_immune)

annotation_colors <- list(
  Age    = c("4 months"   = "#E2F0D9", "16 months"  = "#A9D18E", "28 months"  ="#385700"),
  Time = c(T0 = "white",  T3 = "#CD5C5C", T6 = "#B22222", T12 = "darkred", T24 = "red" ))

GO_classified_full <- read_excel("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/GO_classified_full.xlsx")

GO_terms <- data.frame(GO_term = rownames(mat_immune))
GO_map   <- GO_classified_full %>%  dplyr::select(GO_term, Category)

annotation_row <- GO_terms %>%  left_join(GO_map, by = "GO_term") %>%
                  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>%
                  column_to_rownames("GO_term")

annotation_row <- annotation_row[rownames(mat_immune), , drop = FALSE]

category_order <- c("Immune system procces",  "Response to stimulus")


annotation_row$Category <- factor(annotation_row$Category, levels = category_order)
ord                     <- order(annotation_row$Category)

mat_ord <- mat_immune[ord, ]
annotation_row_ord <- annotation_row[ord, , drop = FALSE]
gaps_row           <- cumsum(table(annotation_row_ord$Category))


category_colors <- c(
  "Immune system procces" = "purple",
  "Response to stimulus"  = "darkorange")


annotation_colors$Category <- category_colors

mycolor  <- colorRampPalette(c("blue","black","yellow"))(100)
myBreaks <- c(seq(min(mat_ord), 0, length.out = 50), seq(max(mat_ord)/100, max(mat_ord), length.out = 50))

H2D_transcriptome_category_Infection_heatmap <- pheatmap(
  mat_ord,
  
  cluster_rows = F,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = T,
  
  scale = "none",
  color = mycolor,
  breaks = myBreaks,
  
  border_color = "black",
  
  annotation_col    = annotation_col,
  annotation_row    = annotation_row_ord,
  annotation_colors = annotation_colors,
  gaps_row = gaps_row,
  gaps_col = c(5,10,15),
  fontsize_row = 14,
  
  legend = T,              
  annotation_legend = FALSE)

H2D_transcriptome_category_Infection_heatmap
ggsave("C:/Users/avaldivi/Desktop/BMC_Biology_revision/transcriptomics/H2D_transcriptome_category_Infection_heatmap.tiff", plot = H2D_transcriptome_category_Infection_heatmap$gtable,  width = 8,    height = 6,  dpi = 600)






