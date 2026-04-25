

library(data.table)
library(tidyverse)
library(dplyr)
library(methylKit)
library(readxl)
library(DESeq2)

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


# GENE Body methyaltion ----

# EPIGENETICS F14R ----
load("D:/Decicomp/R/Melthilkit projects/coverage8X/meth_F14R.RData")
meth_F14R  # <- methylKit::unite(normalize.myobj_F14R, destrand=T)

Coordenates_all_cpg_meth_F14R        <- data.frame(chr   = meth_F14R$chr, 
                                                   start = meth_F14R$start, 
                                                   end   = meth_F14R$end)

Methylation_meth_F14R                <- methylKit::percMethylation(meth_F14R, rowids = F)
Betas_all_cpg_meth_F14R              <- cbind(Coordenates_all_cpg_meth_F14R, Methylation_meth_F14R)
setDT(Betas_all_cpg_meth_F14R,  key = c("chr", "start", "end"))

Genes_body_region         <- read.delim("D:/Gestinnov/1) Oyster proyect/R in datarmor/Immuno_genes_list/Genes_body_region.txt")
Genes_body_region         <- Genes_body_region[,c(4,1,2,3)]
setDT(Genes_body_region,  key = c("chr", "start", "end"))

Meth_F14R_in_genes <- data.table::foverlaps(  Betas_all_cpg_meth_F14R,   Genes_body_region,    nomatch = NULL ) %>% 
                      dplyr::select(2, 5, 7:24) %>%   as.data.frame()   %>%   
                      tidyr::pivot_longer(
                        cols = c(3:20),
                      names_to = "sample",
                      values_to = "methylation") %>% 
                      mutate(sample = str_replace(sample, "_T0", ""))    %>%
                      mutate(sample = str_replace(sample, "_DNA", ""))  %>%
                      mutate(Age = case_when(
                          grepl("C3", sample) ~ "4",
                          grepl("C7", sample) ~ "16",
                          grepl("C8", sample) ~ "28",
                          TRUE ~ NA_character_))  %>%
  group_by(Gene, Age, i.start) %>%
  summarise(mean_cpg = mean(methylation, na.rm = TRUE),.groups = "drop")


# EPIGENETICS H2D ----
load("D:/Decicomp/R/Melthilkit projects/coverage8X/meth_H2D.RData")
meth_H2D  # <- methylKit::unite(normalize.myobj_H2D, destrand=T)

Coordenates_all_cpg_meth_H2D        <- data.frame(chr   = meth_H2D$chr, 
                                                  start = meth_H2D$start, 
                                                  end   = meth_H2D$end)

Methylation_meth_H2D                <- methylKit::percMethylation(meth_H2D, rowids = F)
Betas_all_cpg_meth_H2D              <- cbind(Coordenates_all_cpg_meth_H2D, Methylation_meth_H2D)
setDT(Betas_all_cpg_meth_H2D,  key = c("chr", "start", "end"))

Genes_body_region         <- read.delim("D:/Gestinnov/1) Oyster proyect/R in datarmor/Immuno_genes_list/Genes_body_region.txt")
Genes_body_region         <- Genes_body_region[,c(4,1,2,3)]
setDT(Genes_body_region,  key = c("chr", "start", "end"))

Meth_H2D_in_genes <- data.table::foverlaps(  Betas_all_cpg_meth_H2D,   Genes_body_region,    nomatch = NULL ) %>% 
  dplyr::select(2,5, 7:21) %>%   as.data.frame()   %>%   
  tidyr::pivot_longer(
    cols = c(3:17),
    names_to = "sample",
    values_to = "methylation") %>% 
  mutate(sample = str_replace(sample, "_T0", ""))    %>%
  mutate(sample = str_replace(sample, "_DNA", ""))  %>%
  mutate(Age = case_when(
    grepl("C3", sample) ~ "4",
    grepl("C7", sample) ~ "16",
    grepl("C8", sample) ~ "28",
    TRUE ~ NA_character_))  %>%
  group_by(Gene, Age, i.start) %>%
  summarise(mean_cpg = mean(methylation, na.rm = TRUE),.groups = "drop")



Gene_I <- "G22089"

# F14R
Gene_methylation_F14R_CpG  <- Meth_F14R_in_genes %>% filter(Gene==Gene_I) #%>%
                              #filter(Age != 16)
                              #filter(!all(mean_cpg == 0, na.rm = TRUE)) %>%
                              #ungroup()

Lenght_CpG_F14R <-length(Gene_methylation_F14R_CpG$i.start)
Lenght_CpG_F14R

Gene_methylation_F14R_CpG_mean <- Gene_methylation_F14R_CpG %>% 
                                  group_by(Gene, Age) %>%
                                  summarise(mean_gene_methylation = mean(mean_cpg, na.rm = TRUE),
                                  sd_gene_methylation   = sd(mean_cpg, na.rm = TRUE), .groups = "drop")


#Gene_body_methyaltion_F14R <- ggplot(Gene_methylation_F14R_CpG, aes(x = i.start, y = mean_cpg, color = Age, group = Age)) +
                               #geom_line(linewidth = 1.5) +
                               #scale_color_manual(values = c("4"  = "#dae3f3",  "16" = "#8faadc",  "28" = "#2f5597")) +
                               #labs(x = "Genomic position (CpG)",  y = "Methylation (%)",  color = "Age (months)" ) +
                               #scale_x_continuous(expand = c(0, 0)) +
                               #scale_y_continuous(expand = c(0, 0)) + 
                              
#Gene_body_methyaltion_F14R
                      
#Gene_body_methyaltion_F14R_smooth <- ggplot(Gene_methylation_F14R_CpG, aes(x = i.start, y = mean_cpg, color = Age, fill = Age)) +
                                      #geom_point(alpha = 0.85, size = 1.75) +
                                      #geom_smooth(method = "loess", se = F, linewidth = 1.2, span = 1.0, alpha = 0.5) +
                                      # scale_color_manual(values = c("4"  = "#dae3f3",  "16" = "#8faadc",  "28" = "#2f5597")) +
                                      #scale_fill_manual(values = c("4"   = "#dae3f3",  "16" = "#8faadc",  "28" = "#2f5597")) +
                                      #scale_x_continuous(expand = c(0, 0)) +
                                      #scale_y_continuous(expand = c(0, 0)) +
                                      #labs(x = "Genomic position (CpG)",  y = "Methylation (%)",  color = "Age (months)" ) +
                                      #annotate("text", x = Inf, y = Inf,  label = paste0(Gene_I, "\nCpGs: ", Lenght_CpG_F14R),  hjust = 1.1,  vjust = 1.5, size = 4 ) +
                                      #Style_format_theme

#Gene_body_methyaltion_F14R_smooth


x_min    <- min(Gene_methylation_F14R_CpG$i.start, na.rm = TRUE)
x_max    <- max(Gene_methylation_F14R_CpG$i.start, na.rm = TRUE)
x_breaks <- sort(unique(Gene_methylation_F14R_CpG$i.start))

Gene_body_methyaltion_F14R_smooth_2 <- ggplot(Gene_methylation_F14R_CpG, aes(x=i.start, y=mean_cpg, color=Age, fill=Age)) +
                                       # geom_point(alpha = 0.85, size = 1.75) +                
                                       geom_smooth(method="loess", se= F, linewidth=1.2, span=1, alpha=0.1) +
                                       scale_color_manual(values=c("4"="#dae3f3","16"="#8faadc","28"="#2f5597")) +
                                       scale_fill_manual(values=c( "4"="#dae3f3","16"="#8faadc","28"="#2f5597")) +
                                       scale_x_continuous(limits=c(x_min, x_max), breaks=x_breaks, expand=c(0,0)) +
                                       scale_y_continuous(expand=c(0,0)) +
                                       coord_cartesian(ylim = c(0, 100)) +
                                       labs(x="Genomic position (CpG)", y="Methylation (%)") +
                                       #annotate("text", x=Inf, y=Inf, label=paste0(Gene_I,"\nCpGs: ", Lenght_CpG_F14R), hjust=1.1, vjust=1.5, size=4) +
                                       Style_format_theme +
                                       theme(axis.text.x = element_text(color = "white"))

Gene_body_methyaltion_F14R_smooth_2
ggsave("D:/Decicomp/R/Figures/Gene_body_methyaltion_F14R_smooth_2.tiff",  Gene_body_methyaltion_F14R_smooth_2,  width = 8,   height = 5, dpi = 600)


Gene_methylation_F14R_CpG$Age        <- factor(Gene_methylation_F14R_CpG$Age, levels = c("4", "16", "28"))

Gene_body_methyaltion_F14R_mean_plot <- ggplot(Gene_methylation_F14R_CpG,aes(x = Age, y = mean_cpg, fill = Age)) + 
                                        geom_boxplot() +
                                        scale_fill_manual(values = c("4"   = "#dae3f3",  "16" = "#8faadc",  "28" = "#2f5597")) +
                                        scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                                        labs(x = "Age (months)", y = "Mean Methylation (%)") +
                                        Style_format_theme
Gene_body_methyaltion_F14R_mean_plot
ggsave("D:/Decicomp/R/Figures/Gene_body_methyaltion_F14R_mean_plot.tiff",  Gene_body_methyaltion_F14R_mean_plot,  width = 5,   height = 5, dpi = 600)


# H2D
Gene_methylation_H2D_CpG  <- Meth_H2D_in_genes  %>% filter(Gene==Gene_I) #%>%
                             #filter(Age != 16)
                             #filter(!all(mean_cpg == 0, na.rm = TRUE)) %>%
                             #ungroup()
  

Lenght_CpG_H2D <-length(Gene_methylation_H2D_CpG$i.start)
Lenght_CpG_H2D

#Gene_body_methyaltion_H2D <- ggplot(Gene_methylation_H2D_CpG, aes(x = i.start, y = mean_cpg, color = Age, group = Age)) +
                              #geom_line(linewidth = 1.5) +
                              #scale_color_manual(values = c("4"  = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                              #labs(x = "Genomic position (CpG)",  y = "Methylation (%)",  color = "Age (months)" ) +
                              #scale_x_continuous(expand = c(0, 0)) +
                              #scale_y_continuous(expand = c(0, 0)) +
                              #Style_format_theme
#Gene_body_methyaltion_H2D

#Gene_body_methyaltion_H2D_smooth <- ggplot(Gene_methylation_H2D_CpG, aes(x = i.start, y = mean_cpg, color = Age, fill = Age)) +
                                    #geom_point(alpha = 0.5, size = 1.75) +
                                    #geom_smooth(method = "loess", se = F, linewidth = 1.2, span = 1.0, alpha = 0.5) +
                                    #scale_color_manual(values = c("4" = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                                    #scale_fill_manual(values  = c("4" = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                                    #scale_x_continuous(expand = c(0, 0)) +
                                    #scale_y_continuous(expand = c(0, 0)) +
                                    #labs(x = "Genomic position (CpG)",  y = "Methylation (%)",  color = "Age (months)" ) +
                                    #annotate("text", x = Inf, y = Inf,  label = paste0(Gene_I, "\nCpGs: ", Lenght_CpG_H2D),  hjust = 1.1,  vjust = 1.5, size = 4 ) +
                                    #Style_format_theme

#Gene_body_methyaltion_H2D_smooth

x_min    <- min(Gene_methylation_H2D_CpG$i.start, na.rm = TRUE)
x_max    <- max(Gene_methylation_H2D_CpG$i.start, na.rm = TRUE)
x_breaks <- sort(unique(Gene_methylation_H2D_CpG$i.start))

Gene_body_methyaltion_H2D_smooth_2 <- ggplot(Gene_methylation_H2D_CpG, aes(x=i.start, y=mean_cpg, color=Age, fill=Age)) +
                                      geom_smooth(method="loess", se= F, linewidth=1.2, span=1, alpha=0.1) +
                                      scale_color_manual(values = c("4" = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                                      scale_fill_manual(values  = c("4" = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                                      scale_x_continuous(limits=c(x_min, x_max), breaks=x_breaks, expand=c(0,0)) +
                                      scale_y_continuous(expand=c(0,0)) +
                                      coord_cartesian(ylim = c(0, 100)) +
                                      labs(x="Genomic position (CpG)", y="Methylation (%)") +
                                      #annotate("text", x=Inf, y=Inf, label=paste0(Gene_I,"\nCpGs: ", Lenght_CpG_H2D), hjust=1.1, vjust=1.5, size=4) +
                                      Style_format_theme +
                                      theme(axis.text.x = element_text(color = "white"))

Gene_body_methyaltion_H2D_smooth_2
ggsave("D:/Decicomp/R/Figures/Gene_body_methyaltion_H2D_smooth_2.tiff",  Gene_body_methyaltion_H2D_smooth_2,  width = 8,   height = 5, dpi = 600)



Gene_methylation_H2D_CpG$Age        <- factor(Gene_methylation_H2D_CpG$Age, levels = c("4", "16", "28"))

Gene_body_methyaltion_H2D_mean_plot <- ggplot(Gene_methylation_H2D_CpG, aes(x = Age, y = mean_cpg, fill = Age)) + 
                                       geom_boxplot() +
                                       scale_fill_manual(values  = c("4" = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                                       scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
                                       labs(x = "Age (months)", y = "Mean Methylation (%)") +
                                       Style_format_theme
Gene_body_methyaltion_H2D_mean_plot
ggsave("D:/Decicomp/R/Figures/Gene_body_methyaltion_H2D_mean_plot.tiff",  Gene_body_methyaltion_H2D_mean_plot,  width = 5,   height = 5, dpi = 600)



# RNA seq  ----
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

# TRANSCRIPTOMICS F14R ----
dds_norm      <- vst(dds, blind = F)

F14R_dds_norm <- assay(dds_norm) %>% as.data.frame() %>% 
                 dplyr::select(contains("F14R_C3_T0"),
                               contains("F14R_C7_T0"), 
                               contains("F14R_C8_T0")) %>% 
                as.data.frame()  %>% 
                tibble::rownames_to_column ("Gene")  %>% 
               tidyr::pivot_longer(
                cols = c(2:18),
                names_to = "sample",
                values_to = "expression") %>% 
                mutate(sample = str_replace(sample, "_T0", ""))  %>%
               mutate(Age = case_when(
                  grepl("C3", sample) ~ "4",
                  grepl("C7", sample) ~ "16",
                  grepl("C8", sample) ~ "28",
                  TRUE ~ NA_character_))  %>%
                mutate(Age = factor(Age, levels = c("4", "16", "28")))


Gene_expression_F14R  <- F14R_dds_norm %>% filter(Gene==Gene_I)

Gene_expression_F14R_plot <- ggplot(Gene_expression_F14R, aes(x = Age, y = expression, fill = Age)) +
                             geom_boxplot(alpha = 0.7, outlier.shape = NA) +
                             geom_jitter(shape = 21, width = 0.15, size = 6, alpha = 0.8, stroke = 1) +
                             scale_fill_manual(values = c("4"="#dae3f3","16"="#8faadc","28"="#2f5597")) +
                             coord_cartesian(ylim = c(9, 10.5)) +
                             scale_y_continuous(expand = c(0, 0))+
                             labs(x = "Age (months)",  y = expression("Gene expression ("*Log[2]*")")) +
                             Style_format_theme
Gene_expression_F14R_plot
ggsave("D:/Decicomp/R/Figures/Gene_expression_F14R_plot.tiff",  Gene_expression_F14R_plot,  width = 5,   height = 5, dpi = 600)



H2D_dds_norm <- assay(dds_norm) %>% as.data.frame() %>% 
                dplyr::select(contains("H2D_C3_T0"),
                              contains("H2D_C7_T0"), 
                              contains("H2D_C8_T0")) %>% 
                as.data.frame()  %>% 
                tibble::rownames_to_column ("Gene")  %>% 
                tidyr::pivot_longer(
                  cols = c(2:18),
                  names_to = "sample",
                  values_to = "expression") %>% 
                mutate(sample = str_replace(sample, "_T0", ""))  %>%
                mutate(Age = case_when(
                  grepl("C3", sample) ~ "4",
                  grepl("C7", sample) ~ "16",
                  grepl("C8", sample) ~ "28",
                  TRUE ~ NA_character_))  %>%
                mutate(Age = factor(Age, levels = c("4", "16", "28")))


Gene_expression_H2D  <- H2D_dds_norm %>% filter(Gene==Gene_I)

Gene_expression_H2D_plot <- ggplot(Gene_expression_H2D, aes(x = Age, y = expression, fill = Age)) +
                             geom_boxplot(alpha = 0.7, outlier.shape = NA) +
                             geom_jitter(shape = 21, width = 0.15, size = 6, alpha = 0.8, stroke = 1) +
                             scale_fill_manual(values  = c("4" = "#E2F0D9",  "16" = "#A9D18E",  "28" = "#385700")) +
                            coord_cartesian(ylim = c(9, 10.5)) +
                            scale_y_continuous(expand = c(0, 0))+
                             labs(x = "Age (months)",  y = expression("Gene expression ("*Log[2]*")")) +
                             Style_format_theme
Gene_expression_H2D_plot
ggsave("D:/Decicomp/R/Figures/Gene_expression_H2D_plot.tiff",  Gene_expression_H2D_plot,  width = 5,   height = 5, dpi = 600)
