#  Load viral DECICOMP
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(scales)
library(ggplot2)
library(patchwork)
options(scipen=999)

setwd("D:/Decicomp/R")

Style_format_theme <- theme(
  ## Title of the axis
  axis.title.x       = element_text(color="black", size=26),
  axis.title.y       = element_text(color="black", size=26),
  ## Text in the axis
  axis.text          = element_text(color="black", size=22),
  axis.text.x        = element_text(color="black", size=22),
  axis.text.y        = element_text(color="black", size=22),
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

### DNA viral load ----
Load_viral_decicomp          <- read_excel("D:/Decicomp/Excel Samples DNA to Tremblade Alejandro mars 2023 V2.xlsx", sheet = "Data_viral_load")
Load_viral_decicomp$Copies_2 <- as.numeric(gsub(",", ".", Load_viral_decicomp$Copies))

Load_viral_decicomp <- Load_viral_decicomp  %>%
  mutate(Family = as.factor(Family),
         Age = as.factor(Age),
         Time = as.factor(Time),
         Copies_2 = log10(Copies_2 + 1)) %>%
         dplyr::select(Family, Age, Time, Copies_2)

stats_means_family <- Load_viral_decicomp %>%
  group_by(Family, Age, Time) %>%
  summarise(
    Mean = mean(Copies_2, na.rm = TRUE),
    SD = sd(Copies_2, na.rm = TRUE),
    Num_samples = sum(!is.na(Copies_2)),
    SEM = sd(Copies_2, na.rm = TRUE) / sqrt(sum(!is.na(Copies_2))))

F14R_dna_virus       <- stats_means_family %>% filter(Family=="F14R") %>% dplyr::select(1:4,7)
F14R_dna_virus$Time  <- as.numeric(as.character(F14R_dna_virus$Time))
F14R_dna_virus$Mean  <- as.numeric(F14R_dna_virus$Mean)

F14R_dna_load_graph <- ggplot(F14R_dna_virus, aes(x = Time, y = Mean, fill=factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Age), size=1) +
 # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("4"="#dae3f3", "16"="#8faadc", "28"="#2f5597")) +
  scale_fill_manual(values =  c("4"="#dae3f3",  "16"="#8faadc", "28"="#2f5597")) +
  Style_format_theme
 F14R_dna_load_graph

# F14R_dna_load_graph
# tiff("F14R_dna_load_graph.tiff", units="in", width=10, height=8, res=400)
# F14R_dna_load_graph
# dev.off()

H2D_dna_virus      <- stats_means_family %>% filter(Family=="H2D") %>% dplyr::select(1:4,7)
H2D_dna_virus$Time <- as.numeric(as.character(H2D_dna_virus$Time))
H2D_dna_virus$Mean <- as.numeric(H2D_dna_virus$Mean)

H2D_dna_load_graph <- ggplot(H2D_dna_virus, aes(x = Time, y = Mean,  fill = factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Age), size=1) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  scale_fill_manual(values  = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  Style_format_theme
 H2D_dna_load_graph

# H2D_dna_load_graph
# tiff("H2D_dna_load_graph.tiff", units="in", width=10, height=8, res=400)
# H2D_dna_load_graph
# dev.off()


Age_4_months      <- stats_means_family %>% filter(Age=="4" )
Age_4_months$Time <- as.numeric(as.character(Age_4_months$Time))
Age_4_months$Mean <- as.numeric(Age_4_months$Mean)

Age_4_months_load_graph <- ggplot(Age_4_months, aes(x = Time, y = Mean,  fill = factor(Family), color = factor(Family))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Family), size=1) +
 # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("F14R"="#dae3f3", "H2D" ="#E2F0D9")) +
  scale_fill_manual(values  = c("F14R"="#dae3f3", "H2D" ="#E2F0D9")) +
  Style_format_theme
  Age_4_months_load_graph

# Age_4_months_load_graph
# tiff("Age_4_months_load_graph.tiff", units="in", width=10, height=8, res=400)
# Age_4_months_load_graph
# dev.off()


Age_16_months      <- stats_means_family %>% filter(Age=="16" )
Age_16_months$Time <- as.numeric(as.character(Age_16_months$Time))
Age_16_months$Mean <- as.numeric(Age_16_months$Mean)

Age_16_months_load_graph <- ggplot(Age_16_months, aes(x = Time, y = Mean,  fill = factor(Family), color = factor(Family))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Family), size=1) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("F14R"="#8faadc", "H2D" ="#A9D18E")) +
  scale_fill_manual(values  = c("F14R"="#8faadc", "H2D" ="#A9D18E")) +
  Style_format_theme
  Age_16_months_load_graph

# Age_16_months_load_graph
# tiff("Age_16_months_load_graph.tiff", units="in", width=10, height=8, res=1600)
# Age_16_months_load_graph
# dev.off()


Age_28_months      <- stats_means_family  %>%   filter(Age=="28" )
Age_28_months$Time <- as.numeric(as.character(Age_28_months$Time))
Age_28_months$Mean <- as.numeric(Age_28_months$Mean)

Age_28_months_load_graph <- ggplot(Age_28_months, aes(x = Time, y = Mean,  fill = factor(Family), color = factor(Family))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Family), size=1) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("F14R"="#2f5597", "H2D" ="#385700")) +
  scale_fill_manual(values  = c("F14R"="#2f5597", "H2D" ="#385700")) +
  Style_format_theme
  Age_28_months_load_graph

# Age_28_months_load_graph
# tiff("Age_28_months_load_graph.tiff", units="in", width=10, height=8, res=1600)
# Age_28_months_load_graph
# dev.off()

# --- --- --- --- --- --- --- ---

Load_viral_decicomp          <- read_excel("D:/Decicomp/Excel Samples DNA to Tremblade Alejandro mars 2023 V2.xlsx", sheet = "Data_viral_load")
Load_viral_decicomp$Copies_2 <- as.numeric(gsub(",", ".", Load_viral_decicomp$Copies))

Load_viral_decicomp <- Load_viral_decicomp  %>%
  mutate(Family = as.factor(Family),
         Age = as.factor(Age),
         Time = as.factor(Time),
         Copies_2 = log10(Copies_2 + 1)) %>%
         dplyr::select(Family, Age, Time, Copies_2)

stats_means_all <- Load_viral_decicomp %>%
  group_by(Age, Time) %>%
  summarise(
    Mean = mean(Copies_2, na.rm = TRUE),
    SD = sd(Copies_2, na.rm = TRUE),
    Num_samples = sum(!is.na(Copies_2)),
    SEM = sd(Copies_2, na.rm = TRUE) / sqrt(sum(!is.na(Copies_2)))  )

Age_4_months_all      <- stats_means_all  %>% filter(Age=="4" )
Age_4_months_all$Time <- as.numeric(as.character(Age_4_months_all$Time))
Age_4_months_all$Mean <- as.numeric(Age_4_months_all$Mean)

Age_4_months_all_load_graph <- ggplot(Age_4_months_all, aes(x = Time, y = Mean,  fill = factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1, alpha=1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1, alpha=1) +
  geom_line(aes(group = Age), size=1, alpha=1) +
  #geom_smooth(method = "lm", se = F, linetype = "solid", color = "red", size = 1) +  
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("red")) +
  scale_fill_manual(values  = c("red")) +
  Style_format_theme
  Age_4_months_all_load_graph

# Age_4_months_all_load_graph
# tiff("Age_4_months_all_load_graph.tiff", units="in", width=10, height=8, res=1600)
# Age_4_months_all_load_graph
# dev.off()


Age_16_months_all      <- stats_means_all  %>% filter(Age=="16")
Age_16_months_all$Time <- as.numeric(as.character(Age_16_months_all$Time))
Age_16_months_all$Mean <- as.numeric(Age_16_months_all$Mean)

Age_16_months_all_load_graph <- ggplot(Age_16_months_all, aes(x = Time, y = Mean,  fill = factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1, alpha=1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1, alpha=1) +
  geom_line(aes(group = Age), size=1, alpha=1) +
  # geom_smooth(method = "lm", se = F, linetype = "solid", color = "darkred", size = 1) +  
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("darkred")) +
  scale_fill_manual(values  = c("darkred")) +
  Style_format_theme
Age_16_months_all_load_graph

# Age_16_months_all_load_graph
# tiff("Age_16_months_all_load_graph.tiff", units="in", width=10, height=8, res=1600)
# Age_16_months_all_load_graph
# dev.off()


Age_28_months_all      <- stats_means_all %>% filter(Age=="28" )
Age_28_months_all$Time <- as.numeric(as.character(Age_28_months_all$Time))
Age_28_months_all$Mean <- as.numeric(Age_28_months_all$Mean)

Age_28_months_all_load_graph <- ggplot(Age_28_months_all, aes(x = Time, y = Mean,  fill = factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1, alpha=1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1, alpha=1) +
  geom_line(aes(group = Age), size=1, alpha=1) +
  #geom_smooth(method = "lm", se = F, linetype = "solid", color = "salmon", size = 1) +  
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("salmon")) +
  scale_fill_manual(values  = c("salmon")) +
  Style_format_theme
  Age_28_months_all_load_graph

# Age_28_months_all_load_graph
# tiff("Age_28_months_all_load_graph.tiff", units="in", width=10, height=8, res=1600)
# Age_28_months_all_load_graph
# dev.off()


stats_means_all
stats_means_all$Time <- as.numeric(as.character(stats_means_all$Time))
stats_means_all$Mean <- as.numeric(stats_means_all$Mean)

all_load_graph <- ggplot(stats_means_all, aes(x = Time, y = Mean,  fill = factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1, alpha=1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1, alpha=1) +
  geom_line(aes(group = Age), size=1, alpha=1) +
  #geom_smooth(method = "lm", se = F, linetype = "solid", color = "salmon", size = 1) +  
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar DNA load")) +
  scale_color_manual(values = c("red","darkred","salmon")) +
  scale_fill_manual(values  = c("red","darkred","salmon")) +
  Style_format_theme
  all_load_graph

# all_load_graph
# tiff("all_load_graph.tiff", units="in", width=10, height=8, res=400)
# all_load_graph
# dev.off()

# --- --- --- --- --- --- --- ---

  # Plots DNA based on dots Journal  ----
Load_viral_decicomp          <- read_excel("D:/Decicomp/Excel Samples DNA to Tremblade Alejandro mars 2023 V2.xlsx", sheet = "Data_viral_load")
Load_viral_decicomp$Copies_2 <- as.numeric(gsub(",", ".", Load_viral_decicomp$Copies))
  
Load_viral_decicomp_dots <- Load_viral_decicomp  %>%
                            mutate(Family = as.factor(Family),
                            Age = as.factor(Age),
                            Time = as.factor(Time),
                            Copies_2 = log10(Copies_2 + 1)) %>%
                            dplyr::select(Family, Age, Time, Copies_2)
  
F14R_dna_virus_dots <- Load_viral_decicomp_dots %>%
                       dplyr::filter(Family == "F14R") %>%
                       dplyr::mutate(Time = as.numeric(as.character(Time)),
                       Mean = as.numeric(Copies_2))
  
set.seed(123)
F14R_dna_load_points <- ggplot(F14R_dna_virus_dots, aes(x = Time, y = Mean, group = Age)) +
    # puntos con relleno y borde negro
    geom_point(aes(fill = factor(Age)), shape = 21, size = 8, color = "black", alpha = 0.8) +  
    # línea promedio por grupo
    stat_summary(aes(color = factor(Age)), fun = mean, geom = "line", size = 1.5) +  
    scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0,25), breaks = c(0,3,6,12,24), expand = c(0,0)) +
    scale_fill_manual(values = c("4" =  "#dae3f3", "16" = "#8faadc", "28" = "#2f5597")) +
    scale_color_manual(values = c("4" = "#dae3f3", "16" = "#8faadc", "28" = "#2f5597")) +
    labs(x = "Time (hpc)",  y = expression("Genomic units/ng")) +
    Style_format_theme   +
    theme(axis.title.x = element_text(size = 44),  
          axis.title.y = element_text(size = 38),  
          axis.text.x  = element_text(size = 40),             
          axis.text.y  = element_text(size = 40))
  
F14R_dna_load_points
  
ggsave(filename = "D:/Decicomp/R/Figures/F14R_dna_load_points2.tiff", plot = F14R_dna_load_points,          
       device = "tiff", width = 8, height = 7,  units = "in", dpi = 600)                 
  
  
  
H2D_dna_virus_dots <- Load_viral_decicomp_dots %>%
  dplyr::filter(Family == "H2D") %>%
  dplyr::mutate(Time = as.numeric(as.character(Time)),
                Mean = as.numeric(Copies_2))
  
set.seed(123)
H2D_dna_load_points <- ggplot(H2D_dna_virus_dots,aes(x = Time, y = Mean, group = Age)) +
  # puntos con relleno y borde negro
  geom_point(aes(fill = factor(Age)), shape = 21, size = 8, color = "black", alpha = 0.8) +  
  # línea promedio por grupo
  stat_summary(aes(color = factor(Age)), fun = mean, geom = "line", size = 1.5) +  
  scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0,25), breaks = c(0,3,6,12,24), expand = c(0,0)) +
  scale_color_manual(values = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  scale_fill_manual(values  = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  labs(x = "Time (hpc)",  y = expression("Genomic units/ng")) +
  Style_format_theme   +
  theme(axis.title.x = element_text(size = 44, color = "white"),
        axis.title.y = element_text(size = 38),  
        axis.text.x  = element_text(size = 40),             
        axis.text.y  = element_text(size = 40))

H2D_dna_load_points

ggsave(filename = "D:/Decicomp/R/Figures/H2D_dna_load_points.tiff", plot = H2D_dna_load_points,          
       device = "tiff", width = 8, height = 7,  units = "in", dpi = 600)                 


### RNA viral load ----
  

# Step 1:
family_counts <- read.delim("D:/Decicomp/virus_DNA/family_counts.txt")
ORF_data      <- read.delim("D:/Decicomp/virus_DNA/ORF_bp.txt", header = F) %>% rename("ORF" = V1, "bp"=V2) 
reads_family  <- read.delim("D:/Decicomp/virus_DNA/total_reads_sample.txt", header = F) %>% rename("sample" = V1, "reads"=V2) 


family_counts_rna <- family_counts %>% 
                     dplyr::select(-ORF) %>%
                     summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
                     t() %>%
                     as.data.frame() %>%
                     rownames_to_column("sample") %>% 
                     rename(counts = V1) %>%
                     mutate(counts = if_else(grepl("T0", sample), 0, counts)) %>%
  
                     mutate(original_sample = sample) %>%
                     separate(original_sample, into = c("Family", "Age", "Time", "replicate"), sep = "_") %>%
                     mutate(counts = log10(counts + 1)) %>%
                     mutate(Time = gsub("T", "", Time)) %>%
                     mutate(Time = as.numeric(Time)) %>%
                     mutate(Age = case_when(
                              Age == "C3" ~ 4,
                              Age == "C7" ~ 16,
                              Age == "C8" ~ 28, TRUE ~ as.numeric(Age))) %>%
                     group_by(Family, Age, Time) %>%
                     summarise(
                               Mean = mean(counts, na.rm = TRUE),
                               SD = sd(counts, na.rm = TRUE),
                               Num_samples = n(),  
                               SEM = SD / sqrt(Num_samples),
                               .groups = 'drop'  )


F14R_rna_virus      <- family_counts_rna %>% filter(Family=="F14R") %>% dplyr::select(1:4,7)

F14R_rna_load_graph <- ggplot(F14R_rna_virus, aes(x = Time, y = Mean, fill=factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Age), size=1, linetype= "dashed" ) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar RNA load")) +
  scale_color_manual(values = c("4"="#dae3f3", "16"="#8faadc", "28"="#2f5597")) +
  scale_fill_manual(values =  c("4"="#dae3f3",  "16"="#8faadc", "28"="#2f5597")) +
  Style_format_theme
F14R_rna_load_graph

# F14R_rna_load_graph
# tiff("F14R_rna_load_graph.tiff", units="in", width=10, height=8, res=400)
# F14R_rna_load_graph
# dev.off()

H2D_rna_virus      <- family_counts_rna %>% filter(Family=="H2D") %>% select(1:4,7)

H2D_rna_load_graph <- ggplot(H2D_rna_virus, aes(x = Time, y = Mean,  fill = factor(Age), color = factor(Age))) +
  geom_point(shape = 21, size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.1, size=1) +
  geom_line(aes(group = Age), size=1, linetype= "dashed" ) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
  # labels = trans_format("log10", math_format(10^.x)),limits=c(1, 10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 25), breaks = c(0, 3, 6, 12, 24), expand = c(0, 0)) +
  labs(x = "Time (hours)",
       y = expression("OsHV-1 µvar RNA load")) +
  scale_color_manual(values = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  scale_fill_manual(values  = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  Style_format_theme
H2D_rna_load_graph

# H2D_rna_load_graph
# tiff("H2D_rna_load_graph.tiff", units="in", width=10, height=8, res=400)
# H2D_rna_load_graph
# dev.off()



# Plots RNA based on dots Journal  ----
# Step 1:
family_counts <- read.delim("D:/Decicomp/virus_DNA/family_counts.txt")
ORF_data      <- read.delim("D:/Decicomp/virus_DNA/ORF_bp.txt", header = F) %>% rename("ORF" = V1, "bp"=V2) 
reads_family  <- read.delim("D:/Decicomp/virus_DNA/total_reads_sample.txt", header = F) %>% rename("sample" = V1, "reads"=V2) 

family_counts_rna <- family_counts %>% 
  dplyr::select(-ORF) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>% 
  rename(counts = V1) %>% mutate(counts = if_else(grepl("T0", sample), 0, counts)) %>%
  
  mutate(original_sample = sample) %>%
  separate(original_sample, into = c("Family", "Age", "Time", "replicate"), sep = "_") %>%
  mutate(counts = log10(counts + 1)) %>%
  mutate(Time = gsub("T", "", Time)) %>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(Age = case_when(
    Age == "C3" ~ 4,
    Age == "C7" ~ 16,
    Age == "C8" ~ 28, TRUE ~ as.numeric(Age))) %>%
  group_by(Family, Age, Time) 


F14R_rna_virus      <- family_counts_rna %>% filter(Family=="F14R")

F14R_rna_load_graph <- ggplot(F14R_rna_virus, aes(x = Time, y = counts, fill=factor(Age), color = factor(Age))) +
  geom_point(aes(fill = factor(Age)), shape = 21, size = 8, color = "black", alpha = 0.8) +  
  # línea promedio por grupo
  stat_summary(aes(color = factor(Age)), fun = mean, geom = "line", size = 1.5, linetype= "dashed") +  
  scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0,25), breaks = c(0,3,6,12,24), expand = c(0,0)) +
  scale_fill_manual(values = c("4" = "#dae3f3", "16" = "#8faadc", "28" = "#2f5597")) +
  scale_color_manual(values = c("4" = "#dae3f3", "16" = "#8faadc", "28" = "#2f5597")) +
  labs(x = "Time (hpc)",  y = expression("counts")) +
  Style_format_theme   +
  theme(axis.title.x = element_text(size = 44),  
        axis.title.y = element_text(size = 38),  
        axis.text.x  = element_text(size = 40),             
        axis.text.y  = element_text(size = 40))

F14R_rna_load_graph

ggsave(filename = "D:/Decicomp/R/Figures/F14R_rna_load_graph2.tiff", plot = F14R_rna_load_graph,          
       device = "tiff", width = 8, height = 7,  units = "in", dpi = 600)                 






H2D_rna_virus      <- family_counts_rna %>% filter(Family=="H2D")

H2D_rna_load_graph <- ggplot(H2D_rna_virus, aes(x = Time, y = counts, fill=factor(Age), color = factor(Age))) +
  geom_point(aes(fill = factor(Age)), shape = 21, size = 8, color = "black", alpha = 0.8) +  
  # línea promedio por grupo
  stat_summary(aes(color = factor(Age)), fun = mean, geom = "line", size = 1.5, linetype= "dashed") +  
  scale_y_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0,25), breaks = c(0,3,6,12,24), expand = c(0,0)) +
  scale_color_manual(values = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  scale_fill_manual(values  = c("4"="#E2F0D9", "16"="#A9D18E", "28"="#385700")) +
  labs(x = "Time (hpc)",  y = expression("counts")) +
  Style_format_theme   +
  theme(axis.title.x = element_text(size = 44, color = "white"),  
        axis.title.y = element_text(size = 38),  
        axis.text.x  = element_text(size = 40),             
        axis.text.y  = element_text(size = 40))

H2D_rna_load_graph

ggsave(filename = "D:/Decicomp/R/Figures/H2D_rna_load_graph2.tiff", plot = H2D_rna_load_graph,          
       device = "tiff", width = 8, height = 7,  units = "in", dpi = 600)                 








F14R_dna_load_graph
F14R_rna_load_graph

H2D_dna_load_graph
H2D_rna_load_graph


combined_plot <- (F14R_dna_load_graph |  H2D_dna_load_graph) /
                 (F14R_rna_load_graph  | H2D_rna_load_graph) + 
                 plot_layout(guides = "auto")

combined_plot



#############

ANOVAS

Load_viral_decicomp_anova <- Load_viral_decicomp  %>%
  mutate(Family = as.factor(Family),
         Age = as.factor(Age),
         Time = as.factor(Time),
         Copies_2 = log10(Copies_2 + 1)) %>%
  dplyr::select(Family, Age, Time, Copies_2) 

H2D_at_24H <- Load_viral_decicomp_anova %>% filter(Family=="H2D" & Time ==24)

H2D_at_24H_stat <-   H2D_at_24H %>%   group_by(Age) %>% summarise(p_value = shapiro.test(Copies_2)$p.value)
print(H2D_at_24H_stat)

leveneTest(Copies_2 ~ Age, data = H2D_at_24H)
anova_result_H2D <- aov(Copies_2 ~ Age, data = H2D_at_24H)
summary(anova_result_H2D)
tukey_result_H2D <- TukeyHSD(anova_result_H2D)
print(tukey_result_H2D)
plot(tukey_result_H2D, las = 1, col = "blue")


F14R_at_24H <- Load_viral_decicomp_anova %>% filter(Family=="F14R" & Time ==24)

F14R_at_24H_stat <-   F14R_at_24H %>%   group_by(Age) %>% summarise(p_value = shapiro.test(Copies_2)$p.value)
print(F14R_at_24H_stat)

leveneTest(Copies_2 ~ Age, data = F14R_at_24H)
anova_result_F14R <- aov(Copies_2 ~ Age, data = F14R_at_24H)
summary(anova_result_F14R)
tukey_result_F14R <- TukeyHSD(anova_result_F14R)
print(tukey_result_F14R)
plot(tukey_result_F14R, las = 1, col = "blue")


family_counts_rna <- family_counts %>% 
  select(-ORF) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>% 
  rename(counts = V1) %>%
  mutate(counts = if_else(grepl("T0", sample), 0, counts)) %>%
  
  mutate(original_sample = sample) %>%
  separate(original_sample, into = c("Family", "Age", "Time", "replicate"), sep = "_") %>%
  mutate(counts = log10(counts + 1)) %>%
  mutate(Time = gsub("T", "", Time)) %>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(Age = case_when(
    Age == "C3" ~ 4,
    Age == "C7" ~ 16,
    Age == "C8" ~ 28, TRUE ~ as.numeric(Age))) 


H2D_at_24H_rna <- family_counts_rna %>% filter(Family=="H2D" & Time ==24) %>% 
                  #slice(-13) %>% 
                  mutate(log_count=log(counts+1))
str(H2D_at_24H_rna)

H2D_at_24H_stat <-   H2D_at_24H_rna %>%   group_by(Age) %>% summarise(p_value = shapiro.test(log_count)$p.value)
print(H2D_at_24H_stat)

leveneTest(log_count ~ as.factor(Age), data = H2D_at_24H_rna)

anova_result_H2D <- aov(log_count ~ as.factor(Age), data = H2D_at_24H_rna)
summary(anova_result_H2D)
tukey_result_H2D <- TukeyHSD(anova_result_H2D)
print(tukey_result_H2D)
plot(tukey_result_H2D, las = 1, col = "blue")



F14R_at_24H_rna <- family_counts_rna %>% filter(Family=="F14R" & Time ==24) %>% 
                    mutate(log_count=log(counts+1))
str(F14R_at_24H_rna)

F14R_at_24H_stat <-   F14R_at_24H_rna %>%   group_by(Age) %>% summarise(p_value = shapiro.test(log_count)$p.value)
print(F14R_at_24H_stat)

leveneTest(log_count ~ as.factor(Age), data = F14R_at_24H_rna)

anova_result_F14R <- aov(log_count ~ as.factor(Age), data = F14R_at_24H_rna)
summary(anova_result_F14R)
tukey_result_F14R <- TukeyHSD(anova_result_F14R)
print(tukey_result_F14R)
plot(tukey_result_F14R, las = 1, col = "blue")

