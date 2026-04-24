library(readxl)
library(ggplot2)
library(dplyr)

##### THEME   ----          
Style_format_theme <- theme(
  ## Title of the axis
  axis.title.x       = element_text(color="black", size=24),
  axis.title.y       = element_text(color="black", size=24),
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

Biomass_families <- read_excel("D:/Decicomp/R/Biomas/Biomass_families.xlsx") %>% mutate(Family = factor(Family, levels = c("H2D", "F14R", "F11N", "F14V")))

Biomas_plot <- ggplot(Biomass_families, aes(x = factor(Age), y = mean, fill = Family)) +
  geom_bar(stat = "identity", color= "black", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - 0, ymax = mean + SD),
                position = position_dodge(width = 0.9), width = 0.3) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) +
  scale_fill_manual(values = c("F14R" ="#2f5597","H2D" ="#385700", "F14V"="#ffc000", "F11N"="#7030A0" )) +
  labs(x = "Age (Months)", y = "Biomass (mean ± SD)") +
  Style_format_theme +
  theme(text = element_text(size = 12))

Biomas_plot 


ggsave("D:/Decicomp/R/MOFA_omics/Data_compacted/biomass.tiff", Biomas_plot, width = 8, height = 6, dpi = 600)


