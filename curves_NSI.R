library(readxl)
library(dplyr)
library(tidyr) 
library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(car)
library(scales)
library(plotly)

setwd("D:/Decicomp/R/Figures")

# Theme_format          ---- 
####
Style_format_theme <- theme(
  
  ## Title of the axis
  axis.title.x       = element_text(color="black", size=26,  margin = margin(t = 12, r = 0,  b = 0, l = 0)),
  axis.title.y       = element_text(color="black", size=26,  margin = margin(t = 0,  r = 12, b = 0, l = 0)),
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
  legend.title       = element_text(size=14),
  legend.text        = element_text(size=12),
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
###

NSI       <- read_excel("D:/Decicomp/Paper/Paper Valdi/Draft Age Decicomp V3/DECICOMP mortalité témoins inter EMP.xlsx",    sheet = "NSI")
NSI_long  <- NSI %>% dplyr::select(1, 2,5, 10, 13, 14) 

# Experiment  1 ----

###  Donnors  ----
NSI_long_donors_exp1    <- NSI_long  %>% filter(Experiment==1 & oysters=="donnors") %>%  dplyr::select(3,5,6)
deaths_donors_exp1      <- NSI_long_donors_exp1 %>% filter(!is.na(Dead), Dead > 0) %>%  group_by(time_hour) %>% summarise(total_dead = sum(Dead), .groups = "drop") %>%
                           mutate(num = lapply(total_dead, seq_len)) %>% unnest(num) %>% mutate(cencored = 1) %>% dplyr::select(time_hour, cencored)
max_time_donors_exp1    <- max(NSI_long_donors_exp1$time_hour, na.rm = TRUE)
final_alive_donors_exp1 <- NSI_long_donors_exp1 %>%  filter(time_hour== max_time_donors_exp1) %>% summarise(total_alive = sum(Alive, na.rm = TRUE)) %>% pull(total_alive)
cencored_donors_exp1    <- tibble(time_hour = rep(max(NSI_long_donors_exp1$time_hour), final_alive_donors_exp1), cencored = 0)
long_NSI_donors_exp1    <- bind_rows(deaths_donors_exp1, cencored_donors_exp1) %>%  mutate(numero_oyster = row_number()) %>% 
                           dplyr::select(numero_oyster, time_hour, cencored) %>% mutate(Experiment = 1)
total_muertos_exp1          <- NSI_long_donors_exp1 %>%  filter(!is.na(Dead), Dead > 0) %>% summarise(total_muertos = sum(Dead), .groups = "drop") %>% pull(total_muertos)
long_NSI_curve_donors_exp1  <- Surv(time = long_NSI_donors_exp1$time_hour, event = long_NSI_donors_exp1$cencored)
fit_long_NSI_donors_exp1    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_donors_exp1)

survival_donnors_experiment_1 <- ggsurvplot(fit_long_NSI_donors_exp1, 
                             data = long_NSI_donors_exp1, 
                             #pval = F,
                             #pval.size = 6,
                             surv.median.line = "hv",
                             break.time.by=48, 
                             xlab = "Time (hpc)",
                             #break.y.by = 0.20, 
                             #break.x.by = 50, 
                             palette= c( "red"),
                             conf.int = FALSE,
                             censor=F,
                             #linetype=c(5,1,6),
                             size=2,
                             legend = "none",
                             axes.offset = T,
                             ggtheme = Style_format_theme)

survival_donnors_experiment_1$plot <- survival_donnors_experiment_1$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))               
  
survival_donnors_experiment_1

###  REC  ----
NSI_long_rec_exp1    <- NSI_long  %>% filter(Experiment==1 & oysters=="rec") %>%  dplyr::select(3,5,6)
deaths_rec_exp1      <- NSI_long_rec_exp1 %>% filter(!is.na(Dead), Dead > 0) %>%  group_by(time_hour) %>% summarise(total_dead = sum(Dead), .groups = "drop") %>%
                        mutate(num = lapply(total_dead, seq_len)) %>% unnest(num) %>% mutate(cencored = 1) %>% dplyr::select(time_hour, cencored)
max_time_rec_exp1    <- max(NSI_long_rec_exp1$time_hour, na.rm = TRUE)
final_alive_rec_exp1 <- NSI_long_rec_exp1 %>%  filter(time_hour== max_time_rec_exp1) %>% summarise(total_alive = sum(Alive, na.rm = TRUE)) %>% pull(total_alive)
cencored_rec_exp1    <- tibble(time_hour = rep(max(NSI_long_rec_exp1$time_hour), final_alive_rec_exp1), cencored = 0)
long_NSI_rec_exp1    <- bind_rows(deaths_rec_exp1, cencored_rec_exp1) %>%  mutate(numero_oyster = row_number()) %>% 
                        dplyr::select(numero_oyster, time_hour, cencored) %>% mutate(Experiment = 1)
total_muertos_exp1       <- NSI_long_rec_exp1 %>%  filter(!is.na(Dead), Dead > 0) %>% summarise(total_muertos = sum(Dead), .groups = "drop") %>% pull(total_muertos)
long_NSI_curve_rec_exp1  <- Surv(time = long_NSI_rec_exp1$time_hour, event = long_NSI_rec_exp1$cencored)
fit_long_NSI_rec_exp1    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_rec_exp1)

survival_rec_experiment_1 <- ggsurvplot(fit_long_NSI_rec_exp1, 
                                            data = long_NSI_rec_exp1, 
                                            #pval = F,
                                            #pval.size = 6,
                                            surv.median.line = "hv",
                                            break.time.by=48, 
                                            xlab = "Time (hpc)",
                                            #break.y.by = 0.20, 
                                            #break.x.by = 50, 
                                            palette= c( "green"),
                                            conf.int = FALSE,
                                            censor=F,
                                            #linetype=c(5,1,6),
                                            size=2,
                                            legend = "none",
                                            axes.offset = T,
                                            ggtheme = Style_format_theme)

survival_rec_experiment_1$plot <- survival_rec_experiment_1$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))               

survival_rec_experiment_1


# Experiment  2 ----

###  Donnors  ----
NSI_long_donors_exp2    <- NSI_long  %>% filter(Experiment==2 & oysters=="donnors") %>%  dplyr::select(3,5,6)
deaths_donors_exp2      <- NSI_long_donors_exp2 %>% filter(!is.na(Dead), Dead > 0) %>%  group_by(time_hour) %>% summarise(total_dead = sum(Dead), .groups = "drop") %>%
                           mutate(num = lapply(total_dead, seq_len)) %>% unnest(num) %>% mutate(cencored = 1) %>% dplyr::select(time_hour, cencored)
max_time_donors_exp2    <- max(NSI_long_donors_exp2$time_hour, na.rm = TRUE)

final_alive_donors_exp2 <- NSI_long_donors_exp2 %>%  filter(time_hour   == max_time_donors_exp2) %>% summarise(total_alive = sum(Alive, na.rm = TRUE)) %>% pull(total_alive)
cencored_donors_exp2    <- tibble(time_hour = rep(max(NSI_long_donors_exp2$time_hour), final_alive_donors_exp2), cencored = 0)
long_NSI_donors_exp2    <- bind_rows(deaths_donors_exp2, cencored_donors_exp2) %>%  mutate(numero_oyster = row_number()) %>% 
                           dplyr::select(numero_oyster, time_hour, cencored) %>% mutate(Experiment = 2)
total_muertos_exp2          <- NSI_long_donors_exp2 %>%  filter(!is.na(Dead), Dead > 0) %>% summarise(total_muertos = sum(Dead), .groups = "drop") %>% pull(total_muertos)
long_NSI_curve_donors_exp2  <- Surv(time = long_NSI_donors_exp2$time_hour, event = long_NSI_donors_exp2$cencored)
fit_long_NSI_donors_exp2    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_donors_exp2)

survival_donnors_experiment_2 <- ggsurvplot(fit_long_NSI_donors_exp2, 
                                            data = long_NSI_donors_exp2, 
                                            #pval = F,
                                            #pval.size = 6,
                                            surv.median.line = "hv",
                                            break.time.by=48, 
                                            xlab = "Time (hpc)",
                                            #break.y.by = 0.20, 
                                            #break.x.by = 50, 
                                            palette= c( "darkred"),
                                            conf.int = FALSE,
                                            censor=F,
                                            #linetype=c(5,1,6),
                                            size=2,
                                            legend = "none",
                                            axes.offset = T,
                                            ggtheme = Style_format_theme)

survival_donnors_experiment_2$plot <- survival_donnors_experiment_2$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))               

survival_donnors_experiment_2

###  REC  ----
NSI_long_rec_exp2  <- NSI_long  %>% filter(Experiment==2 & oysters=="rec") %>%  dplyr::select(3,5,6)
deaths_rec_exp2    <- NSI_long_rec_exp2 %>% filter(!is.na(Dead), Dead > 0) %>%  group_by(time_hour) %>% summarise(total_dead = sum(Dead), .groups = "drop") %>%
                      mutate(num = lapply(total_dead, seq_len)) %>% unnest(num) %>% mutate(cencored = 1) %>% dplyr::select(time_hour, cencored)
max_time_rec_exp2    <- max(NSI_long_rec_exp2$time_hour, na.rm = TRUE)
final_alive_rec_exp2 <- NSI_long_rec_exp2 %>%  filter(time_hour   == max_time_rec_exp2) %>% summarise(total_alive = sum(Alive, na.rm = TRUE)) %>% pull(total_alive)
cencored_rec_exp2    <- tibble(time_hour = rep(max(NSI_long_rec_exp2$time_hour), final_alive_rec_exp2), cencored = 0)
long_NSI_rec_exp2    <- bind_rows(deaths_rec_exp2, cencored_rec_exp2) %>%  mutate(numero_oyster = row_number()) %>% 
                        dplyr::select(numero_oyster, time_hour, cencored) %>% mutate(Experiment = 2)
total_muertos_exp2       <- NSI_long_rec_exp2 %>%  filter(!is.na(Dead), Dead > 0) %>% summarise(total_muertos = sum(Dead), .groups = "drop") %>% pull(total_muertos)
long_NSI_curve_rec_exp2  <- Surv(time = long_NSI_rec_exp2$time_hour, event = long_NSI_rec_exp2$cencored)
fit_long_NSI_rec_exp2    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_rec_exp2)

survival_rec_experiment_2 <- ggsurvplot(fit_long_NSI_rec_exp2, 
                                        data = long_NSI_rec_exp2, 
                                        #pval = F,
                                        #pval.size = 6,
                                        surv.median.line = "hv",
                                        break.time.by=48, 
                                        xlab = "Time (hpc)",
                                        #break.y.by = 0.20, 
                                        #break.x.by = 50, 
                                        palette= c( "darkgreen"),
                                        conf.int = FALSE,
                                        censor=F,
                                        #linetype=c(5,1,6),
                                        size=2,
                                        legend = "none",
                                        axes.offset = T,
                                        ggtheme = Style_format_theme)

survival_rec_experiment_2$plot <- survival_rec_experiment_2$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))               

survival_rec_experiment_2

# Experiment  3 ----

###  Donnors  ----
NSI_long_donors_exp3    <- NSI_long  %>% filter(Experiment==3 & oysters=="donnors") %>%  dplyr::select(3,5,6)
deaths_donors_exp3      <- NSI_long_donors_exp3 %>% filter(!is.na(Dead), Dead > 0) %>%  group_by(time_hour) %>% summarise(total_dead = sum(Dead), .groups = "drop") %>%
                            mutate(num = lapply(total_dead, seq_len)) %>% unnest(num) %>% mutate(cencored = 1) %>% dplyr::select(time_hour, cencored)
max_time_donors_exp3    <- max(NSI_long_donors_exp3$time_hour, na.rm = TRUE)
final_alive_donors_exp3 <- NSI_long_donors_exp3 %>%  filter(time_hour== max_time_donors_exp3) %>% summarise(total_alive = sum(Alive, na.rm = TRUE)) %>% pull(total_alive)
cencored_donors_exp3    <- tibble(time_hour = rep(max(NSI_long_donors_exp3$time_hour), final_alive_donors_exp3), cencored = 0)
long_NSI_donors_exp3    <- bind_rows(deaths_donors_exp3, cencored_donors_exp3) %>%  mutate(numero_oyster = row_number()) %>% 
                           dplyr::select(numero_oyster, time_hour, cencored) %>% mutate(Experiment = 3)
total_muertos_exp3          <- NSI_long_donors_exp3 %>%  filter(!is.na(Dead), Dead > 0) %>% summarise(total_muertos = sum(Dead), .groups = "drop") %>% pull(total_muertos)
long_NSI_curve_donors_exp3  <- Surv(time = long_NSI_donors_exp3$time_hour, event = long_NSI_donors_exp3$cencored)
fit_long_NSI_donors_exp3    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_donors_exp3)

survival_donnors_experiment_3 <- ggsurvplot(fit_long_NSI_donors_exp3, 
                                            data = long_NSI_donors_exp3, 
                                            #pval = F,
                                            #pval.size = 6,
                                            surv.median.line = "hv",
                                            break.time.by=48, 
                                            xlab = "Time (hpc)",
                                            #break.y.by = 0.20, 
                                            #break.x.by = 50, 
                                            palette= c( "salmon"),
                                            conf.int = FALSE,
                                            censor=F,
                                            #linetype=c(5,1,6),
                                            size=2,
                                            legend = "none",
                                            axes.offset = T,
                                            ggtheme = Style_format_theme)

survival_donnors_experiment_3$plot <- survival_donnors_experiment_3$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))               

survival_donnors_experiment_3

###  REC  ----
NSI_long_rec_exp3  <- NSI_long  %>% filter(Experiment==3 & oysters=="rec") %>%  dplyr::select(3,5,6)
deaths_rec_exp3    <- NSI_long_rec_exp3 %>% filter(!is.na(Dead), Dead > 0) %>%  group_by(time_hour) %>% summarise(total_dead = sum(Dead), .groups = "drop") %>%
                      mutate(num = lapply(total_dead, seq_len)) %>% unnest(num) %>% mutate(cencored = 1) %>% dplyr::select(time_hour, cencored)
max_time_rec_exp3    <- max(NSI_long_rec_exp3$time_hour, na.rm = TRUE)
final_alive_rec_exp3 <- NSI_long_rec_exp3 %>%  filter(time_hour   == max_time_rec_exp3) %>% summarise(total_alive = sum(Alive, na.rm = TRUE)) %>% pull(total_alive)
cencored_rec_exp3    <- tibble(time_hour = rep(max(NSI_long_rec_exp3$time_hour), final_alive_rec_exp3), cencored = 0)
long_NSI_rec_exp3    <- bind_rows(deaths_rec_exp3, cencored_rec_exp3) %>%  mutate(numero_oyster = row_number()) %>% 
                        dplyr::select(numero_oyster, time_hour, cencored) %>% mutate(Experiment = 3)
total_muertos_exp3       <- NSI_long_rec_exp3 %>%  filter(!is.na(Dead), Dead > 0) %>% summarise(total_muertos = sum(Dead), .groups = "drop") %>% pull(total_muertos)
long_NSI_curve_rec_exp3  <- Surv(time = long_NSI_rec_exp3$time_hour, event = long_NSI_rec_exp3$cencored)
fit_long_NSI_rec_exp3    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_rec_exp3)

survival_rec_experiment_3 <- ggsurvplot(fit_long_NSI_rec_exp3, 
                                        data = long_NSI_rec_exp3, 
                                        #pval = F,
                                        #pval.size = 6,
                                        surv.median.line = "hv",
                                        break.time.by=48, 
                                        xlab = "Time (hpc)",
                                        #break.y.by = 0.20, 
                                        #break.x.by = 50, 
                                        palette= c( "darkgreen"),
                                        conf.int = FALSE,
                                        censor=F,
                                        #linetype=c(5,1,6),
                                        size=2,
                                        legend = "none",
                                        axes.offset = T,
                                        ggtheme = Style_format_theme)

survival_rec_experiment_3$plot <- survival_rec_experiment_3$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))               

survival_rec_experiment_3



# NSI donors all experiments ----
long_NSI_donors_exp1
long_NSI_donors_exp2
long_NSI_donors_exp3

long_NSI_donors_all        <- rbind(long_NSI_donors_exp1, long_NSI_donors_exp2, long_NSI_donors_exp3)
long_NSI_curve_donors_all  <- Surv(time = long_NSI_donors_all$time_hour, event = long_NSI_donors_all$cencored)
fit_long_NSI_donors_all    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_donors_all)

survival_donnors_all <- ggsurvplot(fit_long_NSI_donors_all, 
                                            data = long_NSI_donors_all, 
                                            #pval = F,
                                            #pval.size = 6,
                                            #surv.median.line = "hv",
                                            break.time.by=48, 
                                            xlab = "Time (hpc)",
                                            #break.y.by = 0.20, 
                                            #break.x.by = 50, 
                                            palette= c("darkred", "red", "orange"),
                                            conf.int = FALSE,
                                            censor=F,
                                            #linetype=c(5,1,6),
                                            size=2,
                                            legend = "none",
                                            axes.offset = T,
                                            ggtheme = Style_format_theme)

survival_donnors_all$plot <- survival_donnors_all$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))          

survival_donnors_all

tiff("survival_donnors_all.tiff", units="in", width=9, height=7, res=400)
survival_donnors_all
dev.off()



# NSI rec all experiments ----

long_NSI_rec_exp1
long_NSI_rec_exp2
long_NSI_rec_exp3

long_NSI_rec_all        <- rbind(long_NSI_rec_exp1, long_NSI_rec_exp2, long_NSI_rec_exp3)
long_NSI_curve_rec_all  <- Surv(time = long_NSI_rec_all$time_hour, event = long_NSI_rec_all$cencored)
fit_long_NSI_rec_all    <- survfit(Surv(time_hour, cencored) ~ Experiment, data = long_NSI_rec_all)


survival_rec_all <- ggsurvplot(fit_long_NSI_rec_all, 
                                   data = long_NSI_rec_all, 
                                   #pval = F,
                                   #pval.size = 6,
                                   #surv.median.line = "hv",
                                   break.time.by=48, 
                                   xlab = "Time (hpc)",
                                   #break.y.by = 0.20, 
                                   #break.x.by = 50, 
                                   palette= c("darkgreen", "green", "turquoise"),
                                   conf.int = FALSE,
                                   censor=F,
                                   #linetype=c(5,1,6),
                                   size=2,
                                   legend = "none",
                                   axes.offset = T,
                                   ggtheme = Style_format_theme)

survival_rec_all$plot <- survival_rec_all$plot      +
  scale_y_continuous(breaks = seq(0,1,   by=0.25), 
                     labels = seq(0,100, by=25), 
                     expand=c(0, 0, .05, 0))        +
  scale_x_continuous( expand=c(0, 0))          

survival_rec_all

tiff("survival_rec_all.tiff", units="in", width=9, height=7, res=400)
survival_rec_all
dev.off()


