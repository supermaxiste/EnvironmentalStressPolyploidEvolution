### Script to plot Barplot of DMRs over generations

### Import libraries

library(tidyverse)
library(data.table)
library(patchwork)

## Figure 3

### Import first set of files, progenitors vs synthetic G1 (HM conditions)

setwd("FigS3a/MILD_syn1_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## Filter only significant regions 

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("G1", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import data for progenitors vs synthetic G4 (HM conditions)

setwd("FigS3a/MILD_syn4_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")


## Filter only significant regions 

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("G4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)
## Barplot code

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

hm_conditions <- ggplot(dmrs_context_total, 
                        aes(x=factor(interaction(context, comparison), 
                                     levels = c("CG.G1",
                                                "CG.G4",
                                                "CHG.G1",
                                                "CHG.G4",
                                                "CHH.G1",
                                                "CHH.G4")), 
                            y=value, 
                            fill = comparison, 
                            col = factor(hypo_hyper, levels = c("hyper", "hypo")))) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hyper", "hypo"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=25),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", 
                   labels = c("", "CG", "", "CHG", "", "CHH")) +
  ylim(0, 12000) +
  geom_vline(xintercept = c(2.5, 4.5), linetype = "longdash") +
  ggtitle("Mild Conditions")

## Repeat the procedure for above for hot conditions

### Import first set of files, progenitors vs synthetic G1 (HM conditions)

setwd("FigS3a/STRESS_syn1_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## Filter only significant regions 

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("G1", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)


## Import data for progenitors vs synthetic G4 (HM conditions)

setwd("FigS3a/STRESS_syn4_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")


## Filter only significant regions 

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("G4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)
## Barplot code

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ll_conditions <- ggplot(dmrs_context_total, 
                        aes(x=factor(interaction(context, comparison), 
                                     levels = c("CG.G1",
                                                "CG.G4",
                                                "CHG.G1",
                                                "CHG.G4",
                                                "CHH.G1",
                                                "CHH.G4")), 
                            y=value, 
                            fill = comparison, 
                            col = factor(hypo_hyper, levels = c("hyper", "hypo")))) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hyper", "hypo"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", 
                   labels = c("", "CG", "", "CHG", "", "CHH")) +
  ylim(0, 12000) +
  geom_vline(xintercept = c(2.5, 4.5), linetype = "longdash") +
  ggtitle("Stress Conditions")

hm_conditions + ll_conditions


            