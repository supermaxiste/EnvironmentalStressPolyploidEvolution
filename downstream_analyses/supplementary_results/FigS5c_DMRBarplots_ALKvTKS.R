### Script to plot Barplot of DMRs between ALK G1 and TKS G1

### Import libraries

library(ggplot2)
library(data.table)
library(patchwork)

## Figure X

### Import first set of files, ALK G1 vs TKS G1 (Mild Conditions)
setwd("~/data/FigS3c/MILD_alk1Vtks1/dmrseq")
#setwd("MILD_syn1_v_tks1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/A_v_B_polyploid.txt")
dmrseq_output_CG_lyr <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/A_v_B_polyploid.txt")

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

comparison <- c(rep("Mild", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

hm_conditions <- ggplot(dmrs_context_total, 
                        aes(x=factor(context, 
                                     levels = c("CG",
                                                "CHG",
                                                "CHH")), 
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
        strip.background = element_rect(color="white", fill="white", linewidth=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  #scale_x_discrete(name = "DMR context", 
                   #labels = c("", "CG", "", "CHG", "", "CHH")) +
  ylim(0, 70000) +
  xlab("Context") +
  #geom_vline(xintercept = c(1.5, 2.5), linetype = "longdash") +
  ggtitle("Mild Conditions")

## Repeat the procedure for above for hot conditions

### Import first set of files, progenitors vs synthetic G1 (HM conditions)
setwd("~/data/FigS3c/STRESS_alk1Vtks1/dmrseq")
#setwd("STRESS_syn1_v_tks1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/A_v_B_polyploid.txt")
dmrseq_output_CG_lyr <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/A_v_B_polyploid.txt")

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

comparison <- c(rep("Stress", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)


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
                        aes(x=factor(context, 
                                     levels = c("CG",
                                                "CHG",
                                                "CHH")), 
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
        strip.background = element_rect(color="white", fill="white", linewidth=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  #scale_x_discrete(name = "DMR context", 
  #labels = c("", "CG", "", "CHG", "", "CHH")) +
  ylim(0, 70000) +
  xlab("Context") +
  #geom_vline(xintercept = c(1.5, 2.5), linetype = "longdash") +
  ggtitle("Stress Conditions")

hm_conditions + ll_conditions
