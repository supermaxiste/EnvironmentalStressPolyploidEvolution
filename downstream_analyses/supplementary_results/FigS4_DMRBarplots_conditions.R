### Script to plot Barplot of DMRs between HMvsLL 
### G1 and G4 synthetics against each other

### Import libraries

library(tidyverse)
library(data.table)
library(patchwork)

## Import data

setwd("MILD_v_STRESS_syn1/dmrseq/")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides 

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "r", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "r", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "r", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 9]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 8]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 9]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 8]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 9]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 8]


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

### No significant regions detected!

## set variables to zero

CG_hal <- 0
CHG_hal <- 0
CHH_hal <- 0

CG_lyr <- 0
CHG_lyr <- 0
CHH_lyr <- 0

## First generation

context <- rep(c("CG", "CHG", "CHH"), 2)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(CG_hal, CHG_hal, CHH_hal,
           CG_hal, CHG_hal, CHH_hal,
           CG_lyr, CHG_lyr, CHH_lyr,
           CG_lyr, CHG_lyr, CHH_lyr)

comparison <- c(rep("synthetic_1_HMvLL", 12))

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

dmrs_context_total$hypo_hyper <- factor(dmrs_context_total$hypo_hyper, 
                                        levels = unique(dmrs_context_total$hypo_hyper))

ggplot(dmrs_context_total, aes(x=context, y=value)) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack", 
           aes(col = hypo_hyper)) +
  #scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", breaks = , 
                   labels = c("CG", "CHG", "CHH", "CG", "CHG", "CHH")) +
  ylim(0, 2200) +
  ggtitle("Synthetics G1 - HM vs LL")


## We repeat the same procedure with G4 synthetics

## Import data

setwd("MILD_v_STRESS_syn4/dmrseq/")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides 

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "r", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "r", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "r", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 9]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 8]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 9]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 8]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 9]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 8]


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


## Fourth generation

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_4_HMvLL", 12))

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

dmrs_context_total$hypo_hyper <- factor(dmrs_context_total$hypo_hyper, 
                                        levels = unique(dmrs_context_total$hypo_hyper))

ggplot(dmrs_context_total, aes(x=context, y=value)) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack", 
           aes(col = hypo_hyper)) +
  #scale_fill_manual(values=cbPalette) +
  facet_wrap(~species, labeller = labeller(species = species.labs)) +
  scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", breaks = , 
                   labels = c("CG", "CHG", "CHH", "CG", "CHG", "CHH")) +
  ylim(0, 2200) +
  ggtitle("Synthetics G4 - HM vs LL")

