### Script to plot Barplot of DMRs over generations
### for natural and diploid progenitor species

## Figure S5

### Import libraries

library(tidyverse)
library(data.table)
library(patchwork)
library(GenomicRanges)

### Import first set of files, natural ALK1 vs ALK4 (HM)

setwd("ARPEGGIO_analyses/MILD_alk1_v_alk4")
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

## Filter low coverage scaffolds

HM_hal_lowC_scaffolds  <- read.table("HM_RS7_G4_hal_LowCovRegions.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds  <- read.table("HM_RS7_G4_lyr_LowCovRegions.txt", quote="\"", comment.char="")

## We use GRanges to exclude DMRs overlapping low coverage regions

HM_hal_lowC_scaffolds_G <- GRanges(seqnames = HM_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_hal_lowC_scaffolds$V2,
                                                    end = HM_hal_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                       end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                       end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
#dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
#dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH),]

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr_sig[-c(overlapCHH_lyr),]

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

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We plot the DMRs for both sides

HM_ALK <- ggplot(dmrs_context_total, aes(x=context, 
                                         y=value)) +
  geom_bar(mapping = aes(fill=hypo_hyper), 
           stat = 'identity', 
           show.legend = TRUE, 
           size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("ALK (Cold Conditions)") +
  ylim(0, 12500)

### Now we import the second condition and repeat the same procedure
### Then we'll plot the two one next to the other


### Import first set of files, natural ALK1 vs ALK4 (LL)

setwd("ARPEGGIO_analyses/STRESS_alk1_v_alk4")
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

## Filter low coverage scaffolds

LL_hal_lowC_scaffolds  <- read.table("LL_RS7_G4_hal_LowCovRegions.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("LL_RS7_G4_lyr_LowCovRegions.txt", quote="\"", comment.char="")

## We use GRanges to exclude DMRs overlapping low coverage regions

HM_hal_lowC_scaffolds_G <- GRanges(seqnames = HM_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_hal_lowC_scaffolds$V2,
                                                    end = HM_hal_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH),]

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr_sig[-c(overlapCHH_lyr),]

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

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We plot the DMRs for both sides

LL_ALK <- ggplot(dmrs_context_total, aes(x=context, 
                                         y=value)) +
  geom_bar(mapping = aes(fill=hypo_hyper), 
           stat = 'identity', 
           show.legend = TRUE, 
           size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("ALK (Hot Conditions)") +
  ylim(0, 12500)

HM_ALK | LL_ALK

### We repeat the whole procedure above for TKS

### Import first set of files, natural TKS1 vs TKS5 (HM)

setwd("ARPEGGIO_analyses/MILD_tks1_v_tks5")
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

## Filter low coverage scaffolds
## We use GRanges to exclude DMRs overlapping low coverage regions

HM_hal_lowC_scaffolds_G <- GRanges(seqnames = HM_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_hal_lowC_scaffolds$V2,
                                                    end = HM_hal_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
#dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
#dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr_sig[-c(overlapCHH_lyr),]

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

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We plot the DMRs for both sides

HM_TKS <- ggplot(dmrs_context_total, aes(x=context, 
                                         y=value)) +
  geom_bar(mapping = aes(fill=hypo_hyper), 
           stat = 'identity', 
           show.legend = TRUE, 
           size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("TKS (Cold Conditions)") +
  ylim(0, 12500)

### Now we import the second condition and repeat the same procedure
### Then we'll plot the two one next to the other


### Import first set of files, natural TKS1 vs TKS5 (LL)

setwd("ARPEGGIO_analyses/STRESS_tks1_v_tks5")
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

## Filter low coverage scaffolds
## We use GRanges to exclude DMRs overlapping low coverage regions

HM_hal_lowC_scaffolds_G <- GRanges(seqnames = HM_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_hal_lowC_scaffolds$V2,
                                                    end = HM_hal_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr_sig[-c(overlapCHH_lyr),]

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

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We plot the DMRs for both sides

LL_TKS <- ggplot(dmrs_context_total, aes(x=context, 
                                         y=value)) +
  geom_bar(mapping = aes(fill=hypo_hyper), 
           stat = 'identity', 
           show.legend = TRUE, 
           size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("TKS (Hot Conditions)") +
  ylim(0, 12500)

HM_TKS | LL_TKS

# We repeat the procedure one last time for diploids
### Import first set of files, hal/lyr1 vs hal/lyr4 (HM)

setwd("ARPEGGIO_analyses/MILD_pro1_v_pro4/halleri/")
dmrseq_output_CG_hal <- fread("CG_context/A_v_B_diploid.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/A_v_B_diploid.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/A_v_B_diploid.txt")
setwd("ARPEGGIO_analyses/MILD_pro1_v_pro4/lyrata/")
dmrseq_output_CG_lyr <- fread("CG_context/A_v_B_diploid.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/A_v_B_diploid.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/A_v_B_diploid.txt")

## Filter only significant regions 

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds
## We use GRanges to exclude DMRs overlapping low coverage regions

HM_hal_lowC_scaffolds_G <- GRanges(seqnames = HM_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_hal_lowC_scaffolds$V2,
                                                    end = HM_hal_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr_sig[-c(overlapCHH_lyr),]

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

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We plot the DMRs for both sides

HM_DIP <- ggplot(dmrs_context_total, aes(x=context, 
                                         y=value)) +
  geom_bar(mapping = aes(fill=hypo_hyper), 
           stat = 'identity', 
           show.legend = TRUE, 
           size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("Diploid (Cold Conditions)") +
  ylim(0, 12500)

### Now we import the second condition and repeat the same procedure
### Then we'll plot the two one next to the other

### Import first set of files, hal/lyr1 vs hal/lyr4 (HM)

setwd("ARPEGGIO_analyses/STRESS_pro1_v_pro4/halleri/")
dmrseq_output_CG_hal <- fread("CG_context/A_v_B_diploid.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/A_v_B_diploid.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/A_v_B_diploid.txt")
setwd("ARPEGGIO_analyses/STRESS_pro1_v_pro4/lyrata/")
dmrseq_output_CG_lyr <- fread("CG_context/A_v_B_diploid.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/A_v_B_diploid.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/A_v_B_diploid.txt")


## Filter only significant regions 

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds
## We use GRanges to exclude DMRs overlapping low coverage regions

HM_hal_lowC_scaffolds_G <- GRanges(seqnames = HM_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_hal_lowC_scaffolds$V2,
                                                    end = HM_hal_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, HM_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

HM_lyr_lowC_scaffolds_G <- GRanges(seqnames = HM_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = HM_lyr_lowC_scaffolds$V2,
                                                    end = HM_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, HM_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr_sig[-c(overlapCHH_lyr),]

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

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

## We plot the DMRs for both sides

LL_DIP <- ggplot(dmrs_context_total, aes(x=context, 
                                         y=value)) +
  geom_bar(mapping = aes(fill=hypo_hyper), 
           stat = 'identity', 
           show.legend = TRUE, 
           size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("Diploid (Hot Conditions)") +
  ylim(0, 12500)

HM_DIP | LL_DIP
