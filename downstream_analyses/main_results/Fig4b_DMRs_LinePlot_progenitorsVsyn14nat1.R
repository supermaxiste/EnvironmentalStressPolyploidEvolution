### Script to plot Barplot of DMRs over generations
### Using progenitors as references

## Figure 4

### Import libraries

library(ggplot2)
library(data.table)
library(patchwork)
library(GenomicRanges)

###################################################
# Procedure to compare pro vs SYN1-4 & ALK (HM) ##
##################################################

### Import first set of files, progenitors G1 ALK vs synthetic G1 (HM)

setwd("MILD_syn1_v_pro1/dmrseq")
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

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import first set of files, progenitors G1 ALK vs synthetic G4 (HM)

setwd("MILD_syn4_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

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

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs synthetic generation 4

setwd("MILD_alk1_v_pro1/dmrseq")
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

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
           side = c(rep(c(rep("H-side", 3),
                    rep("L-side", 3)), 3)),
           value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                     dmrs_context_total$value[2] + dmrs_context_total$value[5],
                     dmrs_context_total$value[3] + dmrs_context_total$value[6],
                     dmrs_context_total$value[7] + dmrs_context_total$value[10],
                     dmrs_context_total$value[8] + dmrs_context_total$value[11],
                     dmrs_context_total$value[9] + dmrs_context_total$value[12],
                     dmrs_context_total$value[13] + dmrs_context_total$value[16],
                     dmrs_context_total$value[14] + dmrs_context_total$value[17],
                     dmrs_context_total$value[15] + dmrs_context_total$value[18],
                     dmrs_context_total$value[19] + dmrs_context_total$value[22],
                     dmrs_context_total$value[20] + dmrs_context_total$value[23],
                     dmrs_context_total$value[21] + dmrs_context_total$value[24],
                     dmrs_context_total$value[25] + dmrs_context_total$value[28],
                     dmrs_context_total$value[26] + dmrs_context_total$value[29],
                     dmrs_context_total$value[27] + dmrs_context_total$value[30],
                     dmrs_context_total$value[31] + dmrs_context_total$value[34],
                     dmrs_context_total$value[32] + dmrs_context_total$value[35],
                     dmrs_context_total$value[33] + dmrs_context_total$value[36]),
           comparison = c(rep("PvsS1", 6), rep("PvsS4", 6), rep("PvsN",6)))


HM_ALK <- ggplot(dmrs_context_total_new, aes(x=factor(comparison, levels = c("PvsS1", "PvsS4", "PvsN")), y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, linewidth = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("Progenitors over time - ALK (Mild Conditions)") +
  ylim(0, 36000)

##############################################################################
# Repeat the whole procedure from above to compare pro vs SYN1-4 & ALK (LL) ##
##############################################################################

### Import first set of files, natural TKS vs progenitors (HM)
setwd("STRESS_syn1_v_pro1/dmrseq")
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

## Filter low coverage scaffolds

LL_hal_lowC_scaffolds  <- read.table("LL_RS7_G4_hal_LowCovRegions.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("LL_RS7_G4_lyr_LowCovRegions.txt", quote="\"", comment.char="")

## We use GRanges to exclude DMRs overlapping low coverage regions

LL_hal_lowC_scaffolds_G <- GRanges(seqnames = LL_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_hal_lowC_scaffolds$V2,
                                                    end = LL_hal_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
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

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import first set of files, progenitors G1 ALK vs synthetic G4 (LL)

setwd("STRESS_syn4_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## We use GRanges to exclude DMRs overlapping low coverage regions

LL_hal_lowC_scaffolds_G <- GRanges(seqnames = LL_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_hal_lowC_scaffolds$V2,
                                                    end = LL_hal_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
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

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs progenitor

setwd("STRESS_alk1_v_pro1/dmrseq")
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

## We use GRanges to exclude DMRs overlapping low coverage regions

LL_hal_lowC_scaffolds_G <- GRanges(seqnames = LL_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_hal_lowC_scaffolds$V2,
                                                    end = LL_hal_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
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

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
                                     side = c(rep(c(rep("H-side", 3),
                                                    rep("L-side", 3)), 3)),
                                     value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                                               dmrs_context_total$value[2] + dmrs_context_total$value[5],
                                               dmrs_context_total$value[3] + dmrs_context_total$value[6],
                                               dmrs_context_total$value[7] + dmrs_context_total$value[10],
                                               dmrs_context_total$value[8] + dmrs_context_total$value[11],
                                               dmrs_context_total$value[9] + dmrs_context_total$value[12],
                                               dmrs_context_total$value[13] + dmrs_context_total$value[16],
                                               dmrs_context_total$value[14] + dmrs_context_total$value[17],
                                               dmrs_context_total$value[15] + dmrs_context_total$value[18],
                                               dmrs_context_total$value[19] + dmrs_context_total$value[22],
                                               dmrs_context_total$value[20] + dmrs_context_total$value[23],
                                               dmrs_context_total$value[21] + dmrs_context_total$value[24],
                                               dmrs_context_total$value[25] + dmrs_context_total$value[28],
                                               dmrs_context_total$value[26] + dmrs_context_total$value[29],
                                               dmrs_context_total$value[27] + dmrs_context_total$value[30],
                                               dmrs_context_total$value[31] + dmrs_context_total$value[34],
                                               dmrs_context_total$value[32] + dmrs_context_total$value[35],
                                               dmrs_context_total$value[33] + dmrs_context_total$value[36]),
                                     comparison = c(rep("PvsS1", 6), rep("PvsS4", 6), rep("PvsN",6)))


LL_ALK <- ggplot(dmrs_context_total_new, aes(x=factor(comparison, levels = c("PvsS1", "PvsS4", "PvsN")), y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, linewidth = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("Progenitors over time - ALK (Stress Conditions)") +
  ylim(0, 36000)


###################################################
# Procedure to compare pro vs SYN1-4 & TKS (HM) ##
##################################################

### Import first set of files, progenitors G1 TKS vs synthetic G1 (HM)

setwd("MILD_syn1_v_pro1/dmrseq")
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

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import first set of files, progenitors G1 TKS vs synthetic G4 (HM)

setwd("MILD_syn4_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

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

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs synthetic generation 4

setwd("MILD_tks1_v_pro1/dmrseq")
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

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
                                     side = c(rep(c(rep("H-side", 3),
                                                    rep("L-side", 3)), 3)),
                                     value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                                               dmrs_context_total$value[2] + dmrs_context_total$value[5],
                                               dmrs_context_total$value[3] + dmrs_context_total$value[6],
                                               dmrs_context_total$value[7] + dmrs_context_total$value[10],
                                               dmrs_context_total$value[8] + dmrs_context_total$value[11],
                                               dmrs_context_total$value[9] + dmrs_context_total$value[12],
                                               dmrs_context_total$value[13] + dmrs_context_total$value[16],
                                               dmrs_context_total$value[14] + dmrs_context_total$value[17],
                                               dmrs_context_total$value[15] + dmrs_context_total$value[18],
                                               dmrs_context_total$value[19] + dmrs_context_total$value[22],
                                               dmrs_context_total$value[20] + dmrs_context_total$value[23],
                                               dmrs_context_total$value[21] + dmrs_context_total$value[24],
                                               dmrs_context_total$value[25] + dmrs_context_total$value[28],
                                               dmrs_context_total$value[26] + dmrs_context_total$value[29],
                                               dmrs_context_total$value[27] + dmrs_context_total$value[30],
                                               dmrs_context_total$value[31] + dmrs_context_total$value[34],
                                               dmrs_context_total$value[32] + dmrs_context_total$value[35],
                                               dmrs_context_total$value[33] + dmrs_context_total$value[36]),
                                     comparison = c(rep("PvsS1", 6), rep("PvsS4", 6), rep("PvsN",6)))


HM_TKS <- ggplot(dmrs_context_total_new, aes(x=factor(comparison, levels = c("PvsS1", "PvsS4", "PvsN")), y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, linewidth = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("Progenitors over time - TKS (Mild Conditions)") +
  ylim(0, 36000)

##############################################################################
# Repeat the whole procedure from above to compare pro vs SYN1-4 & TKS (LL) ##
##############################################################################

### Import first set of files, natural TKS vs progenitors (HM)
setwd("STRESS_syn1_v_pro1/dmrseq")
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

## Filter low coverage scaffolds

LL_hal_lowC_scaffolds  <- read.table("LL_RS7_G4_hal_LowCovRegions.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("LL_RS7_G4_lyr_LowCovRegions.txt", quote="\"", comment.char="")

## We use GRanges to exclude DMRs overlapping low coverage regions

LL_hal_lowC_scaffolds_G <- GRanges(seqnames = LL_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_hal_lowC_scaffolds$V2,
                                                    end = LL_hal_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
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

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import first set of files, progenitors G1 TKS vs synthetic G4 (LL)

setwd("STRESS_syn4_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## We use GRanges to exclude DMRs overlapping low coverage regions

LL_hal_lowC_scaffolds_G <- GRanges(seqnames = LL_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_hal_lowC_scaffolds$V2,
                                                    end = LL_hal_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
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

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs progenitor

setwd("STRESS_tks1_v_pro1/dmrseq")
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

## We use GRanges to exclude DMRs overlapping low coverage regions

LL_hal_lowC_scaffolds_G <- GRanges(seqnames = LL_hal_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_hal_lowC_scaffolds$V2,
                                                    end = LL_hal_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CG_hal_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_hal_sig$start,
                                                       end = dmrseq_output_CG_hal_sig$end))

dmrseq_output_CHG_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHG_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_hal_sig$start,
                                                        end = dmrseq_output_CHG_hal_sig$end))

dmrseq_output_CHH_hal_sig_G <- GRanges(seqnames = dmrseq_output_CHH_hal_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_hal_sig$start,
                                                        end = dmrseq_output_CHH_hal_sig$end))

overlapCG_hal <- unique(findOverlaps(dmrseq_output_CG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[-c(overlapCG_hal),]
overlapCHG_hal <- unique(findOverlaps(dmrseq_output_CHG_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
overlapCHH_hal <- unique(findOverlaps(dmrseq_output_CHH_hal_sig_G, LL_hal_lowC_scaffolds_G)@from)
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[-c(overlapCHH_hal),]

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

LL_lyr_lowC_scaffolds_G <- GRanges(seqnames = LL_lyr_lowC_scaffolds$V1,
                                   ranges = IRanges(start = LL_lyr_lowC_scaffolds$V2,
                                                    end = LL_lyr_lowC_scaffolds$V3))

dmrseq_output_CG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CG_lyr_sig$seqnames,
                                      ranges = IRanges(start = dmrseq_output_CG_lyr_sig$start,
                                                       end = dmrseq_output_CG_lyr_sig$end))

dmrseq_output_CHG_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHG_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHG_lyr_sig$start,
                                                        end = dmrseq_output_CHG_lyr_sig$end))

dmrseq_output_CHH_lyr_sig_G <- GRanges(seqnames = dmrseq_output_CHH_lyr_sig$seqnames,
                                       ranges = IRanges(start = dmrseq_output_CHH_lyr_sig$start,
                                                        end = dmrseq_output_CHH_lyr_sig$end))

overlapCG_lyr <- unique(findOverlaps(dmrseq_output_CG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr_sig[-c(overlapCG_lyr),]
overlapCHG_lyr <- unique(findOverlaps(dmrseq_output_CHG_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr_sig[-c(overlapCHG_lyr),]
overlapCHH_lyr <- unique(findOverlaps(dmrseq_output_CHH_lyr_sig_G, LL_lyr_lowC_scaffolds_G)@from)
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

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
                                     side = c(rep(c(rep("H-side", 3),
                                                    rep("L-side", 3)), 3)),
                                     value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                                               dmrs_context_total$value[2] + dmrs_context_total$value[5],
                                               dmrs_context_total$value[3] + dmrs_context_total$value[6],
                                               dmrs_context_total$value[7] + dmrs_context_total$value[10],
                                               dmrs_context_total$value[8] + dmrs_context_total$value[11],
                                               dmrs_context_total$value[9] + dmrs_context_total$value[12],
                                               dmrs_context_total$value[13] + dmrs_context_total$value[16],
                                               dmrs_context_total$value[14] + dmrs_context_total$value[17],
                                               dmrs_context_total$value[15] + dmrs_context_total$value[18],
                                               dmrs_context_total$value[19] + dmrs_context_total$value[22],
                                               dmrs_context_total$value[20] + dmrs_context_total$value[23],
                                               dmrs_context_total$value[21] + dmrs_context_total$value[24],
                                               dmrs_context_total$value[25] + dmrs_context_total$value[28],
                                               dmrs_context_total$value[26] + dmrs_context_total$value[29],
                                               dmrs_context_total$value[27] + dmrs_context_total$value[30],
                                               dmrs_context_total$value[31] + dmrs_context_total$value[34],
                                               dmrs_context_total$value[32] + dmrs_context_total$value[35],
                                               dmrs_context_total$value[33] + dmrs_context_total$value[36]),
                                     comparison = c(rep("PvsS1", 6), rep("PvsS4", 6), rep("PvsN",6)))


LL_TKS <- ggplot(dmrs_context_total_new, aes(x=factor(comparison, levels = c("PvsS1", "PvsS4", "PvsN")), y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, linewidth = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("Progenitors over time - TKS (Stress Conditions)") +
  ylim(0, 36000)



(HM_ALK | LL_ALK) /
  (HM_TKS | LL_TKS)