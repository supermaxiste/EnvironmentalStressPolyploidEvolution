### DMR Spatial Distribution

### Import libraries

library(tidyverse)
library(patchwork)

### Path to DMR files

setwd("~/Library/CloudStorage/OneDrive-Personal/PhD/Project/Chapter_3/new_DMR_results/ARPEGGIO_results_alk1Vsyn4_HM/dmrseq")
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

HM_hal_lowC_scaffolds  <- read.table("~/OneDrive/PhD/Project/Chapter_3/new_coverage_analysis/HM_RS7_G4_hal_LowCovRegions.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds  <- read.table("~/OneDrive/PhD/Project/Chapter_3/new_coverage_analysis/HM_RS7_G4_lyr_LowCovRegions.txt", quote="\"", comment.char="")

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
#dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[-c(overlapCHG_hal),]
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

## Define chromosome sizes

chrom_sizes_hal <- data.frame(Chrom = c(1:8),
                          Start = c(rep(1, 8)),
                          End = c(37868808, 
                                  19340014, 
                                  30847184,
                                  27260850, 
                                  25233939, 
                                  28041738,
                                  30089577, 
                                  26531562),
                          Colors = rep("lightgrey", 
                                       8))

chrom_sizes_lyr <- data.frame(Chrom = c(9:16),
                          Start = c(rep(1, 8)),
                          End = c(28718688, 
                                  21180362, 
                                  27642951,
                                  24725123, 
                                  23209782, 
                                  25845826,
                                  27879834, 
                                  20443964),
                          Colors = rep("lightgrey", 
                                       8))

## Loop through chromosomes and plot relative density distribution of DMRs

plot_list_hal <- list()
plot_list_lyr <- list()

for(i in c(1:8)) {
  chrX_hal <- rbind(dmrseq_output_CG_hal_sig, 
                dmrseq_output_CHG_hal_sig, 
                dmrseq_output_CHH_hal_sig) %>% 
    filter(seqnames==paste0("chr", i))
  
  chr_plot_hal <- ggplot(data = chrX_hal, aes(x = start)) + 
    geom_density() +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14)
    ) +
    ggtitle(paste0("Chromosome ", i, " - H-side"))
  
  plot_list_hal[[i]] <- chr_plot_hal
  
  chrX_lyr <- rbind(dmrseq_output_CG_lyr_sig, 
                    dmrseq_output_CHG_lyr_sig, 
                    dmrseq_output_CHH_lyr_sig) %>% 
    filter(seqnames==paste0("chr", i+8))
  
  chr_plot_lyr <- ggplot(data = chrX_lyr, aes(x = start)) + 
    geom_density() +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14)
    ) +
    ggtitle(paste0("Chromosome ", i, " - L-side"))
  
  plot_list_lyr[[i]] <- chr_plot_lyr
}

plot_list_hal[[1]] /
  plot_list_hal[[2]] /
  plot_list_hal[[3]] /
  plot_list_hal[[4]] /
  plot_list_hal[[5]] /
  plot_list_hal[[6]] /
  plot_list_hal[[7]] /
  plot_list_hal[[8]]

plot_list_lyr[[1]] /
  plot_list_lyr[[2]] /
  plot_list_lyr[[3]] /
  plot_list_lyr[[4]] /
  plot_list_lyr[[5]] /
  plot_list_lyr[[6]] /
  plot_list_lyr[[7]] /
  plot_list_lyr[[8]]
