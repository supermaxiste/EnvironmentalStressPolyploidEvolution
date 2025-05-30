# Visualizing both RNAseq and BSseq results

# Packages

library(data.table)
library(UpSetR)
library(ggplot2)
library(patchwork)
library(tidyverse)

# Methylation files

setwd("~/OneDrive/PhD/Project/Chapter_3/DMR_results/")

# Cold conditions

# halleri-side

HM_hal_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

# lyrata-side

HM_lyr_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

# Hot conditions

# halleri-side

LL_hal_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

# lyrata-side

LL_lyr_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

# Since methylation files do not show the geneID
# We will use the output DM files from ARPEGGIO
# to add the geneID

# Cold conditions

# halleri-side

HM_hal_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

HM_hal_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

HM_hal_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_hal_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

HM_hal_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

HM_hal_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

# lyrata-side

HM_lyr_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

HM_lyr_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

HM_lyr_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_lyr_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

HM_lyr_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

HM_lyr_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

# Hot conditions

# halleri-side

LL_hal_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

LL_hal_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

LL_hal_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_hal_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

LL_hal_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

LL_hal_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

# lyrata-side

LL_lyr_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

LL_lyr_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

LL_lyr_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_lyr_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

LL_lyr_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

LL_lyr_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

# We add the geneID to the dmrseq output files
# First we create a unique identifier for the DMR position

create_ID1 <- function(file){
  new_ID <- paste0(file$seqname,
                   "_",
                   file$start,
                   "_",
                   file$end)
  return(new_ID)
}

create_ID2 <- function(file){
  new_ID <- paste0(file$seqname,
                   "_",
                   file$region_start,
                   "_",
                   file$region_end)
  return(new_ID)
}

# halleri-side 

HM_hal_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(HM_hal_v_synG1_DMG_CG_geneIDs)
HM_hal_v_synG1_DMG_CG_pos <- create_ID1(HM_hal_v_synG1_DMG_CG)
HM_hal_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(HM_hal_v_synG1_DMG_CHG_geneIDs)
HM_hal_v_synG1_DMG_CHG_pos <- create_ID1(HM_hal_v_synG1_DMG_CHG)
HM_hal_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(HM_hal_v_synG1_DMG_CHH_geneIDs)
HM_hal_v_synG1_DMG_CHH_pos <- create_ID1(HM_hal_v_synG1_DMG_CHH)

HM_hal_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(HM_hal_v_synG4_DMG_CG_geneIDs)
HM_hal_v_synG4_DMG_CG_pos <- create_ID1(HM_hal_v_synG4_DMG_CG)
HM_hal_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(HM_hal_v_synG4_DMG_CHG_geneIDs)
HM_hal_v_synG4_DMG_CHG_pos <- create_ID1(HM_hal_v_synG4_DMG_CHG)
HM_hal_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(HM_hal_v_synG4_DMG_CHH_geneIDs)
HM_hal_v_synG4_DMG_CHH_pos <- create_ID1(HM_hal_v_synG4_DMG_CHH)

LL_hal_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(LL_hal_v_synG1_DMG_CG_geneIDs)
LL_hal_v_synG1_DMG_CG_pos <- create_ID1(LL_hal_v_synG1_DMG_CG)
LL_hal_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(LL_hal_v_synG1_DMG_CHG_geneIDs)
LL_hal_v_synG1_DMG_CHG_pos <- create_ID1(LL_hal_v_synG1_DMG_CHG)
LL_hal_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(LL_hal_v_synG1_DMG_CHH_geneIDs)
LL_hal_v_synG1_DMG_CHH_pos <- create_ID1(LL_hal_v_synG1_DMG_CHH)

LL_hal_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(LL_hal_v_synG4_DMG_CG_geneIDs)
LL_hal_v_synG4_DMG_CG_pos <- create_ID1(LL_hal_v_synG4_DMG_CG)
LL_hal_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(LL_hal_v_synG4_DMG_CHG_geneIDs)
LL_hal_v_synG4_DMG_CHG_pos <- create_ID1(LL_hal_v_synG4_DMG_CHG)
LL_hal_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(LL_hal_v_synG4_DMG_CHH_geneIDs)
LL_hal_v_synG4_DMG_CHH_pos <- create_ID1(LL_hal_v_synG4_DMG_CHH)

# lyrata side

HM_lyr_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(HM_lyr_v_synG1_DMG_CG_geneIDs)
HM_lyr_v_synG1_DMG_CG_pos <- create_ID1(HM_lyr_v_synG1_DMG_CG)
HM_lyr_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(HM_lyr_v_synG1_DMG_CHG_geneIDs)
HM_lyr_v_synG1_DMG_CHG_pos <- create_ID1(HM_lyr_v_synG1_DMG_CHG)
HM_lyr_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(HM_lyr_v_synG1_DMG_CHH_geneIDs)
HM_lyr_v_synG1_DMG_CHH_pos <- create_ID1(HM_lyr_v_synG1_DMG_CHH)

HM_lyr_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(HM_lyr_v_synG4_DMG_CG_geneIDs)
HM_lyr_v_synG4_DMG_CG_pos <- create_ID1(HM_lyr_v_synG4_DMG_CG)
HM_lyr_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(HM_lyr_v_synG4_DMG_CHG_geneIDs)
HM_lyr_v_synG4_DMG_CHG_pos <- create_ID1(HM_lyr_v_synG4_DMG_CHG)
HM_lyr_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(HM_lyr_v_synG4_DMG_CHH_geneIDs)
HM_lyr_v_synG4_DMG_CHH_pos <- create_ID1(HM_lyr_v_synG4_DMG_CHH)

LL_lyr_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(LL_lyr_v_synG1_DMG_CG_geneIDs)
LL_lyr_v_synG1_DMG_CG_pos <- create_ID1(LL_lyr_v_synG1_DMG_CG)
LL_lyr_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(LL_lyr_v_synG1_DMG_CHG_geneIDs)
LL_lyr_v_synG1_DMG_CHG_pos <- create_ID1(LL_lyr_v_synG1_DMG_CHG)
LL_lyr_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(LL_lyr_v_synG1_DMG_CHH_geneIDs)
LL_lyr_v_synG1_DMG_CHH_pos <- create_ID1(LL_lyr_v_synG1_DMG_CHH)

LL_lyr_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(LL_lyr_v_synG4_DMG_CG_geneIDs)
LL_lyr_v_synG4_DMG_CG_pos <- create_ID1(LL_lyr_v_synG4_DMG_CG)
LL_lyr_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(LL_lyr_v_synG4_DMG_CHG_geneIDs)
LL_lyr_v_synG4_DMG_CHG_pos <- create_ID1(LL_lyr_v_synG4_DMG_CHG)
LL_lyr_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(LL_lyr_v_synG4_DMG_CHH_geneIDs)
LL_lyr_v_synG4_DMG_CHH_pos <- create_ID1(LL_lyr_v_synG4_DMG_CHH)

# function to add raw methylation difference
# to DM files

add_rawMeth <- function(DM_pos_file, 
                        dmrseq_pos_file,
                        DM_file,
                        dmrseq_file){
  index <- match(DM_pos_file, dmrseq_pos_file)
  new_DM_file <- mutate(DM_file,
                        rawMeth = dmrseq_file$rawDiff[index])
  return(new_DM_file)
}

#halleri side

HM_hal_v_synG1_DMG_CG_final <- add_rawMeth(HM_hal_v_synG1_DMG_CG_geneIDs_pos,
                                           HM_hal_v_synG1_DMG_CG_pos,
                                           HM_hal_v_synG1_DMG_CG_geneIDs,
                                           HM_hal_v_synG1_DMG_CG)

HM_hal_v_synG1_DMG_CHG_final <- add_rawMeth(HM_hal_v_synG1_DMG_CHG_geneIDs_pos,
                                            HM_hal_v_synG1_DMG_CHG_pos,
                                            HM_hal_v_synG1_DMG_CHG_geneIDs,
                                            HM_hal_v_synG1_DMG_CHG)

HM_hal_v_synG1_DMG_CHH_final <- add_rawMeth(HM_hal_v_synG1_DMG_CHH_geneIDs_pos,
                                            HM_hal_v_synG1_DMG_CHH_pos,
                                            HM_hal_v_synG1_DMG_CHH_geneIDs,
                                            HM_hal_v_synG1_DMG_CHH)

HM_hal_v_synG4_DMG_CG_final <- add_rawMeth(HM_hal_v_synG4_DMG_CG_geneIDs_pos,
                                           HM_hal_v_synG4_DMG_CG_pos,
                                           HM_hal_v_synG4_DMG_CG_geneIDs,
                                           HM_hal_v_synG4_DMG_CG)

HM_hal_v_synG4_DMG_CHG_final <- add_rawMeth(HM_hal_v_synG4_DMG_CHG_geneIDs_pos,
                                            HM_hal_v_synG4_DMG_CHG_pos,
                                            HM_hal_v_synG4_DMG_CHG_geneIDs,
                                            HM_hal_v_synG4_DMG_CHG)

HM_hal_v_synG4_DMG_CHH_final <- add_rawMeth(HM_hal_v_synG4_DMG_CHH_geneIDs_pos,
                                            HM_hal_v_synG4_DMG_CHH_pos,
                                            HM_hal_v_synG4_DMG_CHH_geneIDs,
                                            HM_hal_v_synG4_DMG_CHH)

LL_hal_v_synG1_DMG_CG_final <- add_rawMeth(LL_hal_v_synG1_DMG_CG_geneIDs_pos,
                                           LL_hal_v_synG1_DMG_CG_pos,
                                           LL_hal_v_synG1_DMG_CG_geneIDs,
                                           LL_hal_v_synG1_DMG_CG)

LL_hal_v_synG1_DMG_CHG_final <- add_rawMeth(LL_hal_v_synG1_DMG_CHG_geneIDs_pos,
                                            LL_hal_v_synG1_DMG_CHG_pos,
                                            LL_hal_v_synG1_DMG_CHG_geneIDs,
                                            LL_hal_v_synG1_DMG_CHG)

LL_hal_v_synG1_DMG_CHH_final <- add_rawMeth(LL_hal_v_synG1_DMG_CHH_geneIDs_pos,
                                            LL_hal_v_synG1_DMG_CHH_pos,
                                            LL_hal_v_synG1_DMG_CHH_geneIDs,
                                            LL_hal_v_synG1_DMG_CHH)

LL_hal_v_synG4_DMG_CG_final <- add_rawMeth(LL_hal_v_synG4_DMG_CG_geneIDs_pos,
                                           LL_hal_v_synG4_DMG_CG_pos,
                                           LL_hal_v_synG4_DMG_CG_geneIDs,
                                           LL_hal_v_synG4_DMG_CG)

LL_hal_v_synG4_DMG_CHG_final <- add_rawMeth(LL_hal_v_synG4_DMG_CHG_geneIDs_pos,
                                            LL_hal_v_synG4_DMG_CHG_pos,
                                            LL_hal_v_synG4_DMG_CHG_geneIDs,
                                            LL_hal_v_synG4_DMG_CHG)

LL_hal_v_synG4_DMG_CHH_final <- add_rawMeth(LL_hal_v_synG4_DMG_CHH_geneIDs_pos,
                                            LL_hal_v_synG4_DMG_CHH_pos,
                                            LL_hal_v_synG4_DMG_CHH_geneIDs,
                                            LL_hal_v_synG4_DMG_CHH)

#lyrata side

HM_lyr_v_synG1_DMG_CG_final <- add_rawMeth(HM_lyr_v_synG1_DMG_CG_geneIDs_pos,
                                           HM_lyr_v_synG1_DMG_CG_pos,
                                           HM_lyr_v_synG1_DMG_CG_geneIDs,
                                           HM_lyr_v_synG1_DMG_CG)

HM_lyr_v_synG1_DMG_CHG_final <- add_rawMeth(HM_lyr_v_synG1_DMG_CHG_geneIDs_pos,
                                            HM_lyr_v_synG1_DMG_CHG_pos,
                                            HM_lyr_v_synG1_DMG_CHG_geneIDs,
                                            HM_lyr_v_synG1_DMG_CHG)

HM_lyr_v_synG1_DMG_CHH_final <- add_rawMeth(HM_lyr_v_synG1_DMG_CHH_geneIDs_pos,
                                            HM_lyr_v_synG1_DMG_CHH_pos,
                                            HM_lyr_v_synG1_DMG_CHH_geneIDs,
                                            HM_lyr_v_synG1_DMG_CHH)

HM_lyr_v_synG4_DMG_CG_final <- add_rawMeth(HM_lyr_v_synG4_DMG_CG_geneIDs_pos,
                                           HM_lyr_v_synG4_DMG_CG_pos,
                                           HM_lyr_v_synG4_DMG_CG_geneIDs,
                                           HM_lyr_v_synG4_DMG_CG)

HM_lyr_v_synG4_DMG_CHG_final <- add_rawMeth(HM_lyr_v_synG4_DMG_CHG_geneIDs_pos,
                                            HM_lyr_v_synG4_DMG_CHG_pos,
                                            HM_lyr_v_synG4_DMG_CHG_geneIDs,
                                            HM_lyr_v_synG4_DMG_CHG)

HM_lyr_v_synG4_DMG_CHH_final <- add_rawMeth(HM_lyr_v_synG4_DMG_CHH_geneIDs_pos,
                                            HM_lyr_v_synG4_DMG_CHH_pos,
                                            HM_lyr_v_synG4_DMG_CHH_geneIDs,
                                            HM_lyr_v_synG4_DMG_CHH)

LL_lyr_v_synG1_DMG_CG_final <- add_rawMeth(LL_lyr_v_synG1_DMG_CG_geneIDs_pos,
                                           LL_lyr_v_synG1_DMG_CG_pos,
                                           LL_lyr_v_synG1_DMG_CG_geneIDs,
                                           LL_lyr_v_synG1_DMG_CG)

LL_lyr_v_synG1_DMG_CHG_final <- add_rawMeth(LL_lyr_v_synG1_DMG_CHG_geneIDs_pos,
                                            LL_lyr_v_synG1_DMG_CHG_pos,
                                            LL_lyr_v_synG1_DMG_CHG_geneIDs,
                                            LL_lyr_v_synG1_DMG_CHG)

LL_lyr_v_synG1_DMG_CHH_final <- add_rawMeth(LL_lyr_v_synG1_DMG_CHH_geneIDs_pos,
                                            LL_lyr_v_synG1_DMG_CHH_pos,
                                            LL_lyr_v_synG1_DMG_CHH_geneIDs,
                                            LL_lyr_v_synG1_DMG_CHH)

LL_lyr_v_synG4_DMG_CG_final <- add_rawMeth(LL_lyr_v_synG4_DMG_CG_geneIDs_pos,
                                           LL_lyr_v_synG4_DMG_CG_pos,
                                           LL_lyr_v_synG4_DMG_CG_geneIDs,
                                           LL_lyr_v_synG4_DMG_CG)

LL_lyr_v_synG4_DMG_CHG_final <- add_rawMeth(LL_lyr_v_synG4_DMG_CHG_geneIDs_pos,
                                            LL_lyr_v_synG4_DMG_CHG_pos,
                                            LL_lyr_v_synG4_DMG_CHG_geneIDs,
                                            LL_lyr_v_synG4_DMG_CHG)

LL_lyr_v_synG4_DMG_CHH_final <- add_rawMeth(LL_lyr_v_synG4_DMG_CHH_geneIDs_pos,
                                            LL_lyr_v_synG4_DMG_CHH_pos,
                                            LL_lyr_v_synG4_DMG_CHH_geneIDs,
                                            LL_lyr_v_synG4_DMG_CHH)

# Expression files 

setwd("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG")

# Cold conditions

HM_hal_v_synG1_DEG <- read.delim("HM_halG1_v_RS7G1_DEG.txt")
HM_hal_v_synG4_DEG <- read.delim("HM_halG1_v_RS7G4_DEG.txt")
HM_lyr_v_synG1_DEG <- read.delim("HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyr_v_synG4_DEG <- read.delim("HM_lyrG1_v_RS7G4_DEG.txt")

# Hot conditions

LL_hal_v_synG1_DEG <- read.delim("LL_halG1_v_RS7G1_DEG.txt")
LL_hal_v_synG4_DEG <- read.delim("LL_halG1_v_RS7G4_DEG.txt")
LL_lyr_v_synG1_DEG <- read.delim("LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyr_v_synG4_DEG <- read.delim("LL_lyrG1_v_RS7G4_DEG.txt")


# Find genes showing both methylation and expression changes
# and add logFC to DM file

add_logFC <- function(DM_final, DEG_file){
  index <- match(DM_final$geneID, DEG_file$geneID)
  new_DM_file <- mutate(DM_final,
                        logFC = DEG_file$logFC[index])
  return(new_DM_file)
}

# cold conditions 

HM_hal_v_synG1_DMG_CG_all <- add_logFC(HM_hal_v_synG1_DMG_CG_final,
                                       HM_hal_v_synG1_DEG)
HM_hal_v_synG1_DMG_CHG_all <- add_logFC(HM_hal_v_synG1_DMG_CHG_final,
                                        HM_hal_v_synG1_DEG)
HM_hal_v_synG1_DMG_CHH_all <- add_logFC(HM_hal_v_synG1_DMG_CHH_final,
                                        HM_hal_v_synG1_DEG)

HM_hal_v_synG4_DMG_CG_all <- add_logFC(HM_hal_v_synG4_DMG_CG_final,
                                       HM_hal_v_synG4_DEG)
HM_hal_v_synG4_DMG_CHG_all <- add_logFC(HM_hal_v_synG4_DMG_CHG_final,
                                        HM_hal_v_synG4_DEG)
HM_hal_v_synG4_DMG_CHH_all <- add_logFC(HM_hal_v_synG4_DMG_CHH_final,
                                        HM_hal_v_synG4_DEG)

HM_lyr_v_synG1_DMG_CG_all <- add_logFC(HM_lyr_v_synG1_DMG_CG_final,
                                       HM_lyr_v_synG1_DEG)
HM_lyr_v_synG1_DMG_CHG_all <- add_logFC(HM_lyr_v_synG1_DMG_CHG_final,
                                        HM_lyr_v_synG1_DEG)
HM_lyr_v_synG1_DMG_CHH_all <- add_logFC(HM_lyr_v_synG1_DMG_CHH_final,
                                        HM_lyr_v_synG1_DEG)

HM_lyr_v_synG4_DMG_CG_all <- add_logFC(HM_lyr_v_synG4_DMG_CG_final,
                                       HM_lyr_v_synG4_DEG)
HM_lyr_v_synG4_DMG_CHG_all <- add_logFC(HM_lyr_v_synG4_DMG_CHG_final,
                                        HM_lyr_v_synG4_DEG)
HM_lyr_v_synG4_DMG_CHH_all <- add_logFC(HM_lyr_v_synG4_DMG_CHH_final,
                                        HM_lyr_v_synG4_DEG)

# hot conditions 

LL_hal_v_synG1_DMG_CG_all <- add_logFC(LL_hal_v_synG1_DMG_CG_final,
                                       LL_hal_v_synG1_DEG)
LL_hal_v_synG1_DMG_CHG_all <- add_logFC(LL_hal_v_synG1_DMG_CHG_final,
                                        LL_hal_v_synG1_DEG)
LL_hal_v_synG1_DMG_CHH_all <- add_logFC(LL_hal_v_synG1_DMG_CHH_final,
                                        LL_hal_v_synG1_DEG)

LL_hal_v_synG4_DMG_CG_all <- add_logFC(LL_hal_v_synG4_DMG_CG_final,
                                       LL_hal_v_synG4_DEG)
LL_hal_v_synG4_DMG_CHG_all <- add_logFC(LL_hal_v_synG4_DMG_CHG_final,
                                        LL_hal_v_synG4_DEG)
LL_hal_v_synG4_DMG_CHH_all <- add_logFC(LL_hal_v_synG4_DMG_CHH_final,
                                        LL_hal_v_synG4_DEG)

LL_lyr_v_synG1_DMG_CG_all <- add_logFC(LL_lyr_v_synG1_DMG_CG_final,
                                       LL_lyr_v_synG1_DEG)
LL_lyr_v_synG1_DMG_CHG_all <- add_logFC(LL_lyr_v_synG1_DMG_CHG_final,
                                        LL_lyr_v_synG1_DEG)
LL_lyr_v_synG1_DMG_CHH_all <- add_logFC(LL_lyr_v_synG1_DMG_CHH_final,
                                        LL_lyr_v_synG1_DEG)

LL_lyr_v_synG4_DMG_CG_all <- add_logFC(LL_lyr_v_synG4_DMG_CG_final,
                                       LL_lyr_v_synG4_DEG)
LL_lyr_v_synG4_DMG_CHG_all <- add_logFC(LL_lyr_v_synG4_DMG_CHG_final,
                                        LL_lyr_v_synG4_DEG)
LL_lyr_v_synG4_DMG_CHH_all <- add_logFC(LL_lyr_v_synG4_DMG_CHH_final,
                                        LL_lyr_v_synG4_DEG)

# Plot correlation

plot_correlation <- function(file_all){
  file_test <- cor.test(file_all$logFC,
                        file_all$rawMeth)
  file_text <- paste0("Correlation test:\n t = ", 
                      round(file_test$statistic, 2),
                      ", df = ",
                      file_test$parameter,
                      ", CI = [",
                      round(file_test$conf.int[1], 2),
                      ", ",
                      round(file_test$conf.int[2], 2),
                      "], p-value = ",
                      format.pval(file_test$p.value,
                                  digits = 2))
  final_plot <- ggplot(file_all, aes(x=logFC, y=rawMeth)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab("Beta value (methylation change)") +
    ylab("logFC (expression change)") +
    ggtitle(file_text) +
    theme(plot.title = element_text(size=10))
  return(final_plot)
}

# halleri

HM_hal_v_synG1_CG <- plot_correlation(HM_hal_v_synG1_DMG_CG_all)
HM_hal_v_synG1_CHG <- plot_correlation(HM_hal_v_synG1_DMG_CHG_all)
HM_hal_v_synG1_CHH <- plot_correlation(HM_hal_v_synG1_DMG_CHH_all)

LL_hal_v_synG1_CG <- plot_correlation(LL_hal_v_synG1_DMG_CG_all)
LL_hal_v_synG1_CHG <- plot_correlation(LL_hal_v_synG1_DMG_CHG_all)
LL_hal_v_synG1_CHH <- plot_correlation(LL_hal_v_synG1_DMG_CHH_all)

HM_hal_v_synG4_CG <- plot_correlation(HM_hal_v_synG4_DMG_CG_all)
HM_hal_v_synG4_CHG <- plot_correlation(HM_hal_v_synG4_DMG_CHG_all)
HM_hal_v_synG4_CHH <- plot_correlation(HM_hal_v_synG4_DMG_CHH_all)

LL_hal_v_synG4_CG <- plot_correlation(LL_hal_v_synG4_DMG_CG_all)
LL_hal_v_synG4_CHG <- plot_correlation(LL_hal_v_synG4_DMG_CHG_all)
LL_hal_v_synG4_CHH <- plot_correlation(LL_hal_v_synG4_DMG_CHH_all)

# lyrata 

HM_lyr_v_synG1_CG <- plot_correlation(HM_lyr_v_synG1_DMG_CG_all)
HM_lyr_v_synG1_CHG <- plot_correlation(HM_lyr_v_synG1_DMG_CHG_all)
HM_lyr_v_synG1_CHH <- plot_correlation(HM_lyr_v_synG1_DMG_CHH_all)

LL_lyr_v_synG1_CG <- plot_correlation(LL_lyr_v_synG1_DMG_CG_all)
LL_lyr_v_synG1_CHG <- plot_correlation(LL_lyr_v_synG1_DMG_CHG_all)
LL_lyr_v_synG1_CHH <- plot_correlation(LL_lyr_v_synG1_DMG_CHH_all)

HM_lyr_v_synG4_CG <- plot_correlation(HM_lyr_v_synG4_DMG_CG_all)
HM_lyr_v_synG4_CHG <- plot_correlation(HM_lyr_v_synG4_DMG_CHG_all)
HM_lyr_v_synG4_CHH <- plot_correlation(HM_lyr_v_synG4_DMG_CHH_all)

LL_lyr_v_synG4_CG <- plot_correlation(LL_lyr_v_synG4_DMG_CG_all)
LL_lyr_v_synG4_CHG <- plot_correlation(LL_lyr_v_synG4_DMG_CHG_all)
LL_lyr_v_synG4_CHH <- plot_correlation(LL_lyr_v_synG4_DMG_CHH_all)

# Plot all together

(HM_hal_v_synG1_CG | HM_lyr_v_synG1_CG | LL_hal_v_synG1_CG | LL_lyr_v_synG1_CG) /
  (HM_hal_v_synG1_CHG | HM_lyr_v_synG1_CHG | LL_hal_v_synG1_CHG | LL_lyr_v_synG1_CHG) /
  (HM_hal_v_synG1_CHH | HM_lyr_v_synG1_CHH | LL_hal_v_synG1_CHH | LL_lyr_v_synG1_CHH)

(HM_hal_v_synG4_CG | HM_lyr_v_synG4_CG | LL_hal_v_synG4_CG | LL_lyr_v_synG4_CG) /
  (HM_hal_v_synG4_CHG | HM_lyr_v_synG4_CHG | LL_hal_v_synG4_CHG | LL_lyr_v_synG4_CHG) /
  (HM_hal_v_synG4_CHH | HM_lyr_v_synG4_CHH | LL_hal_v_synG4_CHH | LL_lyr_v_synG4_CHH)
