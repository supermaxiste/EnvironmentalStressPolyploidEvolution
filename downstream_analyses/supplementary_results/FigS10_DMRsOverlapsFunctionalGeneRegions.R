## Import libraries

library(grid)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(gridExtra)
library(patchwork)

# Overlap between DMRs and functional genomic regions

DMRs_overlap <- function(hal_CG_all,
                         hal_CHG_all,
                         hal_CHH_all,
                         lyr_CG_all,
                         lyr_CHG_all,
                         lyr_CHH_all,
                         conditions){
  
  ### Import data: annotations and DMRs
  
  ### Paths to annotation files to be modified according to new system!!!!
  
  Ahal_v2_2 <- read.delim("~/OneDrive/PhD/Project/bs_data/Ahal_genome/Ahal_v2_2.gff",
                                        header=FALSE,
                                        col.names = c("scaffold", "tool", "context", "start", "end", "number",
                                                      "strand", "dot", "extra"))
  Alyr_v2_2 <- read.delim("~/OneDrive/PhD/Project/bs_data/Alyr_genome/Alyr_v2_2_renamed.gff", 
                                                header=FALSE,
                                                col.names = c("scaffold", "tool", "context", "start", "end", "number",
                                                              "strand", "dot", "extra"))
  
  Ahal_anno <- select(Ahal_v2_2, scaffold, context, start, end, strand)
  rm(Ahal_v2_2)
  
  
  Alyr_anno <- select(Alyr_v2_2, scaffold, context, start, end, strand)
  rm(Alyr_v2_2)
  
  ### Filter annotations to keep only promoters
  
  Ahal_anno <- filter(Ahal_anno, context!="gene", context !="transcript")
  Alyr_anno <- filter(Alyr_anno, context!="gene", context !="transcript")
  
  ### Filter DMRs to keep only significant ones
  
  hal_CG_all <- filter(hal_CG_all, qval < 0.05)
  hal_CHG_all <- filter(hal_CHG_all, qval < 0.05)
  hal_CHH_all <- filter(hal_CHH_all, qval < 0.05)
  
  lyr_CG_all <- filter(lyr_CG_all, qval < 0.05)
  lyr_CHG_all <- filter(lyr_CHG_all, qval < 0.05)
  lyr_CHH_all <- filter(lyr_CHH_all, qval < 0.05)
  
  ### Convert annotation to GRanges 
  
  Ahal_anno_ranges <- IRanges(start = Ahal_anno$start,
                              end = Ahal_anno$end)
  Ahal_annoG <- GRanges(seqnames = Ahal_anno$scaffold, 
                        ranges = Ahal_anno_ranges,
                        strand = Ahal_anno$strand)
  
  Alyr_anno_ranges <- IRanges(start = Alyr_anno$start,
                              end = Alyr_anno$end)
  Alyr_annoG <- GRanges(seqnames = Alyr_anno$scaffold, 
                        ranges = Alyr_anno_ranges,
                        strand = Alyr_anno$strand)
  
  ### Convert significant DMRs to GRanges
  
  hal_CG_all_ranges <- IRanges(start = hal_CG_all$start,
                               end = hal_CG_all$end)
  hal_CG_allG <- GRanges(seqnames = hal_CG_all$seqnames,
                         ranges = hal_CG_all_ranges)
  
  hal_CHG_all_ranges <- IRanges(start = hal_CHG_all$start,
                                end = hal_CHG_all$end)
  hal_CHG_allG <- GRanges(seqnames = hal_CHG_all$seqnames,
                          ranges = hal_CHG_all_ranges)
  
  hal_CHH_all_ranges <- IRanges(start = hal_CHH_all$start,
                                end = hal_CHH_all$end)
  hal_CHH_allG <- GRanges(seqnames = hal_CHH_all$seqnames,
                          ranges = hal_CHH_all_ranges)
  
  lyr_CG_all_ranges <- IRanges(start = lyr_CG_all$start,
                               end = lyr_CG_all$end)
  lyr_CG_allG <- GRanges(seqnames = lyr_CG_all$seqnames,
                         ranges = lyr_CG_all_ranges)
  
  lyr_CHG_all_ranges <- IRanges(start = lyr_CHG_all$start,
                                end = lyr_CHG_all$end)
  lyr_CHG_allG <- GRanges(seqnames = lyr_CHG_all$seqnames,
                          ranges = lyr_CHG_all_ranges)
  
  lyr_CHH_all_ranges <- IRanges(start = lyr_CHH_all$start,
                                end = lyr_CHH_all$end)
  lyr_CHH_allG <- GRanges(seqnames = lyr_CHH_all$seqnames,
                          ranges = lyr_CHH_all_ranges)
  
  ### Find overlaps for each context between DMRs and annotations
  
  ### CG context - halleri
  
  overlaps_hal_CG <- findOverlaps(hal_CG_allG, Ahal_annoG)
  
  overlaps_hal_CG_final <- cbind(hal_CG_all[overlaps_hal_CG@from,],
                                 Ahal_anno[overlaps_hal_CG@to,])
  
  colnames(overlaps_hal_CG_final) <- c("seqnames", "start", "end", "width", "strand",
                                       "L", "area", "beta", "stat", "pval", "qval", 
                                       "index.start", "index.end", "index.width",
                                       "scaffold", "context", "start2", "end2", "strand2")
  
  ### CHG context - halleri
  
  overlaps_hal_CHG <- findOverlaps(hal_CHG_allG, Ahal_annoG)
  
  overlaps_hal_CHG_final <- cbind(hal_CHG_all[overlaps_hal_CHG@from,],
                                 Ahal_anno[overlaps_hal_CHG@to,])
  
  colnames(overlaps_hal_CHG_final) <- c("seqnames", "start", "end", "width", "strand",
                                       "L", "area", "beta", "stat", "pval", "qval", 
                                       "index.start", "index.end", "index.width",
                                       "scaffold", "context", "start2", "end2", "strand2")
  
  ### CHH context - halleri
  
  overlaps_hal_CHH <- findOverlaps(hal_CHH_allG, Ahal_annoG)
  
  overlaps_hal_CHH_final <- cbind(hal_CHH_all[overlaps_hal_CHH@from,],
                                 Ahal_anno[overlaps_hal_CHH@to,])
  
  colnames(overlaps_hal_CHH_final) <- c("seqnames", "start", "end", "width", "strand",
                                       "L", "area", "beta", "stat", "pval", "qval", 
                                       "index.start", "index.end", "index.width",
                                       "scaffold", "context", "start2", "end2", "strand2")
  
  
  ### CG context - lyrata
  
  overlaps_lyr_CG <- findOverlaps(lyr_CG_allG, Alyr_annoG)
  
  overlaps_lyr_CG_final <- cbind(lyr_CG_all[overlaps_lyr_CG@from,],
                                 Alyr_anno[overlaps_lyr_CG@to,])
  
  colnames(overlaps_lyr_CG_final) <- c("seqnames", "start", "end", "width", "strand",
                                       "L", "area", "beta", "stat", "pval", "qval", 
                                       "index.start", "index.end", "index.width",
                                       "scaffold", "context", "start2", "end2", "strand2")
  
  ### CHG context - lyrata
  
  overlaps_lyr_CHG <- findOverlaps(lyr_CHG_allG, Alyr_annoG)
  
  overlaps_lyr_CHG_final <- cbind(lyr_CHG_all[overlaps_lyr_CHG@from,],
                                  Alyr_anno[overlaps_lyr_CHG@to,])
  
  colnames(overlaps_lyr_CHG_final) <- c("seqnames", "start", "end", "width", "strand",
                                        "L", "area", "beta", "stat", "pval", "qval", 
                                        "index.start", "index.end", "index.width",
                                        "scaffold", "context", "start2", "end2", "strand2")
  
  ### CHH context - lyrata
  
  overlaps_lyr_CHH <- findOverlaps(lyr_CHH_allG, Alyr_annoG)
  
  overlaps_lyr_CHH_final <- cbind(lyr_CHH_all[overlaps_lyr_CHH@from,],
                                  Alyr_anno[overlaps_lyr_CHH@to,])
  
  colnames(overlaps_lyr_CHH_final) <- c("seqnames", "start", "end", "width", "strand",
                                        "L", "area", "beta", "stat", "pval", "qval", 
                                        "index.start", "index.end", "index.width",
                                        "scaffold", "context", "start2", "end2", "strand2")
  
  ### Stopped here
  ### Plotting proportions of DMRs in each gene functional region for each methylation context
  ### halleri
  CG_hal_data <- data.frame(category = c("CDS",
                                         "exon",
                                         "five_prime_UTR",
                                         "three_prime_UTR"),
                            count = c(sum(overlaps_hal_CG_final$context=="CDS"),
                                      sum(overlaps_hal_CG_final$context=="exon"),
                                      sum(overlaps_hal_CG_final$context=="five_prime_UTR"),
                                      sum(overlaps_hal_CG_final$context=="three_prime_UTR")))
  
  CHG_hal_data <- data.frame(category = c("CDS",
                                          "exon",
                                          "five_prime_UTR",
                                          "three_prime_UTR"),
                            count = c(sum(overlaps_hal_CHG_final$context=="CDS"),
                                      sum(overlaps_hal_CHG_final$context=="exon"),
                                      sum(overlaps_hal_CHG_final$context=="five_prime_UTR"),
                                      sum(overlaps_hal_CHG_final$context=="three_prime_UTR")))
  
  CHH_hal_data <- data.frame(category = c("CDS",
                                          "exon",
                                          "five_prime_UTR",
                                          "three_prime_UTR"),
                            count = c(sum(overlaps_hal_CHH_final$context=="CDS"),
                                      sum(overlaps_hal_CHH_final$context=="exon"),
                                      sum(overlaps_hal_CHH_final$context=="five_prime_UTR"),
                                      sum(overlaps_hal_CHH_final$context=="three_prime_UTR")))
  
  ### lyrata
  CG_lyr_data <- data.frame(category = c("CDS",
                                         "exon",
                                         "five_prime_UTR",
                                         "three_prime_UTR"),
                            count = c(sum(overlaps_lyr_CG_final$context=="CDS"),
                                      sum(overlaps_lyr_CG_final$context=="exon"),
                                      sum(overlaps_lyr_CG_final$context=="five_prime_UTR"),
                                      sum(overlaps_lyr_CG_final$context=="three_prime_UTR")))
  
  CHG_lyr_data <- data.frame(category = c("CDS",
                                          "exon",
                                          "five_prime_UTR",
                                          "three_prime_UTR"),
                             count = c(sum(overlaps_lyr_CHG_final$context=="CDS"),
                                       sum(overlaps_lyr_CHG_final$context=="exon"),
                                       sum(overlaps_lyr_CHG_final$context=="five_prime_UTR"),
                                       sum(overlaps_lyr_CHG_final$context=="three_prime_UTR")))
  
  CHH_lyr_data <- data.frame(category = c("CDS",
                                          "exon",
                                          "five_prime_UTR",
                                          "three_prime_UTR"),
                             count = c(sum(overlaps_lyr_CHH_final$context=="CDS"),
                                       sum(overlaps_lyr_CHH_final$context=="exon"),
                                       sum(overlaps_lyr_CHH_final$context=="five_prime_UTR"),
                                       sum(overlaps_lyr_CHH_final$context=="three_prime_UTR")))
  
  # Compute expected counts based on annotation
  
  hal_anno <- data.frame(category = c("CDS",
                                      "exon",
                                      "five_prime_UTR",
                                      "three_prime_UTR"),
                            count = c(sum(Ahal_anno$context=="CDS"),
                                      sum(Ahal_anno$context=="exon"),
                                      sum(Ahal_anno$context=="five_prime_UTR"),
                                      sum(Ahal_anno$context=="three_prime_UTR")))
  
  lyr_anno <- data.frame(category = c("CDS",
                                      "exon",
                                      "five_prime_UTR",
                                      "three_prime_UTR"),
                         count = c(sum(Alyr_anno$context=="CDS"),
                                   sum(Alyr_anno$context=="exon"),
                                   sum(Alyr_anno$context=="five_prime_UTR"),
                                   sum(Alyr_anno$context=="three_prime_UTR")))
  
  # Compute proportions for each data frame
  
  proportions_hal_CG <- CG_hal_data$count / sum(CG_hal_data$count)
  proportions_hal_CHG <- CHG_hal_data$count / sum(CHG_hal_data$count)
  proportions_hal_CHH <- CHH_hal_data$count / sum(CHH_hal_data$count)
  proportions_hal_anno <- hal_anno$count / sum(hal_anno$count)
  
  CG_hal_data <- mutate(CG_hal_data, proportions = proportions_hal_CG)
  CHG_hal_data <- mutate(CHG_hal_data, proportions = proportions_hal_CHG)
  CHH_hal_data <- mutate(CHH_hal_data, proportions = proportions_hal_CHH)
  hal_anno <- mutate(hal_anno, proportions = proportions_hal_anno)
  
  proportions_lyr_CG <- CG_lyr_data$count / sum(CG_lyr_data$count)
  proportions_lyr_CHG <- CHG_lyr_data$count / sum(CHG_lyr_data$count)
  proportions_lyr_CHH <- CHH_lyr_data$count / sum(CHH_lyr_data$count)
  proportions_lyr_anno <- lyr_anno$count / sum(lyr_anno$count)
  
  CG_lyr_data <- mutate(CG_lyr_data, proportions = proportions_lyr_CG)
  CHG_lyr_data <- mutate(CHG_lyr_data, proportions = proportions_lyr_CHG)
  CHH_lyr_data <- mutate(CHH_lyr_data, proportions = proportions_lyr_CHH)
  lyr_anno <- mutate(lyr_anno, proportions = proportions_lyr_anno)
  
  # Combine counts from data with expected counts
  
  origin <- c(rep("CG", 7), rep("CHG", 7), rep("CHH", 7), rep("annotation", 7))
  
  hal_data_all <- rbind(CG_hal_data, CHG_hal_data, CHH_hal_data, hal_anno)
  hal_data_all <- mutate(hal_data_all, origin = origin)
  
  lyr_data_all <- rbind(CG_lyr_data, CHG_lyr_data, CHH_lyr_data, lyr_anno)
  lyr_data_all <- mutate(lyr_data_all, origin = origin)
  
  side <- c(rep("halleri", 28), rep("lyrata", 28))
  data_all <- rbind(hal_data_all, lyr_data_all)
  data_all <- mutate(data_all, side = side)
  
  ## Plotting halleri side
  
  p1 <- ggplot(data = data_all, aes(x = category, y = proportions, fill = origin)) +
    facet_wrap(~side) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Gene region") +
    ylab("Proportions") +
    ggtitle(conditions)
  
  return(p1)
}

DMRs_overlap_special <- function(dmrseq_output_CG_AvB,
                                 dmrseq_output_CHG_AvB,
                                 dmrseq_output_CHH_AvB,
                                 conditions){
  scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
  scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
  scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])
  dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
  dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
  dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
  dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
  dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
  dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]
  DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
               hal_CHG_all = dmrseq_output_CHG_hal,
               hal_CHH_all = dmrseq_output_CHH_hal,
               lyr_CG_all = dmrseq_output_CG_lyr,
               lyr_CHG_all = dmrseq_output_CHG_lyr,
               lyr_CHH_all = dmrseq_output_CHH_lyr,
               conditions = conditions)
}

## Import data

setwd("~/OneDrive/PhD/Project/Chapter_3/DMR_results/ARPEGGIO_results_syn1Vpro1_HM/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## plot

HM_syn1_v_pro1 <- DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
                               hal_CHG_all = dmrseq_output_CHG_hal,
                               hal_CHH_all = dmrseq_output_CHH_hal,
                               lyr_CG_all = dmrseq_output_CG_lyr,
                               lyr_CHG_all = dmrseq_output_CHG_lyr,
                               lyr_CHH_all = dmrseq_output_CHH_lyr,
                               conditions = "Cold conditions")

## Import data

setwd("~/OneDrive/PhD/Project/Chapter_3/DMR_results/ARPEGGIO_results_syn4Vpro1_HM/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## plot

HM_syn4_v_pro1 <- DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
             hal_CHG_all = dmrseq_output_CHG_hal,
             hal_CHH_all = dmrseq_output_CHH_hal,
             lyr_CG_all = dmrseq_output_CG_lyr,
             lyr_CHG_all = dmrseq_output_CHG_lyr,
             lyr_CHH_all = dmrseq_output_CHH_lyr,
             conditions = "Cold conditions")

## Import data

setwd("~/OneDrive/PhD/Project/Chapter_3/DMR_results/ARPEGGIO_results_syn1Vpro1_LL/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## plot

LL_syn1_v_pro1 <-DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
             hal_CHG_all = dmrseq_output_CHG_hal,
             hal_CHH_all = dmrseq_output_CHH_hal,
             lyr_CG_all = dmrseq_output_CG_lyr,
             lyr_CHG_all = dmrseq_output_CHG_lyr,
             lyr_CHH_all = dmrseq_output_CHH_lyr,
             conditions = "Hot conditions")

## Import data

setwd("~/OneDrive/PhD/Project/Chapter_3/DMR_results/ARPEGGIO_results_syn4Vpro1_LL/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## plot

LL_syn4_v_pro1 <- DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
                              hal_CHG_all = dmrseq_output_CHG_hal,
                              hal_CHH_all = dmrseq_output_CHH_hal,
                              lyr_CG_all = dmrseq_output_CG_lyr,
                              lyr_CHG_all = dmrseq_output_CHG_lyr,
                              lyr_CHH_all = dmrseq_output_CHH_lyr,
                              conditions = "Hot conditions")


## combines plots

(HM_syn1_v_pro1 | LL_syn1_v_pro1) /
(HM_syn4_v_pro1 | LL_syn4_v_pro1)
