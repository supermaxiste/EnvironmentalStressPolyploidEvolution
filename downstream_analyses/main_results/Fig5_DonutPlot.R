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
                         extract_legend){
  
  ### Import data: annotations with 500bp flanking regions marked as genes and DMRs
  
  ### Paths to annotation files to be modified according to new system!
  
  Ahal_v2_2_500bp_promoters <- read.delim("/path/to/data/Fig5/hal_r2_500bp_MAKER_32518_8chr.gff",
                                        header=FALSE,
                                        col.names = c("scaffold", "tool", 
                                                      "context", "start", 
                                                      "end", "number",
                                                      "strand", "dot", "extra"))
  Alyr_v2_2_renamed_500bp_promoters <- read.delim("/path/to/data/Fig5/lyr_r2_500bp_MAKER_28737_8chr_renamed.gff", 
                                                header=FALSE,
                                                col.names = c("scaffold", "tool", 
                                                              "context", "start", 
                                                              "end", "number",
                                                              "strand", "dot", "extra"))
  
  Ahal_anno <- select(Ahal_v2_2_500bp_promoters, scaffold, context, start, end, strand)
  rm(Ahal_v2_2_500bp_promoters)
  
  
  Alyr_anno <- select(Alyr_v2_2_renamed_500bp_promoters, scaffold, context, start, end, strand)
  rm(Alyr_v2_2_renamed_500bp_promoters)
  
  ### Filter annotations to keep only promoters
  
  Ahal_anno <- filter(Ahal_anno, context=="gene")
  Alyr_anno <- filter(Alyr_anno, context=="gene")
  
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
  
  overlaps_hal_CG <- findOverlaps(hal_CG_allG, Ahal_annoG)
  
  overlaps_hal_CG_final <- cbind(hal_CG_all[overlaps_hal_CG@from,],
                                 Ahal_anno[overlaps_hal_CG@to,])
  
  colnames(overlaps_hal_CG_final) <- c("seqnames", "start", "end", "width", "strand",
                                       "L", "area", "beta", "stat", "pval", "qval", 
                                       "index.start", "index.end", "index.width",
                                       "scaffold", "context", "start2", "end2", "strand2")
  
  # Some DMRs overlap on promoters twice, so we remove duplicates to count them
  
  overlaps_hal_CG_final_dedup <- overlaps_hal_CG_final[!duplicated(overlaps_hal_CG_final$start),]
  
  # We repeat the same procedure for CHG
  
  overlaps_hal_CHG <- findOverlaps(hal_CHG_allG, Ahal_annoG)
  
  overlaps_hal_CHG_final <- cbind(hal_CHG_all[overlaps_hal_CHG@from,],
                                  Ahal_anno[overlaps_hal_CHG@to,])
  
  colnames(overlaps_hal_CHG_final) <- c("seqnames", "start", "end", "width", "strand",
                                        "L", "area", "beta", "stat", "pval", "qval", 
                                        "index.start", "index.end", "index.width",
                                        "scaffold", "context", "start2", "end2", "strand2")
  
  # Some DMRs overlap on promoters twice, so we remove duplicates to count them
  
  overlaps_hal_CHG_final_dedup <- overlaps_hal_CHG_final[!duplicated(overlaps_hal_CHG_final$start),]
  
  # We repeat the same procedure for CHH
  
  overlaps_hal_CHH <- findOverlaps(hal_CHH_allG, Ahal_annoG)
  
  overlaps_hal_CHH_final <- cbind(hal_CHH_all[overlaps_hal_CHH@from,],
                                  Ahal_anno[overlaps_hal_CHH@to,])
  
  colnames(overlaps_hal_CHH_final) <- c("seqnames", "start", "end", "width", "strand",
                                        "L", "area", "beta", "stat", "pval", "qval", 
                                        "index.start", "index.end", "index.width",
                                        "scaffold", "context", "start2", "end2", "strand2")
  
  # Some DMRs overlap on promoters twice, so we remove duplicates to count them
  
  overlaps_hal_CHH_final_dedup <- overlaps_hal_CHH_final[!duplicated(overlaps_hal_CHH_final$start),]
  
  # We repeat the same procedure for CG, lyrata side
  
  overlaps_lyr_CG <- findOverlaps(lyr_CG_allG, Alyr_annoG)
  
  overlaps_lyr_CG_final <- cbind(lyr_CG_all[overlaps_lyr_CG@from,],
                                 Alyr_anno[overlaps_lyr_CG@to,])
  
  colnames(overlaps_lyr_CG_final) <- c("seqnames", "start", "end", "width", "strand",
                                       "L", "area", "beta", "stat", "pval", "qval", 
                                       "index.start", "index.end", "index.width",
                                       "scaffold", "context", "start2", "end2", "strand2")
  
  # Some DMRs overlap on promoters twice, so we remove duplicates to count them
  
  overlaps_lyr_CG_final_dedup <- overlaps_lyr_CG_final[!duplicated(overlaps_lyr_CG_final$start),]
  
  # We repeat the same procedure for CHG
  
  overlaps_lyr_CHG <- findOverlaps(lyr_CHG_allG, Alyr_annoG)
  
  overlaps_lyr_CHG_final <- cbind(lyr_CHG_all[overlaps_lyr_CHG@from,],
                                  Alyr_anno[overlaps_lyr_CHG@to,])
  
  colnames(overlaps_lyr_CHG_final) <- c("seqnames", "start", "end", "width", "strand",
                                        "L", "area", "beta", "stat", "pval", "qval", 
                                        "index.start", "index.end", "index.width",
                                        "scaffold", "context", "start2", "end2", "strand2")
  
  # Some DMRs overlap on promoters twice, so we remove duplicates to count them
  
  overlaps_lyr_CHG_final_dedup <- overlaps_lyr_CHG_final[!duplicated(overlaps_lyr_CHG_final$start),]
  
  # We repeat the same procedure for CHH
  
  overlaps_lyr_CHH <- findOverlaps(lyr_CHH_allG, Alyr_annoG)
  
  overlaps_lyr_CHH_final <- cbind(lyr_CHH_all[overlaps_lyr_CHH@from,],
                                  Alyr_anno[overlaps_lyr_CHH@to,])
  
  colnames(overlaps_lyr_CHH_final) <- c("seqnames", "start", "end", "width", "strand",
                                        "L", "area", "beta", "stat", "pval", "qval", 
                                        "index.start", "index.end", "index.width",
                                        "scaffold", "context", "start2", "end2", "strand2")
  
  # Some DMRs overlap on promoters twice, so we remove duplicates to count them
  
  overlaps_lyr_CHH_final_dedup <- overlaps_lyr_CHH_final[!duplicated(overlaps_lyr_CHH_final$start),]
  
  
  ### Now that we've got all overlaps, there's still the possibility that a DMR overlaps with a gene
  ### To check on DMRs overlapping with both genes and flanking regions, we will check also 
  ### 1bp overlaps of the same DMRs to genes and exclude those (to keep the ones overlapping with promoter only)
  
  ### Import default annotation
  
  Ahal_v2_2 <- read.delim("/path/to/data/Fig5/hal_r2_MAKER_32518_8chr.gff",
                          header=FALSE,
                          col.names = c("scaffold", "tool", 
                                        "context", "start", 
                                        "end", "number",
                                        "strand", "dot", "extra"))
  Ahal_annoTrue <- select(Ahal_v2_2, scaffold, context, start, end, strand)
  
  Ahal_annoTrue_genes <- filter(Ahal_annoTrue, context == "gene")
  
  rm(Ahal_v2_2)
  
  Alyr_v2_2_renamed <- read.delim("/path/to/data/Fig5/lyr_r2_MAKER_28737_8chr_renamed.gff", 
                                  header=FALSE,
                                  col.names = c("scaffold", "tool", 
                                                "context", "start", 
                                                "end", "number",
                                                "strand", "dot", "extra"))
  Alyr_annoTrue <- select(Alyr_v2_2_renamed, scaffold, context, start, end, strand)
  
  Alyr_annoTrue_genes <- filter(Alyr_annoTrue, context == "gene")
  
  rm(Alyr_v2_2_renamed)
  
  # Convert annotation to GRanges 
  
  Ahal_annoTrue_ranges <- IRanges(start = Ahal_annoTrue_genes$start,
                                  end = Ahal_annoTrue_genes$end)
  Ahal_annoTrueG <- GRanges(seqnames = Ahal_annoTrue_genes$scaffold, 
                            ranges = Ahal_annoTrue_ranges,
                            strand = Ahal_annoTrue_genes$strand)
  
  Alyr_annoTrue_ranges <- IRanges(start = Alyr_annoTrue_genes$start,
                                  end = Alyr_annoTrue_genes$end)
  Alyr_annoTrueG <- GRanges(seqnames = Alyr_annoTrue_genes$scaffold, 
                            ranges = Alyr_annoTrue_ranges,
                            strand = Alyr_annoTrue_genes$strand)
  
  ## Find overlaps once again for halleri
  
  overlaps_hal_CG_2 <- findOverlaps(hal_CG_allG, Ahal_annoTrueG)
  
  overlaps_hal_CG_final_2 <- cbind(hal_CG_all[overlaps_hal_CG_2@from,],
                                   Ahal_annoTrue[overlaps_hal_CG_2@to,])
  
  colnames(overlaps_hal_CG_final_2) <- c("seqnames", "start", "end", "width", "strand",
                                         "L", "area", "beta", "stat", "pval", "qval", 
                                         "index.start", "index.end", "index.width",
                                         "scaffold", "context", "start2", "end2", "strand2")
  
  # We add an "ID" column to th dataframe to find duplicated regions
  
  overlaps_hal_CG_final <- mutate(overlaps_hal_CG_final, ID = paste0(overlaps_hal_CG_final$seqnames, overlaps_hal_CG_final$start, overlaps_hal_CG_final$end))
  overlaps_hal_CG_final_2 <- mutate(overlaps_hal_CG_final_2, ID = paste0(overlaps_hal_CG_final_2$seqnames, overlaps_hal_CG_final_2$start, overlaps_hal_CG_final_2$end))
  
  # Removing duplicates from overlaps
  
  overlaps_hal_CG_final_2 <- overlaps_hal_CG_final_2[!duplicated(overlaps_hal_CG_final_2$ID),]
  
  # We check how many regions overlap with both promoters and genes
  
  promo_genes_CG_hal <- match(overlaps_hal_CG_final$ID, overlaps_hal_CG_final_2$ID)
  
  # We remove those matching regions from the dataframe with promoter overlaps
  
  overlaps_hal_CG_final_prom <- overlaps_hal_CG_final[is.na(promo_genes_CG_hal),]
  
  # Some DMRs overlap on genes twice, so we remove duplicates to count them
  
  hal_CG_promoter <- overlaps_hal_CG_final_prom[!duplicated(overlaps_hal_CG_final_prom$ID),]
  
  # We repeat the same procedure for CHG
  
  overlaps_hal_CHG_2 <- findOverlaps(hal_CHG_allG, Ahal_annoTrueG)
  
  overlaps_hal_CHG_final_2 <- cbind(hal_CHG_all[overlaps_hal_CHG_2@from,],
                                    Ahal_annoTrue[overlaps_hal_CHG_2@to,])
  
  colnames(overlaps_hal_CHG_final_2) <- c("seqnames", "start", "end", "width", "strand",
                                          "L", "area", "beta", "stat", "pval", "qval", 
                                          "index.start", "index.end", "index.width",
                                          "scaffold", "context", "start2", "end2", "strand2")
  
  # We add an "ID" column to th dataframe to find duplicated regions
  
  overlaps_hal_CHG_final <- mutate(overlaps_hal_CHG_final, ID = paste0(overlaps_hal_CHG_final$seqnames, overlaps_hal_CHG_final$start, overlaps_hal_CHG_final$end))
  overlaps_hal_CHG_final_2 <- mutate(overlaps_hal_CHG_final_2, ID = paste0(overlaps_hal_CHG_final_2$seqnames, overlaps_hal_CHG_final_2$start, overlaps_hal_CHG_final_2$end))
  
  # Removing duplicates from overlaps
  
  overlaps_hal_CHG_final_2 <- overlaps_hal_CHG_final_2[!duplicated(overlaps_hal_CHG_final_2$ID),]
  
  # We check how many regions overlap with both promoters and genes
  
  promo_genes_CHG_hal <- match(overlaps_hal_CHG_final$ID, overlaps_hal_CHG_final_2$ID)
  
  # We remove those matching regions from the dataframe with promoter overlaps
  
  overlaps_hal_CHG_final_prom <- overlaps_hal_CHG_final[is.na(promo_genes_CHG_hal),]
  
  # Some DMRs overlap on genes twice, so we remove duplicates to count them
  
  hal_CHG_promoter <- overlaps_hal_CHG_final_prom[!duplicated(overlaps_hal_CHG_final_prom$ID),]
  
  # We repeat the same procedure for CHH
  
  overlaps_hal_CHH_2 <- findOverlaps(hal_CHH_allG, Ahal_annoTrueG)
  
  overlaps_hal_CHH_final_2 <- cbind(hal_CHH_all[overlaps_hal_CHH_2@from,],
                                    Ahal_annoTrue[overlaps_hal_CHH_2@to,])
  
  colnames(overlaps_hal_CHH_final_2) <- c("seqnames", "start", "end", "width", "strand",
                                          "L", "area", "beta", "stat", "pval", "qval", 
                                          "index.start", "index.end", "index.width",
                                          "scaffold", "context", "start2", "end2", "strand2")
  
  # We add an "ID" column to th dataframe to find duplicated regions
  
  overlaps_hal_CHH_final <- mutate(overlaps_hal_CHH_final, ID = paste0(overlaps_hal_CHH_final$seqnames, overlaps_hal_CHH_final$start, overlaps_hal_CHH_final$end))
  overlaps_hal_CHH_final_2 <- mutate(overlaps_hal_CHH_final_2, ID = paste0(overlaps_hal_CHH_final_2$seqnames, overlaps_hal_CHH_final_2$start, overlaps_hal_CHH_final_2$end))
  
  # Removing duplicates from overlaps
  
  overlaps_hal_CHH_final_2 <- overlaps_hal_CHH_final_2[!duplicated(overlaps_hal_CHH_final_2$ID),]
  
  # We check how many regions overlap with both promoters and genes
  
  promo_genes_CHH_hal <- match(overlaps_hal_CHH_final$ID, overlaps_hal_CHH_final_2$ID)
  
  # We remove those matching regions from the dataframe with promoter overlaps
  
  overlaps_hal_CHH_final_prom <- overlaps_hal_CHH_final[is.na(promo_genes_CHH_hal),]
  
  # Some DMRs overlap on genes twice, so we remove duplicates to count them
  
  hal_CHH_promoter <- overlaps_hal_CHH_final_prom[!duplicated(overlaps_hal_CHH_final_prom$ID),]
  
  ## Find overlaps once again for lyrata
  
  overlaps_lyr_CG_2 <- findOverlaps(lyr_CG_allG, Alyr_annoTrueG)
  
  overlaps_lyr_CG_final_2 <- cbind(lyr_CG_all[overlaps_lyr_CG_2@from,],
                                   Alyr_annoTrue[overlaps_lyr_CG_2@to,])
  
  colnames(overlaps_lyr_CG_final_2) <- c("seqnames", "start", "end", "width", "strand",
                                         "L", "area", "beta", "stat", "pval", "qval", 
                                         "index.start", "index.end", "index.width",
                                         "scaffold", "context", "start2", "end2", "strand2")
  
  # We add an "ID" column to th dataframe to find duplicated regions
  
  overlaps_lyr_CG_final <- mutate(overlaps_lyr_CG_final, ID = paste0(overlaps_lyr_CG_final$seqnames, overlaps_lyr_CG_final$start, overlaps_lyr_CG_final$end))
  overlaps_lyr_CG_final_2 <- mutate(overlaps_lyr_CG_final_2, ID = paste0(overlaps_lyr_CG_final_2$seqnames, overlaps_lyr_CG_final_2$start, overlaps_lyr_CG_final_2$end))
  
  # Removing duplicates from overlaps
  
  overlaps_lyr_CG_final_2 <- overlaps_lyr_CG_final_2[!duplicated(overlaps_lyr_CG_final_2$ID),]
  
  # We check how many regions overlap with both promoters and genes
  
  promo_genes_CG_lyr <- match(overlaps_lyr_CG_final$ID, overlaps_lyr_CG_final_2$ID)
  
  # We remove those matching regions from the dataframe with promoter overlaps
  
  overlaps_lyr_CG_final_prom <- overlaps_lyr_CG_final[is.na(promo_genes_CG_lyr),]
  
  # Some DMRs overlap on genes twice, so we remove duplicates to count them
  
  lyr_CG_promoter <- overlaps_lyr_CG_final_prom[!duplicated(overlaps_lyr_CG_final_prom$ID),]
  
  # We repeat the same procedure for CHG
  
  overlaps_lyr_CHG_2 <- findOverlaps(lyr_CHG_allG, Alyr_annoTrueG)
  
  overlaps_lyr_CHG_final_2 <- cbind(lyr_CHG_all[overlaps_lyr_CHG_2@from,],
                                    Alyr_annoTrue[overlaps_lyr_CHG_2@to,])
  
  colnames(overlaps_lyr_CHG_final_2) <- c("seqnames", "start", "end", "width", "strand",
                                          "L", "area", "beta", "stat", "pval", "qval", 
                                          "index.start", "index.end", "index.width",
                                          "scaffold", "context", "start2", "end2", "strand2")
  
  # We add an "ID" column to th dataframe to find duplicated regions
  
  overlaps_lyr_CHG_final <- mutate(overlaps_lyr_CHG_final, ID = paste0(overlaps_lyr_CHG_final$seqnames, overlaps_lyr_CHG_final$start, overlaps_lyr_CHG_final$end))
  overlaps_lyr_CHG_final_2 <- mutate(overlaps_lyr_CHG_final_2, ID = paste0(overlaps_lyr_CHG_final_2$seqnames, overlaps_lyr_CHG_final_2$start, overlaps_lyr_CHG_final_2$end))
  
  # Removing duplicates from overlaps
  
  overlaps_lyr_CHG_final_2 <- overlaps_lyr_CHG_final_2[!duplicated(overlaps_lyr_CHG_final_2$ID),]
  
  # We check how many regions overlap with both promoters and genes
  
  promo_genes_CHG_lyr <- match(overlaps_lyr_CHG_final$ID, overlaps_lyr_CHG_final_2$ID)
  
  # We remove those matching regions from the dataframe with promoter overlaps
  
  overlaps_lyr_CHG_final_prom <- overlaps_lyr_CHG_final[is.na(promo_genes_CHG_lyr),]
  
  # Some DMRs overlap on genes twice, so we remove duplicates to count them
  
  lyr_CHG_promoter <- overlaps_lyr_CHG_final_prom[!duplicated(overlaps_lyr_CHG_final_prom$ID),]
  
  # We repeat the same procedure for CHH
  
  overlaps_lyr_CHH_2 <- findOverlaps(lyr_CHH_allG, Alyr_annoTrueG)
  
  overlaps_lyr_CHH_final_2 <- cbind(lyr_CHH_all[overlaps_lyr_CHH_2@from,],
                                    Alyr_annoTrue[overlaps_lyr_CHH_2@to,])
  
  colnames(overlaps_lyr_CHH_final_2) <- c("seqnames", "start", "end", "width", "strand",
                                          "L", "area", "beta", "stat", "pval", "qval", 
                                          "index.start", "index.end", "index.width",
                                          "scaffold", "context", "start2", "end2", "strand2")
  
  # We add an "ID" column to th dataframe to find duplicated regions
  
  overlaps_lyr_CHH_final <- mutate(overlaps_lyr_CHH_final, ID = paste0(overlaps_lyr_CHH_final$seqnames, overlaps_lyr_CHH_final$start, overlaps_lyr_CHH_final$end))
  overlaps_lyr_CHH_final_2 <- mutate(overlaps_lyr_CHH_final_2, ID = paste0(overlaps_lyr_CHH_final_2$seqnames, overlaps_lyr_CHH_final_2$start, overlaps_lyr_CHH_final_2$end))
  
  # Removing duplicates from overlaps
  
  overlaps_lyr_CHH_final_2 <- overlaps_lyr_CHH_final_2[!duplicated(overlaps_lyr_CHH_final_2$ID),]
  
  # We check how many regions overlap with both promoters and genes
  
  promo_genes_CHH_lyr <- match(overlaps_lyr_CHH_final$ID, overlaps_lyr_CHH_final_2$ID)
  
  # We remove those matching regions from the dataframe with promoter overlaps
  
  overlaps_lyr_CHH_final_prom <- overlaps_lyr_CHH_final[is.na(promo_genes_CHH_lyr),]
  
  # Some DMRs overlap on genes twice, so we remove duplicates to count them
  
  lyr_CHH_promoter <- overlaps_lyr_CHH_final_prom[!duplicated(overlaps_lyr_CHH_final_prom$ID),]
  
  
  ### Plotting proportions of DMRs in gene bodies, flanking regions and genic regions for each methylation context
  
  CG_hal_data <- data.frame(category = c("gene body", "flanking region", "intergenic region"),
                            count = c(nrow(overlaps_hal_CG_final_2),
                                      nrow(hal_CG_promoter), 
                                      nrow(hal_CG_all) - nrow(overlaps_hal_CG_final_2) - nrow(hal_CG_promoter)))
  
  CHG_hal_data <- data.frame(category = c("gene body", "flanking region", "intergenic region"),
                             count = c(nrow(overlaps_hal_CHG_final_2),
                                       nrow(hal_CHG_promoter), 
                                       nrow(hal_CHG_all) - nrow(overlaps_hal_CHG_final_2) - nrow(hal_CHG_promoter)))
  
  CHH_hal_data <- data.frame(category = c("gene body", "flanking region", "intergenic region"),
                             count = c(nrow(overlaps_hal_CHH_final_2),
                                       nrow(hal_CHH_promoter), 
                                       nrow(hal_CHH_all) - nrow(overlaps_hal_CHH_final_2) - nrow(hal_CHH_promoter)))
  
  CG_lyr_data <- data.frame(category = c("gene body", "flanking region", "intergenic region"),
                            count = c(nrow(overlaps_lyr_CG_final_2),
                                      nrow(lyr_CG_promoter), 
                                      nrow(lyr_CG_all) - nrow(overlaps_lyr_CG_final_2) - nrow(lyr_CG_promoter)))
  
  CHG_lyr_data <- data.frame(category = c("gene body", "flanking region", "intergenic region"),
                             count = c(nrow(overlaps_lyr_CHG_final_2),
                                       nrow(lyr_CHG_promoter), 
                                       nrow(lyr_CHG_all) - nrow(overlaps_lyr_CHG_final_2) - nrow(lyr_CHG_promoter)))
  
  CHH_lyr_data <- data.frame(category = c("gene body", "flanking region", "intergenic region"),
                             count = c(nrow(overlaps_lyr_CHH_final_2),
                                       nrow(lyr_CHH_promoter), 
                                       nrow(lyr_CHH_all) - nrow(overlaps_lyr_CHH_final_2) - nrow(lyr_CHH_promoter)))
  
  ## Plotting halleri side
  
  # Compute percentages
  CG_hal_data$fraction <- CG_hal_data$count / sum(CG_hal_data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  CG_hal_data$ymax <- cumsum(CG_hal_data$fraction)
  
  # Compute the bottom of each rectangle
  CG_hal_data$ymin <- c(0, head(CG_hal_data$ymax, n=-1))
  
  # Compute label position
  CG_hal_data$labelPosition <- (CG_hal_data$ymax + CG_hal_data$ymin) / 2
  
  # Compute a good label
  CG_hal_data$label <- paste0(CG_hal_data$count, "\n", round(100*CG_hal_data$fraction, 1), "%")
  
  # Make the plot
  p1 <- ggplot(CG_hal_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=3, show.legend = FALSE) +
    annotate(geom = 'text', x = 2, y = 0, label = "CG") +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") #,
    #legend.text = element_text(size = 16)) +
    #ggtitle("CG")
  
  # Compute percentages
  CHG_hal_data$fraction <- CHG_hal_data$count / sum(CHG_hal_data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  CHG_hal_data$ymax <- cumsum(CHG_hal_data$fraction)
  
  # Compute the bottom of each rectangle
  CHG_hal_data$ymin <- c(0, head(CHG_hal_data$ymax, n=-1))
  
  # Compute label position
  CHG_hal_data$labelPosition <- (CHG_hal_data$ymax + CHG_hal_data$ymin) / 2
  
  # Compute a good label
  CHG_hal_data$label <- paste0(CHG_hal_data$count, "\n", round(100*CHG_hal_data$fraction, 1), "%")
  
  # Make the plot
  p2 <- ggplot(CHG_hal_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=3.5, aes(y=labelPosition, label=label), size=3, show.legend = FALSE) +
    annotate(geom = 'text', x = 2, y = 0, label = "CHG") +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") #,
    #legend.text = element_text(size = 16)) +
    #ggtitle("CHG")
  
  # Compute percentages
  CHH_hal_data$fraction <- CHH_hal_data$count / sum(CHH_hal_data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  CHH_hal_data$ymax <- cumsum(CHH_hal_data$fraction)
  
  # Compute the bottom of each rectangle
  CHH_hal_data$ymin <- c(0, head(CHH_hal_data$ymax, n=-1))
  
  # Compute label position
  CHH_hal_data$labelPosition <- (CHH_hal_data$ymax + CHH_hal_data$ymin) / 2
  
  # Compute a good label
  CHH_hal_data$label <- paste0(CHH_hal_data$count, "\n", round(100*CHH_hal_data$fraction, 1), "%")
  
  # Make the plot
  p3 <- ggplot(CHH_hal_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=3.5, aes(y=labelPosition, label=label), size=3, show.legend = FALSE) +
    annotate(geom = 'text', x = 2, y = 0, label = "CHH") +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") #,
    #legend.text = element_text(size = 16)) +
    #ggtitle("CHH")
    
  # Create general legend
  p1_legend <- ggplot(CHH_hal_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    theme(legend.position = "bottom")
  
  # Create user-defined function, which extracts legends from ggplots
  extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  
  # Apply user-defined function to extract legend
  shared_legend <- extract_legend(p1_legend)
  
  # p_toth <- grid.arrange(arrangeGrob(p1, p2, p3, ncol = 3),
  #              shared_legend, nrow = 2, heights = c(1, 0.1),
  #              top=textGrob("halleri-side", gp=gpar(fontsize=15,font=8)))
  
  ## Plot lyrata side
  
  # Compute percentages
  CG_lyr_data$fraction <- CG_lyr_data$count / sum(CG_lyr_data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  CG_lyr_data$ymax <- cumsum(CG_lyr_data$fraction)
  
  # Compute the bottom of each rectangle
  CG_lyr_data$ymin <- c(0, head(CG_lyr_data$ymax, n=-1))
  
  # Compute label position
  CG_lyr_data$labelPosition <- (CG_lyr_data$ymax + CG_lyr_data$ymin) / 2
  
  # Compute a good label
  CG_lyr_data$label <- paste0(CG_lyr_data$count, "\n", round(100*CG_lyr_data$fraction, 1), "%")
  
  # Make the plot
  p1l <- ggplot(CG_lyr_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=3.5, aes(y=labelPosition, label=label), size=3, show.legend = FALSE) +
    annotate(geom = 'text', x = 2, y = 0, label = "CG") +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") #,
    #legend.text = element_text(size = 16)) +
    #ggtitle("CG")
  
  # Compute percentages
  CHG_lyr_data$fraction <- CHG_lyr_data$count / sum(CHG_lyr_data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  CHG_lyr_data$ymax <- cumsum(CHG_lyr_data$fraction)
  
  # Compute the bottom of each rectangle
  CHG_lyr_data$ymin <- c(0, head(CHG_lyr_data$ymax, n=-1))
  
  # Compute label position
  CHG_lyr_data$labelPosition <- (CHG_lyr_data$ymax + CHG_lyr_data$ymin) / 2
  
  # Compute a good label
  CHG_lyr_data$label <- paste0(CHG_lyr_data$count, "\n", round(100*CHG_lyr_data$fraction, 1), "%")
  
  # Make the plot
  p2l <- ggplot(CHG_lyr_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=3.5, aes(y=labelPosition, label=label), size=3, show.legend = FALSE) +
    annotate(geom = 'text', x = 2, y = 0, label = "CHG") +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") #,
    #legend.text = element_text(size = 16)) +
    #ggtitle("CHG")
  
  # Compute percentages
  CHH_lyr_data$fraction <- CHH_lyr_data$count / sum(CHH_lyr_data$count)
  
  # Compute the cumulative percentages (top of each rectangle)
  CHH_lyr_data$ymax <- cumsum(CHH_lyr_data$fraction)
  
  # Compute the bottom of each rectangle
  CHH_lyr_data$ymin <- c(0, head(CHH_lyr_data$ymax, n=-1))
  
  # Compute label position
  CHH_lyr_data$labelPosition <- (CHH_lyr_data$ymax + CHH_lyr_data$ymin) / 2
  
  # Compute a good label
  CHH_lyr_data$label <- paste0(CHH_lyr_data$count, "\n", round(100*CHH_lyr_data$fraction, 1), "%")
  
  # Make the plot
  p3l <- ggplot(CHH_lyr_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=3.5, aes(y=labelPosition, label=label), size=3, show.legend = FALSE) +
    annotate(geom = 'text', x = 2, y = 0, label = "CHH") +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
    theme_void() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "none") #,
    #legend.text = element_text(size = 16)) +
    #ggtitle("CHH")
  
  # Create general legend
  
  p1_legend <- ggplot(CHH_lyr_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    theme(legend.position = "bottom")
  
  # Create user-defined function, which extracts legends from ggplots
  extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  
  # Apply user-defined function to extract legend
  shared_legend <- extract_legend(p1_legend)
  
  #p_totl <- grid.arrange(arrangeGrob(p1l, p2l, p3l, ncol = 3),
  #             shared_legend, nrow = 2, heights = c(1, 0.1),
  #             top=textGrob("lyrata-side", gp=gpar(fontsize=15,font=8), vjust=0.5))
  
  all <- p1 + p2 + p3 + p1l + p2l + p3l + plot_layout(ncol = 6)
  
  if(isTRUE(extract_legend)){
    return(shared_legend)
  } else{
  return(all)
  }
}

DMRs_overlap_special <- function(dmrseq_output_CG_AvB,
                                 dmrseq_output_CHG_AvB,
                                 dmrseq_output_CHH_AvB,
                                 extract_legend){
  scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "r", 2)[,2])
  scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "r", 2)[,2])
  scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "r", 2)[,2])
  dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 9]
  dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 8]
  dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 9]
  dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 8]
  dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 9]
  dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 8]
  DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
               hal_CHG_all = dmrseq_output_CHG_hal,
               hal_CHH_all = dmrseq_output_CHH_hal,
               lyr_CG_all = dmrseq_output_CG_lyr,
               lyr_CHG_all = dmrseq_output_CHG_lyr,
               lyr_CHH_all = dmrseq_output_CHH_lyr,
               extract_legend = extract_legend)
}

## Import data

setwd("/path/to/data/Fig6a/MILD_syn1_v_pro1/dmrseq")
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
             extract_legend = F)

## Import data

setwd("/path/to/data/Fig6a/MILD_syn4_v_pro1/dmrseq")
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
             extract_legend = F)

legendAll <- DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
                          hal_CHG_all = dmrseq_output_CHG_hal,
                          hal_CHH_all = dmrseq_output_CHH_hal,
                          lyr_CG_all = dmrseq_output_CG_lyr,
                          lyr_CHG_all = dmrseq_output_CHG_lyr,
                          lyr_CHH_all = dmrseq_output_CHH_lyr,
                          extract_legend = T)

## All plots together

(plot_spacer() + plot_spacer() + plot_spacer() + wrap_elements(grid::textGrob('halleri-side')) + plot_spacer() + wrap_elements(grid::textGrob('lyrata-side')) + plot_spacer() + plot_spacer() + plot_layout(ncol = 8, heights = 0.5)) /
(wrap_elements(grid::textGrob('SYN1 vs PRO1')) + HM_syn1_v_pro1 + plot_layout(widths = c(0.2, 1))) / 
  (wrap_elements(grid::textGrob('SYN4 vs PRO1')) + HM_syn4_v_pro1 + plot_layout(widths = c(0.2, 1)))


## We repeat the same plot for LL conditions

## Import data

setwd("/path/to/data/Fig6a/STRESS_syn1_v_pro1/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

## plot

LL_syn1_v_pro1 <- DMRs_overlap(hal_CG_all = dmrseq_output_CG_hal,
                               hal_CHG_all = dmrseq_output_CHG_hal,
                               hal_CHH_all = dmrseq_output_CHH_hal,
                               lyr_CG_all = dmrseq_output_CG_lyr,
                               lyr_CHG_all = dmrseq_output_CHG_lyr,
                               lyr_CHH_all = dmrseq_output_CHH_lyr,
                               extract_legend = F)

## Import data

setwd("/path/to/data/Fig6a/STRESS_syn4_v_pro1/dmrseq")
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
                               extract_legend = F)


## All plots together

(plot_spacer() + plot_spacer() + plot_spacer() + wrap_elements(grid::textGrob('halleri-side')) + plot_spacer() + wrap_elements(grid::textGrob('lyrata-side')) + plot_spacer() + plot_spacer() + plot_layout(ncol = 8, heights = 0.5)) /
  (wrap_elements(grid::textGrob('SYN1 vs PRO1')) + LL_syn1_v_pro1 + plot_layout(widths = c(0.2, 1))) / 
  (wrap_elements(grid::textGrob('SYN4 vs PRO1')) + LL_syn4_v_pro1 + plot_layout(widths = c(0.2, 1)))
