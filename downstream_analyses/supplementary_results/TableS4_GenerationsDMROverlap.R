### Check overlap between DMRs in different generations

### Import libraries

library(GenomicRanges)
library(tidyverse)
library(gridExtra)
library(data.table)


### Import files for all contexts

data <- "~/Library/CloudStorage/OneDrive-Personal/PostDoc/Chapter2_3_manuscript/DMR_output/"

HM_G1_hal_CG_all <- read.csv(paste0(data, "MILD_syn1_v_pro1/dmrseq/CG_context/parent1_v_allo.txt"))
HM_G1_hal_CHG_all <- read.csv(paste0(data, "MILD_syn1_v_pro1/dmrseq/CHG_context/parent1_v_allo.txt"))
HM_G1_hal_CHH_all <- read.csv(paste0(data, "MILD_syn1_v_pro1/dmrseq/CHH_context/parent1_v_allo.txt"))

HM_G1_lyr_CG_all <- read.csv(paste0(data, "MILD_syn1_v_pro1/dmrseq/CG_context/parent2_v_allo.txt"))
HM_G1_lyr_CHG_all <- read.csv(paste0(data, "MILD_syn1_v_pro1/dmrseq/CHG_context/parent2_v_allo.txt"))
HM_G1_lyr_CHH_all <- read.csv(paste0(data, "MILD_syn1_v_pro1/dmrseq/CHH_context/parent2_v_allo.txt"))

HM_G4_hal_CG_all <- read.csv(paste0(data, "MILD_syn4_v_pro1/dmrseq/CG_context/parent1_v_allo.txt"))
HM_G4_hal_CHG_all <- read.csv(paste0(data, "MILD_syn4_v_pro1/dmrseq/CHG_context/parent1_v_allo.txt"))
HM_G4_hal_CHH_all <- read.csv(paste0(data, "MILD_syn4_v_pro1/dmrseq/CHH_context/parent1_v_allo.txt"))

HM_G4_lyr_CG_all <- read.csv(paste0(data, "MILD_syn4_v_pro1/dmrseq/CG_context/parent2_v_allo.txt"))
HM_G4_lyr_CHG_all <- read.csv(paste0(data, "MILD_syn4_v_pro1/dmrseq/CHG_context/parent2_v_allo.txt"))
HM_G4_lyr_CHH_all <- read.csv(paste0(data, "MILD_syn4_v_pro1/dmrseq/CHH_context/parent2_v_allo.txt"))


LL_G1_hal_CG_all <- read.csv(paste0(data, "STRESS_syn1_v_pro1/dmrseq/CG_context/parent1_v_allo.txt"))
LL_G1_hal_CHG_all <- read.csv(paste0(data, "STRESS_syn1_v_pro1/dmrseq/CHG_context/parent1_v_allo.txt"))
LL_G1_hal_CHH_all <- read.csv(paste0(data, "STRESS_syn1_v_pro1/dmrseq/CHH_context/parent1_v_allo.txt"))

LL_G1_lyr_CG_all <- read.csv(paste0(data, "STRESS_syn1_v_pro1/dmrseq/CG_context/parent2_v_allo.txt"))
LL_G1_lyr_CHG_all <- read.csv(paste0(data, "STRESS_syn1_v_pro1/dmrseq/CHG_context/parent2_v_allo.txt"))
LL_G1_lyr_CHH_all <- read.csv(paste0(data, "STRESS_syn1_v_pro1/dmrseq/CHH_context/parent2_v_allo.txt"))

LL_G4_hal_CG_all <- read.csv(paste0(data, "STRESS_syn4_v_pro1/dmrseq/CG_context/parent1_v_allo.txt"))
LL_G4_hal_CHG_all <- read.csv(paste0(data, "STRESS_syn4_v_pro1/dmrseq/CHG_context/parent1_v_allo.txt"))
LL_G4_hal_CHH_all <- read.csv(paste0(data, "STRESS_syn4_v_pro1/dmrseq/CHH_context/parent1_v_allo.txt"))

LL_G4_lyr_CG_all <- read.csv(paste0(data, "STRESS_syn4_v_pro1/dmrseq/CG_context/parent2_v_allo.txt"))
LL_G4_lyr_CHG_all <- read.csv(paste0(data, "STRESS_syn4_v_pro1/dmrseq/CHG_context/parent2_v_allo.txt"))
LL_G4_lyr_CHH_all <- read.csv(paste0(data, "STRESS_syn4_v_pro1/dmrseq/CHH_context/parent2_v_allo.txt"))

### Filter DMRs to keep only significant ones

HM_G1_hal_CG_all <- filter(HM_G1_hal_CG_all, qval < 0.05)
HM_G1_hal_CHG_all <- filter(HM_G1_hal_CHG_all, qval < 0.05)
HM_G1_hal_CHH_all <- filter(HM_G1_hal_CHH_all, qval < 0.05)

HM_G1_lyr_CG_all <- filter(HM_G1_lyr_CG_all, qval < 0.05)
HM_G1_lyr_CHG_all <- filter(HM_G1_lyr_CHG_all, qval < 0.05)
HM_G1_lyr_CHH_all <- filter(HM_G1_lyr_CHH_all, qval < 0.05)

HM_G4_hal_CG_all <- filter(HM_G4_hal_CG_all, qval < 0.05)
HM_G4_hal_CHG_all <- filter(HM_G4_hal_CHG_all, qval < 0.05)
HM_G4_hal_CHH_all <- filter(HM_G4_hal_CHH_all, qval < 0.05)

HM_G4_lyr_CG_all <- filter(HM_G4_lyr_CG_all, qval < 0.05)
HM_G4_lyr_CHG_all <- filter(HM_G4_lyr_CHG_all, qval < 0.05)
HM_G4_lyr_CHH_all <- filter(HM_G4_lyr_CHH_all, qval < 0.05)


LL_G1_hal_CG_all <- filter(LL_G1_hal_CG_all, qval < 0.05)
LL_G1_hal_CHG_all <- filter(LL_G1_hal_CHG_all, qval < 0.05)
LL_G1_hal_CHH_all <- filter(LL_G1_hal_CHH_all, qval < 0.05)

LL_G1_lyr_CG_all <- filter(LL_G1_lyr_CG_all, qval < 0.05)
LL_G1_lyr_CHG_all <- filter(LL_G1_lyr_CHG_all, qval < 0.05)
LL_G1_lyr_CHH_all <- filter(LL_G1_lyr_CHH_all, qval < 0.05)

LL_G4_hal_CG_all <- filter(LL_G4_hal_CG_all, qval < 0.05)
LL_G4_hal_CHG_all <- filter(LL_G4_hal_CHG_all, qval < 0.05)
LL_G4_hal_CHH_all <- filter(LL_G4_hal_CHH_all, qval < 0.05)

LL_G4_lyr_CG_all <- filter(LL_G4_lyr_CG_all, qval < 0.05)
LL_G4_lyr_CHG_all <- filter(LL_G4_lyr_CHG_all, qval < 0.05)
LL_G4_lyr_CHH_all <- filter(LL_G4_lyr_CHH_all, qval < 0.05)


### Start from CG context, halleri side

# Convert DMRs to GRanges objects 

HM_G1_hal_CG_all_ranges <- IRanges(start = HM_G1_hal_CG_all$start,
                             end = HM_G1_hal_CG_all$end)
HM_G1_hal_CG_allG <- GRanges(seqnames = HM_G1_hal_CG_all$seqnames,
                       ranges = HM_G1_hal_CG_all_ranges)

HM_G4_hal_CG_all_ranges <- IRanges(start = HM_G4_hal_CG_all$start,
                                end = HM_G4_hal_CG_all$end)
HM_G4_hal_CG_allG <- GRanges(seqnames = HM_G4_hal_CG_all$seqnames,
                          ranges = HM_G4_hal_CG_all_ranges)


# Find overlaps

HM_overlaps_hal_CG <- findOverlaps(HM_G1_hal_CG_allG, HM_G4_hal_CG_allG)

length(unique(HM_overlaps_hal_CG@from))

length(unique(HM_overlaps_hal_CG@from)) / nrow(HM_G1_hal_CG_all)

# 1352 our of 1962 DMRs are conserved in the 4th gen from the 1st (~68.9%)

# Convert DMRs to GRanges objects 

LL_G1_hal_CG_all_ranges <- IRanges(start = LL_G1_hal_CG_all$start,
                                   end = LL_G1_hal_CG_all$end)
LL_G1_hal_CG_allG <- GRanges(seqnames = LL_G1_hal_CG_all$seqnames,
                             ranges = LL_G1_hal_CG_all_ranges)

LL_G4_hal_CG_all_ranges <- IRanges(start = LL_G4_hal_CG_all$start,
                                   end = LL_G4_hal_CG_all$end)
LL_G4_hal_CG_allG <- GRanges(seqnames = LL_G4_hal_CG_all$seqnames,
                             ranges = LL_G4_hal_CG_all_ranges)

# Find overlaps

LL_overlaps_hal_CG <- findOverlaps(LL_G1_hal_CG_allG, LL_G4_hal_CG_allG)

length(unique(LL_overlaps_hal_CG@from))

length(unique(LL_overlaps_hal_CG@from)) / nrow(LL_G1_hal_CG_all)

# 1544 our of 2518 DMRs are conserved in the 4th gen from the 1st (~61.3%)
# for hot conditions

### CHG context, halleri side

# Convert DMRs to GRanges objects 

HM_G1_hal_CHG_all_ranges <- IRanges(start = HM_G1_hal_CHG_all$start,
                                end = HM_G1_hal_CHG_all$end)
HM_G1_hal_CHG_allG <- GRanges(seqnames = HM_G1_hal_CHG_all$seqnames,
                          ranges = HM_G1_hal_CHG_all_ranges)

HM_G4_hal_CHG_all_ranges <- IRanges(start = HM_G4_hal_CHG_all$start,
                                end = HM_G4_hal_CHG_all$end)
HM_G4_hal_CHG_allG <- GRanges(seqnames = HM_G4_hal_CHG_all$seqnames,
                          ranges = HM_G4_hal_CHG_all_ranges)

# Find overlaps

HM_overlaps_hal_CHG <- findOverlaps(HM_G1_hal_CHG_allG, HM_G4_hal_CHG_allG)

length(unique(HM_overlaps_hal_CHG@from))

length(unique(HM_overlaps_hal_CHG@from)) / nrow(HM_G1_hal_CHG_all)

# 1413 our of 2439 DMRs are conserved in the 4th gen from the 1st (~57.9%)

### CHG context, halleri side

# Convert DMRs to GRanges objects 

LL_G1_hal_CHG_all_ranges <- IRanges(start = LL_G1_hal_CHG_all$start,
                                    end = LL_G1_hal_CHG_all$end)
LL_G1_hal_CHG_allG <- GRanges(seqnames = LL_G1_hal_CHG_all$seqnames,
                              ranges = LL_G1_hal_CHG_all_ranges)

LL_G4_hal_CHG_all_ranges <- IRanges(start = LL_G4_hal_CHG_all$start,
                                    end = LL_G4_hal_CHG_all$end)
LL_G4_hal_CHG_allG <- GRanges(seqnames = LL_G4_hal_CHG_all$seqnames,
                              ranges = LL_G4_hal_CHG_all_ranges)

# Find overlaps

LL_overlaps_hal_CHG <- findOverlaps(LL_G1_hal_CHG_allG, LL_G4_hal_CHG_allG)

length(unique(LL_overlaps_hal_CHG@from))

length(unique(LL_overlaps_hal_CHG@from)) / nrow(LL_G1_hal_CHG_all)

# 3757 our of 8057 DMRs are conserved in the 4th gen from the 1st (~46.6 %)
# in hot conditions


### CHH context, halleri side

# Convert DMRs to GRanges objects 

HM_G1_hal_CHH_all_ranges <- IRanges(start = HM_G1_hal_CHH_all$start,
                                 end = HM_G1_hal_CHH_all$end)
HM_G1_hal_CHH_allG <- GRanges(seqnames = HM_G1_hal_CHH_all$seqnames,
                           ranges = HM_G1_hal_CHH_all_ranges)

HM_G4_hal_CHH_all_ranges <- IRanges(start = HM_G4_hal_CHH_all$start,
                                 end = HM_G4_hal_CHH_all$end)
HM_G4_hal_CHH_allG <- GRanges(seqnames = HM_G4_hal_CHH_all$seqnames,
                           ranges = HM_G4_hal_CHH_all_ranges)

# Find overlaps

HM_overlaps_hal_CHH <- findOverlaps(HM_G1_hal_CHH_allG, HM_G4_hal_CHH_allG)

length(unique(HM_overlaps_hal_CHH@from))

length(unique(HM_overlaps_hal_CHH@from)) / nrow(HM_G1_hal_CHH_all)

# 471 our of 1282 DMRs are conserved in the 4th gen from the 1st (~36.7 %)

# Convert DMRs to GRanges objects 

LL_G1_hal_CHH_all_ranges <- IRanges(start = LL_G1_hal_CHH_all$start,
                                    end = LL_G1_hal_CHH_all$end)
LL_G1_hal_CHH_allG <- GRanges(seqnames = LL_G1_hal_CHH_all$seqnames,
                              ranges = LL_G1_hal_CHH_all_ranges)

LL_G4_hal_CHH_all_ranges <- IRanges(start = LL_G4_hal_CHH_all$start,
                                    end = LL_G4_hal_CHH_all$end)
LL_G4_hal_CHH_allG <- GRanges(seqnames = LL_G4_hal_CHH_all$seqnames,
                              ranges = LL_G4_hal_CHH_all_ranges)

# Find overlaps

LL_overlaps_hal_CHH <- findOverlaps(LL_G1_hal_CHH_allG, LL_G4_hal_CHH_allG)

length(unique(LL_overlaps_hal_CHH@from))

length(unique(LL_overlaps_hal_CHH@from)) / nrow(LL_G1_hal_CHH_all)

# 1797 our of 2775 DMRs are conserved in the 4th gen from the 1st (~64.8 %)
# in hot conditions

### Start from CG context, lyrata side now

# Convert DMRs to GRanges objects 

HM_G1_lyr_CG_all_ranges <- IRanges(start = HM_G1_lyr_CG_all$start,
                                end = HM_G1_lyr_CG_all$end)
HM_G1_lyr_CG_allG <- GRanges(seqnames = HM_G1_lyr_CG_all$seqnames,
                          ranges = HM_G1_lyr_CG_all_ranges)

HM_G4_lyr_CG_all_ranges <- IRanges(start = HM_G4_lyr_CG_all$start,
                                end = HM_G4_lyr_CG_all$end)
HM_G4_lyr_CG_allG <- GRanges(seqnames = HM_G4_lyr_CG_all$seqnames,
                          ranges = HM_G4_lyr_CG_all_ranges)

# Find overlaps

HM_overlaps_lyr_CG <- findOverlaps(HM_G1_lyr_CG_allG, HM_G4_lyr_CG_allG)

length(unique(HM_overlaps_lyr_CG@from))

length(unique(HM_overlaps_lyr_CG@from)) / nrow(HM_G1_lyr_CG_all)

# 2025 our of 2724 DMRs are conserved in the 4th gen from the 1st (~74.3 %)

# Convert DMRs to GRanges objects 

LL_G1_lyr_CG_all_ranges <- IRanges(start = LL_G1_lyr_CG_all$start,
                                   end = LL_G1_lyr_CG_all$end)
LL_G1_lyr_CG_allG <- GRanges(seqnames = LL_G1_lyr_CG_all$seqnames,
                             ranges = LL_G1_lyr_CG_all_ranges)

LL_G4_lyr_CG_all_ranges <- IRanges(start = LL_G4_lyr_CG_all$start,
                                   end = LL_G4_lyr_CG_all$end)
LL_G4_lyr_CG_allG <- GRanges(seqnames = LL_G4_lyr_CG_all$seqnames,
                             ranges = LL_G4_lyr_CG_all_ranges)

# Find overlaps

LL_overlaps_lyr_CG <- findOverlaps(LL_G1_lyr_CG_allG, LL_G4_lyr_CG_allG)

length(unique(LL_overlaps_lyr_CG@from))

length(unique(LL_overlaps_lyr_CG@from)) / nrow(LL_G1_lyr_CG_all)

# 1886 our of 2611 DMRs are conserved in the 4th gen from the 1st (~69.5 %)

### CHG context, lyrata side

# Convert DMRs to GRanges objects 

HM_G1_lyr_CHG_all_ranges <- IRanges(start = HM_G1_lyr_CHG_all$start,
                                 end = HM_G1_lyr_CHG_all$end)
HM_G1_lyr_CHG_allG <- GRanges(seqnames = HM_G1_lyr_CHG_all$seqnames,
                           ranges = HM_G1_lyr_CHG_all_ranges)

HM_G4_lyr_CHG_all_ranges <- IRanges(start = HM_G4_lyr_CHG_all$start,
                                 end = HM_G4_lyr_CHG_all$end)
HM_G4_lyr_CHG_allG <- GRanges(seqnames = HM_G4_lyr_CHG_all$seqnames,
                           ranges = HM_G4_lyr_CHG_all_ranges)

# Find overlaps

HM_overlaps_lyr_CHG <- findOverlaps(HM_G1_lyr_CHG_allG, HM_G4_lyr_CHG_allG)

length(unique(HM_overlaps_lyr_CHG@from))

length(unique(HM_overlaps_lyr_CHG@from)) / nrow(HM_G1_lyr_CHG_all)

# 2047 our of 3567 DMRs are conserved in the 4th gen from the 1st (~57.4 %)

# Convert DMRs to GRanges objects 

LL_G1_lyr_CHG_all_ranges <- IRanges(start = LL_G1_lyr_CHG_all$start,
                                    end = LL_G1_lyr_CHG_all$end)
LL_G1_lyr_CHG_allG <- GRanges(seqnames = LL_G1_lyr_CHG_all$seqnames,
                              ranges = LL_G1_lyr_CHG_all_ranges)

LL_G4_lyr_CHG_all_ranges <- IRanges(start = LL_G4_lyr_CHG_all$start,
                                    end = LL_G4_lyr_CHG_all$end)
LL_G4_lyr_CHG_allG <- GRanges(seqnames = LL_G4_lyr_CHG_all$seqnames,
                              ranges = LL_G4_lyr_CHG_all_ranges)

# Find overlaps

LL_overlaps_lyr_CHG <- findOverlaps(LL_G1_lyr_CHG_allG, LL_G4_lyr_CHG_allG)

length(unique(LL_overlaps_lyr_CHG@from))

length(unique(LL_overlaps_lyr_CHG@from)) / nrow(LL_G1_lyr_CHG_all)

# 2108 our of 3435 DMRs are conserved in the 4th gen from the 1st (~61.4 %)
# in hot conditions

### CHH context, lyrata side

# Convert DMRs to GRanges objects 

HM_G1_lyr_CHH_all_ranges <- IRanges(start = HM_G1_lyr_CHH_all$start,
                                 end = HM_G1_lyr_CHH_all$end)
HM_G1_lyr_CHH_allG <- GRanges(seqnames = HM_G1_lyr_CHH_all$seqnames,
                           ranges = HM_G1_lyr_CHH_all_ranges)

HM_G4_lyr_CHH_all_ranges <- IRanges(start = HM_G4_lyr_CHH_all$start,
                                 end = HM_G4_lyr_CHH_all$end)
HM_G4_lyr_CHH_allG <- GRanges(seqnames = HM_G4_lyr_CHH_all$seqnames,
                           ranges = HM_G4_lyr_CHH_all_ranges)

# Find overlaps

HM_overlaps_lyr_CHH <- findOverlaps(HM_G1_lyr_CHH_allG, HM_G4_lyr_CHH_allG)

length(unique(HM_overlaps_lyr_CHH@from))

length(unique(HM_overlaps_lyr_CHH@from)) / nrow(HM_G1_lyr_CHH_all)

# 1569 our of 2362 DMRs are conserved in the 4th gen from the 1st (~62.8 %)

# Convert DMRs to GRanges objects 

LL_G1_lyr_CHH_all_ranges <- IRanges(start = LL_G1_lyr_CHH_all$start,
                                    end = LL_G1_lyr_CHH_all$end)
LL_G1_lyr_CHH_allG <- GRanges(seqnames = LL_G1_lyr_CHH_all$seqnames,
                              ranges = LL_G1_lyr_CHH_all_ranges)

LL_G4_lyr_CHH_all_ranges <- IRanges(start = LL_G4_lyr_CHH_all$start,
                                    end = LL_G4_lyr_CHH_all$end)
LL_G4_lyr_CHH_allG <- GRanges(seqnames = LL_G4_lyr_CHH_all$seqnames,
                              ranges = LL_G4_lyr_CHH_all_ranges)

# Find overlaps

LL_overlaps_lyr_CHH <- findOverlaps(LL_G1_lyr_CHH_allG, LL_G4_lyr_CHH_allG)

length(unique(LL_overlaps_lyr_CHH@from))

length(unique(LL_overlaps_lyr_CHH@from)) / nrow(LL_G1_lyr_CHH_all)

# 1429 our of 5044 DMRs are conserved in the 4th gen from the 1st (~28.3 %)
# in hot conditions