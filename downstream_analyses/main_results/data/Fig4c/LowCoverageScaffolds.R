### This script has the goal to produce a list
### of scaffolds showing a low coverage after
### sequencing

### Parameters

### Coverage threshold

cov_lim <- 2

### Import libraries

library(data.table)
library(tidyverse)
library(GenomicRanges)

### Import all density files

setwd("/path/to/dens.txt")

## coverage files

HM_RS7_G4_1_hal_dens <- read.csv("HM_RS7_G4_1_hal_dens.txt")
HM_RS7_G4_2_hal_dens <- read.csv("HM_RS7_G4_2_hal_dens.txt")
HM_RS7_G4_3_hal_dens <- read.csv("HM_RS7_G4_3_hal_dens.txt")

HM_RS7_G4_1_lyr_dens <- read.csv("HM_RS7_G4_1_lyr_dens.txt")
HM_RS7_G4_2_lyr_dens <- read.csv("HM_RS7_G4_2_lyr_dens.txt")
HM_RS7_G4_3_lyr_dens <- read.csv("HM_RS7_G4_3_lyr_dens.txt")

LL_RS7_G4_1_hal_dens <- read.csv("LL_RS7_G4_1_hal_dens.txt")
LL_RS7_G4_2_hal_dens <- read.csv("LL_RS7_G4_2_hal_dens.txt")
LL_RS7_G4_3_hal_dens <- read.csv("LL_RS7_G4_3_hal_dens.txt")

LL_RS7_G4_1_lyr_dens <- read.csv("LL_RS7_G4_1_lyr_dens.txt")
LL_RS7_G4_2_lyr_dens <- read.csv("LL_RS7_G4_2_lyr_dens.txt")
LL_RS7_G4_3_lyr_dens <- read.csv("LL_RS7_G4_3_lyr_dens.txt")

## total density

hal_totC_dens <- read.delim("hal_totC_dens.txt")
lyr_totC_dens <- read.delim("lyr_totC_dens.txt")

### Compute real coverage

HM_RS7_G4_1_hal_real_coverage <- HM_RS7_G4_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_2_hal_real_coverage <- HM_RS7_G4_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_3_hal_real_coverage <- HM_RS7_G4_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7_G4_1_lyr_real_coverage <- HM_RS7_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_2_lyr_real_coverage <- HM_RS7_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_3_lyr_real_coverage <- HM_RS7_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7_G4_1_hal_real_coverage <- LL_RS7_G4_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G4_2_hal_real_coverage <- LL_RS7_G4_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G4_3_hal_real_coverage <- LL_RS7_G4_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7_G4_1_lyr_real_coverage <- LL_RS7_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G4_2_lyr_real_coverage <- LL_RS7_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G4_3_lyr_real_coverage <- LL_RS7_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens


### Compute average coverage

HM_RS7_G4_hal_avg_coverage <- c((HM_RS7_G4_1_hal_real_coverage +
                                   HM_RS7_G4_2_hal_real_coverage +
                                   HM_RS7_G4_3_hal_real_coverage) / 3)

HM_RS7_G4_lyr_avg_coverage <- c((HM_RS7_G4_1_lyr_real_coverage +
                                   HM_RS7_G4_2_lyr_real_coverage +
                                   HM_RS7_G4_3_lyr_real_coverage) / 3)

LL_RS7_G4_hal_avg_coverage <- c((LL_RS7_G4_1_hal_real_coverage +
                                    LL_RS7_G4_2_hal_real_coverage +
                                    LL_RS7_G4_3_hal_real_coverage) / 3)

LL_RS7_G4_lyr_avg_coverage <- c((LL_RS7_G4_1_lyr_real_coverage +
                                    LL_RS7_G4_2_lyr_real_coverage +
                                    LL_RS7_G4_3_lyr_real_coverage) / 3)

### Check regions with low coverage

HM_RS7_G4_hal_low_coverage <- which(HM_RS7_G4_hal_avg_coverage < cov_lim)
HM_RS7_G4_lyr_low_coverage <- which(HM_RS7_G4_lyr_avg_coverage < cov_lim)

LL_RS7_G4_hal_low_coverage <- which(LL_RS7_G4_hal_avg_coverage < cov_lim)
LL_RS7_G4_lyr_low_coverage <- which(LL_RS7_G4_lyr_avg_coverage < cov_lim)

### Check regions with high coverage

HM_RS7_G4_hal_high_coverage <- which(HM_RS7_G4_hal_avg_coverage > cov_lim_high)
HM_RS7_G4_lyr_high_coverage <- which(HM_RS7_G4_lyr_avg_coverage > cov_lim_high)

LL_RS7_G4_hal_high_coverage <- which(LL_RS7_G4_hal_avg_coverage > cov_lim_high)
LL_RS7_G4_lyr_high_coverage <- which(LL_RS7_G4_lyr_avg_coverage > cov_lim_high)

### Match regions to original coordinates

HM_RS7_G4_hal_low_coord <- hal_totC_dens[HM_RS7_G4_hal_low_coverage, 1:3]
HM_RS7_G4_lyr_low_coord <- lyr_totC_dens[HM_RS7_G4_lyr_low_coverage, 1:3]

LL_RS7_G4_hal_low_coord <- hal_totC_dens[LL_RS7_G4_hal_low_coverage, 1:3]
LL_RS7_G4_lyr_low_coord <- lyr_totC_dens[LL_RS7_G4_lyr_low_coverage, 1:3]


HM_RS7_G4_hal_high_coord <- hal_totC_dens[HM_RS7_G4_hal_high_coverage, 1:3]
HM_RS7_G4_lyr_high_coord <- lyr_totC_dens[HM_RS7_G4_lyr_high_coverage, 1:3]

LL_RS7_G4_hal_high_coord <- hal_totC_dens[LL_RS7_G4_hal_high_coverage, 1:3]
LL_RS7_G4_lyr_high_coord <- lyr_totC_dens[LL_RS7_G4_lyr_high_coverage, 1:3]

### Export files with low coverage scaffolds

setwd("path/to/output")

write.table(HM_RS7_G4_hal_low_coord, 
            file = "HM_RS7_G4_hal_LowCovRegions.txt",
            quote = F,
            col.names = F,
            row.names = F)

write.table(HM_RS7_G4_lyr_low_coord, 
            file = "HM_RS7_G4_lyr_LowCovRegions.txt",
            quote = F,
            col.names = F,
            row.names = F)

write.table(LL_RS7_G4_hal_low_coord, 
            file = "LL_RS7_G4_hal_LowCovRegions.txt",
            quote = F,
            col.names = F,
            row.names = F)

write.table(LL_RS7_G4_lyr_low_coord, 
            file = "LL_RS7_G4_lyr_LowCovRegions.txt",
            quote = F,
            col.names = F,
            row.names = F)
