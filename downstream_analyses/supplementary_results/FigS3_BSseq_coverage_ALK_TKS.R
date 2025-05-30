### This script aims at plotting cytosine coverage
### over the whole halleri and lyrata genome.
### The input data is binned cytosine coverage

### import libraries

library(ggplot2)
library(karyoploteR)
library(tidyverse)
library(patchwork)

### import data

## progenitors data
folder_path <- c("~/path/to/data/Fig4c")

HM_hal_G1_1_dens <- read.csv(paste0(folder_path, "HM_hal_G1_1_dens.txt"))
HM_hal_G1_2_dens <- read.csv(paste0(folder_path, "HM_hal_G1_2_dens.txt"))
HM_hal_G1_3_dens <- read.csv(paste0(folder_path, "HM_hal_G1_3_dens.txt"))

HM_hal_G4_1_dens <- read.csv(paste0(folder_path, "HM_hal_G4_1_dens.txt"))
HM_hal_G4_2_dens <- read.csv(paste0(folder_path, "HM_hal_G4_2_dens.txt"))
HM_hal_G4_3_dens <- read.csv(paste0(folder_path, "HM_hal_G4_3_dens.txt"))

HM_lyr_G1_1_dens <- read.csv(paste0(folder_path, "HM_lyr_G1_1_dens.txt"))
HM_lyr_G1_2_dens <- read.csv(paste0(folder_path, "HM_lyr_G1_2_dens.txt"))
HM_lyr_G1_3_dens <- read.csv(paste0(folder_path, "HM_lyr_G1_3_dens.txt"))

HM_lyr_G4_1_dens <- read.csv(paste0(folder_path, "HM_lyr_G4_1_dens.txt"))
HM_lyr_G4_2_dens <- read.csv(paste0(folder_path, "HM_lyr_G4_2_dens.txt"))
HM_lyr_G4_3_dens <- read.csv(paste0(folder_path, "HM_lyr_G4_3_dens.txt"))

LL_hal_G1_1_dens <- read.csv(paste0(folder_path, "LL_hal_G1_1_dens.txt"))
LL_hal_G1_2_dens <- read.csv(paste0(folder_path, "LL_hal_G1_2_dens.txt"))
LL_hal_G1_3_dens <- read.csv(paste0(folder_path, "LL_hal_G1_3_dens.txt"))

LL_hal_G4_1_dens <- read.csv(paste0(folder_path, "LL_hal_G4_1_dens.txt"))
LL_hal_G4_2_dens <- read.csv(paste0(folder_path, "LL_hal_G4_2_dens.txt"))

LL_lyr_G1_1_dens <- read.csv(paste0(folder_path, "LL_lyr_G1_1_dens.txt"))
LL_lyr_G1_2_dens <- read.csv(paste0(folder_path, "LL_lyr_G1_2_dens.txt"))
LL_lyr_G1_3_dens <- read.csv(paste0(folder_path, "LL_lyr_G1_3_dens.txt"))

LL_lyr_G4_1_dens <- read.csv(paste0(folder_path, "LL_lyr_G4_1_dens.txt"))
LL_lyr_G4_2_dens <- read.csv(paste0(folder_path, "LL_lyr_G4_2_dens.txt"))

## polyploid data

HM_RS7_G1_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7K_G1_1_hal_dens.txt"))
HM_RS7_G1_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7K_G1_2_hal_dens.txt"))
HM_RS7_G1_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7K_G1_3_hal_dens.txt"))

HM_RS7_G1_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7K_G1_1_lyr_dens.txt"))
HM_RS7_G1_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7K_G1_2_lyr_dens.txt"))
HM_RS7_G1_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7K_G1_3_lyr_dens.txt"))

HM_RS7_G4_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_hal_dens.txt"))
HM_RS7_G4_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_hal_dens.txt"))
HM_RS7_G4_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_hal_dens.txt"))

HM_RS7_G4_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_lyr_dens.txt"))
HM_RS7_G4_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_lyr_dens.txt"))
HM_RS7_G4_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_lyr_dens.txt"))

LL_RS7_G1_1_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_1_hal_dens.txt"))
LL_RS7_G1_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_2_hal_dens.txt"))
LL_RS7_G1_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_3_hal_dens.txt"))

LL_RS7_G1_1_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_1_lyr_dens.txt"))
LL_RS7_G1_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_2_lyr_dens.txt"))
LL_RS7_G1_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_3_lyr_dens.txt"))

LL_RS7K_G4_1_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_1_hal_dens.txt"))
LL_RS7K_G4_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_2_hal_dens.txt"))
LL_RS7K_G4_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_3_hal_dens.txt"))

LL_RS7K_G4_1_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_1_lyr_dens.txt"))
LL_RS7K_G4_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_2_lyr_dens.txt"))
LL_RS7K_G4_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_3_lyr_dens.txt"))

HM_ALK_G1_1_hal_dens <- read.csv(paste0(folder_path, "HM_ALK_G1_1_hal_dens.txt"))
HM_ALK_G1_2_hal_dens <- read.csv(paste0(folder_path, "HM_ALK_G1_2_hal_dens.txt"))
HM_ALK_G1_3_hal_dens <- read.csv(paste0(folder_path, "HM_ALK_G1_3_hal_dens.txt"))

HM_ALK_G1_1_lyr_dens <- read.csv(paste0(folder_path, "HM_ALK_G1_1_lyr_dens.txt"))
HM_ALK_G1_2_lyr_dens <- read.csv(paste0(folder_path, "HM_ALK_G1_2_lyr_dens.txt"))
HM_ALK_G1_3_lyr_dens <- read.csv(paste0(folder_path, "HM_ALK_G1_3_lyr_dens.txt"))

HM_ALK_G4_1_hal_dens <- read.csv(paste0(folder_path, "HM_ALK_G4_1_hal_dens.txt"))
HM_ALK_G4_2_hal_dens <- read.csv(paste0(folder_path, "HM_ALK_G4_2_hal_dens.txt"))
HM_ALK_G4_3_hal_dens <- read.csv(paste0(folder_path, "HM_ALK_G4_3_hal_dens.txt"))

HM_ALK_G4_1_lyr_dens <- read.csv(paste0(folder_path, "HM_ALK_G4_1_lyr_dens.txt"))
HM_ALK_G4_2_lyr_dens <- read.csv(paste0(folder_path, "HM_ALK_G4_2_lyr_dens.txt"))
HM_ALK_G4_3_lyr_dens <- read.csv(paste0(folder_path, "HM_ALK_G4_3_lyr_dens.txt"))

LL_ALK_G1_1_hal_dens <- read.csv(paste0(folder_path, "LL_ALK_G1_1_hal_dens.txt"))
LL_ALK_G1_2_hal_dens <- read.csv(paste0(folder_path, "LL_ALK_G1_2_hal_dens.txt"))
LL_ALK_G1_3_hal_dens <- read.csv(paste0(folder_path, "LL_ALK_G1_3_hal_dens.txt"))

LL_ALK_G1_1_lyr_dens <- read.csv(paste0(folder_path, "LL_ALK_G1_1_lyr_dens.txt"))
LL_ALK_G1_2_lyr_dens <- read.csv(paste0(folder_path, "LL_ALK_G1_2_lyr_dens.txt"))
LL_ALK_G1_3_lyr_dens <- read.csv(paste0(folder_path, "LL_ALK_G1_3_lyr_dens.txt"))

LL_ALK_G4_1_hal_dens <- read.csv(paste0(folder_path, "LL_ALK_G4_1_hal_dens.txt"))
LL_ALK_G4_2_hal_dens <- read.csv(paste0(folder_path, "LL_ALK_G4_2_hal_dens.txt"))
LL_ALK_G4_3_hal_dens <- read.csv(paste0(folder_path, "LL_ALK_G4_3_hal_dens.txt"))

LL_ALK_G4_1_lyr_dens <- read.csv(paste0(folder_path, "LL_ALK_G4_1_lyr_dens.txt"))
LL_ALK_G4_2_lyr_dens <- read.csv(paste0(folder_path, "LL_ALK_G4_2_lyr_dens.txt"))
LL_ALK_G4_3_lyr_dens <- read.csv(paste0(folder_path, "LL_ALK_G4_3_lyr_dens.txt"))

HM_TKS_G1_1_hal_dens <- read.csv(paste0(folder_path, "HM_TKS_G1_1_hal_dens.txt"))
HM_TKS_G1_2_hal_dens <- read.csv(paste0(folder_path, "HM_TKS_G1_2_hal_dens.txt"))
HM_TKS_G1_3_hal_dens <- read.csv(paste0(folder_path, "HM_TKS_G1_3_hal_dens.txt"))

HM_TKS_G1_1_lyr_dens <- read.csv(paste0(folder_path, "HM_TKS_G1_1_lyr_dens.txt"))
HM_TKS_G1_2_lyr_dens <- read.csv(paste0(folder_path, "HM_TKS_G1_2_lyr_dens.txt"))
HM_TKS_G1_3_lyr_dens <- read.csv(paste0(folder_path, "HM_TKS_G1_3_lyr_dens.txt"))

HM_TKS_G5_1_hal_dens <- read.csv(paste0(folder_path, "HM_TKS_G5_1_hal_dens.txt"))
HM_TKS_G5_2_hal_dens <- read.csv(paste0(folder_path, "HM_TKS_G5_2_hal_dens.txt"))
HM_TKS_G5_3_hal_dens <- read.csv(paste0(folder_path, "HM_TKS_G5_3_hal_dens.txt"))

HM_TKS_G5_1_lyr_dens <- read.csv(paste0(folder_path, "HM_TKS_G5_1_lyr_dens.txt"))
HM_TKS_G5_2_lyr_dens <- read.csv(paste0(folder_path, "HM_TKS_G5_2_lyr_dens.txt"))
HM_TKS_G5_3_lyr_dens <- read.csv(paste0(folder_path, "HM_TKS_G5_3_lyr_dens.txt"))

LL_TKS_G1_1_hal_dens <- read.csv(paste0(folder_path, "LL_TKS_G1_1_hal_dens.txt"))
LL_TKS_G1_2_hal_dens <- read.csv(paste0(folder_path, "LL_TKS_G1_2_hal_dens.txt"))
LL_TKS_G1_3_hal_dens <- read.csv(paste0(folder_path, "LL_TKS_G1_3_hal_dens.txt"))

LL_TKS_G1_1_lyr_dens <- read.csv(paste0(folder_path, "LL_TKS_G1_1_lyr_dens.txt"))
LL_TKS_G1_2_lyr_dens <- read.csv(paste0(folder_path, "LL_TKS_G1_2_lyr_dens.txt"))
LL_TKS_G1_3_lyr_dens <- read.csv(paste0(folder_path, "LL_TKS_G1_3_lyr_dens.txt"))

LL_TKS_G5_1_hal_dens <- read.csv(paste0(folder_path, "LL_TKS_G5_1_hal_dens.txt"))
LL_TKS_G5_2_hal_dens <- read.csv(paste0(folder_path, "LL_TKS_G5_2_hal_dens.txt"))
LL_TKS_G5_3_hal_dens <- read.csv(paste0(folder_path, "LL_TKS_G5_3_hal_dens.txt"))

LL_TKS_G5_1_lyr_dens <- read.csv(paste0(folder_path, "LL_TKS_G5_1_lyr_dens.txt"))
LL_TKS_G5_2_lyr_dens <- read.csv(paste0(folder_path, "LL_TKS_G5_2_lyr_dens.txt"))
LL_TKS_G5_3_lyr_dens <- read.csv(paste0(folder_path, "LL_TKS_G5_3_lyr_dens.txt"))

## total density

hal_totC_dens <- read.delim(paste0(folder_path, "hal_totC_dens_Dario.txt"))
lyr_totC_dens <- read.delim(paste0(folder_path, "lyr_totC_dens_Dario.txt"))

### We first compute the real cytosine coverage per window for all samples

## progenitors data

HM_hal_G1_1_real_coverage <- HM_hal_G1_1_dens$window_scores / hal_totC_dens$dens
HM_hal_G1_2_real_coverage <- HM_hal_G1_2_dens$window_scores / hal_totC_dens$dens
HM_hal_G1_3_real_coverage <- HM_hal_G1_3_dens$window_scores / hal_totC_dens$dens

HM_lyr_G1_1_real_coverage <- HM_lyr_G1_1_dens$window_scores / lyr_totC_dens$dens
HM_lyr_G1_2_real_coverage <- HM_lyr_G1_2_dens$window_scores / lyr_totC_dens$dens
HM_lyr_G1_3_real_coverage <- HM_lyr_G1_3_dens$window_scores / lyr_totC_dens$dens

HM_hal_G4_1_real_coverage <- HM_hal_G4_1_dens$window_scores / hal_totC_dens$dens
HM_hal_G4_2_real_coverage <- HM_hal_G4_2_dens$window_scores / hal_totC_dens$dens
HM_hal_G4_3_real_coverage <- HM_hal_G4_3_dens$window_scores / hal_totC_dens$dens

HM_lyr_G4_1_real_coverage <- HM_lyr_G4_1_dens$window_scores / lyr_totC_dens$dens
HM_lyr_G4_2_real_coverage <- HM_lyr_G4_2_dens$window_scores / lyr_totC_dens$dens
HM_lyr_G4_3_real_coverage <- HM_lyr_G4_3_dens$window_scores / lyr_totC_dens$dens

LL_hal_G1_1_real_coverage <- LL_hal_G1_1_dens$window_scores / hal_totC_dens$dens
LL_hal_G1_2_real_coverage <- LL_hal_G1_2_dens$window_scores / hal_totC_dens$dens
LL_hal_G1_3_real_coverage <- LL_hal_G1_3_dens$window_scores / hal_totC_dens$dens

LL_lyr_G1_1_real_coverage <- LL_lyr_G1_1_dens$window_scores / lyr_totC_dens$dens
LL_lyr_G1_2_real_coverage <- LL_lyr_G1_2_dens$window_scores / lyr_totC_dens$dens
LL_lyr_G1_3_real_coverage <- LL_lyr_G1_3_dens$window_scores / lyr_totC_dens$dens

LL_hal_G4_1_real_coverage <- LL_hal_G4_1_dens$window_scores / hal_totC_dens$dens
LL_hal_G4_2_real_coverage <- LL_hal_G4_2_dens$window_scores / hal_totC_dens$dens

LL_lyr_G4_1_real_coverage <- LL_lyr_G4_1_dens$window_scores / lyr_totC_dens$dens
LL_lyr_G4_2_real_coverage <- LL_lyr_G4_2_dens$window_scores / lyr_totC_dens$dens


## polyploid data

HM_RS7_G1_1_hal_real_coverage <- HM_RS7_G1_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G1_2_hal_real_coverage <- HM_RS7_G1_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G1_3_hal_real_coverage <- HM_RS7_G1_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7_G1_1_lyr_real_coverage <- HM_RS7_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G1_2_lyr_real_coverage <- HM_RS7_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G1_3_lyr_real_coverage <- HM_RS7_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_RS7_G4_1_hal_real_coverage <- HM_RS7_G4_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_2_hal_real_coverage <- HM_RS7_G4_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_3_hal_real_coverage <- HM_RS7_G4_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7_G4_1_lyr_real_coverage <- HM_RS7_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_2_lyr_real_coverage <- HM_RS7_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_3_lyr_real_coverage <- HM_RS7_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7_G1_1_hal_real_coverage <- LL_RS7_G1_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G1_2_hal_real_coverage <- LL_RS7_G1_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G1_3_hal_real_coverage <- LL_RS7_G1_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7_G1_1_lyr_real_coverage <- LL_RS7_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G1_2_lyr_real_coverage <- LL_RS7_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G1_3_lyr_real_coverage <- LL_RS7_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7K_G4_1_hal_real_coverage <- LL_RS7K_G4_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_2_hal_real_coverage <- LL_RS7K_G4_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_3_hal_real_coverage <- LL_RS7K_G4_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7K_G4_1_lyr_real_coverage <- LL_RS7K_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_2_lyr_real_coverage <- LL_RS7K_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_3_lyr_real_coverage <- LL_RS7K_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_ALK_G1_1_hal_real_coverage <- HM_ALK_G1_1_hal_dens$window_scores / hal_totC_dens$dens
HM_ALK_G1_2_hal_real_coverage <- HM_ALK_G1_2_hal_dens$window_scores / hal_totC_dens$dens
HM_ALK_G1_3_hal_real_coverage <- HM_ALK_G1_3_hal_dens$window_scores / hal_totC_dens$dens

HM_ALK_G1_1_lyr_real_coverage <- HM_ALK_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_ALK_G1_2_lyr_real_coverage <- HM_ALK_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_ALK_G1_3_lyr_real_coverage <- HM_ALK_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_ALK_G4_1_hal_real_coverage <- HM_ALK_G4_1_hal_dens$window_scores / hal_totC_dens$dens
HM_ALK_G4_2_hal_real_coverage <- HM_ALK_G4_2_hal_dens$window_scores / hal_totC_dens$dens
HM_ALK_G4_3_hal_real_coverage <- HM_ALK_G4_3_hal_dens$window_scores / hal_totC_dens$dens

HM_ALK_G4_1_lyr_real_coverage <- HM_ALK_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_ALK_G4_2_lyr_real_coverage <- HM_ALK_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_ALK_G4_3_lyr_real_coverage <- HM_ALK_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_ALK_G1_1_hal_real_coverage <- LL_ALK_G1_1_hal_dens$window_scores / hal_totC_dens$dens
LL_ALK_G1_2_hal_real_coverage <- LL_ALK_G1_2_hal_dens$window_scores / hal_totC_dens$dens
LL_ALK_G1_3_hal_real_coverage <- LL_ALK_G1_3_hal_dens$window_scores / hal_totC_dens$dens

LL_ALK_G1_1_lyr_real_coverage <- LL_ALK_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_ALK_G1_2_lyr_real_coverage <- LL_ALK_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_ALK_G1_3_lyr_real_coverage <- LL_ALK_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_ALK_G4_1_hal_real_coverage <- LL_ALK_G4_1_hal_dens$window_scores / hal_totC_dens$dens
LL_ALK_G4_2_hal_real_coverage <- LL_ALK_G4_2_hal_dens$window_scores / hal_totC_dens$dens
LL_ALK_G4_3_hal_real_coverage <- LL_ALK_G4_3_hal_dens$window_scores / hal_totC_dens$dens

LL_ALK_G4_1_lyr_real_coverage <- LL_ALK_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_ALK_G4_2_lyr_real_coverage <- LL_ALK_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_ALK_G4_3_lyr_real_coverage <- LL_ALK_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_TKS_G1_1_hal_real_coverage <- HM_TKS_G1_1_hal_dens$window_scores / hal_totC_dens$dens
HM_TKS_G1_2_hal_real_coverage <- HM_TKS_G1_2_hal_dens$window_scores / hal_totC_dens$dens
HM_TKS_G1_3_hal_real_coverage <- HM_TKS_G1_3_hal_dens$window_scores / hal_totC_dens$dens

HM_TKS_G1_1_lyr_real_coverage <- HM_TKS_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_TKS_G1_2_lyr_real_coverage <- HM_TKS_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_TKS_G1_3_lyr_real_coverage <- HM_TKS_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_TKS_G5_1_hal_real_coverage <- HM_TKS_G5_1_hal_dens$window_scores / hal_totC_dens$dens
HM_TKS_G5_2_hal_real_coverage <- HM_TKS_G5_2_hal_dens$window_scores / hal_totC_dens$dens
HM_TKS_G5_3_hal_real_coverage <- HM_TKS_G5_3_hal_dens$window_scores / hal_totC_dens$dens

HM_TKS_G5_1_lyr_real_coverage <- HM_TKS_G5_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_TKS_G5_2_lyr_real_coverage <- HM_TKS_G5_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_TKS_G5_3_lyr_real_coverage <- HM_TKS_G5_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_TKS_G1_1_hal_real_coverage <- LL_TKS_G1_1_hal_dens$window_scores / hal_totC_dens$dens
LL_TKS_G1_2_hal_real_coverage <- LL_TKS_G1_2_hal_dens$window_scores / hal_totC_dens$dens
LL_TKS_G1_3_hal_real_coverage <- LL_TKS_G1_3_hal_dens$window_scores / hal_totC_dens$dens

LL_TKS_G1_1_lyr_real_coverage <- LL_TKS_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_TKS_G1_2_lyr_real_coverage <- LL_TKS_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_TKS_G1_3_lyr_real_coverage <- LL_TKS_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_TKS_G5_1_hal_real_coverage <- LL_TKS_G5_1_hal_dens$window_scores / hal_totC_dens$dens
LL_TKS_G5_2_hal_real_coverage <- LL_TKS_G5_2_hal_dens$window_scores / hal_totC_dens$dens
LL_TKS_G5_3_hal_real_coverage <- LL_TKS_G5_3_hal_dens$window_scores / hal_totC_dens$dens

LL_TKS_G5_1_lyr_real_coverage <- LL_TKS_G5_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_TKS_G5_2_lyr_real_coverage <- LL_TKS_G5_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_TKS_G5_3_lyr_real_coverage <- LL_TKS_G5_3_lyr_dens$window_scores / lyr_totC_dens$dens


### Compute the average for all samples

HM_hal_G1_avg_coverage <- c((HM_hal_G1_1_real_coverage +
                            HM_hal_G1_2_real_coverage +
                            HM_hal_G1_3_real_coverage) / 3)

HM_lyr_G1_avg_coverage <- c((HM_lyr_G1_1_real_coverage +
                            HM_lyr_G1_2_real_coverage +
                            HM_lyr_G1_3_real_coverage) / 3)

HM_hal_G4_avg_coverage <- c((HM_hal_G4_1_real_coverage +
                               HM_hal_G4_2_real_coverage +
                               HM_hal_G4_3_real_coverage) / 3)

HM_lyr_G4_avg_coverage <- c((HM_lyr_G4_1_real_coverage +
                               HM_lyr_G4_2_real_coverage +
                               HM_lyr_G4_3_real_coverage) / 3)

LL_hal_G1_avg_coverage <- c((LL_hal_G1_1_real_coverage +
                               LL_hal_G1_2_real_coverage +
                               LL_hal_G1_3_real_coverage) / 3)

LL_lyr_G1_avg_coverage <- c((LL_lyr_G1_1_real_coverage +
                               LL_lyr_G1_2_real_coverage +
                               LL_lyr_G1_3_real_coverage) / 3)

LL_hal_G4_avg_coverage <- c((LL_hal_G4_1_real_coverage +
                               LL_hal_G4_2_real_coverage) / 2)

LL_lyr_G4_avg_coverage <- c((LL_lyr_G4_1_real_coverage +
                               LL_lyr_G4_2_real_coverage) / 2)

HM_RS7_G1_hal_avg_coverage <- c((HM_RS7_G1_1_hal_real_coverage +
                    HM_RS7_G1_2_hal_real_coverage +
                    HM_RS7_G1_3_hal_real_coverage) / 3)

HM_RS7_G1_lyr_avg_coverage <- c((HM_RS7_G1_1_lyr_real_coverage +
                    HM_RS7_G1_2_lyr_real_coverage +
                    HM_RS7_G1_3_lyr_real_coverage) / 3)

HM_RS7_G4_hal_avg_coverage <- c((HM_RS7_G4_1_hal_real_coverage +
                                HM_RS7_G4_2_hal_real_coverage +
                                HM_RS7_G4_3_hal_real_coverage) / 3)

HM_RS7_G4_lyr_avg_coverage <- c((HM_RS7_G4_1_lyr_real_coverage +
                                HM_RS7_G4_2_lyr_real_coverage +
                                HM_RS7_G4_3_lyr_real_coverage) / 3)

LL_RS7_G1_hal_avg_coverage <- c((LL_RS7_G1_1_hal_real_coverage +
                                    LL_RS7_G1_2_hal_real_coverage +
                                    LL_RS7_G1_3_hal_real_coverage) / 3)

LL_RS7_G1_lyr_avg_coverage <- c((LL_RS7_G1_1_lyr_real_coverage +
                                    LL_RS7_G1_2_lyr_real_coverage +
                                    LL_RS7_G1_3_lyr_real_coverage) / 3)

LL_RS7K_G4_hal_avg_coverage <- c((LL_RS7K_G4_1_hal_real_coverage +
                                    LL_RS7K_G4_2_hal_real_coverage +
                                    LL_RS7K_G4_3_hal_real_coverage) / 3)

LL_RS7K_G4_lyr_avg_coverage <- c((LL_RS7K_G4_1_lyr_real_coverage +
                                    LL_RS7K_G4_2_lyr_real_coverage +
                                    LL_RS7K_G4_3_lyr_real_coverage) / 3)

HM_ALK_G1_hal_avg_coverage <- c((HM_ALK_G1_1_hal_real_coverage +
                                   HM_ALK_G1_2_hal_real_coverage +
                                   HM_ALK_G1_3_hal_real_coverage) / 3)

HM_ALK_G1_lyr_avg_coverage <- c((HM_ALK_G1_1_lyr_real_coverage +
                                   HM_ALK_G1_2_lyr_real_coverage +
                                   HM_ALK_G1_3_lyr_real_coverage) / 3)

HM_ALK_G4_hal_avg_coverage <- c((HM_ALK_G4_1_hal_real_coverage +
                                 HM_ALK_G4_2_hal_real_coverage +
                                 HM_ALK_G4_3_hal_real_coverage) / 3)

HM_ALK_G4_lyr_avg_coverage <- c((HM_ALK_G4_1_lyr_real_coverage +
                                 HM_ALK_G4_2_lyr_real_coverage +
                                 HM_ALK_G4_3_lyr_real_coverage) / 3)

LL_ALK_G1_hal_avg_coverage <- c((LL_ALK_G1_1_hal_real_coverage +
                                   LL_ALK_G1_2_hal_real_coverage +
                                   LL_ALK_G1_3_hal_real_coverage) / 3)

LL_ALK_G1_lyr_avg_coverage <- c((LL_ALK_G1_1_lyr_real_coverage +
                                   LL_ALK_G1_2_lyr_real_coverage +
                                   LL_ALK_G1_3_lyr_real_coverage) / 3)

LL_ALK_G4_hal_avg_coverage <- c((LL_ALK_G4_1_hal_real_coverage +
                                   LL_ALK_G4_2_hal_real_coverage +
                                   LL_ALK_G4_3_hal_real_coverage) / 3)

LL_ALK_G4_lyr_avg_coverage <- c((LL_ALK_G4_1_lyr_real_coverage +
                                   LL_ALK_G4_2_lyr_real_coverage +
                                   LL_ALK_G4_3_lyr_real_coverage) / 3)

HM_TKS_G1_hal_avg_coverage <- c((HM_TKS_G1_1_hal_real_coverage +
                                   HM_TKS_G1_2_hal_real_coverage +
                                   HM_TKS_G1_3_hal_real_coverage) / 3)

HM_TKS_G1_lyr_avg_coverage <- c((HM_TKS_G1_1_lyr_real_coverage +
                                   HM_TKS_G1_2_lyr_real_coverage +
                                   HM_TKS_G1_3_lyr_real_coverage) / 3)

HM_TKS_G5_hal_avg_coverage <- c((HM_TKS_G5_1_hal_real_coverage +
                                   HM_TKS_G5_2_hal_real_coverage +
                                   HM_TKS_G5_3_hal_real_coverage) / 3)

HM_TKS_G5_lyr_avg_coverage <- c((HM_TKS_G5_1_lyr_real_coverage +
                                   HM_TKS_G5_2_lyr_real_coverage +
                                  HM_TKS_G5_3_lyr_real_coverage) / 3)

LL_TKS_G1_hal_avg_coverage <- c((LL_TKS_G1_1_hal_real_coverage +
                                   LL_TKS_G1_2_hal_real_coverage +
                                   LL_TKS_G1_3_hal_real_coverage) / 3)

LL_TKS_G1_lyr_avg_coverage <- c((LL_TKS_G1_1_lyr_real_coverage +
                                   LL_TKS_G1_2_lyr_real_coverage +
                                   LL_TKS_G1_3_lyr_real_coverage) / 3)

LL_TKS_G5_hal_avg_coverage <- c((LL_TKS_G5_1_hal_real_coverage +
                                   LL_TKS_G5_2_hal_real_coverage +
                                   LL_TKS_G5_3_hal_real_coverage) / 3)

LL_TKS_G5_lyr_avg_coverage <- c((LL_TKS_G5_1_lyr_real_coverage +
                                   LL_TKS_G5_2_lyr_real_coverage +
                                   LL_TKS_G5_3_lyr_real_coverage) / 3)

### Create new data frame for all the samples

HM_hal_G1 <- HM_hal_G1_1_dens
HM_hal_G1$coverage <- HM_hal_G1_avg_coverage

HM_lyr_G1 <- HM_lyr_G1_1_dens
HM_lyr_G1$coverage <- HM_lyr_G1_avg_coverage

HM_hal_G4 <- HM_hal_G1_1_dens
HM_hal_G4$coverage <- HM_hal_G4_avg_coverage

HM_lyr_G4 <- HM_lyr_G4_1_dens
HM_lyr_G4$coverage <- HM_lyr_G4_avg_coverage

LL_hal_G1 <- LL_hal_G1_1_dens
LL_hal_G1$coverage <- LL_hal_G1_avg_coverage

LL_lyr_G1 <- LL_lyr_G1_1_dens
LL_lyr_G1$coverage <- LL_lyr_G1_avg_coverage

LL_hal_G4 <- LL_hal_G1_1_dens
LL_hal_G4$coverage <- LL_hal_G4_avg_coverage

LL_lyr_G4 <- LL_lyr_G1_1_dens
LL_lyr_G4$coverage <- LL_lyr_G4_avg_coverage

HM_RS7_G1_hal <- HM_RS7_G1_1_hal_dens
HM_RS7_G1_hal$coverage <- HM_RS7_G1_hal_avg_coverage

HM_RS7_G1_lyr <- HM_RS7_G1_1_lyr_dens
HM_RS7_G1_lyr$coverage <- HM_RS7_G1_lyr_avg_coverage

HM_RS7_G4_hal <- HM_RS7_G4_1_hal_dens
HM_RS7_G4_hal$coverage <- HM_RS7_G4_hal_avg_coverage

HM_RS7_G4_lyr <- HM_RS7_G4_1_lyr_dens
HM_RS7_G4_lyr$coverage <- HM_RS7_G4_lyr_avg_coverage

LL_RS7_G1_lyr <- LL_RS7_G1_1_lyr_dens
LL_RS7_G1_lyr$coverage <- LL_RS7_G1_lyr_avg_coverage

LL_RS7_G1_hal <- LL_RS7_G1_1_hal_dens
LL_RS7_G1_hal$coverage <- LL_RS7_G1_hal_avg_coverage

LL_RS7K_G4_lyr <- LL_RS7K_G4_1_lyr_dens
LL_RS7K_G4_lyr$coverage <- LL_RS7K_G4_lyr_avg_coverage

LL_RS7K_G4_hal <- LL_RS7K_G4_1_hal_dens
LL_RS7K_G4_hal$coverage <- LL_RS7K_G4_hal_avg_coverage

HM_ALK_G1_hal <- HM_ALK_G1_1_hal_dens
HM_ALK_G1_hal$coverage <- HM_ALK_G1_hal_avg_coverage

HM_ALK_G1_lyr <- HM_ALK_G1_1_lyr_dens
HM_ALK_G1_lyr$coverage <- HM_ALK_G1_lyr_avg_coverage

HM_ALK_G4_hal <- HM_ALK_G4_1_hal_dens
HM_ALK_G4_hal$coverage <- HM_ALK_G4_hal_avg_coverage

HM_ALK_G4_lyr <- HM_ALK_G4_1_lyr_dens
HM_ALK_G4_lyr$coverage <- HM_ALK_G4_lyr_avg_coverage

LL_ALK_G1_hal <- LL_ALK_G1_1_hal_dens
LL_ALK_G1_hal$coverage <- LL_ALK_G1_hal_avg_coverage

LL_ALK_G1_lyr <- LL_ALK_G1_1_lyr_dens
LL_ALK_G1_lyr$coverage <- LL_ALK_G1_lyr_avg_coverage

LL_ALK_G4_hal <- LL_ALK_G4_1_hal_dens
LL_ALK_G4_hal$coverage <- LL_ALK_G4_hal_avg_coverage

LL_ALK_G4_lyr <- LL_ALK_G4_1_lyr_dens
LL_ALK_G4_lyr$coverage <- LL_ALK_G4_lyr_avg_coverage

HM_TKS_G1_hal <- HM_TKS_G1_1_hal_dens
HM_TKS_G1_hal$coverage <- HM_TKS_G1_hal_avg_coverage

HM_TKS_G1_lyr <- HM_TKS_G1_1_lyr_dens
HM_TKS_G1_lyr$coverage <- HM_TKS_G1_lyr_avg_coverage

HM_TKS_G5_hal <- HM_TKS_G5_1_hal_dens
HM_TKS_G5_hal$coverage <- HM_TKS_G5_hal_avg_coverage

HM_TKS_G5_lyr <- HM_TKS_G5_1_lyr_dens
HM_TKS_G5_lyr$coverage <- HM_TKS_G5_lyr_avg_coverage

LL_TKS_G1_hal <- LL_TKS_G1_1_hal_dens
LL_TKS_G1_hal$coverage <- LL_TKS_G1_hal_avg_coverage

LL_TKS_G1_lyr <- LL_TKS_G1_1_lyr_dens
LL_TKS_G1_lyr$coverage <- LL_TKS_G1_lyr_avg_coverage

LL_TKS_G5_hal <- LL_TKS_G5_1_hal_dens
LL_TKS_G5_hal$coverage <- LL_TKS_G5_hal_avg_coverage

LL_TKS_G5_lyr <- LL_TKS_G5_1_lyr_dens
LL_TKS_G5_lyr$coverage <- LL_TKS_G5_lyr_avg_coverage


### Compute total windows per chromosome

chrom_length_hal <- rep(0, 8)
chrom_length_lyr <- rep(0, 8)
counter <- 1
for (i in c("chr1", "chr2", "chr3",
            "chr4", "chr5", "chr6",
            "chr7", "chr8")){
  chrom_length_hal[counter] <- sum(LL_hal_G1_1_dens$seqnames == i)
  counter <- counter + 1
}

counter <- 1
for (i in c("chr9", "chr10", "chr11",
            "chr12", "chr13", "chr14",
            "chr15", "chr16")){
  chrom_length_lyr[counter] <- sum(LL_lyr_G1_1_dens$seqnames == i)
  counter <- counter + 1
}

### Plot coverage over each chromosome and each side

HM_hal_G1 <- mutate(HM_hal_G1, sample = "HM_hal_G1")
HM_hal_G4 <- mutate(HM_hal_G4, sample = "HM_hal_G4")
LL_hal_G1 <- mutate(LL_hal_G1, sample = "LL_hal_G1")
LL_hal_G4 <- mutate(LL_hal_G4, sample = "LL_hal_G4")
HM_RS7_G1_hal <- mutate(HM_RS7_G1_hal, sample = "HM_RS7_G1_hal")
HM_RS7_G4_hal <- mutate(HM_RS7_G4_hal, sample = "HM_RS7_G4_hal")
LL_RS7_G1_hal <- mutate(LL_RS7_G1_hal, sample = "LL_RS7_G1_hal")
LL_RS7K_G4_hal <- mutate(LL_RS7K_G4_hal, sample = "LL_RS7K_G4_hal")
HM_ALK_G1_hal <- mutate(HM_ALK_G1_hal, sample = "HM_ALK_G1_hal")
HM_ALK_G4_hal <- mutate(HM_ALK_G4_hal, sample = "HM_ALK_G4_hal")
LL_ALK_G1_hal <- mutate(LL_ALK_G1_hal, sample = "LL_ALK_G1_hal")
LL_ALK_G4_hal <- mutate(LL_ALK_G4_hal, sample = "LL_ALK_G4_hal")
HM_TKS_G1_hal <- mutate(HM_TKS_G1_hal, sample = "HM_TKS_G1_hal")
HM_TKS_G5_hal <- mutate(HM_TKS_G5_hal, sample = "HM_TKS_G5_hal")
LL_TKS_G1_hal <- mutate(LL_TKS_G1_hal, sample = "LL_TKS_G1_hal")
LL_TKS_G5_hal <- mutate(LL_TKS_G5_hal, sample = "LL_TKS_G5_hal")

HM_lyr_G1 <- mutate(HM_lyr_G1, sample = "HM_lyr_G1")
HM_lyr_G4 <- mutate(HM_lyr_G4, sample = "HM_lyr_G4")
LL_lyr_G1 <- mutate(LL_lyr_G1, sample = "LL_lyr_G1")
LL_lyr_G4 <- mutate(LL_lyr_G4, sample = "LL_lyr_G4")
HM_RS7_G1_lyr <- mutate(HM_RS7_G1_lyr, sample = "HM_RS7_G1_lyr")
HM_RS7_G4_lyr <- mutate(HM_RS7_G4_lyr, sample = "HM_RS7_G4_lyr")
LL_RS7_G1_lyr <- mutate(LL_RS7_G1_lyr, sample = "LL_RS7_G1_lyr")
LL_RS7K_G4_lyr <- mutate(LL_RS7K_G4_lyr, sample = "LL_RS7K_G4_lyr")
HM_ALK_G1_lyr <- mutate(HM_ALK_G1_lyr, sample = "HM_ALK_G1_lyr")
HM_ALK_G4_lyr <- mutate(HM_ALK_G4_lyr, sample = "HM_ALK_G4_lyr")
LL_ALK_G1_lyr <- mutate(LL_ALK_G1_lyr, sample = "LL_ALK_G1_lyr")
LL_ALK_G4_lyr <- mutate(LL_ALK_G4_lyr, sample = "LL_ALK_G4_lyr")
HM_TKS_G1_lyr <- mutate(HM_TKS_G1_lyr, sample = "HM_TKS_G1_lyr")
HM_TKS_G5_lyr <- mutate(HM_TKS_G5_lyr, sample = "HM_TKS_G5_lyr")
LL_TKS_G1_lyr <- mutate(LL_TKS_G1_lyr, sample = "LL_TKS_G1_lyr")
LL_TKS_G5_lyr <- mutate(LL_TKS_G5_lyr, sample = "LL_TKS_G5_lyr")

all_hal <- rbind(HM_hal_G1,
                 HM_hal_G4,
                 HM_RS7_G1_hal,
                 HM_RS7_G4_hal,
                 HM_ALK_G1_hal,
                 HM_ALK_G4_hal,
                 HM_TKS_G1_hal,
                 HM_TKS_G5_hal)
all_lyr <- rbind(HM_lyr_G1,
                 HM_lyr_G4,
                 HM_RS7_G1_lyr,
                 HM_RS7_G4_lyr, 
                 HM_ALK_G1_lyr,
                 HM_ALK_G4_lyr,
                 HM_TKS_G1_lyr,
                 HM_TKS_G5_lyr)

LL_all_hal <- rbind(LL_hal_G1,
                    LL_hal_G4,
                    LL_RS7_G1_hal,
                    LL_RS7K_G4_hal,
                    LL_ALK_G1_hal,
                    LL_ALK_G4_hal,
                    LL_TKS_G1_hal,
                    LL_TKS_G5_hal)
LL_all_lyr <- rbind(LL_lyr_G1,
                    LL_lyr_G4,
                    LL_RS7_G1_lyr,
                    LL_RS7K_G4_lyr,
                    LL_ALK_G1_lyr,
                    LL_ALK_G4_lyr,
                    LL_TKS_G1_lyr,
                    LL_TKS_G5_lyr)

## We plot the coverage across 8 chromosomes for RS7 G4

## Plot selected areas

### HM conditions

HM_both_chr1 <- ggplot() +
  #geom_rect(data = nonhomeo_HM_overlap_continuous %>% filter(start_scaffold == "chr1"), aes(xmin = start_coordinate/1e5, xmax = end_coordinate/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "chartreuse4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr1", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr9", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr1", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr9", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr1", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr9", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr1", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr9", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  xlab("Chromosome 1") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr2 <- ggplot() +
  #geom_rect(data = homeo_HM_hal_overlap_continuous %>% filter(seqnames == "chr2"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 25), alpha = 0.5, fill = "brown4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr2", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr10", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr2", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr10", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr2", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr10", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr2", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr10", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  xlab("Chromosome 2") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr3 <- ggplot() +
  #geom_rect(data = homeo_HM_hal_overlap_continuous %>% filter(seqnames == "chr3"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "brown4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr3", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr11", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr3", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr11", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr3", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr11", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr3", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr11", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  xlab("Chromosome 3") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr4 <- ggplot() +
  #geom_rect(data = homeo_HM_hal_overlap_continuous %>% filter(seqnames == "chr4"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "brown4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr4", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr12", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr4", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr12", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr4", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr12", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr4", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr12", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  xlab("Chromosome 4") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr5 <- ggplot() +
  #geom_rect(data = homeo_HM_hal_overlap_continuous %>% filter(seqnames == "chr5"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax =13), alpha = 0.5, fill = "brown4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr5", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr13", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr5", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr13", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr5", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr13", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr5", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr13", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  xlab("Chromosome 5") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr6 <- ggplot() +
  #geom_rect(data = nonhomeo_HM_overlap_continuous %>% filter(start_scaffold == "chr6"), aes(xmin = start_coordinate/1e5, xmax = end_coordinate/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "chartreuse4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr6", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr14", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr6", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr14", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr6", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr14", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr6", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr14", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  xlab("Chromosome 6") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr7 <- ggplot() +
  #geom_rect(data = homeo_HM_hal %>% filter(seqnames == "scaffold_7_RagTag"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 25), alpha = 0.5) +
  geom_line(data = all_hal %>% filter(seqnames == "chr7", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr15", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr7", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr15", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr7", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr15", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr7", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr15", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  xlab("Chromosome 7") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr8 <- ggplot() +
  #geom_rect(data = nonhomeo_HM_overlap_continuous %>% filter(start_scaffold == "scaffold_8_RagTag"), aes(xmin = start_coordinate/1e5, xmax = end_coordinate/1e5, ymin = 0, ymax = 25), alpha = 0.5, fill = "chartreuse4") +
  geom_line(data = all_hal %>% filter(seqnames == "chr8", sample == "HM_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr16", sample == "HM_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr8", sample == "HM_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr16", sample == "HM_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr8", sample == "HM_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr16", sample == "HM_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  geom_line(data = all_hal %>% filter(seqnames == "chr8", sample == "HM_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = all_lyr %>% filter(seqnames == "chr16", sample == "HM_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  xlab("Chromosome 8") +
  theme_linedraw() +
  #theme(legend.position = "none") +
  ylim(0,13)

HM_both_chr1 + HM_both_chr2 +
  HM_both_chr3 + HM_both_chr4 +
  HM_both_chr5 + HM_both_chr6 + 
  HM_both_chr7 + HM_both_chr8 +
  plot_layout(ncol = 3)

### LL conditions

LL_both_chr1 <- ggplot() +
  #geom_rect(data = homeo_LL_hal_overlap_continuous %>% filter(seqnames == "chr1"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "brown4") +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr1", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr9", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr1", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr9", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr1", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr9", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr1", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[1]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr9", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[1]), y = coverage, color = sample)) +
  xlab("Chromosome 1") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr2 <- ggplot() +
  #geom_rect(data = homeo_LL_hal_overlap_continuous %>% filter(seqnames == "chr2"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "brown4") +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr2", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr10", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr2", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr10", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr2", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr10", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr2", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[2]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr10", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[2]), y = coverage, color = sample)) +
  xlab("Chromosome 2") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr3 <- ggplot() +
  #geom_rect(data = homeo_LL_hal_overlap_continuous %>% filter(seqnames == "chr3"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "brown4") +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr3", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr11", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr3", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr11", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr3", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr11", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr3", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[3]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr11", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[3]), y = coverage, color = sample)) +
  xlab("Chromosome 3") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr4 <- ggplot() +
  #geom_rect(data = homeo_LL_hal %>% filter(seqnames == "scaffold_4_RagTag"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 25), alpha = 0.5) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr4", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr12", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr4", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr12", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr4", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr12", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr4", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[4]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr12", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[4]), y = coverage, color = sample)) +
  xlab("Chromosome 4") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr5 <- ggplot() +
  #geom_rect(data = homeo_LL_hal %>% filter(seqnames == "scaffold_5_RagTag"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 25), alpha = 0.5) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr5", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr13", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr5", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr13", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr5", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr13", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr5", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[5]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr13", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[5]), y = coverage, color = sample)) +
  xlab("Chromosome 5") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr6 <- ggplot() +
  #geom_rect(data = nonhomeo_LL_overlap_continuous %>% filter(start_scaffold == "chr6"), aes(xmin = start_coordinate/1e5, xmax = end_coordinate/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "chartreuse4") +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr6", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr14", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr6", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr14", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr6", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr14", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr6", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[6]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr14", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[6]), y = coverage, color = sample)) +
  xlab("Chromosome 6") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr7 <- ggplot() +
  #geom_rect(data = homeo_LL_hal %>% filter(seqnames == "scaffold_7_RagTag"), aes(xmin = start/1e5, xmax = end/1e5, ymin = 0, ymax = 25), alpha = 0.5) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr7", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr15", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr7", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr15", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr7", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr15", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr7", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[7]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr15", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[7]), y = coverage, color = sample)) +
  xlab("Chromosome 7") +
  theme_linedraw() +
  theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr8 <- ggplot() +
  #geom_rect(data = nonhomeo_LL_overlap_continuous %>% filter(start_scaffold == "chr8"), aes(xmin = start_coordinate/1e5, xmax = end_coordinate/1e5, ymin = 0, ymax = 13), alpha = 0.5, fill = "chartreuse4") +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr8", sample == "LL_ALK_G1_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr16", sample == "LL_ALK_G1_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr8", sample == "LL_ALK_G4_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr16", sample == "LL_ALK_G4_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr8", sample == "LL_TKS_G1_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr16", sample == "LL_TKS_G1_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_hal %>% filter(seqnames == "chr8", sample == "LL_TKS_G5_hal"), aes(x = c(1:chrom_length_hal[8]), y = coverage, color = sample)) +
  geom_line(data = LL_all_lyr %>% filter(seqnames == "chr16", sample == "LL_TKS_G5_lyr"), aes(x = c(1:chrom_length_lyr[8]), y = coverage, color = sample)) +
  xlab("Chromosome 8") +
  theme_linedraw() +
  #theme(legend.position = "none") +
  ylim(0,13)

LL_both_chr1 + LL_both_chr2 +
  LL_both_chr3 + LL_both_chr4 +
  LL_both_chr5 + LL_both_chr6 + 
  LL_both_chr7 + LL_both_chr8 +
  plot_layout(ncol = 3)

