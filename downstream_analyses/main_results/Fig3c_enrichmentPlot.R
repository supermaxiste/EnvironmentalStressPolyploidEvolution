### This script will import the gene enrichment output from
### METHImpute and will plot the average methylation level
### within gene bodies and flanking regions (+-500bp)

### Import libraries

library(data.table)
library(tidyverse)
library(scales)
library(patchwork)

### Input data, METHImpute output
setwd("/path/to/data")

### halleri side

HM_hal_G1_1 <- fread("HM_hal_G1_1_dataframe_rcmeth.txt")
HM_hal_G1_2 <- fread("HM_hal_G1_2_dataframe_rcmeth.txt")
HM_hal_G1_3 <- fread("HM_hal_G1_3_dataframe_rcmeth.txt")

HM_hal_G4_1 <- fread("HM_hal_G4_1_dataframe_rcmeth.txt")
HM_hal_G4_2 <- fread("HM_hal_G4_2_dataframe_rcmeth.txt")
HM_hal_G4_3 <- fread("HM_hal_G4_3_dataframe_rcmeth.txt")

HM_RS7_G1_1_hal <- read.delim("HM_RS7K_G1_1_hal_dataframe_rcmeth.txt")
HM_RS7_G1_2_hal <- read.delim("HM_RS7K_G1_2_hal_dataframe_rcmeth.txt")
HM_RS7_G1_3_hal <- read.delim("HM_RS7K_G1_3_hal_dataframe_rcmeth.txt")

HM_RS7_G4_1_hal <- fread("HM_RS7_G4_1_hal_dataframe_rcmeth.txt")
HM_RS7_G4_2_hal <- fread("HM_RS7_G4_2_hal_dataframe_rcmeth.txt")
HM_RS7_G4_3_hal <- fread("HM_RS7_G4_3_hal_dataframe_rcmeth.txt")

HM_ALK_G1_1_hal <- fread("HM_ALK_G1_1_hal_dataframe_rcmeth.txt")
HM_ALK_G1_2_hal <- fread("HM_ALK_G1_2_hal_dataframe_rcmeth.txt")
HM_ALK_G1_3_hal <- fread("HM_ALK_G1_3_hal_dataframe_rcmeth.txt")

HM_ALK_G4_1_hal <- fread("HM_ALK_G4_1_hal_dataframe_rcmeth.txt")
HM_ALK_G4_2_hal <- fread("HM_ALK_G4_2_hal_dataframe_rcmeth.txt")
HM_ALK_G4_3_hal <- fread("HM_ALK_G4_3_hal_dataframe_rcmeth.txt")

HM_TKS_G1_1_hal <- fread("HM_TKS_G1_1_hal_dataframe_rcmeth.txt")
HM_TKS_G1_2_hal <- fread("HM_TKS_G1_2_hal_dataframe_rcmeth.txt")
HM_TKS_G1_3_hal <- fread("HM_TKS_G1_3_hal_dataframe_rcmeth.txt")

HM_TKS_G5_1_hal <- fread("HM_TKS_G5_1_hal_dataframe_rcmeth.txt")
HM_TKS_G5_2_hal <- fread("HM_TKS_G5_2_hal_dataframe_rcmeth.txt")
HM_TKS_G5_3_hal <- fread("HM_TKS_G5_3_hal_dataframe_rcmeth.txt")

LL_hal_G1_1 <- fread("LL_hal_G1_1_dataframe_rcmeth.txt")
LL_hal_G1_2 <- fread("LL_hal_G1_2_dataframe_rcmeth.txt")
LL_hal_G1_3 <- fread("LL_hal_G1_3_dataframe_rcmeth.txt")

LL_hal_G4_1 <- fread("LL_hal_G4_1_dataframe_rcmeth.txt")
LL_hal_G4_2 <- fread("LL_hal_G4_2_dataframe_rcmeth.txt")

LL_RS7_G1_1_hal <- fread("LL_RS7_G1_1_hal_dataframe_rcmeth.txt")
LL_RS7_G1_2_hal <- fread("LL_RS7_G1_2_hal_dataframe_rcmeth.txt")
LL_RS7_G1_3_hal <- fread("LL_RS7_G1_3_hal_dataframe_rcmeth.txt")

LL_RS7_G4_1_hal <- fread("LL_RS7_G4_1_hal_dataframe_rcmeth.txt")
LL_RS7_G4_2_hal <- fread("LL_RS7_G4_2_hal_dataframe_rcmeth.txt")
LL_RS7_G4_3_hal <- fread("LL_RS7_G4_3_hal_dataframe_rcmeth.txt")

LL_ALK_G1_1_hal <- fread("LL_ALK_G1_1_hal_dataframe_rcmeth.txt")
LL_ALK_G1_2_hal <- fread("LL_ALK_G1_2_hal_dataframe_rcmeth.txt")
LL_ALK_G1_3_hal <- fread("LL_ALK_G1_3_hal_dataframe_rcmeth.txt")

LL_ALK_G4_1_hal <- fread("LL_ALK_G4_1_hal_dataframe_rcmeth.txt")
LL_ALK_G4_2_hal <- fread("LL_ALK_G4_2_hal_dataframe_rcmeth.txt")
LL_ALK_G4_3_hal <- fread("LL_ALK_G4_3_hal_dataframe_rcmeth.txt")

LL_TKS_G1_1_hal <- fread("LL_TKS_G1_1_hal_dataframe_rcmeth.txt")
LL_TKS_G1_2_hal <- fread("LL_TKS_G1_2_hal_dataframe_rcmeth.txt")
LL_TKS_G1_3_hal <- fread("LL_TKS_G1_3_hal_dataframe_rcmeth.txt")

LL_TKS_G5_1_hal <- fread("LL_TKS_G5_1_hal_dataframe_rcmeth.txt")
LL_TKS_G5_2_hal <- fread("LL_TKS_G5_2_hal_dataframe_rcmeth.txt")
LL_TKS_G5_3_hal <- fread("LL_TKS_G5_3_hal_dataframe_rcmeth.txt")

### lyrata side

HM_lyr_G1_1 <- fread("HM_lyr_G1_1_dataframe_rcmeth.txt")
HM_lyr_G1_2 <- fread("HM_lyr_G1_2_dataframe_rcmeth.txt")
HM_lyr_G1_3 <- fread("HM_lyr_G1_3_dataframe_rcmeth.txt")

HM_lyr_G4_1 <- fread("HM_lyr_G4_1_dataframe_rcmeth.txt")
HM_lyr_G4_2 <- fread("HM_lyr_G4_2_dataframe_rcmeth.txt")
HM_lyr_G4_3 <- fread("HM_lyr_G4_3_dataframe_rcmeth.txt")

HM_RS7_G1_1_lyr <- fread("HM_RS7K_G1_1_lyr_dataframe_rcmeth.txt")
HM_RS7_G1_2_lyr <- fread("HM_RS7K_G1_2_lyr_dataframe_rcmeth.txt")
HM_RS7_G1_3_lyr <- fread("HM_RS7K_G1_3_lyr_dataframe_rcmeth.txt")

HM_RS7_G4_1_lyr <- fread("HM_RS7_G4_1_lyr_dataframe_rcmeth.txt")
HM_RS7_G4_2_lyr <- fread("HM_RS7_G4_2_lyr_dataframe_rcmeth.txt")
HM_RS7_G4_3_lyr <- fread("HM_RS7_G4_3_lyr_dataframe_rcmeth.txt")

HM_ALK_G1_1_lyr <- fread("HM_ALK_G1_1_lyr_dataframe_rcmeth.txt")
HM_ALK_G1_2_lyr <- fread("HM_ALK_G1_2_lyr_dataframe_rcmeth.txt")
HM_ALK_G1_3_lyr <- fread("HM_ALK_G1_3_lyr_dataframe_rcmeth.txt")

HM_ALK_G4_1_lyr <- fread("HM_ALK_G4_1_lyr_dataframe_rcmeth.txt")
HM_ALK_G4_2_lyr <- fread("HM_ALK_G4_2_lyr_dataframe_rcmeth.txt")
HM_ALK_G4_3_lyr <- fread("HM_ALK_G4_3_lyr_dataframe_rcmeth.txt")

HM_TKS_G1_1_lyr <- fread("HM_TKS_G1_1_lyr_dataframe_rcmeth.txt")
HM_TKS_G1_2_lyr <- fread("HM_TKS_G1_2_lyr_dataframe_rcmeth.txt")
HM_TKS_G1_3_lyr <- fread("HM_TKS_G1_3_lyr_dataframe_rcmeth.txt")

HM_TKS_G5_1_lyr <- fread("HM_TKS_G5_1_lyr_dataframe_rcmeth.txt")
HM_TKS_G5_2_lyr <- fread("HM_TKS_G5_2_lyr_dataframe_rcmeth.txt")
HM_TKS_G5_3_lyr <- fread("HM_TKS_G5_3_lyr_dataframe_rcmeth.txt")

LL_lyr_G1_1 <- fread("LL_lyr_G1_1_dataframe_rcmeth.txt")
LL_lyr_G1_2 <- fread("LL_lyr_G1_2_dataframe_rcmeth.txt")
LL_lyr_G1_3 <- fread("LL_lyr_G1_3_dataframe_rcmeth.txt")

LL_lyr_G4_1 <- fread("LL_lyr_G4_1_dataframe_rcmeth.txt")
LL_lyr_G4_2 <- fread("LL_lyr_G4_2_dataframe_rcmeth.txt")

LL_RS7_G1_1_lyr <- read.delim("LL_RS7_G1_1_lyr_dataframe_rcmeth.txt")
LL_RS7_G1_2_lyr <- read.delim("LL_RS7_G1_2_lyr_dataframe_rcmeth.txt")
LL_RS7_G1_3_lyr <- read.delim("LL_RS7_G1_3_lyr_dataframe_rcmeth.txt")

LL_RS7_G4_1_lyr <- fread("LL_RS7_G4_1_lyr_dataframe_rcmeth.txt")
LL_RS7_G4_2_lyr <- fread("LL_RS7_G4_2_lyr_dataframe_rcmeth.txt")
LL_RS7_G4_3_lyr <- fread("LL_RS7_G4_3_lyr_dataframe_rcmeth.txt")

LL_ALK_G1_1_lyr <- fread("LL_ALK_G1_1_lyr_dataframe_rcmeth.txt")
LL_ALK_G1_2_lyr <- fread("LL_ALK_G1_2_lyr_dataframe_rcmeth.txt")
LL_ALK_G1_3_lyr <- fread("LL_ALK_G1_3_lyr_dataframe_rcmeth.txt")

LL_ALK_G4_1_lyr <- fread("LL_ALK_G4_1_lyr_dataframe_rcmeth.txt")
LL_ALK_G4_2_lyr <- fread("LL_ALK_G4_2_lyr_dataframe_rcmeth.txt")
LL_ALK_G4_3_lyr <- fread("LL_ALK_G4_3_lyr_dataframe_rcmeth.txt")

LL_TKS_G1_1_lyr <- fread("LL_TKS_G1_1_lyr_dataframe_rcmeth.txt")
LL_TKS_G1_2_lyr <- fread("LL_TKS_G1_2_lyr_dataframe_rcmeth.txt")
LL_TKS_G1_3_lyr <- fread("LL_TKS_G1_3_lyr_dataframe_rcmeth.txt")

LL_TKS_G5_1_lyr <- fread("LL_TKS_G5_1_lyr_dataframe_rcmeth.txt")
LL_TKS_G5_2_lyr <- fread("LL_TKS_G5_2_lyr_dataframe_rcmeth.txt")
LL_TKS_G5_3_lyr <- fread("LL_TKS_G5_3_lyr_dataframe_rcmeth.txt")

### Add one column to all imported data with their own origin

HM_hal_G1_1 <- mutate(HM_hal_G1_1, sample = rep("HM_hal_G1_1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("halleri", nrow(HM_hal_G1_1)), generation = "G1")
HM_hal_G1_2 <- mutate(HM_hal_G1_2, sample = rep("HM_hal_G1_2", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("halleri", nrow(HM_hal_G1_1)), generation = "G1")
HM_hal_G1_3 <- mutate(HM_hal_G1_3, sample = rep("HM_hal_G1_3", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("halleri", nrow(HM_hal_G1_1)), generation = "G1")

HM_hal_G4_1 <- mutate(HM_hal_G4_1, sample = rep("HM_hal_G4_1", nrow(HM_hal_G4_1)), side = rep("halleri", nrow(HM_hal_G4_1)), species = rep("halleri", nrow(HM_hal_G4_1)), generation = "G4")
HM_hal_G4_2 <- mutate(HM_hal_G4_2, sample = rep("HM_hal_G4_2", nrow(HM_hal_G4_1)), side = rep("halleri", nrow(HM_hal_G4_1)), species = rep("halleri", nrow(HM_hal_G4_1)), generation = "G4")
HM_hal_G4_3 <- mutate(HM_hal_G4_3, sample = rep("HM_hal_G4_3", nrow(HM_hal_G4_1)), side = rep("halleri", nrow(HM_hal_G4_1)), species = rep("halleri", nrow(HM_hal_G4_1)), generation = "G4")

HM_RS7_G1_1_hal <- mutate(HM_RS7_G1_1_hal, sample = rep("HM_kam_hal_1_g1", nrow(HM_RS7_G1_1_hal)), side = rep("halleri", nrow(HM_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_hal)), generation = "G1")
HM_RS7_G1_2_hal <- mutate(HM_RS7_G1_2_hal, sample = rep("HM_kam_hal_2_g1", nrow(HM_RS7_G1_1_hal)), side = rep("halleri", nrow(HM_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_hal)), generation = "G1")
HM_RS7_G1_3_hal <- mutate(HM_RS7_G1_3_hal, sample = rep("HM_kam_hal_3_g1", nrow(HM_RS7_G1_1_hal)), side = rep("halleri", nrow(HM_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_hal)), generation = "G1")

HM_RS7_G4_1_hal <- mutate(HM_RS7_G4_1_hal, sample = rep("HM_kam_hal_1_g4", nrow(HM_RS7_G1_1_hal)), side = rep("halleri", nrow(HM_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_hal)), generation = "G4")
HM_RS7_G4_2_hal <- mutate(HM_RS7_G4_2_hal, sample = rep("HM_kam_hal_2_g4", nrow(HM_RS7_G1_1_hal)), side = rep("halleri", nrow(HM_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_hal)), generation = "G4")
HM_RS7_G4_3_hal <- mutate(HM_RS7_G4_3_hal, sample = rep("HM_kam_hal_3_g4", nrow(HM_RS7_G1_1_hal)), side = rep("halleri", nrow(HM_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_hal)), generation = "G4")

HM_ALK_G1_1_hal <- mutate(HM_ALK_G1_1_hal, sample = rep("HM_ALK_hal_1_g1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G1")
HM_ALK_G1_2_hal <- mutate(HM_ALK_G1_2_hal, sample = rep("HM_ALK_hal_2_g1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G1")
HM_ALK_G1_3_hal <- mutate(HM_ALK_G1_3_hal, sample = rep("HM_ALK_hal_3_g1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G1")

HM_ALK_G4_1_hal <- mutate(HM_ALK_G4_1_hal, sample = rep("HM_ALK_hal_1_g4", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G4")
HM_ALK_G4_2_hal <- mutate(HM_ALK_G4_2_hal, sample = rep("HM_ALK_hal_2_g4", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G4")
HM_ALK_G4_3_hal <- mutate(HM_ALK_G4_3_hal, sample = rep("HM_ALK_hal_3_g4", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G4")

HM_TKS_G1_1_hal <- mutate(HM_TKS_G1_1_hal, sample = rep("HM_TKS_hal_1_g1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G1")
HM_TKS_G1_2_hal <- mutate(HM_TKS_G1_2_hal, sample = rep("HM_TKS_hal_2_g1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G1")
HM_TKS_G1_3_hal <- mutate(HM_TKS_G1_3_hal, sample = rep("HM_TKS_hal_3_g1", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G1")

HM_TKS_G5_1_hal <- mutate(HM_TKS_G5_1_hal, sample = rep("HM_TKS_hal_1_g5", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G4")
HM_TKS_G5_2_hal <- mutate(HM_TKS_G5_2_hal, sample = rep("HM_TKS_hal_2_g5", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G4")
HM_TKS_G5_3_hal <- mutate(HM_TKS_G5_3_hal, sample = rep("HM_TKS_hal_3_g5", nrow(HM_hal_G1_1)), side = rep("halleri", nrow(HM_hal_G1_1)), species = rep("kamchatica-natural", nrow(HM_hal_G1_1)), generation = "G4")

LL_hal_G1_1 <- mutate(LL_hal_G1_1, sample = rep("LL_hal_G1_1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("halleri", nrow(LL_hal_G1_1)), generation = "G1")
LL_hal_G1_2 <- mutate(LL_hal_G1_2, sample = rep("LL_hal_G1_2", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("halleri", nrow(LL_hal_G1_1)), generation = "G1")
LL_hal_G1_3 <- mutate(LL_hal_G1_3, sample = rep("LL_hal_G1_3", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("halleri", nrow(LL_hal_G1_1)), generation = "G1")

LL_hal_G4_1 <- mutate(LL_hal_G4_1, sample = rep("LL_hal_G4_1", nrow(LL_hal_G4_1)), side = rep("halleri", nrow(LL_hal_G4_1)), species = rep("halleri", nrow(LL_hal_G4_1)), generation = "G4")
LL_hal_G4_2 <- mutate(LL_hal_G4_2, sample = rep("LL_hal_G4_2", nrow(LL_hal_G4_1)), side = rep("halleri", nrow(LL_hal_G4_1)), species = rep("halleri", nrow(LL_hal_G4_1)), generation = "G4")

LL_RS7_G1_1_hal <- mutate(LL_RS7_G1_1_hal, sample = rep("LL_kam_hal_1_g1", nrow(LL_RS7_G1_1_hal)), side = rep("halleri", nrow(LL_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_hal)), generation = "G1")
LL_RS7_G1_2_hal <- mutate(LL_RS7_G1_2_hal, sample = rep("LL_kam_hal_2_g1", nrow(LL_RS7_G1_1_hal)), side = rep("halleri", nrow(LL_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_hal)), generation = "G1")
LL_RS7_G1_3_hal <- mutate(LL_RS7_G1_3_hal, sample = rep("LL_kam_hal_3_g1", nrow(LL_RS7_G1_1_hal)), side = rep("halleri", nrow(LL_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_hal)), generation = "G1")

LL_RS7_G4_1_hal <- mutate(LL_RS7_G4_1_hal, sample = rep("LL_kam_hal_1_g4", nrow(LL_RS7_G1_1_hal)), side = rep("halleri", nrow(LL_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_hal)), generation = "G4")
LL_RS7_G4_2_hal <- mutate(LL_RS7_G4_2_hal, sample = rep("LL_kam_hal_2_g4", nrow(LL_RS7_G1_1_hal)), side = rep("halleri", nrow(LL_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_hal)), generation = "G4")
LL_RS7_G4_3_hal <- mutate(LL_RS7_G4_3_hal, sample = rep("LL_kam_hal_3_g4", nrow(LL_RS7_G1_1_hal)), side = rep("halleri", nrow(LL_RS7_G1_1_hal)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_hal)), generation = "G4")

LL_ALK_G1_1_hal <- mutate(LL_ALK_G1_1_hal, sample = rep("LL_ALK_hal_1_g1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G1")
LL_ALK_G1_2_hal <- mutate(LL_ALK_G1_2_hal, sample = rep("LL_ALK_hal_2_g1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G1")
LL_ALK_G1_3_hal <- mutate(LL_ALK_G1_3_hal, sample = rep("LL_ALK_hal_3_g1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G1")

LL_ALK_G4_1_hal <- mutate(LL_ALK_G4_1_hal, sample = rep("LL_ALK_hal_1_g4", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G4")
LL_ALK_G4_2_hal <- mutate(LL_ALK_G4_2_hal, sample = rep("LL_ALK_hal_2_g4", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G4")
LL_ALK_G4_3_hal <- mutate(LL_ALK_G4_3_hal, sample = rep("LL_ALK_hal_3_g4", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G4")

LL_TKS_G1_1_hal <- mutate(LL_TKS_G1_1_hal, sample = rep("LL_TKS_hal_1_g1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G1")
LL_TKS_G1_2_hal <- mutate(LL_TKS_G1_2_hal, sample = rep("LL_TKS_hal_2_g1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G1")
LL_TKS_G1_3_hal <- mutate(LL_TKS_G1_3_hal, sample = rep("LL_TKS_hal_3_g1", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G1")

LL_TKS_G5_1_hal <- mutate(LL_TKS_G5_1_hal, sample = rep("LL_TKS_hal_1_g5", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G4")
LL_TKS_G5_2_hal <- mutate(LL_TKS_G5_2_hal, sample = rep("LL_TKS_hal_2_g5", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G4")
LL_TKS_G5_3_hal <- mutate(LL_TKS_G5_3_hal, sample = rep("LL_TKS_hal_3_g5", nrow(LL_hal_G1_1)), side = rep("halleri", nrow(LL_hal_G1_1)), species = rep("kamchatica-natural", nrow(LL_hal_G1_1)), generation = "G4")


HM_lyr_G1_1 <- mutate(HM_lyr_G1_1, sample = rep("HM_lyr_G1_1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("lyrata", nrow(HM_lyr_G1_1)), generation = "G1")
HM_lyr_G1_2 <- mutate(HM_lyr_G1_2, sample = rep("HM_lyr_G1_2", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("lyrata", nrow(HM_lyr_G1_1)), generation = "G1")
HM_lyr_G1_3 <- mutate(HM_lyr_G1_3, sample = rep("HM_lyr_G1_3", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("lyrata", nrow(HM_lyr_G1_1)), generation = "G1")

HM_lyr_G4_1 <- mutate(HM_lyr_G4_1, sample = rep("HM_lyr_G4_1", nrow(HM_lyr_G4_1)), side = rep("lyrata", nrow(HM_lyr_G4_1)), species = rep("lyrata", nrow(HM_lyr_G4_1)), generation = "G4")
HM_lyr_G4_2 <- mutate(HM_lyr_G4_2, sample = rep("HM_lyr_G4_2", nrow(HM_lyr_G4_1)), side = rep("lyrata", nrow(HM_lyr_G4_1)), species = rep("lyrata", nrow(HM_lyr_G4_1)), generation = "G4")
HM_lyr_G4_3 <- mutate(HM_lyr_G4_3, sample = rep("HM_lyr_G4_3", nrow(HM_lyr_G4_1)), side = rep("lyrata", nrow(HM_lyr_G4_1)), species = rep("lyrata", nrow(HM_lyr_G4_1)), generation = "G4")

HM_RS7_G1_1_lyr <- mutate(HM_RS7_G1_1_lyr, sample = rep("HM_kam_lyr_1_g1", nrow(HM_RS7_G1_1_lyr)), side = rep("lyrata", nrow(HM_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_lyr)), generation = "G1")
HM_RS7_G1_2_lyr <- mutate(HM_RS7_G1_2_lyr, sample = rep("HM_kam_lyr_2_g1", nrow(HM_RS7_G1_1_lyr)), side = rep("lyrata", nrow(HM_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_lyr)), generation = "G1")
HM_RS7_G1_3_lyr <- mutate(HM_RS7_G1_3_lyr, sample = rep("HM_kam_lyr_3_g1", nrow(HM_RS7_G1_1_lyr)), side = rep("lyrata", nrow(HM_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_lyr)), generation = "G1")

HM_RS7_G4_1_lyr <- mutate(HM_RS7_G4_1_lyr, sample = rep("HM_kam_lyr_1_g4", nrow(HM_RS7_G1_1_lyr)), side = rep("lyrata", nrow(HM_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_lyr)), generation = "G4")
HM_RS7_G4_2_lyr <- mutate(HM_RS7_G4_2_lyr, sample = rep("HM_kam_lyr_2_g4", nrow(HM_RS7_G1_1_lyr)), side = rep("lyrata", nrow(HM_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_lyr)), generation = "G4")
HM_RS7_G4_3_lyr <- mutate(HM_RS7_G4_3_lyr, sample = rep("HM_kam_lyr_3_g4", nrow(HM_RS7_G1_1_lyr)), side = rep("lyrata", nrow(HM_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(HM_RS7_G1_1_lyr)), generation = "G4")

HM_ALK_G1_1_lyr <- mutate(HM_ALK_G1_1_lyr, sample = rep("HM_ALK_lyr_1_g1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G1")
HM_ALK_G1_2_lyr <- mutate(HM_ALK_G1_2_lyr, sample = rep("HM_ALK_lyr_2_g1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G1")
HM_ALK_G1_3_lyr <- mutate(HM_ALK_G1_3_lyr, sample = rep("HM_ALK_lyr_3_g1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G1")

HM_ALK_G4_1_lyr <- mutate(HM_ALK_G4_1_lyr, sample = rep("HM_ALK_lyr_1_g4", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G4")
HM_ALK_G4_2_lyr <- mutate(HM_ALK_G4_2_lyr, sample = rep("HM_ALK_lyr_2_g4", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G4")
HM_ALK_G4_3_lyr <- mutate(HM_ALK_G4_3_lyr, sample = rep("HM_ALK_lyr_3_g4", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G4")

HM_TKS_G1_1_lyr <- mutate(HM_TKS_G1_1_lyr, sample = rep("HM_TKS_lyr_1_g1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G1")
HM_TKS_G1_2_lyr <- mutate(HM_TKS_G1_2_lyr, sample = rep("HM_TKS_lyr_2_g1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G1")
HM_TKS_G1_3_lyr <- mutate(HM_TKS_G1_3_lyr, sample = rep("HM_TKS_lyr_3_g1", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G1")

HM_TKS_G5_1_lyr <- mutate(HM_TKS_G5_1_lyr, sample = rep("HM_TKS_lyr_1_g5", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G4")
HM_TKS_G5_2_lyr <- mutate(HM_TKS_G5_2_lyr, sample = rep("HM_TKS_lyr_2_g5", nrow(HM_lyr_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G4")
HM_TKS_G5_3_lyr <- mutate(HM_TKS_G5_3_lyr, sample = rep("HM_TKS_lyr_3_g5", nrow(HM_hal_G1_1)), side = rep("lyrata", nrow(HM_lyr_G1_1)), species = rep("kamchatica-natural", nrow(HM_lyr_G1_1)), generation = "G4")

LL_lyr_G1_1 <- mutate(LL_lyr_G1_1, sample = rep("LL_lyr_G1_1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("lyrata", nrow(LL_lyr_G1_1)), generation = "G1")
LL_lyr_G1_2 <- mutate(LL_lyr_G1_2, sample = rep("LL_lyr_G1_2", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("lyrata", nrow(LL_lyr_G1_1)), generation = "G1")
LL_lyr_G1_3 <- mutate(LL_lyr_G1_3, sample = rep("LL_lyr_G1_3", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("lyrata", nrow(LL_lyr_G1_1)), generation = "G1")

LL_lyr_G4_1 <- mutate(LL_lyr_G4_1, sample = rep("LL_lyr_G4_1", nrow(LL_lyr_G4_1)), side = rep("lyrata", nrow(LL_lyr_G4_1)), species = rep("lyrata", nrow(LL_lyr_G4_1)), generation = "G4")
LL_lyr_G4_2 <- mutate(LL_lyr_G4_2, sample = rep("LL_lyr_G4_2", nrow(LL_lyr_G4_1)), side = rep("lyrata", nrow(LL_lyr_G4_1)), species = rep("lyrata", nrow(LL_lyr_G4_1)), generation = "G4")

LL_RS7_G1_1_lyr <- mutate(LL_RS7_G1_1_lyr, sample = rep("LL_kam_lyr_1_g1", nrow(LL_RS7_G1_1_lyr)), side = rep("lyrata", nrow(LL_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_lyr)), generation = "G1")
LL_RS7_G1_2_lyr <- mutate(LL_RS7_G1_2_lyr, sample = rep("LL_kam_lyr_2_g1", nrow(LL_RS7_G1_1_lyr)), side = rep("lyrata", nrow(LL_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_lyr)), generation = "G1")
LL_RS7_G1_3_lyr <- mutate(LL_RS7_G1_3_lyr, sample = rep("LL_kam_lyr_3_g1", nrow(LL_RS7_G1_1_lyr)), side = rep("lyrata", nrow(LL_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_lyr)), generation = "G1")

LL_RS7_G4_1_lyr <- mutate(LL_RS7_G4_1_lyr, sample = rep("LL_kam_lyr_1_g4", nrow(LL_RS7_G1_1_lyr)), side = rep("lyrata", nrow(LL_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_lyr)), generation = "G4")
LL_RS7_G4_2_lyr <- mutate(LL_RS7_G4_2_lyr, sample = rep("LL_kam_lyr_2_g4", nrow(LL_RS7_G1_1_lyr)), side = rep("lyrata", nrow(LL_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_lyr)), generation = "G4")
LL_RS7_G4_3_lyr <- mutate(LL_RS7_G4_3_lyr, sample = rep("LL_kam_lyr_3_g4", nrow(LL_RS7_G1_1_lyr)), side = rep("lyrata", nrow(LL_RS7_G1_1_lyr)), species = rep("kamchatica-syn", nrow(LL_RS7_G1_1_lyr)), generation = "G4")

LL_ALK_G1_1_lyr <- mutate(LL_ALK_G1_1_lyr, sample = rep("LL_ALK_lyr_1_g1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G1")
LL_ALK_G1_2_lyr <- mutate(LL_ALK_G1_2_lyr, sample = rep("LL_ALK_lyr_2_g1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G1")
LL_ALK_G1_3_lyr <- mutate(LL_ALK_G1_3_lyr, sample = rep("LL_ALK_lyr_3_g1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G1")

LL_ALK_G4_1_lyr <- mutate(LL_ALK_G4_1_lyr, sample = rep("LL_ALK_lyr_1_g4", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G4")
LL_ALK_G4_2_lyr <- mutate(LL_ALK_G4_2_lyr, sample = rep("LL_ALK_lyr_2_g4", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G4")
LL_ALK_G4_3_lyr <- mutate(LL_ALK_G4_3_lyr, sample = rep("LL_ALK_lyr_3_g4", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G4")

LL_TKS_G1_1_lyr <- mutate(LL_TKS_G1_1_lyr, sample = rep("LL_TKS_lyr_1_g1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G1")
LL_TKS_G1_2_lyr <- mutate(LL_TKS_G1_2_lyr, sample = rep("LL_TKS_lyr_2_g1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G1")
LL_TKS_G1_3_lyr <- mutate(LL_TKS_G1_3_lyr, sample = rep("LL_TKS_lyr_3_g1", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G1")

LL_TKS_G5_1_lyr <- mutate(LL_TKS_G5_1_lyr, sample = rep("LL_TKS_lyr_1_g5", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G4")
LL_TKS_G5_2_lyr <- mutate(LL_TKS_G5_2_lyr, sample = rep("LL_TKS_lyr_2_g5", nrow(LL_lyr_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G4")
LL_TKS_G5_3_lyr <- mutate(LL_TKS_G5_3_lyr, sample = rep("LL_TKS_lyr_3_g5", nrow(LL_hal_G1_1)), side = rep("lyrata", nrow(LL_lyr_G1_1)), species = rep("kamchatica-natural", nrow(LL_lyr_G1_1)), generation = "G4")

### Create dataframes with only one context for each of the two sides.
### We use the part of the METHImpute output that's NOT imputed (i.e. covered)

HM_CG_hal_all <- rbind(HM_hal_G1_1[HM_hal_G1_1$context=="CG",],
                HM_hal_G1_2[HM_hal_G1_2$context=="CG",],
                HM_hal_G1_3[HM_hal_G1_3$context=="CG",],
                HM_hal_G4_1[HM_hal_G4_1$context=="CG",],
                HM_hal_G4_2[HM_hal_G4_2$context=="CG",],
                HM_hal_G4_3[HM_hal_G4_3$context=="CG",],
                HM_RS7_G1_1_hal[HM_RS7_G1_1_hal$context=="CG",],
                HM_RS7_G1_2_hal[HM_RS7_G1_2_hal$context=="CG",],
                HM_RS7_G1_3_hal[HM_RS7_G1_3_hal$context=="CG",],
                HM_RS7_G4_1_hal[HM_RS7_G4_1_hal$context=="CG",],
                HM_RS7_G4_2_hal[HM_RS7_G4_2_hal$context=="CG",],
                HM_RS7_G4_3_hal[HM_RS7_G4_3_hal$context=="CG",],
                HM_ALK_G1_1_hal[HM_ALK_G1_1_hal$context=="CG",],
                HM_ALK_G1_2_hal[HM_ALK_G1_2_hal$context=="CG",],
                HM_ALK_G1_3_hal[HM_ALK_G1_3_hal$context=="CG",],
                HM_ALK_G4_1_hal[HM_ALK_G4_1_hal$context=="CG",],
                HM_ALK_G4_2_hal[HM_ALK_G4_2_hal$context=="CG",],
                HM_ALK_G4_3_hal[HM_ALK_G4_3_hal$context=="CG",],
                HM_TKS_G1_1_hal[HM_TKS_G1_1_hal$context=="CG",],
                HM_TKS_G1_2_hal[HM_TKS_G1_2_hal$context=="CG",],
                HM_TKS_G1_3_hal[HM_TKS_G1_3_hal$context=="CG",],
                HM_TKS_G5_1_hal[HM_TKS_G5_1_hal$context=="CG",],
                HM_TKS_G5_2_hal[HM_TKS_G5_2_hal$context=="CG",],
                HM_TKS_G5_3_hal[HM_TKS_G5_3_hal$context=="CG",])

HM_CHG_hal_all <- rbind(HM_hal_G1_1[HM_hal_G1_1$context=="CHG",],
                     HM_hal_G1_2[HM_hal_G1_2$context=="CHG",],
                     HM_hal_G1_3[HM_hal_G1_3$context=="CHG",],
                     HM_hal_G4_1[HM_hal_G4_1$context=="CHG",],
                     HM_hal_G4_2[HM_hal_G4_2$context=="CHG",],
                     HM_hal_G4_3[HM_hal_G4_3$context=="CHG",],
                     HM_RS7_G1_1_hal[HM_RS7_G1_1_hal$context=="CHG",],
                     HM_RS7_G1_2_hal[HM_RS7_G1_2_hal$context=="CHG",],
                     HM_RS7_G1_3_hal[HM_RS7_G1_3_hal$context=="CHG",],
                     HM_RS7_G4_1_hal[HM_RS7_G4_1_hal$context=="CHG",],
                     HM_RS7_G4_2_hal[HM_RS7_G4_2_hal$context=="CHG",],
                     HM_RS7_G4_3_hal[HM_RS7_G4_3_hal$context=="CHG",],
                     HM_ALK_G1_1_hal[HM_ALK_G1_1_hal$context=="CHG",],
                     HM_ALK_G1_2_hal[HM_ALK_G1_2_hal$context=="CHG",],
                     HM_ALK_G1_3_hal[HM_ALK_G1_3_hal$context=="CHG",],
                     HM_ALK_G4_1_hal[HM_ALK_G4_1_hal$context=="CHG",],
                     HM_ALK_G4_2_hal[HM_ALK_G4_2_hal$context=="CHG",],
                     HM_ALK_G4_3_hal[HM_ALK_G4_3_hal$context=="CHG",],
                     HM_TKS_G1_1_hal[HM_TKS_G1_1_hal$context=="CHG",],
                     HM_TKS_G1_2_hal[HM_TKS_G1_2_hal$context=="CHG",],
                     HM_TKS_G1_3_hal[HM_TKS_G1_3_hal$context=="CHG",],
                     HM_TKS_G5_1_hal[HM_TKS_G5_1_hal$context=="CHG",],
                     HM_TKS_G5_2_hal[HM_TKS_G5_2_hal$context=="CHG",],
                     HM_TKS_G5_3_hal[HM_TKS_G5_3_hal$context=="CHG",])

HM_CHH_hal_all <- rbind(HM_hal_G1_1[HM_hal_G1_1$context=="CHH",],
                     HM_hal_G1_2[HM_hal_G1_2$context=="CHH",],
                     HM_hal_G1_3[HM_hal_G1_3$context=="CHH",],
                     HM_hal_G4_1[HM_hal_G4_1$context=="CHH",],
                     HM_hal_G4_2[HM_hal_G4_2$context=="CHH",],
                     HM_hal_G4_3[HM_hal_G4_3$context=="CHH",],
                     HM_RS7_G1_1_hal[HM_RS7_G1_1_hal$context=="CHH",],
                     HM_RS7_G1_2_hal[HM_RS7_G1_2_hal$context=="CHH",],
                     HM_RS7_G1_3_hal[HM_RS7_G1_3_hal$context=="CHH",],
                     HM_RS7_G4_1_hal[HM_RS7_G4_1_hal$context=="CHH",],
                     HM_RS7_G4_2_hal[HM_RS7_G4_2_hal$context=="CHH",],
                     HM_RS7_G4_3_hal[HM_RS7_G4_3_hal$context=="CHH",],
                     HM_ALK_G1_1_hal[HM_ALK_G1_1_hal$context=="CHH",],
                     HM_ALK_G1_2_hal[HM_ALK_G1_2_hal$context=="CHH",],
                     HM_ALK_G1_3_hal[HM_ALK_G1_3_hal$context=="CHH",],
                     HM_ALK_G4_1_hal[HM_ALK_G4_1_hal$context=="CHH",],
                     HM_ALK_G4_2_hal[HM_ALK_G4_2_hal$context=="CHH",],
                     HM_ALK_G4_3_hal[HM_ALK_G4_3_hal$context=="CHH",],
                     HM_TKS_G1_1_hal[HM_TKS_G1_1_hal$context=="CHH",],
                     HM_TKS_G1_2_hal[HM_TKS_G1_2_hal$context=="CHH",],
                     HM_TKS_G1_3_hal[HM_TKS_G1_3_hal$context=="CHH",],
                     HM_TKS_G5_1_hal[HM_TKS_G5_1_hal$context=="CHH",],
                     HM_TKS_G5_2_hal[HM_TKS_G5_2_hal$context=="CHH",],
                     HM_TKS_G5_3_hal[HM_TKS_G5_3_hal$context=="CHH",])

LL_CG_hal_all <- rbind(LL_hal_G1_1[LL_hal_G1_1$context=="CG",],
                       LL_hal_G1_2[LL_hal_G1_2$context=="CG",],
                       LL_hal_G1_3[LL_hal_G1_3$context=="CG",],
                       LL_hal_G4_1[LL_hal_G4_1$context=="CG",],
                       LL_hal_G4_2[LL_hal_G4_2$context=="CG",],
                       LL_RS7_G1_1_hal[LL_RS7_G1_1_hal$context=="CG",],
                       LL_RS7_G1_2_hal[LL_RS7_G1_2_hal$context=="CG",],
                       LL_RS7_G1_3_hal[LL_RS7_G1_3_hal$context=="CG",],
                       LL_RS7_G4_1_hal[LL_RS7_G4_1_hal$context=="CG",],
                       LL_RS7_G4_2_hal[LL_RS7_G4_2_hal$context=="CG",],
                       LL_RS7_G4_3_hal[LL_RS7_G4_3_hal$context=="CG",],
                       LL_ALK_G1_1_hal[LL_ALK_G1_1_hal$context=="CG",],
                       LL_ALK_G1_2_hal[LL_ALK_G1_2_hal$context=="CG",],
                       LL_ALK_G1_3_hal[LL_ALK_G1_3_hal$context=="CG",],
                       LL_ALK_G4_1_hal[LL_ALK_G4_1_hal$context=="CG",],
                       LL_ALK_G4_2_hal[LL_ALK_G4_2_hal$context=="CG",],
                       LL_ALK_G4_3_hal[LL_ALK_G4_3_hal$context=="CG",],
                       LL_TKS_G1_1_hal[LL_TKS_G1_1_hal$context=="CG",],
                       LL_TKS_G1_2_hal[LL_TKS_G1_2_hal$context=="CG",],
                       LL_TKS_G1_3_hal[LL_TKS_G1_3_hal$context=="CG",],
                       LL_TKS_G5_1_hal[LL_TKS_G5_1_hal$context=="CG",],
                       LL_TKS_G5_2_hal[LL_TKS_G5_2_hal$context=="CG",],
                       LL_TKS_G5_3_hal[LL_TKS_G5_3_hal$context=="CG",])

LL_CHG_hal_all <- rbind(LL_hal_G1_1[LL_hal_G1_1$context=="CHG",],
                        LL_hal_G1_2[LL_hal_G1_2$context=="CHG",],
                        LL_hal_G1_3[LL_hal_G1_3$context=="CHG",],
                        LL_hal_G4_1[LL_hal_G4_1$context=="CHG",],
                        LL_hal_G4_2[LL_hal_G4_2$context=="CHG",],
                        LL_RS7_G1_1_hal[LL_RS7_G1_1_hal$context=="CHG",],
                        LL_RS7_G1_2_hal[LL_RS7_G1_2_hal$context=="CHG",],
                        LL_RS7_G1_3_hal[LL_RS7_G1_3_hal$context=="CHG",],
                        LL_RS7_G4_1_hal[LL_RS7_G4_1_hal$context=="CHG",],
                        LL_RS7_G4_2_hal[LL_RS7_G4_2_hal$context=="CHG",],
                        LL_RS7_G4_3_hal[LL_RS7_G4_3_hal$context=="CHG",],
                        LL_ALK_G1_1_hal[LL_ALK_G1_1_hal$context=="CHG",],
                        LL_ALK_G1_2_hal[LL_ALK_G1_2_hal$context=="CHG",],
                        LL_ALK_G1_3_hal[LL_ALK_G1_3_hal$context=="CHG",],
                        LL_ALK_G4_1_hal[LL_ALK_G4_1_hal$context=="CHG",],
                        LL_ALK_G4_2_hal[LL_ALK_G4_2_hal$context=="CHG",],
                        LL_ALK_G4_3_hal[LL_ALK_G4_3_hal$context=="CHG",],
                        LL_TKS_G1_1_hal[LL_TKS_G1_1_hal$context=="CHG",],
                        LL_TKS_G1_2_hal[LL_TKS_G1_2_hal$context=="CHG",],
                        LL_TKS_G1_3_hal[LL_TKS_G1_3_hal$context=="CHG",],
                        LL_TKS_G5_1_hal[LL_TKS_G5_1_hal$context=="CHG",],
                        LL_TKS_G5_2_hal[LL_TKS_G5_2_hal$context=="CHG",],
                        LL_TKS_G5_3_hal[LL_TKS_G5_3_hal$context=="CHG",])

LL_CHH_hal_all <- rbind(LL_hal_G1_1[LL_hal_G1_1$context=="CHH",],
                        LL_hal_G1_2[LL_hal_G1_2$context=="CHH",],
                        LL_hal_G1_3[LL_hal_G1_3$context=="CHH",],
                        LL_hal_G4_1[LL_hal_G4_1$context=="CHH",],
                        LL_hal_G4_2[LL_hal_G4_2$context=="CHH",],
                        LL_RS7_G1_1_hal[LL_RS7_G1_1_hal$context=="CHH",],
                        LL_RS7_G1_2_hal[LL_RS7_G1_2_hal$context=="CHH",],
                        LL_RS7_G1_3_hal[LL_RS7_G1_3_hal$context=="CHH",],
                        LL_RS7_G4_1_hal[LL_RS7_G4_1_hal$context=="CHH",],
                        LL_RS7_G4_2_hal[LL_RS7_G4_2_hal$context=="CHH",],
                        LL_RS7_G4_3_hal[LL_RS7_G4_3_hal$context=="CHH",],
                        LL_ALK_G1_1_hal[LL_ALK_G1_1_hal$context=="CHH",],
                        LL_ALK_G1_2_hal[LL_ALK_G1_2_hal$context=="CHH",],
                        LL_ALK_G1_3_hal[LL_ALK_G1_3_hal$context=="CHH",],
                        LL_ALK_G4_1_hal[LL_ALK_G4_1_hal$context=="CHH",],
                        LL_ALK_G4_2_hal[LL_ALK_G4_2_hal$context=="CHH",],
                        LL_ALK_G4_3_hal[LL_ALK_G4_3_hal$context=="CHH",],
                        LL_TKS_G1_1_hal[LL_TKS_G1_1_hal$context=="CHH",],
                        LL_TKS_G1_2_hal[LL_TKS_G1_2_hal$context=="CHH",],
                        LL_TKS_G1_3_hal[LL_TKS_G1_3_hal$context=="CHH",],
                        LL_TKS_G5_1_hal[LL_TKS_G5_1_hal$context=="CHH",],
                        LL_TKS_G5_2_hal[LL_TKS_G5_2_hal$context=="CHH",],
                        LL_TKS_G5_3_hal[LL_TKS_G5_3_hal$context=="CHH",])


HM_CG_lyr_all <- rbind(HM_lyr_G1_1[HM_lyr_G1_1$context=="CG",],
                    HM_lyr_G1_2[HM_lyr_G1_2$context=="CG",],
                    HM_lyr_G1_3[HM_lyr_G1_3$context=="CG",],
                    HM_lyr_G4_1[HM_lyr_G4_1$context=="CG",],
                    HM_lyr_G4_2[HM_lyr_G4_2$context=="CG",],
                    HM_lyr_G4_3[HM_lyr_G4_3$context=="CG",],
                    HM_RS7_G1_1_lyr[HM_RS7_G1_1_lyr$context=="CG",],
                    HM_RS7_G1_2_lyr[HM_RS7_G1_2_lyr$context=="CG",],
                    HM_RS7_G1_3_lyr[HM_RS7_G1_3_lyr$context=="CG",],
                    HM_RS7_G4_1_lyr[HM_RS7_G4_1_lyr$context=="CG",],
                    HM_RS7_G4_2_lyr[HM_RS7_G4_2_lyr$context=="CG",],
                    HM_RS7_G4_3_lyr[HM_RS7_G4_3_lyr$context=="CG",],
                    HM_ALK_G1_1_lyr[HM_ALK_G1_1_lyr$context=="CG",],
                    HM_ALK_G1_2_lyr[HM_ALK_G1_2_lyr$context=="CG",],
                    HM_ALK_G1_3_lyr[HM_ALK_G1_3_lyr$context=="CG",],
                    HM_ALK_G4_1_lyr[HM_ALK_G4_1_lyr$context=="CG",],
                    HM_ALK_G4_2_lyr[HM_ALK_G4_2_lyr$context=="CG",],
                    HM_ALK_G4_3_lyr[HM_ALK_G4_3_lyr$context=="CG",],
                    HM_TKS_G1_1_lyr[HM_TKS_G1_1_lyr$context=="CG",],
                    HM_TKS_G1_2_lyr[HM_TKS_G1_2_lyr$context=="CG",],
                    HM_TKS_G1_3_lyr[HM_TKS_G1_3_lyr$context=="CG",],
                    HM_TKS_G5_1_lyr[HM_TKS_G5_1_lyr$context=="CG",],
                    HM_TKS_G5_2_lyr[HM_TKS_G5_2_lyr$context=="CG",],
                    HM_TKS_G5_3_lyr[HM_TKS_G5_3_lyr$context=="CG",])

HM_CHG_lyr_all <- rbind(HM_lyr_G1_1[HM_lyr_G1_1$context=="CHG",],
                     HM_lyr_G1_2[HM_lyr_G1_2$context=="CHG",],
                     HM_lyr_G1_3[HM_lyr_G1_3$context=="CHG",],
                     HM_lyr_G4_1[HM_lyr_G4_1$context=="CHG",],
                     HM_lyr_G4_2[HM_lyr_G4_2$context=="CHG",],
                     HM_lyr_G4_3[HM_lyr_G4_3$context=="CHG",],
                     HM_RS7_G1_1_lyr[HM_RS7_G1_1_lyr$context=="CHG",],
                     HM_RS7_G1_2_lyr[HM_RS7_G1_2_lyr$context=="CHG",],
                     HM_RS7_G1_3_lyr[HM_RS7_G1_3_lyr$context=="CHG",],
                     HM_RS7_G4_1_lyr[HM_RS7_G4_1_lyr$context=="CHG",],
                     HM_RS7_G4_2_lyr[HM_RS7_G4_2_lyr$context=="CHG",],
                     HM_RS7_G4_3_lyr[HM_RS7_G4_3_lyr$context=="CHG",],
                     HM_ALK_G1_1_lyr[HM_ALK_G1_1_lyr$context=="CHG",],
                     HM_ALK_G1_2_lyr[HM_ALK_G1_2_lyr$context=="CHG",],
                     HM_ALK_G1_3_lyr[HM_ALK_G1_3_lyr$context=="CHG",],
                     HM_ALK_G4_1_lyr[HM_ALK_G4_1_lyr$context=="CHG",],
                     HM_ALK_G4_2_lyr[HM_ALK_G4_2_lyr$context=="CHG",],
                     HM_ALK_G4_3_lyr[HM_ALK_G4_3_lyr$context=="CHG",],
                     HM_TKS_G1_1_lyr[HM_TKS_G1_1_lyr$context=="CHG",],
                     HM_TKS_G1_2_lyr[HM_TKS_G1_2_lyr$context=="CHG",],
                     HM_TKS_G1_3_lyr[HM_TKS_G1_3_lyr$context=="CHG",],
                     HM_TKS_G5_1_lyr[HM_TKS_G5_1_lyr$context=="CHG",],
                     HM_TKS_G5_2_lyr[HM_TKS_G5_2_lyr$context=="CHG",],
                     HM_TKS_G5_3_lyr[HM_TKS_G5_3_lyr$context=="CHG",])

HM_CHH_lyr_all <- rbind(HM_lyr_G1_1[HM_lyr_G1_1$context=="CHH",],
                     HM_lyr_G1_2[HM_lyr_G1_2$context=="CHH",],
                     HM_lyr_G1_3[HM_lyr_G1_3$context=="CHH",],
                     HM_lyr_G4_1[HM_lyr_G4_1$context=="CHH",],
                     HM_lyr_G4_2[HM_lyr_G4_2$context=="CHH",],
                     HM_lyr_G4_3[HM_lyr_G4_3$context=="CHH",],
                     HM_RS7_G1_1_lyr[HM_RS7_G1_1_lyr$context=="CHH",],
                     HM_RS7_G1_2_lyr[HM_RS7_G1_2_lyr$context=="CHH",],
                     HM_RS7_G1_3_lyr[HM_RS7_G1_3_lyr$context=="CHH",],
                     HM_RS7_G4_1_lyr[HM_RS7_G4_1_lyr$context=="CHH",],
                     HM_RS7_G4_2_lyr[HM_RS7_G4_2_lyr$context=="CHH",],
                     HM_RS7_G4_3_lyr[HM_RS7_G4_3_lyr$context=="CHH",],
                     HM_ALK_G1_1_lyr[HM_ALK_G1_1_lyr$context=="CHH",],
                     HM_ALK_G1_2_lyr[HM_ALK_G1_2_lyr$context=="CHH",],
                     HM_ALK_G1_3_lyr[HM_ALK_G1_3_lyr$context=="CHH",],
                     HM_ALK_G4_1_lyr[HM_ALK_G4_1_lyr$context=="CHH",],
                     HM_ALK_G4_2_lyr[HM_ALK_G4_2_lyr$context=="CHH",],
                     HM_ALK_G4_3_lyr[HM_ALK_G4_3_lyr$context=="CHH",],
                     HM_TKS_G1_1_lyr[HM_TKS_G1_1_lyr$context=="CHH",],
                     HM_TKS_G1_2_lyr[HM_TKS_G1_2_lyr$context=="CHH",],
                     HM_TKS_G1_3_lyr[HM_TKS_G1_3_lyr$context=="CHH",],
                     HM_TKS_G5_1_lyr[HM_TKS_G5_1_lyr$context=="CHH",],
                     HM_TKS_G5_2_lyr[HM_TKS_G5_2_lyr$context=="CHH",],
                     HM_TKS_G5_3_lyr[HM_TKS_G5_3_lyr$context=="CHH",])

LL_CG_lyr_all <- rbind(LL_lyr_G1_1[LL_lyr_G1_1$context=="CG",],
                       LL_lyr_G1_2[LL_lyr_G1_2$context=="CG",],
                       LL_lyr_G1_3[LL_lyr_G1_3$context=="CG",],
                       LL_lyr_G4_1[LL_lyr_G4_1$context=="CG",],
                       LL_lyr_G4_2[LL_lyr_G4_2$context=="CG",],
                       LL_RS7_G1_1_lyr[LL_RS7_G1_1_lyr$context=="CG",],
                       LL_RS7_G1_2_lyr[LL_RS7_G1_2_lyr$context=="CG",],
                       LL_RS7_G1_3_lyr[LL_RS7_G1_3_lyr$context=="CG",],
                       LL_RS7_G4_1_lyr[LL_RS7_G4_1_lyr$context=="CG",],
                       LL_RS7_G4_2_lyr[LL_RS7_G4_2_lyr$context=="CG",],
                       LL_RS7_G4_3_lyr[LL_RS7_G4_3_lyr$context=="CG",],
                       LL_ALK_G1_1_lyr[LL_ALK_G1_1_lyr$context=="CG",],
                       LL_ALK_G1_2_lyr[LL_ALK_G1_2_lyr$context=="CG",],
                       LL_ALK_G1_3_lyr[LL_ALK_G1_3_lyr$context=="CG",],
                       LL_ALK_G4_1_lyr[LL_ALK_G4_1_lyr$context=="CG",],
                       LL_ALK_G4_2_lyr[LL_ALK_G4_2_lyr$context=="CG",],
                       LL_ALK_G4_3_lyr[LL_ALK_G4_3_lyr$context=="CG",],
                       LL_TKS_G1_1_lyr[LL_TKS_G1_1_lyr$context=="CG",],
                       LL_TKS_G1_2_lyr[LL_TKS_G1_2_lyr$context=="CG",],
                       LL_TKS_G1_3_lyr[LL_TKS_G1_3_lyr$context=="CG",],
                       LL_TKS_G5_1_lyr[LL_TKS_G5_1_lyr$context=="CG",],
                       LL_TKS_G5_2_lyr[LL_TKS_G5_2_lyr$context=="CG",],
                       LL_TKS_G5_3_lyr[LL_TKS_G5_3_lyr$context=="CG",])

LL_CHG_lyr_all <- rbind(LL_lyr_G1_1[LL_lyr_G1_1$context=="CHG",],
                        LL_lyr_G1_2[LL_lyr_G1_2$context=="CHG",],
                        LL_lyr_G1_3[LL_lyr_G1_3$context=="CHG",],
                        LL_lyr_G4_1[LL_lyr_G4_1$context=="CHG",],
                        LL_lyr_G4_2[LL_lyr_G4_2$context=="CHG",],
                        LL_RS7_G1_1_lyr[LL_RS7_G1_1_lyr$context=="CHG",],
                        LL_RS7_G1_2_lyr[LL_RS7_G1_2_lyr$context=="CHG",],
                        LL_RS7_G1_3_lyr[LL_RS7_G1_3_lyr$context=="CHG",],
                        LL_RS7_G4_1_lyr[LL_RS7_G4_1_lyr$context=="CHG",],
                        LL_RS7_G4_2_lyr[LL_RS7_G4_2_lyr$context=="CHG",],
                        LL_RS7_G4_3_lyr[LL_RS7_G4_3_lyr$context=="CHG",],
                        LL_ALK_G1_1_lyr[LL_ALK_G1_1_lyr$context=="CHG",],
                        LL_ALK_G1_2_lyr[LL_ALK_G1_2_lyr$context=="CHG",],
                        LL_ALK_G1_3_lyr[LL_ALK_G1_3_lyr$context=="CHG",],
                        LL_ALK_G4_1_lyr[LL_ALK_G4_1_lyr$context=="CHG",],
                        LL_ALK_G4_2_lyr[LL_ALK_G4_2_lyr$context=="CHG",],
                        LL_ALK_G4_3_lyr[LL_ALK_G4_3_lyr$context=="CHG",],
                        LL_TKS_G1_1_lyr[LL_TKS_G1_1_lyr$context=="CHG",],
                        LL_TKS_G1_2_lyr[LL_TKS_G1_2_lyr$context=="CHG",],
                        LL_TKS_G1_3_lyr[LL_TKS_G1_3_lyr$context=="CHG",],
                        LL_TKS_G5_1_lyr[LL_TKS_G5_1_lyr$context=="CHG",],
                        LL_TKS_G5_2_lyr[LL_TKS_G5_2_lyr$context=="CHG",],
                        LL_TKS_G5_3_lyr[LL_TKS_G5_3_lyr$context=="CHG",])

LL_CHH_lyr_all <- rbind(LL_lyr_G1_1[LL_lyr_G1_1$context=="CHH",],
                        LL_lyr_G1_2[LL_lyr_G1_2$context=="CHH",],
                        LL_lyr_G1_3[LL_lyr_G1_3$context=="CHH",],
                        LL_lyr_G4_1[LL_lyr_G4_1$context=="CHH",],
                        LL_lyr_G4_2[LL_lyr_G4_2$context=="CHH",],
                        LL_RS7_G1_1_lyr[LL_RS7_G1_1_lyr$context=="CHH",],
                        LL_RS7_G1_2_lyr[LL_RS7_G1_2_lyr$context=="CHH",],
                        LL_RS7_G1_3_lyr[LL_RS7_G1_3_lyr$context=="CHH",],
                        LL_RS7_G4_1_lyr[LL_RS7_G4_1_lyr$context=="CHH",],
                        LL_RS7_G4_2_lyr[LL_RS7_G4_2_lyr$context=="CHH",],
                        LL_RS7_G4_3_lyr[LL_RS7_G4_3_lyr$context=="CHH",],
                        LL_ALK_G1_1_lyr[LL_ALK_G1_1_lyr$context=="CHH",],
                        LL_ALK_G1_2_lyr[LL_ALK_G1_2_lyr$context=="CHH",],
                        LL_ALK_G1_3_lyr[LL_ALK_G1_3_lyr$context=="CHH",],
                        LL_ALK_G4_1_lyr[LL_ALK_G4_1_lyr$context=="CHH",],
                        LL_ALK_G4_2_lyr[LL_ALK_G4_2_lyr$context=="CHH",],
                        LL_ALK_G4_3_lyr[LL_ALK_G4_3_lyr$context=="CHH",],
                        LL_TKS_G1_1_lyr[LL_TKS_G1_1_lyr$context=="CHH",],
                        LL_TKS_G1_2_lyr[LL_TKS_G1_2_lyr$context=="CHH",],
                        LL_TKS_G1_3_lyr[LL_TKS_G1_3_lyr$context=="CHH",],
                        LL_TKS_G5_1_lyr[LL_TKS_G5_1_lyr$context=="CHH",],
                        LL_TKS_G5_2_lyr[LL_TKS_G5_2_lyr$context=="CHH",],
                        LL_TKS_G5_3_lyr[LL_TKS_G5_3_lyr$context=="CHH",])

## We compute the average methylation over all three samples (for two sides)

## First generate elements to create new dataframe

species_hal <- c(rep("halleri", 60), 
             rep("kamchatica-syn", 60), 
             rep("kamchatica-natural-alk", 60),
             rep("kamchatica-natural-tks", 60))
species_lyr <- c(rep("lyrata", 60), 
                 rep("kamchatica-syn", 60), 
                 rep("kamchatica-natural-alk", 60),
                 rep("kamchatica-natural-tks", 60))
side <- c(rep("halleri", 240))
cont_CG <- rep("CG", 240)
cont_CHG <- rep("CHG", 240)
cont_CHH <- rep("CHH", 240)
category <- c(rep("covered", 240))
distance <- rep(HM_CG_hal_all$distance[1:30], 8)
what <- rep(HM_CG_hal_all$what[1:30], 8)
where <- rep(HM_CG_hal_all$where[1:30], 8)
mean <- rep(0, 240)
generation <- c(rep("G1", 30),
                rep("G4", 30),
                rep("G1", 30),
                rep("G4", 30),
                rep("G1", 30),
                rep("G4", 30),
                rep("G1", 30),
                rep("G4", 30))

## halleri side
  
HM_CG_hal <- data.frame(context = cont_CG,
                     category = category,
                     distance = distance,
                     what = what,
                     where = where,
                     mean = mean,
                     side = side,
                     species = species_hal,
                     generation = generation)

HM_CHG_hal <- data.frame(context = cont_CHG,
                     category = category,
                     distance = distance,
                     what = what,
                     where = where,
                     mean = mean,
                     side = side,
                     species = species_hal,
                     generation = generation)

HM_CHH_hal <- data.frame(context = cont_CHH,
                     category = category,
                     distance = distance,
                     what = what,
                     where = where,
                     mean = mean,
                     side = side,
                     species = species_hal,
                     generation = generation)

LL_CG_hal <- data.frame(context = cont_CG,
                        category = category,
                        distance = distance,
                        what = what,
                        where = where,
                        mean = mean,
                        side = side,
                        species = species_hal,
                        generation = generation)

LL_CHG_hal <- data.frame(context = cont_CHG,
                         category = category,
                         distance = distance,
                         what = what,
                         where = where,
                         mean = mean,
                         side = side,
                         species = species_hal,
                         generation = generation)

LL_CHH_hal <- data.frame(context = cont_CHH,
                         category = category,
                         distance = distance,
                         what = what,
                         where = where,
                         mean = mean,
                         side = side,
                         species = species_hal,
                         generation = generation)


## lyrata side 

HM_CG_lyr <- data.frame(context = cont_CG,
                     category = category,
                     distance = distance,
                     what = what,
                     where = where,
                     mean = mean,
                     side = side,
                     species = species_lyr,
                     generation = generation)

HM_CHG_lyr <- data.frame(context = cont_CHG,
                      category = category,
                      distance = distance,
                      what = what,
                      where = where,
                      mean = mean,
                      side = side,
                      species = species_lyr,
                      generation = generation)

HM_CHH_lyr <- data.frame(context = cont_CHH,
                      category = category,
                      distance = distance,
                      what = what,
                      where = where,
                      mean = mean,
                      side = side,
                      species = species_lyr,
                      generation = generation)

LL_CG_lyr <- data.frame(context = cont_CG,
                        category = category,
                        distance = distance,
                        what = what,
                        where = where,
                        mean = mean,
                        side = side,
                        species = species_lyr,
                        generation = generation)

LL_CHG_lyr <- data.frame(context = cont_CHG,
                         category = category,
                         distance = distance,
                         what = what,
                         where = where,
                         mean = mean,
                         side = side,
                         species = species_lyr,
                         generation = generation)

LL_CHH_lyr <- data.frame(context = cont_CHH,
                         category = category,
                         distance = distance,
                         what = what,
                         where = where,
                         mean = mean,
                         side = side,
                         species = species_lyr,
                         generation = generation)

## We compute the new average values and insert them in the new dataframes

## halleri side

## CG context

HM_sample_average_hal1 <- c(rep(0,30))
HM_sample_average_hal4 <- c(rep(0,30))
HM_sample_average_syn1 <- c(rep(0,30))
HM_sample_average_syn4 <- c(rep(0,30))
HM_sample_average_alk1 <- c(rep(0,30))
HM_sample_average_alk4 <- c(rep(0,30))
HM_sample_average_tks1 <- c(rep(0,30))
HM_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  HM_sample_average_hal1[value] <- mean(c(HM_CG_hal_all$mean[value], HM_CG_hal_all$mean[value+30], HM_CG_hal_all$mean[value+60]))
  HM_sample_average_hal4[value] <- mean(c(HM_CG_hal_all$mean[value+90], HM_CG_hal_all$mean[value+120], HM_CG_hal_all$mean[value+150]))
  HM_sample_average_syn1[value] <- mean(c(HM_CG_hal_all$mean[value+180], HM_CG_hal_all$mean[value+210], HM_CG_hal_all$mean[value+240]))
  HM_sample_average_syn4[value] <- mean(c(HM_CG_hal_all$mean[value+270], HM_CG_hal_all$mean[value+300], HM_CG_hal_all$mean[value+330]))
  HM_sample_average_alk1[value] <- mean(c(HM_CG_hal_all$mean[value+360], HM_CG_hal_all$mean[value+390], HM_CG_hal_all$mean[value+420]))
  HM_sample_average_alk4[value] <- mean(c(HM_CG_hal_all$mean[value+450], HM_CG_hal_all$mean[value+480], HM_CG_hal_all$mean[value+510]))
  HM_sample_average_tks1[value] <- mean(c(HM_CG_hal_all$mean[value+540], HM_CG_hal_all$mean[value+570], HM_CG_hal_all$mean[value+600]))
  HM_sample_average_tks5[value] <- mean(c(HM_CG_hal_all$mean[value+630], HM_CG_hal_all$mean[value+660], HM_CG_hal_all$mean[value+690]))
}

HM_CG_hal$mean[1:30] <- HM_sample_average_hal1
HM_CG_hal$mean[31:60] <- HM_sample_average_hal4
HM_CG_hal$mean[61:90] <- HM_sample_average_syn1
HM_CG_hal$mean[91:120] <- HM_sample_average_syn4
HM_CG_hal$mean[121:150] <- HM_sample_average_alk1
HM_CG_hal$mean[151:180] <- HM_sample_average_alk4
HM_CG_hal$mean[181:210] <- HM_sample_average_tks1
HM_CG_hal$mean[211:240] <- HM_sample_average_tks5

## CHG context

HM_sample_average_hal1 <- c(rep(0,30))
HM_sample_average_hal4 <- c(rep(0,30))
HM_sample_average_syn1 <- c(rep(0,30))
HM_sample_average_syn4 <- c(rep(0,30))
HM_sample_average_alk1 <- c(rep(0,30))
HM_sample_average_alk4 <- c(rep(0,30))
HM_sample_average_tks1 <- c(rep(0,30))
HM_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  HM_sample_average_hal1[value] <- mean(c(HM_CHG_hal_all$mean[value], HM_CHG_hal_all$mean[value+30], HM_CHG_hal_all$mean[value+60]))
  HM_sample_average_hal4[value] <- mean(c(HM_CHG_hal_all$mean[value+90], HM_CHG_hal_all$mean[value+120], HM_CHG_hal_all$mean[value+150]))
  HM_sample_average_syn1[value] <- mean(c(HM_CHG_hal_all$mean[value+180], HM_CHG_hal_all$mean[value+210], HM_CHG_hal_all$mean[value+240]))
  HM_sample_average_syn4[value] <- mean(c(HM_CHG_hal_all$mean[value+270], HM_CHG_hal_all$mean[value+300], HM_CHG_hal_all$mean[value+330]))
  HM_sample_average_alk1[value] <- mean(c(HM_CHG_hal_all$mean[value+360], HM_CHG_hal_all$mean[value+390], HM_CHG_hal_all$mean[value+420]))
  HM_sample_average_alk4[value] <- mean(c(HM_CHG_hal_all$mean[value+450], HM_CHG_hal_all$mean[value+480], HM_CHG_hal_all$mean[value+510]))
  HM_sample_average_tks1[value] <- mean(c(HM_CHG_hal_all$mean[value+540], HM_CHG_hal_all$mean[value+570], HM_CHG_hal_all$mean[value+600]))
  HM_sample_average_tks5[value] <- mean(c(HM_CHG_hal_all$mean[value+630], HM_CHG_hal_all$mean[value+660], HM_CHG_hal_all$mean[value+690]))
}

HM_CHG_hal$mean[1:30] <- HM_sample_average_hal1
HM_CHG_hal$mean[31:60] <- HM_sample_average_hal4
HM_CHG_hal$mean[61:90] <- HM_sample_average_syn1
HM_CHG_hal$mean[91:120] <- HM_sample_average_syn4
HM_CHG_hal$mean[121:150] <- HM_sample_average_alk1
HM_CHG_hal$mean[151:180] <- HM_sample_average_alk4
HM_CHG_hal$mean[181:210] <- HM_sample_average_tks1
HM_CHG_hal$mean[211:240] <- HM_sample_average_tks5

## CHH context

HM_sample_average_hal1 <- c(rep(0,30))
HM_sample_average_hal4 <- c(rep(0,30))
HM_sample_average_syn1 <- c(rep(0,30))
HM_sample_average_syn4 <- c(rep(0,30))
HM_sample_average_alk1 <- c(rep(0,30))
HM_sample_average_alk4 <- c(rep(0,30))
HM_sample_average_tks1 <- c(rep(0,30))
HM_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  HM_sample_average_hal1[value] <- mean(c(HM_CHH_hal_all$mean[value], HM_CHH_hal_all$mean[value+30], HM_CHH_hal_all$mean[value+60]))
  HM_sample_average_hal4[value] <- mean(c(HM_CHH_hal_all$mean[value+90], HM_CHH_hal_all$mean[value+120], HM_CHH_hal_all$mean[value+150]))
  HM_sample_average_syn1[value] <- mean(c(HM_CHH_hal_all$mean[value+180], HM_CHH_hal_all$mean[value+210], HM_CHH_hal_all$mean[value+240]))
  HM_sample_average_syn4[value] <- mean(c(HM_CHH_hal_all$mean[value+270], HM_CHH_hal_all$mean[value+300], HM_CHH_hal_all$mean[value+330]))
  HM_sample_average_alk1[value] <- mean(c(HM_CHH_hal_all$mean[value+360], HM_CHH_hal_all$mean[value+390], HM_CHH_hal_all$mean[value+420]))
  HM_sample_average_alk4[value] <- mean(c(HM_CHH_hal_all$mean[value+450], HM_CHH_hal_all$mean[value+480], HM_CHH_hal_all$mean[value+510]))
  HM_sample_average_tks1[value] <- mean(c(HM_CHH_hal_all$mean[value+540], HM_CHH_hal_all$mean[value+570], HM_CHH_hal_all$mean[value+600]))
  HM_sample_average_tks5[value] <- mean(c(HM_CHH_hal_all$mean[value+630], HM_CHH_hal_all$mean[value+660], HM_CHH_hal_all$mean[value+690]))
}

HM_CHH_hal$mean[1:30] <- HM_sample_average_hal1
HM_CHH_hal$mean[31:60] <- HM_sample_average_hal4
HM_CHH_hal$mean[61:90] <- HM_sample_average_syn1
HM_CHH_hal$mean[91:120] <- HM_sample_average_syn4
HM_CHH_hal$mean[121:150] <- HM_sample_average_alk1
HM_CHH_hal$mean[151:180] <- HM_sample_average_alk4
HM_CHH_hal$mean[181:210] <- HM_sample_average_tks1
HM_CHH_hal$mean[211:240] <- HM_sample_average_tks5

## CG context

LL_sample_average_hal1 <- c(rep(0,30))
LL_sample_average_hal4 <- c(rep(0,30))
LL_sample_average_syn1 <- c(rep(0,30))
LL_sample_average_syn4 <- c(rep(0,30))
LL_sample_average_alk1 <- c(rep(0,30))
LL_sample_average_alk4 <- c(rep(0,30))
LL_sample_average_tks1 <- c(rep(0,30))
LL_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  LL_sample_average_hal1[value] <- mean(c(LL_CG_hal_all$mean[value], LL_CG_hal_all$mean[value+30], LL_CG_hal_all$mean[value+60]))
  LL_sample_average_hal4[value] <- mean(c(LL_CG_hal_all$mean[value+90], LL_CG_hal_all$mean[value+120]))
  LL_sample_average_syn1[value] <- mean(c(LL_CG_hal_all$mean[value+150], LL_CG_hal_all$mean[value+180], LL_CG_hal_all$mean[value+210]))
  LL_sample_average_syn4[value] <- mean(c(LL_CG_hal_all$mean[value+240], LL_CG_hal_all$mean[value+270], LL_CG_hal_all$mean[value+300]))
  LL_sample_average_alk1[value] <- mean(c(LL_CG_hal_all$mean[value+330], LL_CG_hal_all$mean[value+360], LL_CG_hal_all$mean[value+390]))
  LL_sample_average_alk4[value] <- mean(c(LL_CG_hal_all$mean[value+420], LL_CG_hal_all$mean[value+450], LL_CG_hal_all$mean[value+480]))
  LL_sample_average_tks1[value] <- mean(c(LL_CG_hal_all$mean[value+510], LL_CG_hal_all$mean[value+540], LL_CG_hal_all$mean[value+570]))
  LL_sample_average_tks5[value] <- mean(c(LL_CG_hal_all$mean[value+600], LL_CG_hal_all$mean[value+630], LL_CG_hal_all$mean[value+660]))
}

LL_CG_hal$mean[1:30] <- LL_sample_average_hal1
LL_CG_hal$mean[31:60] <- LL_sample_average_hal4
LL_CG_hal$mean[61:90] <- LL_sample_average_syn1
LL_CG_hal$mean[91:120] <- LL_sample_average_syn4
LL_CG_hal$mean[121:150] <- LL_sample_average_alk1
LL_CG_hal$mean[151:180] <- LL_sample_average_alk4
LL_CG_hal$mean[181:210] <- LL_sample_average_tks1
LL_CG_hal$mean[211:240] <- LL_sample_average_tks5

## CHG context

LL_sample_average_hal1 <- c(rep(0,30))
LL_sample_average_hal4 <- c(rep(0,30))
LL_sample_average_syn1 <- c(rep(0,30))
LL_sample_average_syn4 <- c(rep(0,30))
LL_sample_average_alk1 <- c(rep(0,30))
LL_sample_average_alk4 <- c(rep(0,30))
LL_sample_average_tks1 <- c(rep(0,30))
LL_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  LL_sample_average_hal1[value] <- mean(c(LL_CHG_hal_all$mean[value], LL_CHG_hal_all$mean[value+30], LL_CHG_hal_all$mean[value+60]))
  LL_sample_average_hal4[value] <- mean(c(LL_CHG_hal_all$mean[value+90], LL_CHG_hal_all$mean[value+120]))
  LL_sample_average_syn1[value] <- mean(c(LL_CHG_hal_all$mean[value+150], LL_CHG_hal_all$mean[value+180], LL_CHG_hal_all$mean[value+210]))
  LL_sample_average_syn4[value] <- mean(c(LL_CHG_hal_all$mean[value+240], LL_CHG_hal_all$mean[value+270], LL_CHG_hal_all$mean[value+300]))
  LL_sample_average_alk1[value] <- mean(c(LL_CHG_hal_all$mean[value+330], LL_CHG_hal_all$mean[value+360], LL_CHG_hal_all$mean[value+390]))
  LL_sample_average_alk4[value] <- mean(c(LL_CHG_hal_all$mean[value+420], LL_CHG_hal_all$mean[value+450], LL_CHG_hal_all$mean[value+480]))
  LL_sample_average_tks1[value] <- mean(c(LL_CHG_hal_all$mean[value+510], LL_CHG_hal_all$mean[value+540], LL_CHG_hal_all$mean[value+570]))
  LL_sample_average_tks5[value] <- mean(c(LL_CHG_hal_all$mean[value+600], LL_CHG_hal_all$mean[value+630], LL_CHG_hal_all$mean[value+660]))
}

LL_CHG_hal$mean[1:30] <- LL_sample_average_hal1
LL_CHG_hal$mean[31:60] <- LL_sample_average_hal4
LL_CHG_hal$mean[61:90] <- LL_sample_average_syn1
LL_CHG_hal$mean[91:120] <- LL_sample_average_syn4
LL_CHG_hal$mean[121:150] <- LL_sample_average_alk1
LL_CHG_hal$mean[151:180] <- LL_sample_average_alk4
LL_CHG_hal$mean[181:210] <- LL_sample_average_tks1
LL_CHG_hal$mean[211:240] <- LL_sample_average_tks5

## CHH context

LL_sample_average_hal1 <- c(rep(0,30))
LL_sample_average_hal4 <- c(rep(0,30))
LL_sample_average_syn1 <- c(rep(0,30))
LL_sample_average_syn4 <- c(rep(0,30))
LL_sample_average_alk1 <- c(rep(0,30))
LL_sample_average_alk4 <- c(rep(0,30))
LL_sample_average_tks1 <- c(rep(0,30))
LL_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  LL_sample_average_hal1[value] <- mean(c(LL_CHH_hal_all$mean[value], LL_CHH_hal_all$mean[value+30], LL_CHH_hal_all$mean[value+60]))
  LL_sample_average_hal4[value] <- mean(c(LL_CHH_hal_all$mean[value+90], LL_CHH_hal_all$mean[value+120]))
  LL_sample_average_syn1[value] <- mean(c(LL_CHH_hal_all$mean[value+150], LL_CHH_hal_all$mean[value+180], LL_CHH_hal_all$mean[value+210]))
  LL_sample_average_syn4[value] <- mean(c(LL_CHH_hal_all$mean[value+240], LL_CHH_hal_all$mean[value+270], LL_CHH_hal_all$mean[value+300]))
  LL_sample_average_alk1[value] <- mean(c(LL_CHH_hal_all$mean[value+330], LL_CHH_hal_all$mean[value+360], LL_CHH_hal_all$mean[value+390]))
  LL_sample_average_alk4[value] <- mean(c(LL_CHH_hal_all$mean[value+420], LL_CHH_hal_all$mean[value+450], LL_CHH_hal_all$mean[value+480]))
  LL_sample_average_tks1[value] <- mean(c(LL_CHH_hal_all$mean[value+510], LL_CHH_hal_all$mean[value+540], LL_CHH_hal_all$mean[value+570]))
  LL_sample_average_tks5[value] <- mean(c(LL_CHH_hal_all$mean[value+600], LL_CHH_hal_all$mean[value+630], LL_CHH_hal_all$mean[value+660]))
}

LL_CHH_hal$mean[1:30] <- LL_sample_average_hal1
LL_CHH_hal$mean[31:60] <- LL_sample_average_hal4
LL_CHH_hal$mean[61:90] <- LL_sample_average_syn1
LL_CHH_hal$mean[91:120] <- LL_sample_average_syn4
LL_CHH_hal$mean[121:150] <- LL_sample_average_alk1
LL_CHH_hal$mean[151:180] <- LL_sample_average_alk4
LL_CHH_hal$mean[181:210] <- LL_sample_average_tks1
LL_CHH_hal$mean[211:240] <- LL_sample_average_tks5

## lyrata side 

## CG context

HM_sample_average_lyr1 <- c(rep(0,30))
HM_sample_average_lyr4 <- c(rep(0,30))
HM_sample_average_syn1 <- c(rep(0,30))
HM_sample_average_syn4 <- c(rep(0,30))
HM_sample_average_alk1 <- c(rep(0,30))
HM_sample_average_alk4 <- c(rep(0,30))
HM_sample_average_tks1 <- c(rep(0,30))
HM_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  HM_sample_average_lyr1[value] <- mean(c(HM_CG_lyr_all$mean[value], HM_CG_lyr_all$mean[value+30], HM_CG_lyr_all$mean[value+60]))
  HM_sample_average_lyr4[value] <- mean(c(HM_CG_lyr_all$mean[value+90], HM_CG_lyr_all$mean[value+120], HM_CG_lyr_all$mean[value+150]))
  HM_sample_average_syn1[value] <- mean(c(HM_CG_lyr_all$mean[value+180], HM_CG_lyr_all$mean[value+210], HM_CG_lyr_all$mean[value+240]))
  HM_sample_average_syn4[value] <- mean(c(HM_CG_lyr_all$mean[value+270], HM_CG_lyr_all$mean[value+300], HM_CG_lyr_all$mean[value+330]))
  HM_sample_average_alk1[value] <- mean(c(HM_CG_lyr_all$mean[value+360], HM_CG_lyr_all$mean[value+390], HM_CG_lyr_all$mean[value+420]))
  HM_sample_average_alk4[value] <- mean(c(HM_CG_lyr_all$mean[value+450], HM_CG_lyr_all$mean[value+480], HM_CG_lyr_all$mean[value+510]))
  HM_sample_average_tks1[value] <- mean(c(HM_CG_lyr_all$mean[value+540], HM_CG_lyr_all$mean[value+570], HM_CG_lyr_all$mean[value+600]))
  HM_sample_average_tks5[value] <- mean(c(HM_CG_lyr_all$mean[value+630], HM_CG_lyr_all$mean[value+660], HM_CG_lyr_all$mean[value+690]))
}

HM_CG_lyr$mean[1:30] <- HM_sample_average_lyr1
HM_CG_lyr$mean[31:60] <- HM_sample_average_lyr4
HM_CG_lyr$mean[61:90] <- HM_sample_average_syn1
HM_CG_lyr$mean[91:120] <- HM_sample_average_syn4
HM_CG_lyr$mean[121:150] <- HM_sample_average_alk1
HM_CG_lyr$mean[151:180] <- HM_sample_average_alk4
HM_CG_lyr$mean[181:210] <- HM_sample_average_tks1
HM_CG_lyr$mean[211:240] <- HM_sample_average_tks5

## CHG context

HM_sample_average_lyr1 <- c(rep(0,30))
HM_sample_average_lyr4 <- c(rep(0,30))
HM_sample_average_syn1 <- c(rep(0,30))
HM_sample_average_syn4 <- c(rep(0,30))
HM_sample_average_alk1 <- c(rep(0,30))
HM_sample_average_alk4 <- c(rep(0,30))
HM_sample_average_tks1 <- c(rep(0,30))
HM_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  HM_sample_average_lyr1[value] <- mean(c(HM_CHG_lyr_all$mean[value], HM_CHG_lyr_all$mean[value+30], HM_CHG_lyr_all$mean[value+60]))
  HM_sample_average_lyr4[value] <- mean(c(HM_CHG_lyr_all$mean[value+90], HM_CHG_lyr_all$mean[value+120], HM_CHG_lyr_all$mean[value+150]))
  HM_sample_average_syn1[value] <- mean(c(HM_CHG_lyr_all$mean[value+180], HM_CHG_lyr_all$mean[value+210], HM_CHG_lyr_all$mean[value+240]))
  HM_sample_average_syn4[value] <- mean(c(HM_CHG_lyr_all$mean[value+270], HM_CHG_lyr_all$mean[value+300], HM_CHG_lyr_all$mean[value+330]))
  HM_sample_average_alk1[value] <- mean(c(HM_CHG_lyr_all$mean[value+360], HM_CHG_lyr_all$mean[value+390], HM_CHG_lyr_all$mean[value+420]))
  HM_sample_average_alk4[value] <- mean(c(HM_CHG_lyr_all$mean[value+450], HM_CHG_lyr_all$mean[value+480], HM_CHG_lyr_all$mean[value+510]))
  HM_sample_average_tks1[value] <- mean(c(HM_CHG_lyr_all$mean[value+540], HM_CHG_lyr_all$mean[value+570], HM_CHG_lyr_all$mean[value+600]))
  HM_sample_average_tks5[value] <- mean(c(HM_CHG_lyr_all$mean[value+630], HM_CHG_lyr_all$mean[value+660], HM_CHG_lyr_all$mean[value+690]))
}

HM_CHG_lyr$mean[1:30] <- HM_sample_average_lyr1
HM_CHG_lyr$mean[31:60] <- HM_sample_average_lyr4
HM_CHG_lyr$mean[61:90] <- HM_sample_average_syn1
HM_CHG_lyr$mean[91:120] <- HM_sample_average_syn4
HM_CHG_lyr$mean[121:150] <- HM_sample_average_alk1
HM_CHG_lyr$mean[151:180] <- HM_sample_average_alk4
HM_CHG_lyr$mean[181:210] <- HM_sample_average_tks1
HM_CHG_lyr$mean[211:240] <- HM_sample_average_tks5

## CHH context

HM_sample_average_lyr1 <- c(rep(0,30))
HM_sample_average_lyr4 <- c(rep(0,30))
HM_sample_average_syn1 <- c(rep(0,30))
HM_sample_average_syn4 <- c(rep(0,30))
HM_sample_average_alk1 <- c(rep(0,30))
HM_sample_average_alk4 <- c(rep(0,30))
HM_sample_average_tks1 <- c(rep(0,30))
HM_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  HM_sample_average_lyr1[value] <- mean(c(HM_CHH_lyr_all$mean[value], HM_CHH_lyr_all$mean[value+30], HM_CHH_lyr_all$mean[value+60]))
  HM_sample_average_lyr4[value] <- mean(c(HM_CHH_lyr_all$mean[value+90], HM_CHH_lyr_all$mean[value+120], HM_CHH_lyr_all$mean[value+150]))
  HM_sample_average_syn1[value] <- mean(c(HM_CHH_lyr_all$mean[value+180], HM_CHH_lyr_all$mean[value+210], HM_CHH_lyr_all$mean[value+240]))
  HM_sample_average_syn4[value] <- mean(c(HM_CHH_lyr_all$mean[value+270], HM_CHH_lyr_all$mean[value+300], HM_CHH_lyr_all$mean[value+330]))
  HM_sample_average_alk1[value] <- mean(c(HM_CHH_lyr_all$mean[value+360], HM_CHH_lyr_all$mean[value+390], HM_CHH_lyr_all$mean[value+420]))
  HM_sample_average_alk4[value] <- mean(c(HM_CHH_lyr_all$mean[value+450], HM_CHH_lyr_all$mean[value+480], HM_CHH_lyr_all$mean[value+510]))
  HM_sample_average_tks1[value] <- mean(c(HM_CHH_lyr_all$mean[value+540], HM_CHH_lyr_all$mean[value+570], HM_CHH_lyr_all$mean[value+600]))
  HM_sample_average_tks5[value] <- mean(c(HM_CHH_lyr_all$mean[value+630], HM_CHH_lyr_all$mean[value+660], HM_CHH_lyr_all$mean[value+690]))
}

HM_CHH_lyr$mean[1:30] <- HM_sample_average_lyr1
HM_CHH_lyr$mean[31:60] <- HM_sample_average_lyr4
HM_CHH_lyr$mean[61:90] <- HM_sample_average_syn1
HM_CHH_lyr$mean[91:120] <- HM_sample_average_syn4
HM_CHH_lyr$mean[121:150] <- HM_sample_average_alk1
HM_CHH_lyr$mean[151:180] <- HM_sample_average_alk4
HM_CHH_lyr$mean[181:210] <- HM_sample_average_tks1
HM_CHH_lyr$mean[211:240] <- HM_sample_average_tks5

## CG context

LL_sample_average_lyr1 <- c(rep(0,30))
LL_sample_average_lyr4 <- c(rep(0,30))
LL_sample_average_syn1 <- c(rep(0,30))
LL_sample_average_syn4 <- c(rep(0,30))
LL_sample_average_alk1 <- c(rep(0,30))
LL_sample_average_alk4 <- c(rep(0,30))
LL_sample_average_tks1 <- c(rep(0,30))
LL_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  LL_sample_average_lyr1[value] <- mean(c(LL_CG_lyr_all$mean[value], LL_CG_lyr_all$mean[value+30], LL_CG_lyr_all$mean[value+60]))
  LL_sample_average_lyr4[value] <- mean(c(LL_CG_lyr_all$mean[value+90], LL_CG_lyr_all$mean[value+120]))
  LL_sample_average_syn1[value] <- mean(c(LL_CG_lyr_all$mean[value+150], LL_CG_lyr_all$mean[value+180], LL_CG_lyr_all$mean[value+210]))
  LL_sample_average_syn4[value] <- mean(c(LL_CG_lyr_all$mean[value+240], LL_CG_lyr_all$mean[value+270], LL_CG_lyr_all$mean[value+300]))
  LL_sample_average_alk1[value] <- mean(c(LL_CG_lyr_all$mean[value+330], LL_CG_lyr_all$mean[value+360], LL_CG_lyr_all$mean[value+390]))
  LL_sample_average_alk4[value] <- mean(c(LL_CG_lyr_all$mean[value+420], LL_CG_lyr_all$mean[value+450], LL_CG_lyr_all$mean[value+480]))
  LL_sample_average_tks1[value] <- mean(c(LL_CG_lyr_all$mean[value+510], LL_CG_lyr_all$mean[value+540], LL_CG_lyr_all$mean[value+570]))
  LL_sample_average_tks5[value] <- mean(c(LL_CG_lyr_all$mean[value+600], LL_CG_lyr_all$mean[value+630], LL_CG_lyr_all$mean[value+660]))
}

LL_CG_lyr$mean[1:30] <- LL_sample_average_lyr1
LL_CG_lyr$mean[31:60] <- LL_sample_average_lyr4
LL_CG_lyr$mean[61:90] <- LL_sample_average_syn1
LL_CG_lyr$mean[91:120] <- LL_sample_average_syn4
LL_CG_lyr$mean[121:150] <- LL_sample_average_alk1
LL_CG_lyr$mean[151:180] <- LL_sample_average_alk4
LL_CG_lyr$mean[181:210] <- LL_sample_average_tks1
LL_CG_lyr$mean[211:240] <- LL_sample_average_tks5

## CHG context

LL_sample_average_lyr1 <- c(rep(0,30))
LL_sample_average_lyr4 <- c(rep(0,30))
LL_sample_average_syn1 <- c(rep(0,30))
LL_sample_average_syn4 <- c(rep(0,30))
LL_sample_average_alk1 <- c(rep(0,30))
LL_sample_average_alk4 <- c(rep(0,30))
LL_sample_average_tks1 <- c(rep(0,30))
LL_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  LL_sample_average_lyr1[value] <- mean(c(LL_CHG_lyr_all$mean[value], LL_CHG_lyr_all$mean[value+30], LL_CHG_lyr_all$mean[value+60]))
  LL_sample_average_lyr4[value] <- mean(c(LL_CHG_lyr_all$mean[value+90], LL_CHG_lyr_all$mean[value+120]))
  LL_sample_average_syn1[value] <- mean(c(LL_CHG_lyr_all$mean[value+150], LL_CHG_lyr_all$mean[value+180], LL_CHG_lyr_all$mean[value+210]))
  LL_sample_average_syn4[value] <- mean(c(LL_CHG_lyr_all$mean[value+240], LL_CHG_lyr_all$mean[value+270], LL_CHG_lyr_all$mean[value+300]))
  LL_sample_average_alk1[value] <- mean(c(LL_CHG_lyr_all$mean[value+330], LL_CHG_lyr_all$mean[value+360], LL_CHG_lyr_all$mean[value+390]))
  LL_sample_average_alk4[value] <- mean(c(LL_CHG_lyr_all$mean[value+420], LL_CHG_lyr_all$mean[value+450], LL_CHG_lyr_all$mean[value+480]))
  LL_sample_average_tks1[value] <- mean(c(LL_CHG_lyr_all$mean[value+510], LL_CHG_lyr_all$mean[value+540], LL_CHG_lyr_all$mean[value+570]))
  LL_sample_average_tks5[value] <- mean(c(LL_CHG_lyr_all$mean[value+600], LL_CHG_lyr_all$mean[value+630], LL_CHG_lyr_all$mean[value+660]))
}

LL_CHG_lyr$mean[1:30] <- LL_sample_average_lyr1
LL_CHG_lyr$mean[31:60] <- LL_sample_average_lyr4
LL_CHG_lyr$mean[61:90] <- LL_sample_average_syn1
LL_CHG_lyr$mean[91:120] <- LL_sample_average_syn4
LL_CHG_lyr$mean[121:150] <- LL_sample_average_alk1
LL_CHG_lyr$mean[151:180] <- LL_sample_average_alk4
LL_CHG_lyr$mean[181:210] <- LL_sample_average_tks1
LL_CHG_lyr$mean[211:240] <- LL_sample_average_tks5

## CHH context

LL_sample_average_lyr1 <- c(rep(0,30))
LL_sample_average_lyr4 <- c(rep(0,30))
LL_sample_average_syn1 <- c(rep(0,30))
LL_sample_average_syn4 <- c(rep(0,30))
LL_sample_average_alk1 <- c(rep(0,30))
LL_sample_average_alk4 <- c(rep(0,30))
LL_sample_average_tks1 <- c(rep(0,30))
LL_sample_average_tks5 <- c(rep(0,30))


for (value in c(1:30)){
  LL_sample_average_lyr1[value] <- mean(c(LL_CHH_lyr_all$mean[value], LL_CHH_lyr_all$mean[value+30], LL_CHH_lyr_all$mean[value+60]))
  LL_sample_average_lyr4[value] <- mean(c(LL_CHH_lyr_all$mean[value+90], LL_CHH_lyr_all$mean[value+120]))
  LL_sample_average_syn1[value] <- mean(c(LL_CHH_lyr_all$mean[value+150], LL_CHH_lyr_all$mean[value+180], LL_CHH_lyr_all$mean[value+210]))
  LL_sample_average_syn4[value] <- mean(c(LL_CHH_lyr_all$mean[value+240], LL_CHH_lyr_all$mean[value+270], LL_CHH_lyr_all$mean[value+300]))
  LL_sample_average_alk1[value] <- mean(c(LL_CHH_lyr_all$mean[value+330], LL_CHH_lyr_all$mean[value+360], LL_CHH_lyr_all$mean[value+390]))
  LL_sample_average_alk4[value] <- mean(c(LL_CHH_lyr_all$mean[value+420], LL_CHH_lyr_all$mean[value+450], LL_CHH_lyr_all$mean[value+480]))
  LL_sample_average_tks1[value] <- mean(c(LL_CHH_lyr_all$mean[value+510], LL_CHH_lyr_all$mean[value+540], LL_CHH_lyr_all$mean[value+570]))
  LL_sample_average_tks5[value] <- mean(c(LL_CHH_lyr_all$mean[value+600], LL_CHH_lyr_all$mean[value+630], LL_CHH_lyr_all$mean[value+660]))
}

LL_CHH_lyr$mean[1:30] <- LL_sample_average_lyr1
LL_CHH_lyr$mean[31:60] <- LL_sample_average_lyr4
LL_CHH_lyr$mean[61:90] <- LL_sample_average_syn1
LL_CHH_lyr$mean[91:120] <- LL_sample_average_syn4
LL_CHH_lyr$mean[121:150] <- LL_sample_average_alk1
LL_CHH_lyr$mean[151:180] <- LL_sample_average_alk4
LL_CHH_lyr$mean[181:210] <- LL_sample_average_tks1
LL_CHH_lyr$mean[211:240] <- LL_sample_average_tks5


### METHimpute plotting function


df1 <- HM_CG_hal
df2 <- HM_CHG_hal
df3 <- HM_CHH_hal
df4 <- HM_CG_lyr
df5 <- HM_CHG_lyr
df6 <- HM_CHH_lyr

range <- 500

breaks <- c(c(-1, -0.5, 0, 0.5, 1, 1.5, 2) * range)
labels <- c(-range, -range/2, '0%', '50%', '100%', range/2, range)
ylabs <- c('Mean methylation\nlevel', 'Mean recalibrated\nmethylation level')


colours <- hue_pal()(10)

ggplt1 <- ggplot(df1, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1)
ggplt1 <- ggplt1 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt1 <- ggplt1 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt1 <- ggplt1 + scale_y_continuous(limits=c(0,0.60))
ggplt1 <- ggplt1 + xlab(NULL)
ggplt1 <- ggplt1 + ylab(ylabs[2])
ggplt1 <- ggplt1 + theme_bw()
ggplt1 <- ggplt1 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt1 <- ggplt1 + ggtitle(expression(paste("Cold conditions (CG) ",italic("\nhalleri-side"))))
ggplt1

# Enrichment profile halleri CHG

ggplt2 <- ggplot(df2, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt2 <- ggplt2 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt2 <- ggplt2 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt2 <- ggplt2 + scale_y_continuous(limits=c(0,0.25))
ggplt2 <- ggplt2 + xlab(NULL)
ggplt2 <- ggplt2 + ylab(ylabs[2])
ggplt2 <- ggplt2 + theme_bw()
ggplt2 <- ggplt2 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         #legend.title = element_blank(),
                         legend.position = "none")
ggplt2 <- ggplt2 + ggtitle(expression(paste("Cold conditions (CHG) ",italic("\nhalleri-side"))))
ggplt2

# Enrichment profile halleri CHH

ggplt3 <- ggplot(df3, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt3 <- ggplt3 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt3 <- ggplt3 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt3 <- ggplt3 + scale_y_continuous(limits=c(0,0.15))
ggplt3 <- ggplt3 + xlab('Distance from annotation')
ggplt3 <- ggplt3 + ylab(ylabs[2])
ggplt3 <- ggplt3 + theme_bw()
ggplt3 <- ggplt3 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt3 <- ggplt3 + ggtitle(expression(paste("Cold conditions (CHH) ",italic("\nhalleri-side"))))
ggplt3

# Enrichment profile lyrata CG

ggplt4 <- ggplot(df4, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt4 <- ggplt4 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt4 <- ggplt4 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt4 <- ggplt4 + scale_y_continuous(limits=c(0,0.60))
ggplt4 <- ggplt4 + xlab(NULL)
ggplt4 <- ggplt4 + ylab(NULL)
ggplt4 <- ggplt4 + theme_bw()
ggplt4 <- ggplt4 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt4 <- ggplt4 + ggtitle(expression(italic("lyrata-side")))
ggplt4 <- ggplt4 + ggtitle(expression(paste("Cold conditions (CG) ",italic("\nlyrata-side"))))
ggplt4

# Enrichment profile lyrata CHG

ggplt5 <- ggplot(df5, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt5 <- ggplt5 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt5 <- ggplt5 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt5 <- ggplt5 + scale_y_continuous(limits=c(0,0.25))
ggplt5 <- ggplt5 + xlab(NULL)
ggplt5 <- ggplt5 + ylab(NULL)
ggplt5 <- ggplt5 + theme_bw()
ggplt5 <- ggplt5 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         #legend.title = element_blank(),
                         legend.position = "none")
ggplt5 <- ggplt5 + ggtitle(expression(paste("Cold conditions (CHG) ",italic("\nlyrata-side"))))
ggplt5

# Enrichment profile lyrata CHH

ggplt6 <- ggplot(df6, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt6 <- ggplt6 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt6 <- ggplt6 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt6 <- ggplt6 + scale_y_continuous(limits=c(0,0.15))
ggplt6 <- ggplt6 + xlab('Distance from annotation')
ggplt6 <- ggplt6 + ylab(NULL)
ggplt6 <- ggplt6 + theme_bw()
ggplt6 <- ggplt6 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt6 <- ggplt6 + ggtitle(expression(paste("Cold conditions (CHH) ",italic("\nlyrata-side"))))
ggplt6

### All enrichment profiles together for cold conditions

(ggplt1 + ggplt4) / (ggplt2 + ggplt5) / (ggplt3 + ggplt6)

### Repeat previous steps for hot conditions

df1 <- LL_CG_hal
df2 <- LL_CHG_hal
df3 <- LL_CHH_hal
df4 <- LL_CG_lyr
df5 <- LL_CHG_lyr
df6 <- LL_CHH_lyr

range <- 500

breaks <- c(c(-1, -0.5, 0, 0.5, 1, 1.5, 2) * range)
labels <- c(-range, -range/2, '0%', '50%', '100%', range/2, range)
ylabs <- c('Mean methylation\nlevel', 'Mean recalibrated\nmethylation level')


ggplt1 <- ggplot(df1, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt1 <- ggplt1 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt1 <- ggplt1 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt1 <- ggplt1 + scale_y_continuous(limits=c(0,0.60))
ggplt1 <- ggplt1 + xlab(NULL)
ggplt1 <- ggplt1 + ylab(ylabs[2])
ggplt1 <- ggplt1 + theme_bw()
ggplt1 <- ggplt1 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt1 <- ggplt1 + ggtitle(expression(paste("Hot conditions ",italic("\nhalleri-side"))))
ggplt1

# Enrichment profile halleri CHG

ggplt2 <- ggplot(df2, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt2 <- ggplt2 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt2 <- ggplt2 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt2 <- ggplt2 + scale_y_continuous(limits=c(0,0.25))
ggplt2 <- ggplt2 + xlab(NULL)
ggplt2 <- ggplt2 + ylab(ylabs[2])
ggplt2 <- ggplt2 + theme_bw()
ggplt2 <- ggplt2 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.title = element_blank())
ggplt2

# Enrichment profile halleri CHH

ggplt3 <- ggplot(df3, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt3 <- ggplt3 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt3 <- ggplt3 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt3 <- ggplt3 + scale_y_continuous(limits=c(0,0.15))
ggplt3 <- ggplt3 + xlab('Distance from annotation')
ggplt3 <- ggplt3 + ylab(ylabs[2])
ggplt3 <- ggplt3 + theme_bw()
ggplt3 <- ggplt3 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt3

# Enrichment profile lyrata CG

ggplt4 <- ggplot(df4, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt4 <- ggplt4 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt4 <- ggplt4 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt4 <- ggplt4 + scale_y_continuous(limits=c(0,0.6))
ggplt4 <- ggplt4 + xlab(NULL)
ggplt4 <- ggplt4 + ylab(NULL)
ggplt4 <- ggplt4 + theme_bw()
ggplt4 <- ggplt4 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt4 <- ggplt4 + ggtitle(expression(italic("lyrata-side")))
ggplt4

# Enrichment profile lyrata CHG

ggplt5 <- ggplot(df5, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt5 <- ggplt5 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt5 <- ggplt5 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt5 <- ggplt5 + scale_y_continuous(limits=c(0,0.25))
ggplt5 <- ggplt5 + xlab(NULL)
ggplt5 <- ggplt5 + ylab(NULL)
ggplt5 <- ggplt5 + theme_bw()
ggplt5 <- ggplt5 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.title = element_blank())
ggplt5

# Enrichment profile lyrata CHH

ggplt6 <- ggplot(df6, aes(x = distance, y = mean, group = interaction(species, generation), col = interaction(species, generation))) + geom_line(lwd = 1.5)
ggplt6 <- ggplt6 + scale_color_manual(values = c("kamchatica-syn.G1" = colours[4], 
                                                 "kamchatica-syn.G4" = colours[5], 
                                                 "halleri.G1" = colours[1],
                                                 "halleri.G4" = colours[2], 
                                                 "lyrata.G1" = colours[3], 
                                                 "lyrata.G4" = colours[8],
                                                 "kamchatica-natural-alk.G1" = colours[6],
                                                 "kamchatica-natural-alk.G4" = colours[7],
                                                 "kamchatica-natural-tks.G1" = colours[9],
                                                 "kamchatica-natural-tks.G4" = colours[10]))
ggplt6 <- ggplt6 + scale_x_continuous(breaks=breaks, labels=labels)
ggplt6 <- ggplt6 + scale_y_continuous(limits=c(0,0.15))
ggplt6 <- ggplt6 + xlab('Distance from annotation')
ggplt6 <- ggplt6 + ylab(NULL)
ggplt6 <- ggplt6 + theme_bw()
ggplt6 <- ggplt6 + theme(text = element_text(size=15),
                         strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
                         axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 10),
                         legend.position = "none")
ggplt6

### All enrichment profiles together for LL conditions

(ggplt1 + ggplt4) / (ggplt2 + ggplt5) / (ggplt3 + ggplt6)