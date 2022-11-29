### This script plots global methylation levels
### with their standard deviations with a dotplot
### and error bars. This is an alternative version
### to the first script to improve the visualization
### for publication purposes

## Import libraries

library(data.table)
library(tidyverse)
library(patchwork)

## Import data

# Set path to input data

data_path <- c("../data/Fig3a/")

### Progenitor's data

HM_hal_G1_1_MU <- read.csv(paste0(data_path, "HM_hal_G1_1_MU.txt"), sep="")
HM_hal_G1_2_MU <- read.csv(paste0(data_path, "HM_hal_G1_2_MU.txt"), sep="")
HM_hal_G1_3_MU <- read.csv(paste0(data_path, "HM_hal_G1_3_MU.txt"), sep="")

LL_hal_G1_1_MU <- read.csv(paste0(data_path, "LL_hal_G1_1_MU.txt"), sep="")
LL_hal_G1_2_MU <- read.csv(paste0(data_path, "LL_hal_G1_2_MU.txt"), sep="")
LL_hal_G1_3_MU <- read.csv(paste0(data_path, "LL_hal_G1_3_MU.txt"), sep="")

HM_hal_G4_1_MU <- read.csv(paste0(data_path, "HM_hal_G4_1_MU.txt"), sep="")
HM_hal_G4_2_MU <- read.csv(paste0(data_path, "HM_hal_G4_2_MU.txt"), sep="")
HM_hal_G4_3_MU <- read.csv(paste0(data_path, "HM_hal_G4_3_MU.txt"), sep="")

LL_hal_G4_1_MU <- read.csv(paste0(data_path, "LL_hal_G4_1_MU.txt"), sep="")
LL_hal_G4_2_MU <- read.csv(paste0(data_path, "LL_hal_G4_2_MU.txt"), sep="")

HM_lyr_G1_1_MU <- read.csv(paste0(data_path, "HM_lyr_G1_1_MU.txt"), sep="")
HM_lyr_G1_2_MU <- read.csv(paste0(data_path, "HM_lyr_G1_2_MU.txt"), sep="")
HM_lyr_G1_3_MU <- read.csv(paste0(data_path, "HM_lyr_G1_3_MU.txt"), sep="")

LL_lyr_G1_1_MU <- read.csv(paste0(data_path, "LL_lyr_G1_1_MU.txt"), sep="")
LL_lyr_G1_2_MU <- read.csv(paste0(data_path, "LL_lyr_G1_2_MU.txt"), sep="")
LL_lyr_G1_3_MU <- read.csv(paste0(data_path, "LL_lyr_G1_3_MU.txt"), sep="")

HM_lyr_G4_1_MU <- read.csv(paste0(data_path, "HM_lyr_G4_1_MU.txt"), sep="")
HM_lyr_G4_2_MU <- read.csv(paste0(data_path, "HM_lyr_G4_2_MU.txt"), sep="")
HM_lyr_G4_3_MU <- read.csv(paste0(data_path, "HM_lyr_G4_3_MU.txt"), sep="")

LL_lyr_G4_1_MU <- read.csv(paste0(data_path, "LL_lyr_G4_1_MU.txt"), sep="")
LL_lyr_G4_2_MU <- read.csv(paste0(data_path, "LL_lyr_G4_2_MU.txt"), sep="")
LL_lyr_G4_3_MU <- read.csv(paste0(data_path, "LL_lyr_G4_3_MU.txt"), sep="")

### kamchatica data

HM_kamhal_G1_1_MU <- read.csv(paste0(data_path, "HM_RS7K_G1_1_hal_MU.txt"), sep="")
HM_kamhal_G1_2_MU <- read.csv(paste0(data_path, "HM_RS7K_G1_2_hal_MU.txt"), sep="")
HM_kamhal_G1_3_MU <- read.csv(paste0(data_path, "HM_RS7K_G1_3_hal_MU.txt"), sep="")
HM_kamlyr_G1_1_MU <- read.csv(paste0(data_path, "HM_RS7K_G1_1_lyr_MU.txt"), sep="")
HM_kamlyr_G1_2_MU <- read.csv(paste0(data_path, "HM_RS7K_G1_2_lyr_MU.txt"), sep="")
HM_kamlyr_G1_3_MU <- read.csv(paste0(data_path, "HM_RS7K_G1_3_lyr_MU.txt"), sep="")

HM_kamhal_G4_1_MU <- read.csv(paste0(data_path, "HM_RS7_G4_1_hal_MU.txt"), sep="")
HM_kamhal_G4_2_MU <- read.csv(paste0(data_path, "HM_RS7_G4_2_hal_MU.txt"), sep="")
HM_kamhal_G4_3_MU <- read.csv(paste0(data_path, "HM_RS7_G4_3_hal_MU.txt"), sep="")
HM_kamlyr_G4_1_MU <- read.csv(paste0(data_path, "HM_RS7_G4_1_lyr_MU.txt"), sep="")
HM_kamlyr_G4_2_MU <- read.csv(paste0(data_path, "HM_RS7_G4_2_lyr_MU.txt"), sep="")
HM_kamlyr_G4_3_MU <- read.csv(paste0(data_path, "HM_RS7_G4_3_lyr_MU.txt"), sep="")

LL_kamhal_G1_1_MU <- read.csv(paste0(data_path, "LL_RS7_G1_1_hal_MU.txt"), sep="")
LL_kamhal_G1_2_MU <- read.csv(paste0(data_path, "LL_RS7_G1_2_hal_MU.txt"), sep="")
LL_kamhal_G1_3_MU <- read.csv(paste0(data_path, "LL_RS7_G1_3_hal_MU.txt"), sep="")
LL_kamlyr_G1_1_MU <- read.csv(paste0(data_path, "LL_RS7_G1_1_lyr_MU.txt"), sep="")
LL_kamlyr_G1_2_MU <- read.csv(paste0(data_path, "LL_RS7_G1_2_lyr_MU.txt"), sep="")
LL_kamlyr_G1_3_MU <- read.csv(paste0(data_path, "LL_RS7_G1_3_lyr_MU.txt"), sep="")

LL_kamhal_G4_1_MU <- read.csv(paste0(data_path, "LL_RS7_G4_1_hal_MU.txt"), sep="")
LL_kamhal_G4_2_MU <- read.csv(paste0(data_path, "LL_RS7_G4_2_hal_MU.txt"), sep="")
LL_kamhal_G4_3_MU <- read.csv(paste0(data_path, "LL_RS7_G4_3_hal_MU.txt"), sep="")
LL_kamlyr_G4_1_MU <- read.csv(paste0(data_path, "LL_RS7_G4_1_lyr_MU.txt"), sep="")
LL_kamlyr_G4_2_MU <- read.csv(paste0(data_path, "LL_RS7_G4_2_lyr_MU.txt"), sep="")
LL_kamlyr_G4_3_MU <- read.csv(paste0(data_path, "LL_RS7_G4_3_lyr_MU.txt"), sep="")

HM_ALKhal_G1_1_MU <- read.csv(paste0(data_path, "HM_ALK_G1_1_hal_MU.txt"), sep="")
HM_ALKhal_G1_2_MU <- read.csv(paste0(data_path, "HM_ALK_G1_2_hal_MU.txt"), sep="")
HM_ALKhal_G1_3_MU <- read.csv(paste0(data_path, "HM_ALK_G1_3_hal_MU.txt"), sep="")
HM_ALKlyr_G1_1_MU <- read.csv(paste0(data_path, "HM_ALK_G1_1_lyr_MU.txt"), sep="")
HM_ALKlyr_G1_2_MU <- read.csv(paste0(data_path, "HM_ALK_G1_2_lyr_MU.txt"), sep="")
HM_ALKlyr_G1_3_MU <- read.csv(paste0(data_path, "HM_ALK_G1_3_lyr_MU.txt"), sep="")

LL_ALKhal_G1_1_MU <- read.csv(paste0(data_path, "LL_ALK_G1_1_hal_MU.txt"), sep="")
LL_ALKhal_G1_2_MU <- read.csv(paste0(data_path, "LL_ALK_G1_2_hal_MU.txt"), sep="")
LL_ALKhal_G1_3_MU <- read.csv(paste0(data_path, "LL_ALK_G1_3_hal_MU.txt"), sep="")
LL_ALKlyr_G1_1_MU <- read.csv(paste0(data_path, "LL_ALK_G1_1_lyr_MU.txt"), sep="")
LL_ALKlyr_G1_2_MU <- read.csv(paste0(data_path, "LL_ALK_G1_2_lyr_MU.txt"), sep="")
LL_ALKlyr_G1_3_MU <- read.csv(paste0(data_path, "LL_ALK_G1_3_lyr_MU.txt"), sep="")


HM_ALKhal_G4_1_MU <- read.csv(paste0(data_path, "HM_ALK_G4_1_hal_MU.txt"), sep="")
HM_ALKhal_G4_2_MU <- read.csv(paste0(data_path, "HM_ALK_G4_2_hal_MU.txt"), sep="")
HM_ALKhal_G4_3_MU <- read.csv(paste0(data_path, "HM_ALK_G4_3_hal_MU.txt"), sep="")
HM_ALKlyr_G4_1_MU <- read.csv(paste0(data_path, "HM_ALK_G4_1_lyr_MU.txt"), sep="")
HM_ALKlyr_G4_2_MU <- read.csv(paste0(data_path, "HM_ALK_G4_2_lyr_MU.txt"), sep="")
HM_ALKlyr_G4_3_MU <- read.csv(paste0(data_path, "HM_ALK_G4_3_lyr_MU.txt"), sep="")

LL_ALKhal_G4_1_MU <- read.csv(paste0(data_path, "LL_ALK_G4_1_hal_MU.txt"), sep="")
LL_ALKhal_G4_2_MU <- read.csv(paste0(data_path, "LL_ALK_G4_2_hal_MU.txt"), sep="")
LL_ALKhal_G4_3_MU <- read.csv(paste0(data_path, "LL_ALK_G4_3_hal_MU.txt"), sep="")
LL_ALKlyr_G4_1_MU <- read.csv(paste0(data_path, "LL_ALK_G4_1_lyr_MU.txt"), sep="")
LL_ALKlyr_G4_2_MU <- read.csv(paste0(data_path, "LL_ALK_G4_2_lyr_MU.txt"), sep="")
LL_ALKlyr_G4_3_MU <- read.csv(paste0(data_path, "LL_ALK_G4_3_lyr_MU.txt"), sep="")

HM_TKShal_G1_1_MU <- read.csv(paste0(data_path, "HM_TKS_G1_1_hal_MU.txt"), sep="")
HM_TKShal_G1_2_MU <- read.csv(paste0(data_path, "HM_TKS_G1_2_hal_MU.txt"), sep="")
HM_TKShal_G1_3_MU <- read.csv(paste0(data_path, "HM_TKS_G1_3_hal_MU.txt"), sep="")
HM_TKSlyr_G1_1_MU <- read.csv(paste0(data_path, "HM_TKS_G1_1_lyr_MU.txt"), sep="")
HM_TKSlyr_G1_2_MU <- read.csv(paste0(data_path, "HM_TKS_G1_2_lyr_MU.txt"), sep="")
HM_TKSlyr_G1_3_MU <- read.csv(paste0(data_path, "HM_TKS_G1_3_lyr_MU.txt"), sep="")

LL_TKShal_G1_1_MU <- read.csv(paste0(data_path, "LL_TKS_G1_1_hal_MU.txt"), sep="")
LL_TKShal_G1_2_MU <- read.csv(paste0(data_path, "LL_TKS_G1_2_hal_MU.txt"), sep="")
LL_TKShal_G1_3_MU <- read.csv(paste0(data_path, "LL_TKS_G1_3_hal_MU.txt"), sep="")
LL_TKSlyr_G1_1_MU <- read.csv(paste0(data_path, "LL_TKS_G1_1_lyr_MU.txt"), sep="")
LL_TKSlyr_G1_2_MU <- read.csv(paste0(data_path, "LL_TKS_G1_2_lyr_MU.txt"), sep="")
LL_TKSlyr_G1_3_MU <- read.csv(paste0(data_path, "LL_TKS_G1_3_lyr_MU.txt"), sep="")


HM_TKShal_G5_1_MU <- read.csv(paste0(data_path, "HM_TKS_G5_1_hal_MU.txt"), sep="")
HM_TKShal_G5_2_MU <- read.csv(paste0(data_path, "HM_TKS_G5_2_hal_MU.txt"), sep="")
HM_TKShal_G5_3_MU <- read.csv(paste0(data_path, "HM_TKS_G5_3_hal_MU.txt"), sep="")
HM_TKSlyr_G5_1_MU <- read.csv(paste0(data_path, "HM_TKS_G5_1_lyr_MU.txt"), sep="")
HM_TKSlyr_G5_2_MU <- read.csv(paste0(data_path, "HM_TKS_G5_2_lyr_MU.txt"), sep="")
HM_TKSlyr_G5_3_MU <- read.csv(paste0(data_path, "HM_TKS_G5_3_lyr_MU.txt"), sep="")

LL_TKShal_G5_1_MU <- read.csv(paste0(data_path, "LL_TKS_G5_1_hal_MU.txt"), sep="")
LL_TKShal_G5_2_MU <- read.csv(paste0(data_path, "LL_TKS_G5_2_hal_MU.txt"), sep="")
LL_TKShal_G5_3_MU <- read.csv(paste0(data_path, "LL_TKS_G5_3_hal_MU.txt"), sep="")
LL_TKSlyr_G5_1_MU <- read.csv(paste0(data_path, "LL_TKS_G5_1_lyr_MU.txt"), sep="")
LL_TKSlyr_G5_2_MU <- read.csv(paste0(data_path, "LL_TKS_G5_2_lyr_MU.txt"), sep="")
LL_TKSlyr_G5_3_MU <- read.csv(paste0(data_path, "LL_TKS_G5_3_lyr_MU.txt"), sep="")


###########################
#### halleri, HM (G1)  ####
###########################
# Obtain amount of methylated Cs and total Cs

G1_m_hal_HM <- c(HM_hal_G1_1_MU[2,2], HM_hal_G1_2_MU[2,2], HM_hal_G1_3_MU[2,2])

total_hal <- sum(HM_hal_G1_1_MU[,2])

# Divide methylated Cs by total Cs

G1_gml_hal_HM <- mean(G1_m_hal_HM / total_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_hal_HM <- sd(G1_m_hal_HM / total_hal) * 100

###########################
#### halleri, HM (G4)  ####
###########################

# Obtain amount of methylated Cs and total Cs

G4_m_hal_HM <- c(HM_hal_G4_1_MU[2,2], HM_hal_G4_2_MU[2,2], HM_hal_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_hal_HM <- mean(G4_m_hal_HM / total_hal) * 100

# Find sd to estimate uncertainty

G4_gml_sd_hal_HM <- sd(G4_m_hal_HM / total_hal) * 100

###########################
#### halleri, LL (G1)  ####
###########################
# Obtain amount of methylated Cs and total Cs

G1_m_hal_LL <- c(LL_hal_G1_1_MU[2,2], LL_hal_G1_2_MU[2,2], LL_hal_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_hal_LL <- mean(G1_m_hal_LL / total_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_hal_LL <- sd(G1_m_hal_LL / total_hal) * 100

###########################
####  halleri, LL (G4) ####
###########################

# Obtain amount of methylated Cs and total Cs

G4_m_hal_LL <- c(LL_hal_G4_1_MU[2,2], LL_hal_G4_2_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_hal_LL <- mean(G4_m_hal_LL / total_hal) * 100

# Find sd to estimate uncertainty

G4_gml_sd_hal_LL <- sd(G4_m_hal_LL / total_hal) * 100

###########################
####  lyrata, HM (G1)  ####
###########################
# Obtain amount of methylated Cs and total Cs

G1_m_lyr_HM <- c(HM_lyr_G1_1_MU[2,2], HM_lyr_G1_2_MU[2,2], HM_lyr_G1_3_MU[2,2])

total_lyr <- sum(HM_lyr_G1_1_MU[,2])

# Divide methylated Cs by total Cs

G1_gml_lyr_HM <- mean(G1_m_lyr_HM / total_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_lyr_HM <- sd(G1_m_lyr_HM / total_lyr) * 100

###########################
####  lyrata, HM (G4)  ####
###########################
# Obtain amount of methylated Cs and total Cs

G4_m_lyr_HM <- c(HM_lyr_G4_1_MU[2,2], HM_lyr_G4_2_MU[2,2], HM_lyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_lyr_HM <- mean(G4_m_lyr_HM / total_lyr) * 100

# Find sd to estimate uncertainty

G4_gml_sd_lyr_HM <- sd(G4_m_lyr_HM / total_lyr) * 100

###########################
####  lyrata, LL (G1)  ####
###########################

# Obtain amount of methylated Cs and total Cs

G1_m_lyr_LL <- c(LL_lyr_G1_1_MU[2,2], LL_lyr_G1_2_MU[2,2], LL_lyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_lyr_LL <- mean(G1_m_lyr_LL / total_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_lyr_LL <- sd(G1_m_lyr_LL / total_lyr) * 100

###########################
####  lyrata, LL (G4) #####
###########################

# Obtain amount of methylated Cs and total Cs

G4_m_lyr_LL <- c(LL_lyr_G4_1_MU[2,2], LL_lyr_G4_2_MU[2,2], LL_lyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_lyr_LL <- mean(G4_m_lyr_LL / total_lyr) * 100

# Find sd to estimate uncertainty

G4_gml_sd_lyr_LL <- sd(G4_m_lyr_LL / total_lyr) * 100

######################################################
#####  synthetic kamchatica HM: halleri-side (G1) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_kamhal_HM <- c(HM_kamhal_G1_1_MU[2,2], HM_kamhal_G1_2_MU[2,2], HM_kamhal_G1_3_MU[2,2])

# Total masked Cs halleri side

total_masked_hal <- sum(HM_kamhal_G1_1_MU[,2])

# Divide methylated Cs by total Cs

G1_gml_kamhal_HM <- mean(G1_m_kamhal_HM / total_masked_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_kamhal_HM <- sd(G1_m_kamhal_HM / total_masked_hal) * 100

######################################################
#### synthetic kamchatica HM: halleri-side (G4) ######
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_kamhal_HM <- c(HM_kamhal_G4_1_MU[2,2], HM_kamhal_G4_2_MU[2,2], HM_kamhal_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_kamhal_HM <- mean(G4_m_kamhal_HM / total_masked_hal) * 100

# Find sd to estimate uncertainty

G4_gml_sd_kamhal_HM <- sd(G4_m_kamhal_HM / total_masked_hal) * 100

######################################################
#### synthetic kamchatica LL: halleri-side (G1)  #####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_kamhal_LL <- c(LL_kamhal_G1_1_MU[2,2], LL_kamhal_G1_2_MU[2,2], LL_kamhal_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_kamhal_LL <- mean(G1_m_kamhal_LL / total_masked_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_kamhal_LL <- sd(G1_m_kamhal_LL / total_masked_hal) * 100

######################################################
#### synthetic kamchatica LL: halleri-side (G4)  #####
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_kamhal_LL <- c(LL_kamhal_G4_1_MU[2,2], LL_kamhal_G4_2_MU[2,2], LL_kamhal_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_kamhal_LL <- mean(G4_m_kamhal_LL / total_masked_hal) * 100

# Find sd to estimate uncertainty

G4_gml_sd_kamhal_LL <- sd(G4_m_kamhal_LL / total_masked_hal) * 100

######################################################
##  natural kamchatica (ALK) HM: halleri-side (G1) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_ALKhal_HM <- c(HM_ALKhal_G1_1_MU[2,2], HM_ALKhal_G1_2_MU[2,2], HM_ALKhal_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_ALKhal_HM <- mean(G1_m_ALKhal_HM / total_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_ALKhal_HM <- sd(G1_m_ALKhal_HM / total_hal) * 100

######################################################
##  natural kamchatica (ALK) HM: halleri-side (G4) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_ALKhal_HM <- c(HM_ALKhal_G4_1_MU[2,2], HM_ALKhal_G4_2_MU[2,2], HM_ALKhal_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_ALKhal_HM <- mean(G4_m_ALKhal_HM / total_hal) * 100

# Find sd to estimate uncertainty

G4_gml_sd_ALKhal_HM <- sd(G4_m_ALKhal_HM / total_hal) * 100

######################################################
##  natural kamchatica (TKS) HM: halleri-side (G1) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_TKShal_HM <- c(HM_TKShal_G1_1_MU[2,2], HM_TKShal_G1_2_MU[2,2], HM_TKShal_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_TKShal_HM <- mean(G1_m_TKShal_HM / total_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_TKShal_HM <- sd(G1_m_TKShal_HM / total_hal) * 100

######################################################
##  natural kamchatica (TKS) HM: halleri-side (G5) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G5_m_TKShal_HM <- c(HM_TKShal_G5_1_MU[2,2], HM_TKShal_G5_2_MU[2,2], HM_TKShal_G5_3_MU[2,2])

# Divide methylated Cs by total Cs

G5_gml_TKShal_HM <- mean(G5_m_TKShal_HM / total_hal) * 100

# Find sd to estimate uncertainty

G5_gml_sd_TKShal_HM <- sd(G5_m_TKShal_HM / total_hal) * 100

######################################################
##  natural kamchatica (ALK) LL: halleri-side (G1) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_ALKhal_LL <- c(LL_ALKhal_G1_1_MU[2,2], LL_ALKhal_G1_2_MU[2,2], LL_ALKhal_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_ALKhal_LL <- mean(G1_m_ALKhal_LL / total_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_ALKhal_LL <- sd(G1_m_ALKhal_LL / total_hal) * 100

######################################################
##  natural kamchatica (ALK) LL: halleri-side (G4) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_ALKhal_LL <- c(LL_ALKhal_G4_1_MU[2,2], LL_ALKhal_G4_2_MU[2,2], LL_ALKhal_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_ALKhal_LL <- mean(G4_m_ALKhal_LL / total_hal) * 100

# Find sd to estimate uncertainty

G4_gml_sd_ALKhal_LL <- sd(G4_m_ALKhal_LL / total_hal) * 100

######################################################
##  natural kamchatica (TKS) LL: halleri-side (G1) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_TKShal_LL <- c(LL_TKShal_G1_1_MU[2,2], LL_TKShal_G1_2_MU[2,2], LL_TKShal_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_TKShal_LL <- mean(G1_m_TKShal_LL / total_hal) * 100

# Find sd to estimate uncertainty

G1_gml_sd_TKShal_LL <- sd(G1_m_TKShal_LL / total_hal) * 100

######################################################
##  natural kamchatica (TKS) LL: halleri-side (G5) ###
######################################################

# Obtain amount of methylated Cs and total Cs

G5_m_TKShal_LL <- c(LL_TKShal_G5_1_MU[2,2], LL_TKShal_G5_2_MU[2,2], LL_TKShal_G5_3_MU[2,2])

# Divide methylated Cs by total Cs

G5_gml_TKShal_LL <- mean(G5_m_TKShal_LL / total_hal) * 100

# Find sd to estimate uncertainty

G5_gml_sd_TKShal_LL <- sd(G5_m_TKShal_LL / total_hal) * 100

######################################################
####   synthetic kamchatica HM: lyrata-side (G1) #####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_kamlyr_HM <- c(HM_kamlyr_G1_1_MU[2,2], HM_kamlyr_G1_2_MU[2,2], HM_kamlyr_G1_3_MU[2,2])

# Total masked Cs halleri side

total_masked_lyr <- sum(HM_kamlyr_G1_1_MU[,2])

# Divide methylated Cs by total Cs

G1_gml_kamlyr_HM <- mean(G1_m_kamlyr_HM / total_masked_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_kamlyr_HM <- sd(G1_m_kamlyr_HM / total_masked_lyr) * 100

######################################################
##### synthetic kamchatica HM: lyrata-side (G4) ######
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_kamlyr_HM <- c(HM_kamlyr_G4_1_MU[2,2], HM_kamlyr_G4_2_MU[2,2], HM_kamlyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_kamlyr_HM <- mean(G4_m_kamlyr_HM / total_masked_lyr) * 100

# Find sd to estimate uncertainty

G4_gml_sd_kamlyr_HM <- sd(G4_m_kamlyr_HM / total_masked_lyr) * 100

######################################################
####  synthetic kamchatica LL: lyrata-side (G1)  #####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_kamlyr_LL <- c(LL_kamlyr_G1_1_MU[2,2], LL_kamlyr_G1_2_MU[2,2], LL_kamlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_kamlyr_LL <- mean(G1_m_kamlyr_LL / total_masked_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_kamlyr_LL <- sd(G1_m_kamlyr_LL / total_masked_lyr) * 100

######################################################
####  synthetic kamchatica LL: lyrata-side (G4)  #####
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_kamlyr_LL <- c(LL_kamlyr_G4_1_MU[2,2], LL_kamlyr_G4_2_MU[2,2], LL_kamlyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_kamlyr_LL <- mean(G4_m_kamlyr_LL / total_masked_lyr) * 100

# Find sd to estimate uncertainty

G4_gml_sd_kamlyr_LL <- sd(G4_m_kamlyr_LL / total_masked_lyr) * 100

######################################################
### natural kamchatica (ALK) HM: lyrata-side (G1) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_ALKlyr_HM <- c(HM_ALKlyr_G1_1_MU[2,2], HM_ALKlyr_G1_2_MU[2,2], HM_ALKlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_ALKlyr_HM <- mean(G1_m_ALKlyr_HM / total_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_ALKlyr_HM <- sd(G1_m_ALKlyr_HM / total_lyr) * 100

######################################################
### natural kamchatica (ALK) HM: lyrata-side (G4) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_ALKlyr_HM <- c(HM_ALKlyr_G4_1_MU[2,2], HM_ALKlyr_G4_2_MU[2,2], HM_ALKlyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_ALKlyr_HM <- mean(G4_m_ALKlyr_HM / total_lyr) * 100

# Find sd to estimate uncertainty

G4_gml_sd_ALKlyr_HM <- sd(G4_m_ALKlyr_HM / total_lyr) * 100

######################################################
### natural kamchatica (TKS) HM: lyrata-side (G1) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_TKSlyr_HM <- c(HM_TKSlyr_G1_1_MU[2,2], HM_TKSlyr_G1_2_MU[2,2], HM_TKSlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_TKSlyr_HM <- mean(G1_m_TKSlyr_HM / total_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_TKSlyr_HM <- sd(G1_m_TKSlyr_HM / total_lyr) * 100

######################################################
### natural kamchatica (TKS) HM: lyrata-side (G5) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G5_m_TKSlyr_HM <- c(HM_TKSlyr_G5_1_MU[2,2], HM_TKSlyr_G5_2_MU[2,2], HM_TKSlyr_G5_3_MU[2,2])

# Divide methylated Cs by total Cs

G5_gml_TKSlyr_HM <- mean(G5_m_TKSlyr_HM / total_lyr) * 100

# Find sd to estimate uncertainty

G5_gml_sd_TKSlyr_HM <- sd(G5_m_TKSlyr_HM / total_lyr) * 100

######################################################
### natural kamchatica (ALK) LL: lyrata-side (G1) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_ALKlyr_LL <- c(LL_ALKlyr_G1_1_MU[2,2], LL_ALKlyr_G1_2_MU[2,2], LL_ALKlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_ALKlyr_LL <- mean(G1_m_ALKlyr_LL / total_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_ALKlyr_LL <- sd(G1_m_ALKlyr_LL / total_lyr) * 100

######################################################
### natural kamchatica (ALK) LL: lyrata-side (G4) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_ALKlyr_LL <- c(LL_ALKlyr_G4_1_MU[2,2], LL_ALKlyr_G4_2_MU[2,2], LL_ALKlyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_ALKlyr_LL <- mean(G4_m_ALKlyr_LL / total_lyr) * 100

# Find sd to estimate uncertainty

G4_gml_sd_ALKlyr_LL <- sd(G4_m_ALKlyr_LL / total_lyr) * 100

######################################################
### natural kamchatica (TKS) LL: lyrata-side (G1) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_TKSlyr_LL <- c(LL_TKSlyr_G1_1_MU[2,2], LL_TKSlyr_G1_2_MU[2,2], LL_TKSlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_TKSlyr_LL <- mean(G1_m_TKSlyr_LL / total_lyr) * 100

# Find sd to estimate uncertainty

G1_gml_sd_TKSlyr_LL <- sd(G1_m_TKSlyr_LL / total_lyr) * 100

######################################################
### natural kamchatica (TKS) LL: lyrata-side (G5) ####
######################################################

# Obtain amount of methylated Cs and total Cs

G5_m_TKSlyr_LL <- c(LL_TKSlyr_G5_1_MU[2,2], LL_TKSlyr_G5_2_MU[2,2], LL_TKSlyr_G5_3_MU[2,2])

# Divide methylated Cs by total Cs

G5_gml_TKSlyr_LL <- mean(G5_m_TKSlyr_LL / total_lyr) * 100

# Find sd to estimate uncertainty

G5_gml_sd_TKSlyr_LL <- sd(G5_m_TKSlyr_LL / total_lyr) * 100

######################################################
##### synthetic kamchatica HM: overall (G1)   ########
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_kam_HM <- c(HM_kamhal_G1_1_MU[2,2]+HM_kamlyr_G1_1_MU[2,2], 
                 HM_kamhal_G1_2_MU[2,2]+HM_kamlyr_G1_2_MU[2,2], 
                 HM_kamhal_G1_3_MU[2,2]+HM_kamlyr_G1_3_MU[2,2])

total_kam <- total_lyr + total_hal

# Divide methylated Cs by total Cs

G1_gml_kam_HM <- mean(G1_m_kam_HM / total_kam) * 100

# Find sd to estimate uncertainty

G1_gml_sd_kam_HM <- sd(G1_m_kam_HM / total_kam) * 100

######################################################
######  synthetic kamchatica HM: overall (G4) ########
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_kam_HM <- c(HM_kamhal_G4_1_MU[2,2]+HM_kamlyr_G4_1_MU[2,2], 
                 HM_kamhal_G4_2_MU[2,2]+HM_kamlyr_G4_2_MU[2,2], 
                 HM_kamhal_G4_3_MU[2,2]+HM_kamlyr_G4_3_MU[2,2])

total_kam <- total_lyr + total_hal

# Divide methylated Cs by total Cs

G4_gml_kam_HM <- mean(G4_m_kam_HM / total_kam) * 100

# Find sd to estimate uncertainty

G4_gml_sd_kam_HM <- sd(G4_m_kam_HM / total_kam) * 100

######################################################
###### synthetic kamchatica LL: overall (G1) #########
######################################################

# Obtain amount of methylated Cs and total Cs

G1_m_kam_LL <- c(LL_kamhal_G1_1_MU[2,2]+LL_kamlyr_G1_1_MU[2,2], 
                 LL_kamhal_G1_2_MU[2,2]+LL_kamlyr_G1_2_MU[2,2], 
                 LL_kamhal_G1_3_MU[2,2]+LL_kamlyr_G1_3_MU[2,2])

total_kam <- total_lyr + total_hal

# Divide methylated Cs by total Cs

G1_gml_kam_LL <- mean(G1_m_kam_LL / total_kam) * 100

# Find sd to estimate uncertainty

G1_gml_sd_kam_LL <- sd(G1_m_kam_LL / total_kam) * 100

######################################################
###### synthetic kamchatica LL: overall (G4) #########
######################################################

# Obtain amount of methylated Cs and total Cs

G4_m_kam_LL <- c(LL_kamhal_G4_1_MU[2,2]+LL_kamlyr_G4_1_MU[2,2], 
                 LL_kamhal_G4_2_MU[2,2]+LL_kamlyr_G4_2_MU[2,2], 
                 LL_kamhal_G4_3_MU[2,2]+LL_kamlyr_G4_3_MU[2,2])

total_kam <- total_lyr + total_hal

# Divide methylated Cs by total Cs

G4_gml_kam_LL <- mean(G4_m_kam_LL / total_kam) * 100

# Find sd to estimate uncertainty

G4_gml_sd_kam_LL <- sd(G4_m_kam_LL / total_kam) * 100

######################################################
#### natural kamchatica (ALK) HM: overall (G1) #######
######################################################

G1_m_ALK_HM <- c(HM_ALKhal_G1_1_MU[2,2] + HM_ALKlyr_G1_1_MU[2,2],
                 HM_ALKhal_G1_2_MU[2,2] + HM_ALKlyr_G1_2_MU[2,2], 
                 HM_ALKhal_G1_3_MU[2,2] + HM_ALKlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_ALK_HM <- mean(G1_m_ALK_HM / total_kam) * 100

# Find sd to estimate uncertainty

G1_gml_sd_ALK_HM <- sd(G1_m_ALK_HM / total_kam) * 100

######################################################
#### natural kamchatica (ALK) HM: overall (G4) #######
######################################################

G4_m_ALK_HM <- c(HM_ALKhal_G4_1_MU[2,2] + HM_ALKlyr_G4_1_MU[2,2],
                 HM_ALKhal_G4_2_MU[2,2] + HM_ALKlyr_G4_2_MU[2,2], 
                 HM_ALKhal_G4_3_MU[2,2] + HM_ALKlyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_ALK_HM <- mean(G4_m_ALK_HM / total_kam) * 100

# Find sd to estimate uncertainty

G4_gml_sd_ALK_HM <- sd(G4_m_ALK_HM / total_kam) * 100

######################################################
#### natural kamchatica (TKS) HM: overall (G1) #######
######################################################

G1_m_TKS_HM <- c(HM_TKShal_G1_1_MU[2,2] + HM_TKSlyr_G1_1_MU[2,2],
                 HM_TKShal_G1_2_MU[2,2] + HM_TKSlyr_G1_2_MU[2,2], 
                 HM_TKShal_G1_3_MU[2,2] + HM_TKSlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_TKS_HM <- mean(G1_m_TKS_HM / total_kam) * 100

# Find sd to estimate uncertainty

G1_gml_sd_TKS_HM <- sd(G1_m_TKS_HM / total_kam) * 100

######################################################
#### natural kamchatica (TKS) HM: overall (G5) #######
######################################################

G5_m_TKS_HM <- c(HM_TKShal_G5_1_MU[2,2] + HM_TKSlyr_G5_1_MU[2,2],
                 HM_TKShal_G5_2_MU[2,2] + HM_TKSlyr_G5_2_MU[2,2], 
                 HM_TKShal_G5_3_MU[2,2] + HM_TKSlyr_G5_3_MU[2,2])

# Divide methylated Cs by total Cs

G5_gml_TKS_HM <- mean(G5_m_TKS_HM / total_kam) * 100

# Find sd to estimate uncertainty

G5_gml_sd_TKS_HM <- sd(G5_m_TKS_HM / total_kam) * 100

######################################################
#### natural kamchatica (ALK) LL: overall (G1) #######
######################################################

G1_m_ALK_LL <- c(LL_ALKhal_G1_1_MU[2,2] + LL_ALKlyr_G4_1_MU[2,2], 
                 LL_ALKhal_G1_2_MU[2,2] + LL_ALKlyr_G1_2_MU[2,2], 
                 LL_ALKhal_G1_3_MU[2,2] + LL_ALKlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_ALK_LL <- mean(G1_m_ALK_LL / total_kam) * 100

# Find sd to estimate uncertainty

G1_gml_sd_ALK_LL <- sd(G1_m_ALK_LL / total_kam) * 100

######################################################
#### natural kamchatica (ALK) LL: overall (G4) #######
######################################################

G4_m_ALK_LL <- c(LL_ALKhal_G4_1_MU[2,2] + LL_ALKlyr_G4_1_MU[2,2], 
                 LL_ALKhal_G4_2_MU[2,2] + LL_ALKlyr_G4_2_MU[2,2], 
                 LL_ALKhal_G4_3_MU[2,2] + LL_ALKlyr_G4_3_MU[2,2])

# Divide methylated Cs by total Cs

G4_gml_ALK_LL <- mean(G4_m_ALK_LL / total_kam) * 100

# Find sd to estimate uncertainty

G4_gml_sd_ALK_LL <- sd(G4_m_ALK_LL / total_kam) * 100

######################################################
#### natural kamchatica (TKS) LL: overall (G1) #######
######################################################

G1_m_TKS_LL <- c(LL_TKShal_G1_1_MU[2,2] + LL_TKSlyr_G1_1_MU[2,2], 
                 LL_TKShal_G1_2_MU[2,2] + LL_TKSlyr_G1_2_MU[2,2], 
                 LL_TKShal_G1_3_MU[2,2] + LL_TKSlyr_G1_3_MU[2,2])

# Divide methylated Cs by total Cs

G1_gml_TKS_LL <- mean(G1_m_TKS_LL / total_kam) * 100

# Find sd to estimate uncertainty

G1_gml_sd_TKS_LL <- sd(G1_m_TKS_LL / total_kam) * 100

######################################################
#### natural kamchatica (TKS) LL: overall (G5) #######
######################################################

G5_m_TKS_LL <- c(LL_TKShal_G5_1_MU[2,2] + LL_TKSlyr_G5_1_MU[2,2], 
                 LL_TKShal_G5_2_MU[2,2] + LL_TKSlyr_G5_2_MU[2,2], 
                 LL_TKShal_G5_3_MU[2,2] + LL_TKSlyr_G5_3_MU[2,2])

# Divide methylated Cs by total Cs

G5_gml_TKS_LL <- mean(G5_m_TKS_LL / total_kam) * 100

# Find sd to estimate uncertainty

G5_gml_sd_TKS_LL <- sd(G5_m_TKS_LL / total_kam) * 100

# Adding natural kamchatica populations to the plot

### Create second dataframe for changes over generations, natural kamchatica and different conditions

conditions <- c(rep("HM", 4), rep("LL", 4),
                rep("HM", 6), rep("LL", 6),
                rep("HM", 6), rep("LL", 6),
                rep("HM", 6), rep("LL", 6))

generation <- c(rep("G1", 2), rep("G4", 2),
                rep("G1", 2), rep("G4", 2),
                rep("G1", 3), rep("G4", 3),
                rep("G1", 3), rep("G4", 3),
                rep("G1", 3), rep("G4", 3),
                rep("G1", 3), rep("G4", 3),
                rep("G1", 3), rep("G4", 3),
                rep("G1", 3), rep("G4", 3))

species <- c(rep(c("halleri", "lyrata"), 4),
             rep("kamchatica-syn", 12),
             rep("kamchatica-natural-alk", 12),
             rep("kamchatica-natural-tks", 12))

ploidy <- c(rep("diploid", 8),
            rep("syn", 12),
            rep("nat", 24))

side <- c(rep(c("halleri", "lyrata"), 4),
          rep(c("halleri", "lyrata", "K"), 12))

count <- c(1:44)

GML <- c(G1_gml_hal_HM, G1_gml_lyr_HM, 
         G4_gml_hal_HM, G4_gml_lyr_HM,
         G1_gml_hal_LL, G1_gml_lyr_LL, 
         G4_gml_hal_LL, G4_gml_lyr_LL, 
         G1_gml_kamhal_HM, G1_gml_kamlyr_HM, G1_gml_kam_HM, 
         G4_gml_kamhal_HM, G4_gml_kamlyr_HM, G4_gml_kam_HM,
         G1_gml_kamhal_LL, G1_gml_kamlyr_LL, G1_gml_kam_LL,
         G4_gml_kamhal_LL, G4_gml_kamlyr_LL, G4_gml_kam_LL,
         G1_gml_ALKhal_HM, G1_gml_ALKlyr_HM, G1_gml_ALK_HM,
         G4_gml_ALKhal_HM, G4_gml_ALKlyr_HM, G4_gml_ALK_HM,
         G1_gml_ALKhal_LL, G1_gml_ALKlyr_LL, G1_gml_ALK_LL,
         G4_gml_ALKhal_LL, G4_gml_ALKlyr_LL, G4_gml_ALK_LL,
         G1_gml_TKShal_HM, G1_gml_TKSlyr_HM, G1_gml_TKS_HM,
         G5_gml_TKShal_HM, G5_gml_TKSlyr_HM, G5_gml_TKS_HM,
         G1_gml_TKShal_LL, G1_gml_TKSlyr_LL, G1_gml_TKS_LL,
         G5_gml_TKShal_LL, G5_gml_TKSlyr_LL, G5_gml_TKS_LL)

upper <- c(G1_gml_hal_HM + G1_gml_sd_hal_HM,
           G1_gml_lyr_HM + G1_gml_sd_lyr_HM,
           G4_gml_hal_HM + G4_gml_sd_hal_HM,
           G4_gml_lyr_HM + G4_gml_sd_lyr_HM,
           G1_gml_hal_LL + G1_gml_sd_hal_LL,
           G1_gml_lyr_LL + G1_gml_sd_lyr_LL,
           G4_gml_hal_LL + G4_gml_sd_hal_LL,
           G4_gml_lyr_LL + G4_gml_sd_lyr_LL,
           G1_gml_kamhal_HM + G1_gml_sd_kamhal_HM,
           G1_gml_kamlyr_HM + G1_gml_sd_kamlyr_HM,
           G1_gml_kam_HM + G1_gml_sd_kam_HM,
           G4_gml_kamhal_HM + G4_gml_sd_kamhal_HM,
           G4_gml_kamlyr_HM + G4_gml_sd_kamlyr_HM,
           G4_gml_kam_HM + G4_gml_sd_kam_HM,
           G1_gml_kamhal_LL + G1_gml_sd_kamhal_LL,
           G1_gml_kamlyr_LL + G1_gml_sd_kamlyr_LL,
           G1_gml_kam_LL + G1_gml_sd_kam_LL,
           G4_gml_kamhal_LL + G4_gml_sd_kamhal_LL,
           G4_gml_kamlyr_LL + G4_gml_sd_kamlyr_LL,
           G4_gml_kam_LL + G4_gml_sd_kam_LL,
           G1_gml_ALKhal_HM + G1_gml_sd_ALKhal_HM, 
           G1_gml_ALKlyr_HM + G1_gml_sd_ALKlyr_HM,
           G1_gml_ALK_HM + G1_gml_sd_ALK_HM,
           G4_gml_ALKhal_HM + G4_gml_sd_ALKhal_HM, 
           G4_gml_ALKlyr_HM + G4_gml_sd_ALKlyr_HM,
           G4_gml_ALK_HM + G4_gml_sd_ALK_HM,
           G1_gml_ALKhal_LL + G1_gml_sd_ALKhal_LL, 
           G1_gml_ALKlyr_LL + G1_gml_sd_ALKhal_LL,
           G1_gml_ALK_LL + G1_gml_sd_ALK_LL,
           G4_gml_ALKhal_LL + G4_gml_sd_ALKhal_LL, 
           G4_gml_ALKlyr_LL + G4_gml_sd_ALKhal_LL,
           G4_gml_ALK_LL + G4_gml_sd_ALK_LL,
           G1_gml_TKShal_HM + G1_gml_sd_TKShal_HM, 
           G1_gml_TKSlyr_HM + G1_gml_sd_TKSlyr_HM,
           G1_gml_TKS_HM + G1_gml_sd_TKS_HM,
           G5_gml_TKShal_HM + G5_gml_sd_TKShal_HM, 
           G5_gml_TKSlyr_HM + G5_gml_sd_TKSlyr_HM,
           G5_gml_TKS_HM + G5_gml_sd_TKS_HM,
           G1_gml_TKShal_LL + G1_gml_sd_TKShal_LL, 
           G1_gml_TKSlyr_LL + G1_gml_sd_TKShal_LL,
           G1_gml_TKS_LL + G1_gml_sd_TKS_LL,
           G5_gml_TKShal_LL + G5_gml_sd_TKShal_LL, 
           G5_gml_TKSlyr_LL + G5_gml_sd_TKShal_LL,
           G5_gml_TKS_LL + G5_gml_sd_TKS_LL)

lower <- c(G1_gml_hal_HM - G1_gml_sd_hal_HM,
           G1_gml_lyr_HM - G1_gml_sd_lyr_HM,
           G4_gml_hal_HM - G4_gml_sd_hal_HM,
           G4_gml_lyr_HM - G4_gml_sd_lyr_HM,
           G1_gml_hal_LL - G1_gml_sd_hal_LL,
           G1_gml_lyr_LL - G1_gml_sd_lyr_LL,
           G4_gml_hal_LL - G4_gml_sd_hal_LL,
           G4_gml_lyr_LL - G4_gml_sd_lyr_LL,
           G1_gml_kamhal_HM - G1_gml_sd_kamhal_HM,
           G1_gml_kamlyr_HM - G1_gml_sd_kamlyr_HM,
           G1_gml_kam_HM - G1_gml_sd_kam_HM,
           G4_gml_kamhal_HM - G4_gml_sd_kamhal_HM,
           G4_gml_kamlyr_HM - G4_gml_sd_kamlyr_HM,
           G4_gml_kam_HM - G4_gml_sd_kam_HM,
           G1_gml_kamhal_LL - G1_gml_sd_kamhal_LL,
           G1_gml_kamlyr_LL - G1_gml_sd_kamlyr_LL,
           G1_gml_kam_LL - G1_gml_sd_kam_LL,
           G4_gml_kamhal_LL - G4_gml_sd_kamhal_LL,
           G4_gml_kamlyr_LL - G4_gml_sd_kamlyr_LL,
           G4_gml_kam_LL - G4_gml_sd_kam_LL,
           G1_gml_ALKhal_HM - G1_gml_sd_ALKhal_HM, 
           G1_gml_ALKlyr_HM - G1_gml_sd_ALKlyr_HM,
           G1_gml_ALK_HM - G1_gml_sd_ALK_HM,
           G4_gml_ALKhal_HM - G4_gml_sd_ALKhal_HM, 
           G4_gml_ALKlyr_HM - G4_gml_sd_ALKlyr_HM,
           G4_gml_ALK_HM - G4_gml_sd_ALK_HM,
           G1_gml_ALKhal_LL - G1_gml_sd_ALKhal_LL, 
           G1_gml_ALKlyr_LL - G1_gml_sd_ALKhal_LL,
           G1_gml_ALK_LL - G1_gml_sd_ALK_LL,
           G4_gml_ALKhal_LL - G4_gml_sd_ALKhal_LL, 
           G4_gml_ALKlyr_LL - G4_gml_sd_ALKhal_LL,
           G4_gml_ALK_LL - G4_gml_sd_ALK_LL,
           G1_gml_TKShal_HM - G1_gml_sd_TKShal_HM, 
           G1_gml_TKSlyr_HM - G1_gml_sd_TKSlyr_HM,
           G1_gml_TKS_HM - G1_gml_sd_TKS_HM,
           G5_gml_TKShal_HM - G5_gml_sd_TKShal_HM, 
           G5_gml_TKSlyr_HM - G5_gml_sd_TKSlyr_HM,
           G5_gml_TKS_HM - G5_gml_sd_TKS_HM,
           G1_gml_TKShal_LL - G1_gml_sd_TKShal_LL, 
           G1_gml_TKSlyr_LL - G1_gml_sd_TKShal_LL,
           G1_gml_TKS_LL - G1_gml_sd_TKS_LL,
           G5_gml_TKShal_LL - G5_gml_sd_TKShal_LL, 
           G5_gml_TKSlyr_LL - G5_gml_sd_TKShal_LL,
           G5_gml_TKS_LL - G5_gml_sd_TKS_LL)

df_GML2 <- data.frame(conditions = conditions,
                      generation = generation,
                      species = species,
                      ploidy = ploidy,
                      count = count,
                      side = side,
                      GML = GML,
                      upper = upper,
                      lower = lower)

df_GML2 <- df_GML2 %>%
  arrange(species) %>%
  mutate(species_new = factor(species, levels=c("halleri", 
                                                "lyrata", 
                                                "kamchatica-syn",
                                                "kamchatica-natural-alk", 
                                                "kamchatica-natural-tks")))

df_GML2 <- df_GML2 %>%
  arrange(ploidy) %>%
  mutate(ploidy_new = factor(ploidy, levels=c("diploid", "syn", "nat")))

df_GML2_noK_HM <- filter(df_GML2, side != "K", conditions == "HM")
df_GML2_noK_LL <- filter(df_GML2, side != "K", conditions == "LL")


# Plot with points and error bars for all sides

hm_condition <- ggplot(data = df_GML2_noK_HM, aes(x = interaction(generation, ploidy_new), y = GML, colour = species_new)) +
  geom_pointrange(aes(ymin = lower, ymax = upper, shape = species_new), position = position_dodge(width = 1), size = 1, show.legend = F) +
  facet_wrap(~side) +
  theme_bw() +
  theme(text = element_text(size=25)) +
  theme(legend.title = element_blank(),
        plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(19, 28), expand = FALSE, clip = "off") +
  annotate(geom = "text", x = 1, y = 18.5, label = "G1", size = 5) +
  annotate(geom = "text", x = 2, y = 18.5, label = "G4", size = 5) +
  annotate(geom = "text", x = 3, y = 18.5, label = "G1", size = 5) +
  annotate(geom = "text", x = 4, y = 18.5, label = "G4", size = 5) +
  annotate(geom = "text", x = 5, y = 18.5, label = "G1", size = 5) +
  annotate(geom = "text", x = 6, y = 18.5, label = "G4", size = 5) +
  annotate(geom = "text", x = 1.5, y = 18, label = "diploid", size = 5) +
  annotate(geom = "text", x = 3.5, y = 18, label = "syn", size = 5) +
  annotate(geom = "text", x = 5.5, y = 18, label = "nat", size = 5) +
  ggtitle("Cold conditions")

ll_condition <- ggplot(data = df_GML2_noK_LL, aes(x = interaction(generation, ploidy_new), y = GML, colour = species_new)) +
  geom_pointrange(aes(ymin = lower, ymax = upper, shape = species_new), position = position_dodge(width = 1), size = 1) +
  facet_wrap(~side) +
  theme_bw() +
  theme(text = element_text(size=25)) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 16),
        plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(19, 28), expand = FALSE, clip = "off") +
  annotate(geom = "text", x = 1, y = 18.5, label = "G1", size = 5) +
  annotate(geom = "text", x = 2, y = 18.5, label = "G4", size = 5) +
  annotate(geom = "text", x = 3, y = 18.5, label = "G1", size = 5) +
  annotate(geom = "text", x = 4, y = 18.5, label = "G4", size = 5) +
  annotate(geom = "text", x = 5, y = 18.5, label = "G1", size = 5) +
  annotate(geom = "text", x = 6, y = 18.5, label = "G4", size = 5) +
  annotate(geom = "text", x = 1.5, y = 18, label = "diploid", size = 5) +
  annotate(geom = "text", x = 3.5, y = 18, label = "syn", size = 5) +
  annotate(geom = "text", x = 5.5, y = 18, label = "nat", size = 5) +
  ggtitle("Hot conditions")

hm_condition + ll_condition


