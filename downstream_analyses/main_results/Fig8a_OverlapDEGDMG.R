### This script will compute and plot the overlap
### between differentially methylated and expressed
### genes. Inputs are from dmrseq outputs for
### differential methylation analysis and edgeR for
### differential expression analysis


### Import libraries

library(tidyverse)
library(patchwork)
library(VennDiagram)
library(ggVennDiagram)
library(RColorBrewer)
library(svglite)

### Import methylation data

setwd("downstream_analyses/main_results/data/Fig8a/")

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DMG_CG <- read.delim("MILD_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_halG1_v_RS7G1_DMG_CHG <- read.delim("MILD_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_halG1_v_RS7G1_DMG_CHH <- read.delim("MILD_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_halG1_v_RS7G4_DMG_CG <- read.delim("MILD_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_halG1_v_RS7G4_DMG_CHG <- read.delim("MILD_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_halG1_v_RS7G4_DMG_CHH <- read.delim("MILD_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_RS7G1_v_RS7G4_DMG_CG_hal <- read.delim("MILD_syn1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_RS7G1_v_RS7G4_DMG_CHG_hal <- read.delim("MILD_syn1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_RS7G1_v_RS7G4_DMG_CHH_hal <- read.delim("MILD_syn1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_ALKG1_v_halG1_DMG_CG <- read.delim("MILD_alk1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_ALKG1_v_halG1_DMG_CHG <- read.delim("MILD_alk1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_ALKG1_v_halG1_DMG_CHH <- read.delim("MILD_alk1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_ALKG1_v_RS7G1_DMG_CG_hal <- read.delim("MILD_alk1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_ALKG1_v_RS7G1_DMG_CHG_hal <- read.delim("MILD_alk1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_ALKG1_v_RS7G1_DMG_CHH_hal <- read.delim("MILD_alk1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_ALKG1_v_RS7G4_DMG_CG_hal <- read.delim("MILD_alk1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_ALKG1_v_RS7G4_DMG_CHG_hal <- read.delim("MILD_alk1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_ALKG1_v_RS7G4_DMG_CHH_hal <- read.delim("MILD_alk1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_TKSG1_v_halG1_DMG_CG <- read.delim("MILD_tks1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_TKSG1_v_halG1_DMG_CHG <- read.delim("MILD_tks1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_TKSG1_v_halG1_DMG_CHH <- read.delim("MILD_tks1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_TKSG1_v_RS7G1_DMG_CG_hal <- read.delim("MILD_tks1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_TKSG1_v_RS7G1_DMG_CHG_hal <- read.delim("MILD_tks1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_TKSG1_v_RS7G1_DMG_CHH_hal <- read.delim("MILD_tks1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_TKSG1_v_RS7G4_DMG_CG_hal <- read.delim("MILD_tks1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_TKSG1_v_RS7G4_DMG_CHG_hal <- read.delim("MILD_tks1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_TKSG1_v_RS7G4_DMG_CHH_hal <- read.delim("MILD_tks1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")


### lyrata side

HM_lyrG1_v_RS7G1_DMG_CG <- read.delim("MILD_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_lyrG1_v_RS7G1_DMG_CHG <- read.delim("MILD_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_lyrG1_v_RS7G1_DMG_CHH <- read.delim("MILD_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_lyrG1_v_RS7G4_DMG_CG <- read.delim("MILD_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_lyrG1_v_RS7G4_DMG_CHG <- read.delim("MILD_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_lyrG1_v_RS7G4_DMG_CHH <- read.delim("MILD_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_RS7G1_v_RS7G4_DMG_CG_lyr <- read.delim("MILD_syn1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_RS7G1_v_RS7G4_DMG_CHG_lyr <- read.delim("MILD_syn1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_RS7G1_v_RS7G4_DMG_CHH_lyr <- read.delim("MILD_syn1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_ALKG1_v_lyrG1_DMG_CG <- read.delim("MILD_alk1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_ALKG1_v_lyrG1_DMG_CHG <- read.delim("MILD_alk1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_ALKG1_v_lyrG1_DMG_CHH <- read.delim("MILD_alk1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_ALKG1_v_RS7G1_DMG_CG_lyr <- read.delim("MILD_alk1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_ALKG1_v_RS7G1_DMG_CHG_lyr <- read.delim("MILD_alk1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_ALKG1_v_RS7G1_DMG_CHH_lyr <- read.delim("MILD_alk1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_ALKG1_v_RS7G4_DMG_CG_lyr <- read.delim("MILD_alk1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_ALKG1_v_RS7G4_DMG_CHG_lyr <- read.delim("MILD_alk1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_ALKG1_v_RS7G4_DMG_CHH_lyr <- read.delim("MILD_alk1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_TKSG1_v_lyrG1_DMG_CG <- read.delim("MILD_tks1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_TKSG1_v_lyrG1_DMG_CHG <- read.delim("MILD_tks1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_TKSG1_v_lyrG1_DMG_CHH <- read.delim("MILD_tks1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_TKSG1_v_RS7G1_DMG_CG_lyr <- read.delim("MILD_tks1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_TKSG1_v_RS7G1_DMG_CHG_lyr <- read.delim("MILD_tks1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_TKSG1_v_RS7G1_DMG_CHH_lyr <- read.delim("MILD_tks1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_TKSG1_v_RS7G4_DMG_CG_lyr <- read.delim("MILD_tks1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_TKSG1_v_RS7G4_DMG_CHG_lyr <- read.delim("MILD_tks1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_TKSG1_v_RS7G4_DMG_CHH_lyr <- read.delim("MILD_tks1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")


### Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DMG_CG <- read.delim("STRESS_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_halG1_v_RS7G1_DMG_CHG <- read.delim("STRESS_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_halG1_v_RS7G1_DMG_CHH <- read.delim("STRESS_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_halG1_v_RS7G4_DMG_CG <- read.delim("STRESS_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_halG1_v_RS7G4_DMG_CHG <- read.delim("STRESS_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_halG1_v_RS7G4_DMG_CHH <- read.delim("STRESS_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_RS7G1_v_RS7G4_DMG_CG_hal <- read.delim("STRESS_syn1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_RS7G1_v_RS7G4_DMG_CHG_hal <- read.delim("STRESS_syn1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_RS7G1_v_RS7G4_DMG_CHH_hal <- read.delim("STRESS_syn1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_ALKG1_v_halG1_DMG_CG <- read.delim("STRESS_alk1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_ALKG1_v_halG1_DMG_CHG <- read.delim("STRESS_alk1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_ALKG1_v_halG1_DMG_CHH <- read.delim("STRESS_alk1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_ALKG1_v_RS7G1_DMG_CG_hal <- read.delim("STRESS_alk1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_ALKG1_v_RS7G1_DMG_CHG_hal <- read.delim("STRESS_alk1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_ALKG1_v_RS7G1_DMG_CHH_hal <- read.delim("STRESS_alk1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_ALKG1_v_RS7G4_DMG_CG_hal <- read.delim("STRESS_alk1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_ALKG1_v_RS7G4_DMG_CHG_hal <- read.delim("STRESS_alk1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_ALKG1_v_RS7G4_DMG_CHH_hal <- read.delim("STRESS_alk1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_TKSG1_v_halG1_DMG_CG <- read.delim("STRESS_tks1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_TKSG1_v_halG1_DMG_CHG <- read.delim("STRESS_tks1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_TKSG1_v_halG1_DMG_CHH <- read.delim("STRESS_tks1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_TKSG1_v_RS7G1_DMG_CG_hal <- read.delim("STRESS_tks1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_TKSG1_v_RS7G1_DMG_CHG_hal <- read.delim("STRESS_tks1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_TKSG1_v_RS7G1_DMG_CHH_hal <- read.delim("STRESS_tks1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_TKSG1_v_RS7G4_DMG_CG_hal <- read.delim("STRESS_tks1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_TKSG1_v_RS7G4_DMG_CHG_hal <- read.delim("STRESS_tks1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_TKSG1_v_RS7G4_DMG_CHH_hal <- read.delim("STRESS_tks1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")


### lyrata side

LL_lyrG1_v_RS7G1_DMG_CG <- read.delim("STRESS_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_lyrG1_v_RS7G1_DMG_CHG <- read.delim("STRESS_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_lyrG1_v_RS7G1_DMG_CHH <- read.delim("STRESS_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_lyrG1_v_RS7G4_DMG_CG <- read.delim("STRESS_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_lyrG1_v_RS7G4_DMG_CHG <- read.delim("STRESS_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_lyrG1_v_RS7G4_DMG_CHH <- read.delim("STRESS_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_RS7G1_v_RS7G4_DMG_CG_lyr <- read.delim("STRESS_syn1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_RS7G1_v_RS7G4_DMG_CHG_lyr <- read.delim("STRESS_syn1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_RS7G1_v_RS7G4_DMG_CHH_lyr <- read.delim("STRESS_syn1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_ALKG1_v_lyrG1_DMG_CG <- read.delim("STRESS_alk1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_ALKG1_v_lyrG1_DMG_CHG <- read.delim("STRESS_alk1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_ALKG1_v_lyrG1_DMG_CHH <- read.delim("STRESS_alk1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_ALKG1_v_RS7G1_DMG_CG_lyr <- read.delim("STRESS_alk1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_ALKG1_v_RS7G1_DMG_CHG_lyr <- read.delim("STRESS_alk1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_ALKG1_v_RS7G1_DMG_CHH_lyr <- read.delim("STRESS_alk1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_ALKG1_v_RS7G4_DMG_CG_lyr <- read.delim("STRESS_alk1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_ALKG1_v_RS7G4_DMG_CHG_lyr <- read.delim("STRESS_alk1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_ALKG1_v_RS7G4_DMG_CHH_lyr <- read.delim("STRESS_alk1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_TKSG1_v_lyrG1_DMG_CG <- read.delim("STRESS_tks1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_TKSG1_v_lyrG1_DMG_CHG <- read.delim("STRESS_tks1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_TKSG1_v_lyrG1_DMG_CHH <- read.delim("STRESS_tks1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_TKSG1_v_RS7G1_DMG_CG_lyr <- read.delim("STRESS_tks1_v_syn1/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_TKSG1_v_RS7G1_DMG_CHG_lyr <- read.delim("STRESS_tks1_v_syn1/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_TKSG1_v_RS7G1_DMG_CHH_lyr <- read.delim("STRESS_tks1_v_syn1/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_TKSG1_v_RS7G4_DMG_CG_lyr <- read.delim("STRESS_tks1_v_syn4/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_TKSG1_v_RS7G4_DMG_CHG_lyr <- read.delim("STRESS_tks1_v_syn4/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_TKSG1_v_RS7G4_DMG_CHH_lyr <- read.delim("STRESS_tks1_v_syn4/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

# We filter all DMGs that are part of low coverage regions

# Import file with low coverage genes

setwd("downstream_analyses/main_results/data/Fig5a/")

HM_hal_lowC_genes  <- read.delim("HM_lowC_genes_hal.txt", header = FALSE)
HM_lyr_lowC_genes  <- read.table("HM_lowC_genes_lyr.txt", quote="\"", comment.char="")

LL_hal_lowC_genes  <- read.table("LL_lowC_genes_hal.txt", quote="\"", comment.char="")
LL_lyr_lowC_genes  <- read.table("LL_lowC_genes_lyr.txt", quote="\"", comment.char="")

HM_halG1_v_RS7G1_DMG_CG <- filter(HM_halG1_v_RS7G1_DMG_CG, 
                                  !(geneID %in% HM_hal_lowC_genes$V1))
HM_halG1_v_RS7G1_DMG_CHG <- filter(HM_halG1_v_RS7G1_DMG_CHG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))
HM_halG1_v_RS7G1_DMG_CHH <- filter(HM_halG1_v_RS7G1_DMG_CHH, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))

HM_halG1_v_RS7G4_DMG_CG <- filter(HM_halG1_v_RS7G4_DMG_CG, 
                                  !(geneID %in% HM_hal_lowC_genes$V1))
HM_halG1_v_RS7G4_DMG_CHG <- filter(HM_halG1_v_RS7G4_DMG_CHG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))
HM_halG1_v_RS7G4_DMG_CHH <- filter(HM_halG1_v_RS7G4_DMG_CHH, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))

HM_RS7G1_v_RS7G4_DMG_CG_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CG_hal, 
                                      !(geneID %in% HM_hal_lowC_genes$V1))
HM_RS7G1_v_RS7G4_DMG_CHG_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CHG_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))
HM_RS7G1_v_RS7G4_DMG_CHH_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CHH_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))

HM_ALKG1_v_halG1_DMG_CG <- filter(HM_ALKG1_v_halG1_DMG_CG, 
                                  !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_halG1_DMG_CHG <- filter(HM_ALKG1_v_halG1_DMG_CHG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_halG1_DMG_CHH <- filter(HM_ALKG1_v_halG1_DMG_CHH, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))

HM_ALKG1_v_RS7G1_DMG_CG_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CG_hal, 
                                      !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_RS7G1_DMG_CHG_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CHG_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_RS7G1_DMG_CHH_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CHH_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))

HM_ALKG1_v_RS7G4_DMG_CG_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CG_hal, 
                                      !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_RS7G4_DMG_CHG_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CHG_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_RS7G4_DMG_CHH_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CHH_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))

HM_TKSG1_v_halG1_DMG_CG <- filter(HM_TKSG1_v_halG1_DMG_CG, 
                                  !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_halG1_DMG_CHG <- filter(HM_TKSG1_v_halG1_DMG_CHG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_halG1_DMG_CHH <- filter(HM_TKSG1_v_halG1_DMG_CHH, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))

HM_TKSG1_v_RS7G1_DMG_CG_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CG_hal, 
                                      !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_RS7G1_DMG_CHG_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CHG_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_RS7G1_DMG_CHH_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CHH_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))

HM_TKSG1_v_RS7G4_DMG_CG_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CG_hal, 
                                      !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_RS7G4_DMG_CHG_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CHG_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_RS7G4_DMG_CHH_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CHH_hal, 
                                       !(geneID %in% HM_hal_lowC_genes$V1))


### lyrata side

HM_lyrG1_v_RS7G1_DMG_CG <- filter(HM_lyrG1_v_RS7G1_DMG_CG, 
                                  !(geneID %in% HM_lyr_lowC_genes$V1))
HM_lyrG1_v_RS7G1_DMG_CHG <- filter(HM_lyrG1_v_RS7G1_DMG_CHG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))
HM_lyrG1_v_RS7G1_DMG_CHH <- filter(HM_lyrG1_v_RS7G1_DMG_CHH, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))

HM_lyrG1_v_RS7G4_DMG_CG <- filter(HM_lyrG1_v_RS7G4_DMG_CG, 
                                  !(geneID %in% HM_lyr_lowC_genes$V1))
HM_lyrG1_v_RS7G4_DMG_CHG <- filter(HM_lyrG1_v_RS7G4_DMG_CHG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))
HM_lyrG1_v_RS7G4_DMG_CHH <- filter(HM_lyrG1_v_RS7G4_DMG_CHH, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))

HM_RS7G1_v_RS7G4_DMG_CG_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CG_lyr, 
                                      !(geneID %in% HM_lyr_lowC_genes$V1))
HM_RS7G1_v_RS7G4_DMG_CHG_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CHG_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))
HM_RS7G1_v_RS7G4_DMG_CHH_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CHH_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))

HM_ALKG1_v_lyrG1_DMG_CG <- filter(HM_ALKG1_v_lyrG1_DMG_CG, 
                                  !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_lyrG1_DMG_CHG <- filter(HM_ALKG1_v_lyrG1_DMG_CHG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_lyrG1_DMG_CHH <- filter(HM_ALKG1_v_lyrG1_DMG_CHH, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))

HM_ALKG1_v_RS7G1_DMG_CG_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CG_lyr, 
                                      !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_RS7G1_DMG_CHG_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_RS7G1_DMG_CHH_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))

HM_ALKG1_v_RS7G4_DMG_CG_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CG_lyr, 
                                      !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_RS7G4_DMG_CHG_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_RS7G4_DMG_CHH_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))

HM_TKSG1_v_lyrG1_DMG_CG <- filter(HM_TKSG1_v_lyrG1_DMG_CG, 
                                  !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_lyrG1_DMG_CHG <- filter(HM_TKSG1_v_lyrG1_DMG_CHG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_lyrG1_DMG_CHH <- filter(HM_TKSG1_v_lyrG1_DMG_CHH, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))

HM_TKSG1_v_RS7G1_DMG_CG_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CG_lyr, 
                                      !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_RS7G1_DMG_CHG_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_RS7G1_DMG_CHH_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))

HM_TKSG1_v_RS7G4_DMG_CG_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CG_lyr, 
                                      !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_RS7G4_DMG_CHG_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_RS7G4_DMG_CHH_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(geneID %in% HM_lyr_lowC_genes$V1))


### Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DMG_CG <- filter(LL_halG1_v_RS7G1_DMG_CG, 
                                  !(geneID %in% LL_hal_lowC_genes$V1))
LL_halG1_v_RS7G1_DMG_CHG <- filter(LL_halG1_v_RS7G1_DMG_CHG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))
LL_halG1_v_RS7G1_DMG_CHH <- filter(LL_halG1_v_RS7G1_DMG_CHH, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

LL_halG1_v_RS7G4_DMG_CG <- filter(LL_halG1_v_RS7G4_DMG_CG, 
                                  !(geneID %in% LL_hal_lowC_genes$V1))
LL_halG1_v_RS7G4_DMG_CHG <- filter(LL_halG1_v_RS7G4_DMG_CHG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))
LL_halG1_v_RS7G4_DMG_CHH <- filter(LL_halG1_v_RS7G4_DMG_CHH, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

LL_RS7G1_v_RS7G4_DMG_CG_hal <- filter(LL_RS7G1_v_RS7G4_DMG_CG_hal, 
                                      !(geneID %in% LL_hal_lowC_genes$V1))
LL_RS7G1_v_RS7G4_DMG_CHG_hal <- filter(LL_RS7G1_v_RS7G4_DMG_CHG_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))
LL_RS7G1_v_RS7G4_DMG_CHH_hal <- filter(LL_RS7G1_v_RS7G4_DMG_CHH_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))

LL_ALKG1_v_halG1_DMG_CG <- filter(LL_ALKG1_v_halG1_DMG_CG, 
                                  !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_halG1_DMG_CHG <- filter(LL_ALKG1_v_halG1_DMG_CHG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_halG1_DMG_CHH <- filter(LL_ALKG1_v_halG1_DMG_CHH, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

LL_ALKG1_v_RS7G1_DMG_CG_hal <- filter(LL_ALKG1_v_RS7G1_DMG_CG_hal, 
                                      !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_RS7G1_DMG_CHG_hal <- filter(LL_ALKG1_v_RS7G1_DMG_CHG_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_RS7G1_DMG_CHH_hal <- filter(LL_ALKG1_v_RS7G1_DMG_CHH_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))

LL_ALKG1_v_RS7G4_DMG_CG_hal <- filter(LL_ALKG1_v_RS7G4_DMG_CG_hal, 
                                      !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_RS7G4_DMG_CHG_hal <- filter(LL_ALKG1_v_RS7G4_DMG_CHG_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_RS7G4_DMG_CHH_hal <- filter(LL_ALKG1_v_RS7G4_DMG_CHH_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))

LL_TKSG1_v_halG1_DMG_CG <- filter(LL_TKSG1_v_halG1_DMG_CG, 
                                  !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_halG1_DMG_CHG <- filter(LL_TKSG1_v_halG1_DMG_CHG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_halG1_DMG_CHH <- filter(LL_TKSG1_v_halG1_DMG_CHH, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

LL_TKSG1_v_RS7G1_DMG_CG_hal <- filter(LL_TKSG1_v_RS7G1_DMG_CG_hal, 
                                      !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_RS7G1_DMG_CHG_hal <- filter(LL_TKSG1_v_RS7G1_DMG_CHG_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_RS7G1_DMG_CHH_hal <- filter(LL_TKSG1_v_RS7G1_DMG_CHH_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))

LL_TKSG1_v_RS7G4_DMG_CG_hal <- filter(LL_TKSG1_v_RS7G4_DMG_CG_hal, 
                                      !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_RS7G4_DMG_CHG_hal <- filter(LL_TKSG1_v_RS7G4_DMG_CHG_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_RS7G4_DMG_CHH_hal <- filter(LL_TKSG1_v_RS7G4_DMG_CHH_hal, 
                                       !(geneID %in% LL_hal_lowC_genes$V1))


### lyrata side

LL_lyrG1_v_RS7G1_DMG_CG <- filter(LL_lyrG1_v_RS7G1_DMG_CG, 
                                  !(geneID %in% LL_lyr_lowC_genes$V1))
LL_lyrG1_v_RS7G1_DMG_CHG <- filter(LL_lyrG1_v_RS7G1_DMG_CHG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))
LL_lyrG1_v_RS7G1_DMG_CHH <- filter(LL_lyrG1_v_RS7G1_DMG_CHH, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))

LL_lyrG1_v_RS7G4_DMG_CG <- filter(LL_lyrG1_v_RS7G4_DMG_CG, 
                                  !(geneID %in% LL_lyr_lowC_genes$V1))
LL_lyrG1_v_RS7G4_DMG_CHG <- filter(LL_lyrG1_v_RS7G4_DMG_CHG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))
LL_lyrG1_v_RS7G4_DMG_CHH <- filter(LL_lyrG1_v_RS7G4_DMG_CHH, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))

LL_RS7G1_v_RS7G4_DMG_CG_lyr <- filter(LL_RS7G1_v_RS7G4_DMG_CG_lyr, 
                                      !(geneID %in% LL_lyr_lowC_genes$V1))
LL_RS7G1_v_RS7G4_DMG_CHG_lyr <- filter(LL_RS7G1_v_RS7G4_DMG_CHG_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))
LL_RS7G1_v_RS7G4_DMG_CHH_lyr <- filter(LL_RS7G1_v_RS7G4_DMG_CHH_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))

LL_ALKG1_v_lyrG1_DMG_CG <- filter(LL_ALKG1_v_lyrG1_DMG_CG, 
                                  !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_lyrG1_DMG_CHG <- filter(LL_ALKG1_v_lyrG1_DMG_CHG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_lyrG1_DMG_CHH <- filter(LL_ALKG1_v_lyrG1_DMG_CHH, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))

LL_ALKG1_v_RS7G1_DMG_CG_lyr <- filter(LL_ALKG1_v_RS7G1_DMG_CG_lyr, 
                                      !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_RS7G1_DMG_CHG_lyr <- filter(LL_ALKG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_RS7G1_DMG_CHH_lyr <- filter(LL_ALKG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))

LL_ALKG1_v_RS7G4_DMG_CG_lyr <- filter(LL_ALKG1_v_RS7G4_DMG_CG_lyr, 
                                      !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_RS7G4_DMG_CHG_lyr <- filter(LL_ALKG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_RS7G4_DMG_CHH_lyr <- filter(LL_ALKG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))

LL_TKSG1_v_lyrG1_DMG_CG <- filter(LL_TKSG1_v_lyrG1_DMG_CG, 
                                  !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_lyrG1_DMG_CHG <- filter(LL_TKSG1_v_lyrG1_DMG_CHG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_lyrG1_DMG_CHH <- filter(LL_TKSG1_v_lyrG1_DMG_CHH, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))

LL_TKSG1_v_RS7G1_DMG_CG_lyr <- filter(LL_TKSG1_v_RS7G1_DMG_CG_lyr, 
                                      !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_RS7G1_DMG_CHG_lyr <- filter(LL_TKSG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_RS7G1_DMG_CHH_lyr <- filter(LL_TKSG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))

LL_TKSG1_v_RS7G4_DMG_CG_lyr <- filter(LL_TKSG1_v_RS7G4_DMG_CG_lyr, 
                                      !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_RS7G4_DMG_CHG_lyr <- filter(LL_TKSG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_RS7G4_DMG_CHH_lyr <- filter(LL_TKSG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(geneID %in% LL_lyr_lowC_genes$V1))

### Import expression data

setwd("/Users/ste/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/")

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- read.delim("HM_halG1_v_RS7G1_DEG.txt")
HM_halG1_v_RS7G4_DEG <- read.delim("HM_halG1_v_RS7G4_DEG.txt")

HM_RS7G1_v_RS7G4_hal_DEG <- read.delim("HM_RS7G1_v_G4_hal_DEG.txt")

HM_ALKG1_v_halG1_DEG <- read.delim("HM_ALKG1_v_halG1_DEG.txt")
HM_ALKG1_v_RS7G1_hal_DEG <- read.delim("HM_ALKG1_v_RS7G1_hal_DEG.txt")
HM_ALKG1_v_RS7G4_hal_DEG <- read.delim("HM_ALKG1_v_RS7G4_hal_DEG.txt")

HM_TKSG1_v_halG1_DEG <- read.delim("HM_TKSG1_v_halG1_DEG.txt")
HM_TKSG1_v_RS7G1_hal_DEG <- read.delim("HM_TKSG1_v_RS7G1_hal_DEG.txt")
HM_TKSG1_v_RS7G4_hal_DEG <- read.delim("HM_TKSG1_v_RS7G4_hal_DEG.txt")


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- read.delim("HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyrG1_v_RS7G4_DEG <- read.delim("HM_lyrG1_v_RS7G4_DEG.txt")

HM_RS7G1_v_RS7G4_lyr_DEG <- read.delim("HM_RS7G1_v_G4_lyr_DEG.txt")

HM_ALKG1_v_lyrG1_DEG <- read.delim("HM_ALKG1_v_lyrG1_DEG.txt")
HM_ALKG1_v_RS7G1_lyr_DEG <- read.delim("HM_ALKG1_v_RS7G1_lyr_DEG.txt")
HM_ALKG1_v_RS7G4_lyr_DEG <- read.delim("HM_ALKG1_v_RS7G4_lyr_DEG.txt")

HM_TKSG1_v_lyrG1_DEG <- read.delim("HM_TKSG1_v_lyrG1_DEG.txt")
HM_TKSG1_v_RS7G1_lyr_DEG <- read.delim("HM_TKSG1_v_RS7G1_lyr_DEG.txt")
HM_TKSG1_v_RS7G4_lyr_DEG <- read.delim("HM_TKSG1_v_RS7G4_lyr_DEG.txt")


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- read.delim("LL_halG1_v_RS7G1_DEG.txt")
LL_halG1_v_RS7G4_DEG <- read.delim("LL_halG1_v_RS7G4_DEG.txt")

LL_RS7G1_v_RS7G4_hal_DEG <- read.delim("LL_RS7G1_v_G4_hal_DEG.txt")

LL_ALKG1_v_halG1_DEG <- read.delim("LL_ALKG1_v_halG1_DEG.txt")
LL_ALKG1_v_RS7G1_hal_DEG <- read.delim("LL_ALKG1_v_RS7G1_hal_DEG.txt")
LL_ALKG1_v_RS7G4_hal_DEG <- read.delim("LL_ALKG1_v_RS7G4_hal_DEG.txt")

LL_TKSG1_v_halG1_DEG <- read.delim("LL_TKSG1_v_halG1_DEG.txt")
LL_TKSG1_v_RS7G1_hal_DEG <- read.delim("LL_TKSG1_v_RS7G1_hal_DEG.txt")
LL_TKSG1_v_RS7G4_hal_DEG <- read.delim("LL_TKSG1_v_RS7G4_hal_DEG.txt")

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- read.delim("LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyrG1_v_RS7G4_DEG <- read.delim("LL_lyrG1_v_RS7G4_DEG.txt")

LL_RS7G1_v_RS7G4_lyr_DEG <- read.delim("LL_RS7G1_v_G4_lyr_DEG.txt")

LL_ALKG1_v_lyrG1_DEG <- read.delim("LL_ALKG1_v_lyrG1_DEG.txt")
LL_ALKG1_v_RS7G1_lyr_DEG <- read.delim("LL_ALKG1_v_RS7G1_lyr_DEG.txt")
LL_ALKG1_v_RS7G4_lyr_DEG <- read.delim("LL_ALKG1_v_RS7G4_lyr_DEG.txt")

LL_TKSG1_v_lyrG1_DEG <- read.delim("LL_TKSG1_v_lyrG1_DEG.txt")
LL_TKSG1_v_RS7G1_lyr_DEG <- read.delim("LL_TKSG1_v_RS7G1_lyr_DEG.txt")
LL_TKSG1_v_RS7G4_lyr_DEG <- read.delim("LL_TKSG1_v_RS7G4_lyr_DEG.txt")

# We filter DEG falling in low coverage genes as well

# Filtering

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- filter(HM_halG1_v_RS7G1_DEG, 
                               !(geneID %in% HM_hal_lowC_genes$V1))
HM_halG1_v_RS7G4_DEG <- filter(HM_halG1_v_RS7G4_DEG, 
                               !(geneID %in% HM_hal_lowC_genes$V1))

HM_RS7G1_v_RS7G4_hal_DEG <- filter(HM_RS7G1_v_RS7G4_hal_DEG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))

HM_ALKG1_v_halG1_DEG <- filter(HM_ALKG1_v_halG1_DEG, 
                               !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_RS7G1_hal_DEG <- filter(HM_ALKG1_v_RS7G1_hal_DEG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))
HM_ALKG1_v_RS7G4_hal_DEG <- filter(HM_ALKG1_v_RS7G4_hal_DEG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))

HM_TKSG1_v_halG1_DEG <- filter(HM_TKSG1_v_halG1_DEG, 
                               !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_RS7G1_hal_DEG <- filter(HM_TKSG1_v_RS7G1_hal_DEG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))
HM_TKSG1_v_RS7G4_hal_DEG <- filter(HM_TKSG1_v_RS7G4_hal_DEG, 
                                   !(geneID %in% HM_hal_lowC_genes$V1))


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- filter(HM_lyrG1_v_RS7G1_DEG, 
                               !(geneID %in% HM_lyr_lowC_genes$V1))
HM_lyrG1_v_RS7G4_DEG <- filter(HM_lyrG1_v_RS7G4_DEG, 
                               !(geneID %in% HM_lyr_lowC_genes$V1))

HM_RS7G1_v_RS7G4_lyr_DEG <- filter(HM_RS7G1_v_RS7G4_lyr_DEG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))

HM_ALKG1_v_lyrG1_DEG <- filter(HM_ALKG1_v_lyrG1_DEG, 
                               !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_RS7G1_lyr_DEG <- filter(HM_ALKG1_v_RS7G1_lyr_DEG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))
HM_ALKG1_v_RS7G4_lyr_DEG <- filter(HM_ALKG1_v_RS7G4_lyr_DEG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))

HM_TKSG1_v_lyrG1_DEG <- filter(HM_TKSG1_v_lyrG1_DEG, 
                               !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_RS7G1_lyr_DEG <- filter(HM_TKSG1_v_RS7G1_lyr_DEG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))
HM_TKSG1_v_RS7G4_lyr_DEG <- filter(HM_TKSG1_v_RS7G4_lyr_DEG, 
                                   !(geneID %in% HM_lyr_lowC_genes$V1))


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- filter(LL_halG1_v_RS7G1_DEG, 
                               !(geneID %in% LL_hal_lowC_genes$V1))
LL_halG1_v_RS7G4_DEG <- filter(LL_halG1_v_RS7G4_DEG, 
                               !(geneID %in% LL_hal_lowC_genes$V1))

LL_RS7G1_v_RS7G4_hal_DEG <- filter(LL_RS7G1_v_RS7G4_hal_DEG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

LL_ALKG1_v_halG1_DEG <- filter(LL_ALKG1_v_halG1_DEG, 
                               !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_RS7G1_hal_DEG <- filter(LL_ALKG1_v_RS7G1_hal_DEG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))
LL_ALKG1_v_RS7G4_hal_DEG <- filter(LL_ALKG1_v_RS7G4_hal_DEG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

LL_TKSG1_v_halG1_DEG <- filter(LL_TKSG1_v_halG1_DEG, 
                               !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_RS7G1_hal_DEG <- filter(LL_TKSG1_v_RS7G1_hal_DEG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))
LL_TKSG1_v_RS7G4_hal_DEG <- filter(LL_TKSG1_v_RS7G4_hal_DEG, 
                                   !(geneID %in% LL_hal_lowC_genes$V1))

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- filter(LL_lyrG1_v_RS7G1_DEG, 
                               !(geneID %in% LL_lyr_lowC_genes$V1))
LL_lyrG1_v_RS7G4_DEG <- filter(LL_lyrG1_v_RS7G4_DEG, 
                               !(geneID %in% LL_lyr_lowC_genes$V1))

LL_RS7G1_v_RS7G4_lyr_DEG <- filter(LL_RS7G1_v_RS7G4_lyr_DEG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))

LL_ALKG1_v_lyrG1_DEG <- filter(LL_ALKG1_v_lyrG1_DEG, 
                               !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_RS7G1_lyr_DEG <- filter(LL_ALKG1_v_RS7G1_lyr_DEG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))
LL_ALKG1_v_RS7G4_lyr_DEG <- filter(LL_ALKG1_v_RS7G4_lyr_DEG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))

LL_TKSG1_v_lyrG1_DEG <- filter(LL_TKSG1_v_lyrG1_DEG, 
                               !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_RS7G1_lyr_DEG <- filter(LL_TKSG1_v_RS7G1_lyr_DEG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))
LL_TKSG1_v_RS7G4_lyr_DEG <- filter(LL_TKSG1_v_RS7G4_lyr_DEG, 
                                   !(geneID %in% LL_lyr_lowC_genes$V1))


### Consolidate methylation data
### 1) merge all contexts
### 2) remove duplicated genes

consolidate <- function(DMG1, DMG2, DMG3){
  new_DMG <- rbind(DMG1, DMG2, DMG3)
  new_DMG <- arrange(new_DMG, geneID)
  new_DMG <- new_DMG[!duplicated(new_DMG$geneID),]
  return(new_DMG)
}

### Cold conditions

## halleri side

HM_halG1_v_RS7G1_DMG <- consolidate(HM_halG1_v_RS7G1_DMG_CG,
                                    HM_halG1_v_RS7G1_DMG_CHG,
                                    HM_halG1_v_RS7G1_DMG_CHH)

HM_halG1_v_RS7G4_DMG <- consolidate(HM_halG1_v_RS7G4_DMG_CG,
                                    HM_halG1_v_RS7G4_DMG_CHG,
                                    HM_halG1_v_RS7G4_DMG_CHH)

HM_RS7G1_v_RS7G4_DMG_hal <- consolidate(HM_RS7G1_v_RS7G4_DMG_CG_hal,
                                    HM_RS7G1_v_RS7G4_DMG_CHG_hal,
                                    HM_RS7G1_v_RS7G4_DMG_CHH_hal)

HM_ALKG1_v_halG1_DMG <- consolidate(HM_ALKG1_v_halG1_DMG_CG,
                                    HM_ALKG1_v_halG1_DMG_CHG,
                                    HM_ALKG1_v_halG1_DMG_CHH)

HM_ALKG1_v_RS7G1_DMG_hal <- consolidate(HM_ALKG1_v_RS7G1_DMG_CG_hal,
                                    HM_ALKG1_v_RS7G1_DMG_CHG_hal,
                                    HM_ALKG1_v_RS7G1_DMG_CHH_hal)

HM_ALKG1_v_RS7G4_DMG_hal <- consolidate(HM_ALKG1_v_RS7G4_DMG_CG_hal,
                                    HM_ALKG1_v_RS7G4_DMG_CHG_hal,
                                    HM_ALKG1_v_RS7G4_DMG_CHH_hal)

HM_TKSG1_v_halG1_DMG <- consolidate(HM_TKSG1_v_halG1_DMG_CG,
                                    HM_TKSG1_v_halG1_DMG_CHG,
                                    HM_TKSG1_v_halG1_DMG_CHH)

HM_TKSG1_v_RS7G1_DMG_hal <- consolidate(HM_TKSG1_v_RS7G1_DMG_CG_hal,
                                    HM_TKSG1_v_RS7G1_DMG_CHG_hal,
                                    HM_TKSG1_v_RS7G1_DMG_CHH_hal)

HM_TKSG1_v_RS7G4_DMG_hal <- consolidate(HM_TKSG1_v_RS7G4_DMG_CG_hal,
                                    HM_TKSG1_v_RS7G4_DMG_CHG_hal,
                                    HM_TKSG1_v_RS7G4_DMG_CHH_hal)

## lyrata side

HM_lyrG1_v_RS7G1_DMG <- consolidate(HM_lyrG1_v_RS7G1_DMG_CG,
                                    HM_lyrG1_v_RS7G1_DMG_CHG,
                                    HM_lyrG1_v_RS7G1_DMG_CHH)

HM_lyrG1_v_RS7G4_DMG <- consolidate(HM_lyrG1_v_RS7G4_DMG_CG,
                                    HM_lyrG1_v_RS7G4_DMG_CHG,
                                    HM_lyrG1_v_RS7G4_DMG_CHH)

HM_RS7G1_v_RS7G4_DMG_lyr <- consolidate(HM_RS7G1_v_RS7G4_DMG_CG_lyr,
                                        HM_RS7G1_v_RS7G4_DMG_CHG_lyr,
                                        HM_RS7G1_v_RS7G4_DMG_CHH_lyr)

HM_ALKG1_v_lyrG1_DMG <- consolidate(HM_ALKG1_v_lyrG1_DMG_CG,
                                    HM_ALKG1_v_lyrG1_DMG_CHG,
                                    HM_ALKG1_v_lyrG1_DMG_CHH)

HM_ALKG1_v_RS7G1_DMG_lyr <- consolidate(HM_ALKG1_v_RS7G1_DMG_CG_lyr,
                                        HM_ALKG1_v_RS7G1_DMG_CHG_lyr,
                                        HM_ALKG1_v_RS7G1_DMG_CHH_lyr)

HM_ALKG1_v_RS7G4_DMG_lyr <- consolidate(HM_ALKG1_v_RS7G4_DMG_CG_lyr,
                                        HM_ALKG1_v_RS7G4_DMG_CHG_lyr,
                                        HM_ALKG1_v_RS7G4_DMG_CHH_lyr)

HM_TKSG1_v_lyrG1_DMG <- consolidate(HM_TKSG1_v_lyrG1_DMG_CG,
                                    HM_TKSG1_v_lyrG1_DMG_CHG,
                                    HM_TKSG1_v_lyrG1_DMG_CHH)

HM_TKSG1_v_RS7G1_DMG_lyr <- consolidate(HM_TKSG1_v_RS7G1_DMG_CG_lyr,
                                        HM_TKSG1_v_RS7G1_DMG_CHG_lyr,
                                        HM_TKSG1_v_RS7G1_DMG_CHH_lyr)

HM_TKSG1_v_RS7G4_DMG_lyr <- consolidate(HM_TKSG1_v_RS7G4_DMG_CG_lyr,
                                        HM_TKSG1_v_RS7G4_DMG_CHG_lyr,
                                        HM_TKSG1_v_RS7G4_DMG_CHH_lyr)


### Hot conditions

## halleri side

LL_halG1_v_RS7G1_DMG <- consolidate(LL_halG1_v_RS7G1_DMG_CG,
                                    LL_halG1_v_RS7G1_DMG_CHG,
                                    LL_halG1_v_RS7G1_DMG_CHH)

LL_halG1_v_RS7G4_DMG <- consolidate(LL_halG1_v_RS7G4_DMG_CG,
                                    LL_halG1_v_RS7G4_DMG_CHG,
                                    LL_halG1_v_RS7G4_DMG_CHH)

LL_RS7G1_v_RS7G4_DMG_hal <- consolidate(LL_RS7G1_v_RS7G4_DMG_CG_hal,
                                        LL_RS7G1_v_RS7G4_DMG_CHG_hal,
                                        LL_RS7G1_v_RS7G4_DMG_CHH_hal)

LL_ALKG1_v_halG1_DMG <- consolidate(LL_ALKG1_v_halG1_DMG_CG,
                                    LL_ALKG1_v_halG1_DMG_CHG,
                                    LL_ALKG1_v_halG1_DMG_CHH)

LL_ALKG1_v_RS7G1_DMG_hal <- consolidate(LL_ALKG1_v_RS7G1_DMG_CG_hal,
                                        LL_ALKG1_v_RS7G1_DMG_CHG_hal,
                                        LL_ALKG1_v_RS7G1_DMG_CHH_hal)

LL_ALKG1_v_RS7G4_DMG_hal <- consolidate(LL_ALKG1_v_RS7G4_DMG_CG_hal,
                                        LL_ALKG1_v_RS7G4_DMG_CHG_hal,
                                        LL_ALKG1_v_RS7G4_DMG_CHH_hal)

LL_TKSG1_v_halG1_DMG <- consolidate(LL_TKSG1_v_halG1_DMG_CG,
                                    LL_TKSG1_v_halG1_DMG_CHG,
                                    LL_TKSG1_v_halG1_DMG_CHH)

LL_TKSG1_v_RS7G1_DMG_hal <- consolidate(LL_TKSG1_v_RS7G1_DMG_CG_hal,
                                        LL_TKSG1_v_RS7G1_DMG_CHG_hal,
                                        LL_TKSG1_v_RS7G1_DMG_CHH_hal)

LL_TKSG1_v_RS7G4_DMG_hal <- consolidate(LL_TKSG1_v_RS7G4_DMG_CG_hal,
                                        LL_TKSG1_v_RS7G4_DMG_CHG_hal,
                                        LL_TKSG1_v_RS7G4_DMG_CHH_hal)

## lyrata side

LL_lyrG1_v_RS7G1_DMG <- consolidate(LL_lyrG1_v_RS7G1_DMG_CG,
                                    LL_lyrG1_v_RS7G1_DMG_CHG,
                                    LL_lyrG1_v_RS7G1_DMG_CHH)

LL_lyrG1_v_RS7G4_DMG <- consolidate(LL_lyrG1_v_RS7G4_DMG_CG,
                                    LL_lyrG1_v_RS7G4_DMG_CHG,
                                    LL_lyrG1_v_RS7G4_DMG_CHH)

LL_RS7G1_v_RS7G4_DMG_lyr <- consolidate(LL_RS7G1_v_RS7G4_DMG_CG_lyr,
                                        LL_RS7G1_v_RS7G4_DMG_CHG_lyr,
                                        LL_RS7G1_v_RS7G4_DMG_CHH_lyr)

LL_ALKG1_v_lyrG1_DMG <- consolidate(LL_ALKG1_v_lyrG1_DMG_CG,
                                    LL_ALKG1_v_lyrG1_DMG_CHG,
                                    LL_ALKG1_v_lyrG1_DMG_CHH)

LL_ALKG1_v_RS7G1_DMG_lyr <- consolidate(LL_ALKG1_v_RS7G1_DMG_CG_lyr,
                                        LL_ALKG1_v_RS7G1_DMG_CHG_lyr,
                                        LL_ALKG1_v_RS7G1_DMG_CHH_lyr)

LL_ALKG1_v_RS7G4_DMG_lyr <- consolidate(LL_ALKG1_v_RS7G4_DMG_CG_lyr,
                                        LL_ALKG1_v_RS7G4_DMG_CHG_lyr,
                                        LL_ALKG1_v_RS7G4_DMG_CHH_lyr)

LL_TKSG1_v_lyrG1_DMG <- consolidate(LL_TKSG1_v_lyrG1_DMG_CG,
                                    LL_TKSG1_v_lyrG1_DMG_CHG,
                                    LL_TKSG1_v_lyrG1_DMG_CHH)

LL_TKSG1_v_RS7G1_DMG_lyr <- consolidate(LL_TKSG1_v_RS7G1_DMG_CG_lyr,
                                        LL_TKSG1_v_RS7G1_DMG_CHG_lyr,
                                        LL_TKSG1_v_RS7G1_DMG_CHH_lyr)

LL_TKSG1_v_RS7G4_DMG_lyr <- consolidate(LL_TKSG1_v_RS7G4_DMG_CG_lyr,
                                        LL_TKSG1_v_RS7G4_DMG_CHG_lyr,
                                        LL_TKSG1_v_RS7G4_DMG_CHH_lyr)

### Plot overlaps

setwd("~/OneDrive/PhD/Project/Chapter_3/Pictures/Paper_picsV5")

# Prepare a palette of 4 colors with R colorbrewer:
myCol <- brewer.pal(4, "Paired")

# Function to output venn diagram for two sets
venn <- function(set1, set2, title, filename){
  venn.diagram(
  x = list(set1, set2),
  category.names = c("DEG" , "DMG"),
  cat.cex = 0.7,
  cat.pos = 1,
  filename = filename,
  output=TRUE,
  main = title,
  main.cex = 0.5,
  print.mode = c("raw", "percent"),

  # Output features
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 400,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],

  # Numbers
  cex = 0.5,
  fontface = "bold"
)
}

# Function to output venn diagram for four sets

venn4 <- function(set1, set2, set3, set4, title, filename){
  venn.diagram(
    x = list(set1, set2, set3, set4),
    category.names = c("DEG" , "DMG", "DEG" , "DMG"),
    cat.cex = 0.7,
    cat.pos = 1,
    filename = filename,
    output=TRUE,
    main = title,
    main.cex = 0.5,
    print.mode = c("raw", "percent"),
    
    # Output features
    imagetype="svg" ,
    height = 480 , 
    width = 480 , 
    resolution = 400,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol[1:4],
    
    # Numbers
    cex = 0.5,
    fontface = "bold"
  )
}

# venn_DMG <- function(set1, set2, title, filename){
#   venn.diagram(
#     x = list(set1, set2),
#     category.names = c("DMG lyr" , "DMG hal"),
#     cat.cex = 0.7,
#     cat.pos = 1,
#     filename = filename,
#     output=TRUE,
#     main = title,
#     main.cex = 0.5,
#     
#     # Output features
#     imagetype="png" ,
#     height = 480 , 
#     width = 480 , 
#     resolution = 400,
#     compression = "lzw",
#     
#     # Circles
#     lwd = 2,
#     lty = 'blank',
#     fill = myCol[1:2],
#     
#     # Numbers
#     cex = 0.5,
#     fontface = "bold"
#   )
# }

## HM halleri G1 vs synthetic G1

venn_HM_hal1Vsyn1_DEG_v_DMG <- venn(set1 = HM_halG1_v_RS7G1_DEG$geneID,
                                    set2 = HM_halG1_v_RS7G1_DMG$geneID,
                                    title = "HM_halG1_v_RS7G1",
                                    filename = NULL)

ggsave(venn_HM_hal1Vsyn1_DEG_v_DMG, 
       filename = "HM_hal1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM halleri DMG vs lyrata DEG G1 vs synthetic G1

# venn_HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G1_DEG$geneID,
#                                              set2 = HM_halG1_v_RS7G1_DMG$geneID,
#                                              title = "HM_lyrDEGhalDMGG1_v_RS7G1",
#                                              filename = NULL)
# 
# ggsave(venn_HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG, 
#        filename = "HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM lyrata DMG vs halleri DEG G1 vs synthetic G1

# venn_HM_lyrDMGhalDEG1Vsyn1_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G1_DMG$geneID,
#      set2 = HM_halG1_v_RS7G1_DEG$geneID,
#      title = "HM_lyrDMGhalDEGG1_v_RS7G1",
#      filename = NULL)
# 
# ggsave(venn_HM_lyrDMGhalDEG1Vsyn1_DEG_v_DMG, 
#        filename = "HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM halleri DMG vs lyrata DMG (G1)

# venn_HM_lyrDMGhalDMG1Vsyn1 <- venn_DMG(set1 = HM_lyrG1_v_RS7G1_DMG$geneID,
#      set2 = HM_halG1_v_RS7G1_DMG$geneID,
#      title = "HM_lyrhalDMGG1_v_RS7G1",
#      filename = NULL)
# 
# ggsave(venn_HM_lyrDMGhalDMG1Vsyn1, 
#        filename = "HM_lyrDMGhalDMG1Vsyn1.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM halleri DMG vs lyrata DMG (G4)

# venn_HM_lyrDMGhalDMG1Vsyn4 <- venn_DMG(set1 = HM_lyrG1_v_RS7G4_DMG$geneID,
#      set2 = HM_halG1_v_RS7G4_DMG$geneID,
#      title = "HM_lyrhalDMGG4_v_RS7G1",
#      filename = NULL)
# 
# ggsave(venn_HM_lyrDMGhalDMG1Vsyn4, 
#        filename = "HM_lyrDMGhalDMG1Vsyn4.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL halleri DMG vs lyrata DMG (G1)

# venn_LL_lyrDMGhalDMG1Vsyn1 <- venn_DMG(set1 = LL_lyrG1_v_RS7G1_DMG$geneID,
#          set2 = LL_halG1_v_RS7G1_DMG$geneID,
#          title = "LL_lyrhalDMGG1_v_RS7G1",
#          filename = NULL)
# 
# ggsave(venn_LL_lyrDMGhalDMG1Vsyn1, 
#        filename = "LL_lyrDMGhalDMG1Vsyn1.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL halleri DMG vs lyrata DMG (G4)

# venn_LL_lyrDMGhalDMG1Vsyn4 <- venn_DMG(set1 = LL_lyrG1_v_RS7G4_DMG$geneID,
#          set2 = LL_halG1_v_RS7G4_DMG$geneID,
#          title = "LL_lyrhalDMGG4_v_RS7G1",
#          filename = NULL)
# 
# ggsave(venn_LL_lyrDMGhalDMG1Vsyn4, 
#        filename = "LL_lyrDMGhalDMG1Vsyn4.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM halleri G1 vs synthetic G4

venn_HM_hal1Vsyn4_DEG_v_DMG <- venn(set1 = HM_halG1_v_RS7G4_DEG$geneID,
     set2 = HM_halG1_v_RS7G4_DMG$geneID,
     title = "HM_halG1_v_RS7G4",
     filename = NULL)

ggsave(venn_HM_hal1Vsyn4_DEG_v_DMG, 
       filename = "HM_hal1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM lyrata G1 vs synthetic G1

venn_HM_lyr1Vsyn1_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G1_DEG$geneID,
     set2 = HM_lyrG1_v_RS7G1_DMG$geneID,
     title = "HM_lyrG1_v_RS7G1",
     filename = NULL)

ggsave(venn_HM_lyr1Vsyn1_DEG_v_DMG, 
       filename = "HM_lyr1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM lyrata G1 vs synthetic G4

venn_HM_lyr1Vsyn4_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G4_DEG$geneID,
     set2 = HM_lyrG1_v_RS7G4_DMG$geneID,
     title = "HM_lyrG1_v_RS7G4",
     filename = NULL)

ggsave(venn_HM_lyr1Vsyn4_DEG_v_DMG, 
       filename = "HM_lyr1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL halleri G1 vs synthetic G1

venn_LL_hal1Vsyn1_DEG_v_DMG <- venn(set1 = LL_halG1_v_RS7G1_DEG$geneID,
     set2 = LL_halG1_v_RS7G1_DMG$geneID,
     title = "LL_halG1_v_RS7G1",
     filename = NULL)

ggsave(venn_LL_hal1Vsyn1_DEG_v_DMG, 
       filename = "LL_hal1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL halleri G1 vs synthetic G4

venn_LL_hal1Vsyn4_DEG_v_DMG <- venn(set1 = LL_halG1_v_RS7G4_DEG$geneID,
     set2 = LL_halG1_v_RS7G4_DMG$geneID,
     title = "LL_halG1_v_RS7G4",
     filename = NULL)

ggsave(venn_LL_hal1Vsyn4_DEG_v_DMG, 
       filename = "LL_hal1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL lyrata G1 vs synthetic G1

venn_LL_lyr1Vsyn1_DEG_v_DMG <- venn(set1 = LL_lyrG1_v_RS7G1_DEG$geneID,
     set2 = LL_lyrG1_v_RS7G1_DMG$geneID,
     title = "LL_lyrG1_v_RS7G1",
     filename = NULL)

ggsave(venn_LL_lyr1Vsyn1_DEG_v_DMG, 
       filename = "LL_lyr1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL lyrata G1 vs synthetic G4

venn_LL_lyr1Vsyn4_DEG_v_DMG <- venn(set1 = LL_lyrG1_v_RS7G4_DEG$geneID,
     set2 = LL_lyrG1_v_RS7G4_DMG$geneID,
     title = "LL_lyrG1_v_RS7G4",
     filename = NULL)

ggsave(venn_LL_lyr1Vsyn4_DEG_v_DMG, 
       filename = "LL_lyr1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM synthetic G1 vs synthetic G4

# venn_HM_syn1Vsyn4_DEG_v_DMG <- venn(set1 = HM_RS7G1_v_RS7G4_hal_DEG$geneID,
#      set2 = HM_RS7G1_v_RS7G4_DMG_hal$geneID,
#      title = "HM_RS7G1_v_RS7G4 (hal)",
#      filename = NULL)
# 
# ggsave(venn_HM_syn1Vsyn4_DEG_v_DMG, 
#        filename = "HM_syn1Vsyn4_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM ALK G1 vs pro G1 (hal)

# venn_HM_alk1Vhal1_DEG_v_DMG <- venn(set1 = HM_ALKG1_v_halG1_DEG$geneID,
#      set2 = HM_ALKG1_v_halG1_DMG$geneID,
#      title = "HM_ALKG1_v_halG1",
#      filename = NULL)
# 
# ggsave(venn_HM_alk1Vhal1_DEG_v_DMG, 
#        filename = "HM_alk1Vhal1_DEG_v_DMG.svvg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM ALK G1 vs synthetic G1 (hal)

# venn_HM_alk1Vsyn1_DEG_v_DMG_hal <- venn(set1 = HM_ALKG1_v_RS7G1_hal_DEG$geneID,
#      set2 = HM_ALKG1_v_RS7G1_DMG_hal$geneID,
#      title = "HM_ALKG1_v_RS7G1 (hal)",
#      filename = NULL)
# 
# ggsave(venn_HM_alk1Vsyn1_DEG_v_DMG_hal, 
#        filename = "HM_alk1Vsyn1_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM ALK G1 vs synthetic G4 (hal)

# venn_HM_alk1Vsyn4_DEG_v_DMG_hal <- venn(set1 = HM_ALKG1_v_RS7G4_hal_DEG$geneID,
#      set2 = HM_ALKG1_v_RS7G4_DMG_hal$geneID,
#      title = "HM_ALKG1_v_RS7G4 (hal)",
#      filename = NULL)
# 
# ggsave(venn_HM_alk1Vsyn4_DEG_v_DMG_hal, 
#        filename = "HM_alk1Vsyn4_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM TKS G1 vs pro G1 (hal)

# venn_HM_tks1Vhal1_DEG_v_DMG <- venn(set1 = HM_TKSG1_v_halG1_DEG$geneID,
#      set2 = HM_TKSG1_v_halG1_DMG$geneID,
#      title = "HM_TKSG1_v_halG1",
#      filename = NULL)
# 
# ggsave(venn_HM_tks1Vhal1_DEG_v_DMG, 
#        filename = "HM_tks1Vhal1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM TKS G1 vs synthetic G1 (hal)

# venn_HM_tks1Vsyn1_DEG_v_DMG_hal <- venn(set1 = HM_TKSG1_v_RS7G1_hal_DEG$geneID,
#      set2 = HM_TKSG1_v_RS7G1_DMG_hal$geneID,
#      title = "HM_TKSG1_v_RS7G1 (hal)",
#      filename = NULL)
# 
# ggsave(venn_HM_tks1Vsyn1_DEG_v_DMG_hal, 
#        filename = "HM_tks1Vsyn1_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM TKS G1 vs synthetic G4 (hal)

# venn_HM_tks1Vsyn4_DEG_v_DMG_hal <- venn(set1 = HM_TKSG1_v_RS7G4_hal_DEG$geneID,
#      set2 = HM_TKSG1_v_RS7G4_DMG_hal$geneID,
#      title = "HM_TKSG1_v_RS7G4 (hal)",
#      filename = NULL)
# 
# ggsave(venn_HM_tks1Vsyn4_DEG_v_DMG_hal, 
#        filename = "HM_tks1Vsyn4_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL ALK G1 vs pro G1 (hal)

# venn_LL_alk1Vhal1_DEG_v_DMG <- venn(set1 = LL_ALKG1_v_halG1_DEG$geneID,
#      set2 = LL_ALKG1_v_halG1_DMG$geneID,
#      title = "LL_ALKG1_v_halG1",
#      filename = NULL)
# 
# ggsave(venn_LL_alk1Vhal1_DEG_v_DMG, 
#        filename = "LL_alk1Vhal1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL ALK G1 vs synthetic G1 (hal)

# venn_LL_alk1Vsyn1_DEG_v_DMG_hal <- venn(set1 = LL_ALKG1_v_RS7G1_hal_DEG$geneID,
#      set2 = LL_ALKG1_v_RS7G1_DMG_hal$geneID,
#      title = "LL_ALKG1_v_RS7G1 (hal)",
#      filename = NULL)
# 
# ggsave(venn_LL_alk1Vsyn1_DEG_v_DMG_hal, 
#        filename = "LL_alk1Vsyn1_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL ALK G1 vs synthetic G4 (hal)

# venn_LL_alk1Vsyn4_DEG_v_DMG_hal <- venn(set1 = LL_ALKG1_v_RS7G4_hal_DEG$geneID,
#      set2 = LL_ALKG1_v_RS7G4_DMG_hal$geneID,
#      title = "LL_ALKG1_v_RS7G4 (hal)",
#      filename = NULL)
# 
# ggsave(venn_LL_alk1Vsyn4_DEG_v_DMG_hal, 
#        filename = "LL_alk1Vsyn4_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL TKS G1 vs pro G1 (hal)

# venn_LL_tks1Vhal1_DEG_v_DMG <- venn(set1 = LL_TKSG1_v_halG1_DEG$geneID,
#      set2 = LL_TKSG1_v_halG1_DMG$geneID,
#      title = "LL_TKSG1_v_halG1",
#      filename = NULL)
# 
# ggsave(venn_LL_tks1Vhal1_DEG_v_DMG, 
#        filename = "LL_tks1Vhal1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL TKS G1 vs synthetic G1 (hal)

# venn_LL_tks1Vsyn1_DEG_v_DMG_hal <- venn(set1 = LL_TKSG1_v_RS7G1_hal_DEG$geneID,
#      set2 = LL_TKSG1_v_RS7G1_DMG_hal$geneID,
#      title = "LL_TKSG1_v_RS7G1 (hal)",
#      filename = NULL)
# 
# ggsave(venn_LL_tks1Vsyn1_DEG_v_DMG_hal, 
#        filename = "LL_tks1Vsyn1_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL TKS G1 vs synthetic G4 (hal)

# venn_LL_tks1Vsyn4_DEG_v_DMG_hal <- venn(set1 = LL_TKSG1_v_RS7G4_hal_DEG$geneID,
#      set2 = LL_TKSG1_v_RS7G4_DMG_hal$geneID,
#      title = "LL_TKSG1_v_RS7G4 (hal)",
#      filename = NULL)
# 
# ggsave(venn_LL_tks1Vsyn4_DEG_v_DMG_hal, 
#        filename = "LL_tks1Vsyn4_DEG_v_DMG_hal.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM ALK G1 vs pro G1 (lyr)

# venn_HM_alk1Vlyr1_DEG_v_DMG <- venn(set1 = HM_ALKG1_v_lyrG1_DEG$geneID,
#      set2 = HM_ALKG1_v_lyrG1_DMG$geneID,
#      title = "HM_ALKG1_v_lyrG1",
#      filename = NULL)
# 
# ggsave(venn_HM_alk1Vlyr1_DEG_v_DMG, 
#        filename = "HM_alk1Vlyr1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM ALK G1 vs synthetic G1 (lyr)

# venn_HM_alk1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = HM_ALKG1_v_RS7G1_lyr_DEG$geneID,
#      set2 = HM_ALKG1_v_RS7G1_DMG_lyr$geneID,
#      title = "HM_ALKG1_v_RS7G1 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_HM_alk1Vsyn1_DEG_v_DMG_lyr, 
#        filename = "HM_alk1Vsyn1_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM ALK G1 vs synthetic G4 (lyr)

# venn_HM_alk1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = HM_ALKG1_v_RS7G4_lyr_DEG$geneID,
#      set2 = HM_ALKG1_v_RS7G4_DMG_lyr$geneID,
#      title = "HM_ALKG1_v_RS7G4 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_HM_alk1Vsyn4_DEG_v_DMG_lyr, 
#        filename = "HM_alk1Vsyn4_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM TKS G1 vs pro G1 (lyr)

# venn_HM_tks1Vlyr1_DEG_v_DMG <- venn(set1 = HM_TKSG1_v_lyrG1_DEG$geneID,
#      set2 = HM_TKSG1_v_lyrG1_DMG$geneID,
#      title = "HM_TKSG1_v_lyrG1",
#      filename = NULL)
# 
# ggsave(venn_HM_tks1Vlyr1_DEG_v_DMG, 
#        filename = "HM_tks1Vlyr1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM TKS G1 vs synthetic G1 (lyr)

# venn_HM_tks1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = HM_TKSG1_v_RS7G1_lyr_DEG$geneID,
#      set2 = HM_TKSG1_v_RS7G1_DMG_lyr$geneID,
#      title = "HM_TKSG1_v_RS7G1 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_HM_tks1Vsyn1_DEG_v_DMG_lyr, 
#        filename = "HM_tks1Vsyn1_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## HM TKS G1 vs synthetic G4 (lyr)

# venn_HM_tks1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = HM_TKSG1_v_RS7G4_lyr_DEG$geneID,
#      set2 = HM_TKSG1_v_RS7G4_DMG_lyr$geneID,
#      title = "HM_TKSG1_v_RS7G4 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_HM_tks1Vsyn4_DEG_v_DMG_lyr, 
#        filename = "HM_tks1Vsyn4_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL ALK G1 vs pro G1 (lyr)

# venn_LL_alk1Vlyr1_DEG_v_DMG <- venn(set1 = LL_ALKG1_v_lyrG1_DEG$geneID,
#      set2 = LL_ALKG1_v_lyrG1_DMG$geneID,
#      title = "LL_ALKG1_v_lyrG1",
#      filename = NULL)
# 
# ggsave(venn_LL_alk1Vlyr1_DEG_v_DMG, 
#        filename = "LL_alk1Vlyr1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL ALK G1 vs synthetic G1 (lyr)

# venn_LL_alk1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = LL_ALKG1_v_RS7G1_lyr_DEG$geneID,
#      set2 = LL_ALKG1_v_RS7G1_DMG_lyr$geneID,
#      title = "LL_ALKG1_v_RS7G1 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_LL_alk1Vsyn1_DEG_v_DMG_lyr, 
#        filename = "LL_alk1Vsyn1_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL ALK G1 vs synthetic G4 (lyr)

# venn_LL_alk1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = LL_ALKG1_v_RS7G4_lyr_DEG$geneID,
#      set2 = LL_ALKG1_v_RS7G4_DMG_lyr$geneID,
#      title = "LL_ALKG1_v_RS7G4 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_LL_alk1Vsyn4_DEG_v_DMG_lyr, 
#        filename = "LL_alk1Vsyn4_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL TKS G1 vs pro G1 (lyr)

# venn_LL_tks1Vlyr1_DEG_v_DMG <- venn(set1 = LL_TKSG1_v_lyrG1_DEG$geneID,
#      set2 = LL_TKSG1_v_lyrG1_DMG$geneID,
#      title = "LL_TKSG1_v_lyrG1",
#      filename = NULL)
# 
# ggsave(venn_LL_tks1Vlyr1_DEG_v_DMG, 
#        filename = "LL_tks1Vlyr1_DEG_v_DMG.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL TKS G1 vs synthetic G1 (lyr)

# venn_LL_tks1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = LL_TKSG1_v_RS7G1_lyr_DEG$geneID,
#      set2 = LL_TKSG1_v_RS7G1_DMG_lyr$geneID,
#      title = "LL_TKSG1_v_RS7G1 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_LL_tks1Vsyn1_DEG_v_DMG_lyr, 
#        filename = "LL_tks1Vsyn1_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)

## LL TKS G1 vs synthetic G4 (lyr)

# venn_LL_tks1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = LL_TKSG1_v_RS7G4_lyr_DEG$geneID,
#      set2 = LL_TKSG1_v_RS7G4_DMG_lyr$geneID,
#      title = "LL_TKSG1_v_RS7G4 (lyr)",
#      filename = NULL)
# 
# ggsave(venn_LL_tks1Vsyn4_DEG_v_DMG_lyr, 
#        filename = "LL_tks1Vsyn4_DEG_v_DMG_lyr.svg",
#        device = "svg",
#        width = 1.5,
#        height = 1.5)


## HM halleri G1 vs synthetic G1/4

# venn4(set1 = HM_halG1_v_RS7G1_DEG$geneID,
#      set2 = HM_halG1_v_RS7G1_DMG$geneID,
#      set3 = HM_halG1_v_RS7G4_DEG$geneID,
#      set4 = HM_halG1_v_RS7G4_DMG$geneID,
#      title = "Cold conditions - halleri side",
#      filename = "HM_hal1Vsyn14_DEG_v_DMG.png")



## Other things below - detailed numbers for further insights

## Checking overlap across generations for DEGs and DMGs

## First look at diploid vs synthetic overlap

HM_pro1Vsyn1_halDEG <- sum(!is.na(match(HM_halG1_v_RS7G1_DEG$geneID, HM_halG1_v_RS7G4_DEG$geneID)))
paste0(round(HM_pro1Vsyn1_halDEG/length(HM_halG1_v_RS7G1_DEG$geneID)*100, 2), "%")

HM_pro1Vsyn1_halDMG <- sum(!is.na(match(HM_halG1_v_RS7G1_DMG$geneID, HM_halG1_v_RS7G4_DMG$geneID)))
paste0(round(HM_pro1Vsyn1_halDMG/length(HM_halG1_v_RS7G1_DMG$geneID)*100, 2), "%")

HM_pro1Vsyn1_lyrDEG <- sum(!is.na(match(HM_lyrG1_v_RS7G1_DEG$geneID, HM_lyrG1_v_RS7G4_DEG$geneID)))
paste0(round(HM_pro1Vsyn1_lyrDEG/length(HM_lyrG1_v_RS7G1_DEG$geneID)*100, 2), "%")

HM_pro1Vsyn1_lyrDMG <- sum(!is.na(match(HM_lyrG1_v_RS7G1_DMG$geneID, HM_lyrG1_v_RS7G4_DMG$geneID)))
paste0(round(HM_pro1Vsyn1_lyrDMG/length(HM_lyrG1_v_RS7G1_DMG$geneID)*100, 2), "%")

LL_pro1Vsyn1_halDEG <- sum(!is.na(match(LL_halG1_v_RS7G1_DEG$geneID, LL_halG1_v_RS7G4_DEG$geneID)))
paste0(round(LL_pro1Vsyn1_halDEG/length(LL_halG1_v_RS7G1_DEG$geneID)*100, 2), "%")

LL_pro1Vsyn1_halDMG <- sum(!is.na(match(LL_halG1_v_RS7G1_DMG$geneID, LL_halG1_v_RS7G4_DMG$geneID)))
paste0(round(LL_pro1Vsyn1_halDMG/length(LL_halG1_v_RS7G1_DMG$geneID)*100, 2), "%")

LL_pro1Vsyn1_lyrDEG <- sum(!is.na(match(LL_lyrG1_v_RS7G1_DEG$geneID, LL_lyrG1_v_RS7G4_DEG$geneID)))
paste0(round(LL_pro1Vsyn1_lyrDEG/length(LL_lyrG1_v_RS7G1_DEG$geneID)*100, 2), "%")

LL_pro1Vsyn1_lyrDMG <- sum(!is.na(match(LL_lyrG1_v_RS7G1_DMG$geneID, LL_lyrG1_v_RS7G4_DMG$geneID)))
paste0(round(LL_pro1Vsyn1_lyrDMG/length(LL_lyrG1_v_RS7G1_DMG$geneID)*100, 2), "%")

## Then we look at ALK vs all overlap

HM_ALKvpro1syn1_halDEG <- sum(!is.na(match(HM_ALKG1_v_halG1_DEG$geneID, HM_ALKG1_v_RS7G1_hal_DEG$geneID)))
paste0(round(HM_ALKvpro1syn1_halDEG/length(HM_ALKG1_v_halG1_DEG$geneID)*100, 2), "%")

HM_ALKvpro1syn1_halDMG <- sum(!is.na(match(HM_ALKG1_v_halG1_DMG$geneID, HM_ALKG1_v_RS7G1_DMG_hal$geneID)))
paste0(round(HM_ALKvpro1syn1_halDMG/length(HM_ALKG1_v_halG1_DMG$geneID)*100, 2), "%")

HM_ALKvpro1syn1_lyrDEG <- sum(!is.na(match(HM_ALKG1_v_lyrG1_DEG$geneID, HM_ALKG1_v_RS7G1_lyr_DEG$geneID)))
paste0(round(HM_ALKvpro1syn1_lyrDEG/length(HM_ALKG1_v_lyrG1_DEG$geneID)*100, 2), "%")

HM_ALKvpro1syn1_lyrDMG <- sum(!is.na(match(HM_ALKG1_v_lyrG1_DMG$geneID, HM_ALKG1_v_RS7G1_DMG_lyr$geneID)))
paste0(round(HM_ALKvpro1syn1_lyrDMG/length(HM_ALKG1_v_lyrG1_DMG$geneID)*100, 2), "%")

HM_ALKvsyn1syn4_halDEG <- sum(!is.na(match(HM_ALKG1_v_RS7G1_hal_DEG$geneID, HM_ALKG1_v_RS7G4_hal_DEG$geneID)))
paste0(round(HM_ALKvsyn1syn4_halDEG/length(HM_ALKG1_v_RS7G1_hal_DEG$geneID)*100, 2), "%")

HM_ALKvsyn1syn4_halDMG <- sum(!is.na(match(HM_ALKG1_v_RS7G1_DMG_hal$geneID, HM_ALKG1_v_RS7G4_DMG_hal$geneID)))
paste0(round(HM_ALKvsyn1syn4_halDMG/length(HM_ALKG1_v_RS7G1_DMG_hal$geneID)*100, 2), "%")

HM_ALKvsyn1syn4_lyrDEG <- sum(!is.na(match(HM_ALKG1_v_RS7G1_lyr_DEG$geneID, HM_ALKG1_v_RS7G4_lyr_DEG$geneID)))
paste0(round(HM_ALKvsyn1syn4_lyrDEG/length(HM_ALKG1_v_RS7G1_lyr_DEG$geneID)*100, 2), "%")

HM_ALKvsyn1syn4_lyrDMG <- sum(!is.na(match(HM_ALKG1_v_RS7G1_DMG_lyr$geneID, HM_ALKG1_v_RS7G4_DMG_lyr$geneID)))
paste0(round(HM_ALKvsyn1syn4_lyrDMG/length(HM_ALKG1_v_RS7G1_DMG_lyr$geneID)*100, 2), "%")


LL_ALKvpro1syn1_halDEG <- sum(!is.na(match(LL_ALKG1_v_halG1_DEG$geneID, LL_ALKG1_v_RS7G1_hal_DEG$geneID)))
paste0(round(LL_ALKvpro1syn1_halDEG/length(LL_ALKG1_v_halG1_DEG$geneID)*100, 2), "%")

LL_ALKvpro1syn1_halDMG <- sum(!is.na(match(LL_ALKG1_v_halG1_DMG$geneID, LL_ALKG1_v_RS7G1_DMG_hal$geneID)))
paste0(round(LL_ALKvpro1syn1_halDMG/length(LL_ALKG1_v_halG1_DMG$geneID)*100, 2), "%")

LL_ALKvpro1syn1_lyrDEG <- sum(!is.na(match(LL_ALKG1_v_lyrG1_DEG$geneID, LL_ALKG1_v_RS7G1_lyr_DEG$geneID)))
paste0(round(LL_ALKvpro1syn1_lyrDEG/length(LL_ALKG1_v_lyrG1_DEG$geneID)*100, 2), "%")

LL_ALKvpro1syn1_lyrDMG <- sum(!is.na(match(LL_ALKG1_v_lyrG1_DMG$geneID, LL_ALKG1_v_RS7G1_DMG_lyr$geneID)))
paste0(round(LL_ALKvpro1syn1_lyrDMG/length(LL_ALKG1_v_lyrG1_DMG$geneID)*100, 2), "%")

LL_ALKvsyn1syn4_halDEG <- sum(!is.na(match(LL_ALKG1_v_RS7G1_hal_DEG$geneID, LL_ALKG1_v_RS7G4_hal_DEG$geneID)))
paste0(round(LL_ALKvsyn1syn4_halDEG/length(LL_ALKG1_v_RS7G1_hal_DEG$geneID)*100, 2), "%")

LL_ALKvsyn1syn4_halDMG <- sum(!is.na(match(LL_ALKG1_v_RS7G1_DMG_hal$geneID, LL_ALKG1_v_RS7G4_DMG_hal$geneID)))
paste0(round(LL_ALKvsyn1syn4_halDMG/length(LL_ALKG1_v_RS7G1_DMG_hal$geneID)*100, 2), "%")

LL_ALKvsyn1syn4_lyrDEG <- sum(!is.na(match(LL_ALKG1_v_RS7G1_lyr_DEG$geneID, LL_ALKG1_v_RS7G4_lyr_DEG$geneID)))
paste0(round(LL_ALKvsyn1syn4_lyrDEG/length(LL_ALKG1_v_RS7G1_lyr_DEG$geneID)*100, 2), "%")

LL_ALKvsyn1syn4_lyrDMG <- sum(!is.na(match(LL_ALKG1_v_RS7G1_DMG_lyr$geneID, LL_ALKG1_v_RS7G4_DMG_lyr$geneID)))
paste0(round(LL_ALKvsyn1syn4_lyrDMG/length(LL_ALKG1_v_RS7G1_DMG_lyr$geneID)*100, 2), "%")


## TKS vs all overlap

HM_TKSvpro1syn1_halDEG <- sum(!is.na(match(HM_TKSG1_v_halG1_DEG$geneID, HM_TKSG1_v_RS7G1_hal_DEG$geneID)))
paste0(round(HM_TKSvpro1syn1_halDEG/length(HM_TKSG1_v_halG1_DEG$geneID)*100, 2), "%")

HM_TKSvpro1syn1_halDMG <- sum(!is.na(match(HM_TKSG1_v_halG1_DMG$geneID, HM_TKSG1_v_RS7G1_DMG_hal$geneID)))
paste0(round(HM_TKSvpro1syn1_halDMG/length(HM_TKSG1_v_halG1_DMG$geneID)*100, 2), "%")

HM_TKSvpro1syn1_lyrDEG <- sum(!is.na(match(HM_TKSG1_v_lyrG1_DEG$geneID, HM_TKSG1_v_RS7G1_lyr_DEG$geneID)))
paste0(round(HM_TKSvpro1syn1_lyrDEG/length(HM_TKSG1_v_lyrG1_DEG$geneID)*100, 2), "%")

HM_TKSvpro1syn1_lyrDMG <- sum(!is.na(match(HM_TKSG1_v_lyrG1_DMG$geneID, HM_TKSG1_v_RS7G1_DMG_lyr$geneID)))
paste0(round(HM_TKSvpro1syn1_lyrDMG/length(HM_TKSG1_v_lyrG1_DMG$geneID)*100, 2), "%")

HM_TKSvsyn1syn4_halDEG <- sum(!is.na(match(HM_TKSG1_v_RS7G1_hal_DEG$geneID, HM_TKSG1_v_RS7G4_hal_DEG$geneID)))
paste0(round(HM_TKSvsyn1syn4_halDEG/length(HM_TKSG1_v_RS7G1_hal_DEG$geneID)*100, 2), "%")

HM_TKSvsyn1syn4_halDMG <- sum(!is.na(match(HM_TKSG1_v_RS7G1_DMG_hal$geneID, HM_TKSG1_v_RS7G4_DMG_hal$geneID)))
paste0(round(HM_TKSvsyn1syn4_halDMG/length(HM_TKSG1_v_RS7G1_DMG_hal$geneID)*100, 2), "%")

HM_TKSvsyn1syn4_lyrDEG <- sum(!is.na(match(HM_TKSG1_v_RS7G1_lyr_DEG$geneID, HM_TKSG1_v_RS7G4_lyr_DEG$geneID)))
paste0(round(HM_TKSvsyn1syn4_lyrDEG/length(HM_TKSG1_v_RS7G1_lyr_DEG$geneID)*100, 2), "%")

HM_TKSvsyn1syn4_lyrDMG <- sum(!is.na(match(HM_TKSG1_v_RS7G1_DMG_lyr$geneID, HM_TKSG1_v_RS7G4_DMG_lyr$geneID)))
paste0(round(HM_TKSvsyn1syn4_lyrDMG/length(HM_TKSG1_v_RS7G1_DMG_lyr$geneID)*100, 2), "%")


LL_TKSvpro1syn1_halDEG <- sum(!is.na(match(LL_TKSG1_v_halG1_DEG$geneID, LL_TKSG1_v_RS7G1_hal_DEG$geneID)))
paste0(round(LL_TKSvpro1syn1_halDEG/length(LL_TKSG1_v_halG1_DEG$geneID)*100, 2), "%")

LL_TKSvpro1syn1_halDMG <- sum(!is.na(match(LL_TKSG1_v_halG1_DMG$geneID, LL_TKSG1_v_RS7G1_DMG_hal$geneID)))
paste0(round(LL_TKSvpro1syn1_halDMG/length(LL_TKSG1_v_halG1_DMG$geneID)*100, 2), "%")

LL_TKSvpro1syn1_lyrDEG <- sum(!is.na(match(LL_TKSG1_v_lyrG1_DEG$geneID, LL_TKSG1_v_RS7G1_lyr_DEG$geneID)))
paste0(round(LL_TKSvpro1syn1_lyrDEG/length(LL_TKSG1_v_lyrG1_DEG$geneID)*100, 2), "%")

LL_TKSvpro1syn1_lyrDMG <- sum(!is.na(match(LL_TKSG1_v_lyrG1_DMG$geneID, LL_TKSG1_v_RS7G1_DMG_lyr$geneID)))
paste0(round(LL_TKSvpro1syn1_lyrDMG/length(LL_TKSG1_v_lyrG1_DMG$geneID)*100, 2), "%")

LL_TKSvsyn1syn4_halDEG <- sum(!is.na(match(LL_TKSG1_v_RS7G1_hal_DEG$geneID, LL_TKSG1_v_RS7G4_hal_DEG$geneID)))
paste0(round(LL_TKSvsyn1syn4_halDEG/length(LL_TKSG1_v_RS7G1_hal_DEG$geneID)*100, 2), "%")

LL_TKSvsyn1syn4_halDMG <- sum(!is.na(match(LL_TKSG1_v_RS7G1_DMG_hal$geneID, LL_TKSG1_v_RS7G4_DMG_hal$geneID)))
paste0(round(LL_TKSvsyn1syn4_halDMG/length(LL_TKSG1_v_RS7G1_DMG_hal$geneID)*100, 2), "%")

LL_TKSvsyn1syn4_lyrDEG <- sum(!is.na(match(LL_TKSG1_v_RS7G1_lyr_DEG$geneID, LL_TKSG1_v_RS7G4_lyr_DEG$geneID)))
paste0(round(LL_TKSvsyn1syn4_lyrDEG/length(LL_TKSG1_v_RS7G1_lyr_DEG$geneID)*100, 2), "%")

LL_TKSvsyn1syn4_lyrDMG <- sum(!is.na(match(LL_TKSG1_v_RS7G1_DMG_lyr$geneID, LL_TKSG1_v_RS7G4_DMG_lyr$geneID)))
paste0(round(LL_TKSvsyn1syn4_lyrDMG/length(LL_TKSG1_v_RS7G1_DMG_lyr$geneID)*100, 2), "%")


# output overlaps

setwd("downstream_analyses/main_results/data/Fig8a")

write.table(data.frame(geneID=intersect(HM_halG1_v_RS7G1_DEG$geneID,
                      HM_halG1_v_RS7G1_DMG$geneID)),
            file = "HM_halG1_v_RS7G1_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(HM_halG1_v_RS7G1_DMG$geneID,
                    intersect(HM_halG1_v_RS7G1_DEG$geneID,
                              HM_halG1_v_RS7G1_DMG$geneID))),
            file = "HM_halG1_v_RS7G1_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=intersect(HM_halG1_v_RS7G4_DEG$geneID,
                      HM_halG1_v_RS7G4_DMG$geneID)),
            file = "HM_halG1_v_RS7G4_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(HM_halG1_v_RS7G4_DMG$geneID,
                    intersect(HM_halG1_v_RS7G4_DEG$geneID,
                              HM_halG1_v_RS7G4_DMG$geneID))),
            file = "HM_halG1_v_RS7G4_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)


write.table(data.frame(geneID=intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                      HM_lyrG1_v_RS7G1_DMG$geneID)),
            file = "HM_lyrG1_v_RS7G1_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(HM_lyrG1_v_RS7G1_DMG$geneID,
                    intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                              HM_lyrG1_v_RS7G1_DMG$geneID))),
            file = "HM_lyrG1_v_RS7G1_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                      HM_lyrG1_v_RS7G4_DMG$geneID)),
            file = "HM_lyrG1_v_RS7G4_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(HM_lyrG1_v_RS7G4_DMG$geneID,
                    intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                              HM_lyrG1_v_RS7G4_DMG$geneID))),
            file = "HM_lyrG1_v_RS7G4_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)


write.table(data.frame(geneID=intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                        LL_halG1_v_RS7G1_DMG$geneID)),
            file = "LL_halG1_v_RS7G1_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(LL_halG1_v_RS7G1_DMG$geneID,
                                      intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                LL_halG1_v_RS7G1_DMG$geneID))),
            file = "LL_halG1_v_RS7G1_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                        LL_halG1_v_RS7G4_DMG$geneID)),
            file = "LL_halG1_v_RS7G4_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(LL_halG1_v_RS7G4_DMG$geneID,
                                      intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                LL_halG1_v_RS7G4_DMG$geneID))),
            file = "LL_halG1_v_RS7G4_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)


write.table(data.frame(geneID=intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                        LL_lyrG1_v_RS7G1_DMG$geneID)),
            file = "LL_lyrG1_v_RS7G1_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(LL_lyrG1_v_RS7G1_DMG$geneID,
                                      intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                LL_lyrG1_v_RS7G1_DMG$geneID))),
            file = "LL_lyrG1_v_RS7G1_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                        LL_lyrG1_v_RS7G4_DMG$geneID)),
            file = "LL_lyrG1_v_RS7G4_DEG_DMG_overlap.txt",
            row.names = F, col.names = T, quote = F)

write.table(data.frame(geneID=setdiff(LL_lyrG1_v_RS7G4_DMG$geneID,
                                      intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                LL_lyrG1_v_RS7G4_DMG$geneID))),
            file = "LL_lyrG1_v_RS7G4_DEG_DMG_notOverlap.txt",
            row.names = F, col.names = T, quote = F)
