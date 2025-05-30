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

### Import methylation data

setwd("ARPEGGIO_analyses")

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

# We filter all DMGs that are part of low coverage genes

# Import file with low coverage genes

HM_hal_lowC_scaffolds <- read.table("path/to/data/Fig5a/HM_lowC_genes_hal.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds <- read.table("path/to/data/Fig5a/HM_lowC_genes_lyr.txt", quote="\"", comment.char="")

LL_hal_lowC_scaffolds <- read.table("path/to/data/Fig5a/LL_lowC_genes_hal.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds <- read.table("path/to/data/Fig5a/LL_lowC_genes_lyr.txt", quote="\"", comment.char="")


# Filter out low coverage genes 

HM_halG1_v_RS7G1_DMG_CG <- filter(HM_halG1_v_RS7G1_DMG_CG,
                                  !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G1_DMG_CHG <- filter(HM_halG1_v_RS7G1_DMG_CHG,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G1_DMG_CHH <- filter(HM_halG1_v_RS7G1_DMG_CHH,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_halG1_v_RS7G4_DMG_CG <- filter(HM_halG1_v_RS7G4_DMG_CG,
                                  !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DMG_CHG <- filter(HM_halG1_v_RS7G4_DMG_CHG,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DMG_CHH <- filter(HM_halG1_v_RS7G4_DMG_CHH,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_RS7G1_v_RS7G4_DMG_CG_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CG_hal,
                                      !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHG_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CHG_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHH_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CHH_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_halG1_DMG_CG <- filter(HM_ALKG1_v_halG1_DMG_CG,
                                  !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_halG1_DMG_CHG <- filter(HM_ALKG1_v_halG1_DMG_CHG,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_halG1_DMG_CHH <- filter(HM_ALKG1_v_halG1_DMG_CHH,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G1_DMG_CG_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CG_hal,
                                      !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHG_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CHG_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHH_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CHH_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G4_DMG_CG_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CG_hal,
                                      !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHG_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CHG_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHH_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CHH_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_halG1_DMG_CG <- filter(HM_TKSG1_v_halG1_DMG_CG,
                                  !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_halG1_DMG_CHG <- filter(HM_TKSG1_v_halG1_DMG_CHG,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_halG1_DMG_CHH <- filter(HM_TKSG1_v_halG1_DMG_CHH,
                                   !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G1_DMG_CG_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CG_hal,
                                      !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHG_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CHG_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHH_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CHH_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G4_DMG_CG_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CG_hal,
                                      !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHG_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CHG_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHH_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CHH_hal,
                                       !(geneID %in% HM_hal_lowC_scaffolds$V1))


### lyrata side

HM_lyrG1_v_RS7G1_DMG_CG <- filter(HM_lyrG1_v_RS7G1_DMG_CG,
                                  !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G1_DMG_CHG <- filter(HM_lyrG1_v_RS7G1_DMG_CHG,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G1_DMG_CHH <- filter(HM_lyrG1_v_RS7G1_DMG_CHH,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_lyrG1_v_RS7G4_DMG_CG <- filter(HM_lyrG1_v_RS7G4_DMG_CG,
                                  !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DMG_CHG <- filter(HM_lyrG1_v_RS7G4_DMG_CHG,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DMG_CHH <- filter(HM_lyrG1_v_RS7G4_DMG_CHH,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_RS7G1_v_RS7G4_DMG_CG_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CG_lyr,
                                      !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHG_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CHG_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHH_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CHH_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_lyrG1_DMG_CG <- filter(HM_ALKG1_v_lyrG1_DMG_CG,
                                  !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_lyrG1_DMG_CHG <- filter(HM_ALKG1_v_lyrG1_DMG_CHG,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_lyrG1_DMG_CHH <- filter(HM_ALKG1_v_lyrG1_DMG_CHH,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G1_DMG_CG_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CG_lyr,
                                      !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHG_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CHG_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHH_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CHH_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G4_DMG_CG_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CG_lyr,
                                      !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHG_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CHG_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHH_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CHH_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_lyrG1_DMG_CG <- filter(HM_TKSG1_v_lyrG1_DMG_CG,
                                  !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_lyrG1_DMG_CHG <- filter(HM_TKSG1_v_lyrG1_DMG_CHG,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_lyrG1_DMG_CHH <- filter(HM_TKSG1_v_lyrG1_DMG_CHH,
                                   !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G1_DMG_CG_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CG_lyr,
                                      !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHG_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CHG_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHH_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CHH_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G4_DMG_CG_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CG_lyr,
                                      !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHG_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CHG_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHH_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CHH_lyr,
                                       !(geneID %in% HM_lyr_lowC_scaffolds$V1))


### Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DMG_CG <- filter(LL_halG1_v_RS7G1_DMG_CG,
                                  !(geneID %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G1_DMG_CHG <- filter(LL_halG1_v_RS7G1_DMG_CHG,
                                   !(geneID %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G1_DMG_CHH <- filter(LL_halG1_v_RS7G1_DMG_CHH,
                                   !(geneID %in% LL_hal_lowC_scaffolds$V1))

LL_halG1_v_RS7G4_DMG_CG <- filter(LL_halG1_v_RS7G4_DMG_CG,
                                  !(geneID %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DMG_CHG <- filter(LL_halG1_v_RS7G4_DMG_CHG,
                                   !(geneID %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DMG_CHH <- filter(LL_halG1_v_RS7G4_DMG_CHH,
                                   !(geneID %in% LL_hal_lowC_scaffolds$V1))

### lyrata side

LL_lyrG1_v_RS7G1_DMG_CG <- filter(LL_lyrG1_v_RS7G1_DMG_CG,
                                  !(geneID %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G1_DMG_CHG <- filter(LL_lyrG1_v_RS7G1_DMG_CHG,
                                   !(geneID %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G1_DMG_CHH <- filter(LL_lyrG1_v_RS7G1_DMG_CHH,
                                   !(geneID %in% LL_lyr_lowC_scaffolds$V1))

LL_lyrG1_v_RS7G4_DMG_CG <- filter(LL_lyrG1_v_RS7G4_DMG_CG,
                                  !(geneID %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DMG_CHG <- filter(LL_lyrG1_v_RS7G4_DMG_CHG,
                                   !(geneID %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DMG_CHH <- filter(LL_lyrG1_v_RS7G4_DMG_CHH,
                                   !(geneID %in% LL_lyr_lowC_scaffolds$V1))


### Import expression data

setwd("DEGs/")

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- read.delim("HM_halG1_v_RS7G1_DEG.txt")
HM_halG1_v_RS7G4_DEG <- read.delim("HM_halG1_v_RS7G4_DEG.txt")


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- read.delim("HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyrG1_v_RS7G4_DEG <- read.delim("HM_lyrG1_v_RS7G4_DEG.txt")


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- read.delim("LL_halG1_v_RS7G1_DEG.txt")
LL_halG1_v_RS7G4_DEG <- read.delim("LL_halG1_v_RS7G4_DEG.txt")

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- read.delim("LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyrG1_v_RS7G4_DEG <- read.delim("LL_lyrG1_v_RS7G4_DEG.txt")

# We filter DEG falling in low coverage regions as well

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- filter(HM_halG1_v_RS7G1_DEG,
                               !(geneID %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DEG <- filter(HM_halG1_v_RS7G4_DEG,
                               !(geneID %in% HM_hal_lowC_scaffolds$V1))

### lyrata side

HM_lyrG1_v_RS7G1_DEG <- filter(HM_lyrG1_v_RS7G1_DEG,
                               !(geneID %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DEG <- filter(HM_lyrG1_v_RS7G4_DEG,
                               !(geneID %in% HM_lyr_lowC_scaffolds$V1))

# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- filter(LL_halG1_v_RS7G1_DEG,
                               !(geneID %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DEG <- filter(LL_halG1_v_RS7G4_DEG,
                               !(geneID %in% LL_hal_lowC_scaffolds$V1))

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- filter(LL_lyrG1_v_RS7G1_DEG,
                               !(geneID %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DEG <- filter(LL_lyrG1_v_RS7G4_DEG,
                               !(geneID %in% LL_lyr_lowC_scaffolds$V1))


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


## lyrata side

HM_lyrG1_v_RS7G1_DMG <- consolidate(HM_lyrG1_v_RS7G1_DMG_CG,
                                    HM_lyrG1_v_RS7G1_DMG_CHG,
                                    HM_lyrG1_v_RS7G1_DMG_CHH)

HM_lyrG1_v_RS7G4_DMG <- consolidate(HM_lyrG1_v_RS7G4_DMG_CG,
                                    HM_lyrG1_v_RS7G4_DMG_CHG,
                                    HM_lyrG1_v_RS7G4_DMG_CHH)

### Hot conditions

## halleri side

LL_halG1_v_RS7G1_DMG <- consolidate(LL_halG1_v_RS7G1_DMG_CG,
                                    LL_halG1_v_RS7G1_DMG_CHG,
                                    LL_halG1_v_RS7G1_DMG_CHH)

LL_halG1_v_RS7G4_DMG <- consolidate(LL_halG1_v_RS7G4_DMG_CG,
                                    LL_halG1_v_RS7G4_DMG_CHG,
                                    LL_halG1_v_RS7G4_DMG_CHH)

## lyrata side

LL_lyrG1_v_RS7G1_DMG <- consolidate(LL_lyrG1_v_RS7G1_DMG_CG,
                                    LL_lyrG1_v_RS7G1_DMG_CHG,
                                    LL_lyrG1_v_RS7G1_DMG_CHH)

LL_lyrG1_v_RS7G4_DMG <- consolidate(LL_lyrG1_v_RS7G4_DMG_CG,
                                    LL_lyrG1_v_RS7G4_DMG_CHG,
                                    LL_lyrG1_v_RS7G4_DMG_CHH)


## Checking overlap across generations for DEGs and DMGs

## First look at diploid vs synthetic overlap

# HM_pro1Vsyn1_halDEG <- sum(!is.na(match(HM_halG1_v_RS7G1_DEG$geneID, HM_halG1_v_RS7G4_DEG$geneID)))
# HM_pro1Vsyn1_halDEG
# paste0(round(HM_pro1Vsyn1_halDEG/length(HM_halG1_v_RS7G1_DEG$geneID)*100, 2), "%")
# 
# HM_pro1Vsyn1_halDMG <- sum(!is.na(match(HM_halG1_v_RS7G1_DMG$geneID, HM_halG1_v_RS7G4_DMG$geneID)))
# HM_pro1Vsyn1_halDMG
# paste0(round(HM_pro1Vsyn1_halDMG/length(HM_halG1_v_RS7G1_DMG$geneID)*100, 2), "%")
# 
# HM_pro1Vsyn1_lyrDEG <- sum(!is.na(match(HM_lyrG1_v_RS7G1_DEG$geneID, HM_lyrG1_v_RS7G4_DEG$geneID)))
# HM_pro1Vsyn1_lyrDEG
# paste0(round(HM_pro1Vsyn1_lyrDEG/length(HM_lyrG1_v_RS7G1_DEG$geneID)*100, 2), "%")
# 
# HM_pro1Vsyn1_lyrDMG <- sum(!is.na(match(HM_lyrG1_v_RS7G1_DMG$geneID, HM_lyrG1_v_RS7G4_DMG$geneID)))
# HM_pro1Vsyn1_lyrDMG
# paste0(round(HM_pro1Vsyn1_lyrDMG/length(HM_lyrG1_v_RS7G1_DMG$geneID)*100, 2), "%")
# 
# LL_pro1Vsyn1_halDEG <- sum(!is.na(match(LL_halG1_v_RS7G1_DEG$geneID, LL_halG1_v_RS7G4_DEG$geneID)))
# LL_pro1Vsyn1_halDEG
# paste0(round(LL_pro1Vsyn1_halDEG/length(LL_halG1_v_RS7G1_DEG$geneID)*100, 2), "%")
# 
# LL_pro1Vsyn1_halDMG <- sum(!is.na(match(LL_halG1_v_RS7G1_DMG$geneID, LL_halG1_v_RS7G4_DMG$geneID)))
# LL_pro1Vsyn1_halDMG
# paste0(round(LL_pro1Vsyn1_halDMG/length(LL_halG1_v_RS7G1_DMG$geneID)*100, 2), "%")
# 
# LL_pro1Vsyn1_lyrDEG <- sum(!is.na(match(LL_lyrG1_v_RS7G1_DEG$geneID, LL_lyrG1_v_RS7G4_DEG$geneID)))
# LL_pro1Vsyn1_lyrDEG
# paste0(round(LL_pro1Vsyn1_lyrDEG/length(LL_lyrG1_v_RS7G1_DEG$geneID)*100, 2), "%")
# 
# LL_pro1Vsyn1_lyrDMG <- sum(!is.na(match(LL_lyrG1_v_RS7G1_DMG$geneID, LL_lyrG1_v_RS7G4_DMG$geneID)))
# LL_pro1Vsyn1_lyrDMG
# paste0(round(LL_pro1Vsyn1_lyrDMG/length(LL_lyrG1_v_RS7G1_DMG$geneID)*100, 2), "%")


### We can also do a simple contingency
### table and test if the DMGs and DEGs are independent (H0)
### or dependent (H1)

# create contingency table for HM hal side

## G1
HM_table_halG1_v_RS7G1 <- matrix(c(length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                  HM_halG1_v_RS7G1_DMG$geneID)),
                                 nrow(HM_halG1_v_RS7G1_DEG) - length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                                                 HM_halG1_v_RS7G1_DMG$geneID)),
                                 nrow(HM_halG1_v_RS7G1_DMG) - length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                                               HM_halG1_v_RS7G1_DMG$geneID)),
                                 21518 - nrow(HM_halG1_v_RS7G1_DEG) - nrow(HM_halG1_v_RS7G1_DMG) + length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                                                                                    HM_halG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_halG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(HM_table_halG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(HM_table_halG1_v_RS7G1, correct = F)

## G4
HM_table_halG1_v_RS7G4 <- matrix(c(length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                    HM_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_halG1_v_RS7G4_DEG) - length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                                                 HM_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_halG1_v_RS7G4_DMG) - length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                                                 HM_halG1_v_RS7G4_DMG$geneID)),
                                   21518 - nrow(HM_halG1_v_RS7G4_DEG) - nrow(HM_halG1_v_RS7G4_DMG) + length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                                                                                      HM_halG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_halG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(HM_table_halG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(HM_table_halG1_v_RS7G4, correct = F)

# create contingency table for HM lyr side

## G1
HM_table_lyrG1_v_RS7G1 <- matrix(c(length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                    HM_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G1_DEG) - length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G1_DMG) - length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G1_DMG$geneID)),
                                   19595 - nrow(HM_lyrG1_v_RS7G1_DEG) - nrow(HM_lyrG1_v_RS7G1_DMG) + length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                                                                                      HM_lyrG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_lyrG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(HM_table_lyrG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(HM_table_lyrG1_v_RS7G1)

## G4
HM_table_lyrG1_v_RS7G4 <- matrix(c(length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                    HM_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G4_DEG) - length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G4_DMG) - length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G4_DMG$geneID)),
                                   19595 - nrow(HM_lyrG1_v_RS7G4_DEG) - nrow(HM_lyrG1_v_RS7G4_DMG) + length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                                                                                      HM_lyrG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_lyrG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(HM_table_lyrG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(HM_table_lyrG1_v_RS7G4)

# create contingency table for LL hal side

## G1
LL_table_halG1_v_RS7G1 <- matrix(c(length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                    LL_halG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G1_DEG) - length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                                                 LL_halG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G1_DMG) - length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                                                 LL_halG1_v_RS7G1_DMG$geneID)),
                                   21518 - nrow(LL_halG1_v_RS7G1_DEG) - nrow(LL_halG1_v_RS7G1_DMG) + length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                                                                                      LL_halG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_halG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(LL_table_halG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(LL_table_halG1_v_RS7G1, correct = F)

## G4
LL_table_halG1_v_RS7G4 <- matrix(c(length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                    LL_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G4_DEG) - length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                                                 LL_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G4_DMG) - length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                                                 LL_halG1_v_RS7G4_DMG$geneID)),
                                   21518 - nrow(LL_halG1_v_RS7G4_DEG) - nrow(LL_halG1_v_RS7G4_DMG) + length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                                                                                      LL_halG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_halG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(LL_table_halG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(LL_table_halG1_v_RS7G4, correct = F)

# create contingency table for LL lyr side

## G1
LL_table_lyrG1_v_RS7G1 <- matrix(c(length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                    LL_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G1_DEG) - length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G1_DMG) - length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G1_DMG$geneID)),
                                   19595 - nrow(LL_lyrG1_v_RS7G1_DEG) - nrow(LL_lyrG1_v_RS7G1_DMG) + length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                                                                                      LL_lyrG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_lyrG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(LL_table_lyrG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(LL_table_lyrG1_v_RS7G1)

## G4
LL_table_lyrG1_v_RS7G4 <- matrix(c(length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                    LL_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G4_DEG) - length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G4_DMG) - length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G4_DMG$geneID)),
                                   19595 - nrow(LL_lyrG1_v_RS7G4_DEG) - nrow(LL_lyrG1_v_RS7G4_DMG) + length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                                                                                      LL_lyrG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_lyrG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(LL_table_lyrG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(LL_table_lyrG1_v_RS7G4)
