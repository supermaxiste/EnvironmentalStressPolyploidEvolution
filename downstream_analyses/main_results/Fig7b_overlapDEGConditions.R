### This script will compute and plot the overlap
### between differentially expressed genes. 
### Inputs are results from edgeR


### Import libraries

library(tidyverse)
library(patchwork)
library(VennDiagram)
library(ggVennDiagram)
library(RColorBrewer)
library(svglite)

# Folder with differentially expressed genes
DEG_folder <- c("path/to/deg_files_folder/")
# Folder with list of low coverage genes
lowC <- c("path/to/downstream_analyses/main_results/data/Fig7a/")

# import low coverage genes list

HM_lowC_genes_hal <- read.table(paste0(lowC, "HM_lowC_genes_hal.txt"), quote="\"", comment.char="")
HM_lowC_genes_lyr <- read.table(paste0(lowC, "HM_lowC_genes_lyr.txt"), quote="\"", comment.char="")

LL_lowC_genes_hal <- read.table(paste0(lowC, "LL_lowC_genes_hal.txt"), quote="\"", comment.char="")
LL_lowC_genes_lyr <- read.table(paste0(lowC, "LL_lowC_genes_lyr.txt"), quote="\"", comment.char="")

# import DEGs between conditions

HMvLL_synG1_hal_DEG

# across conditions

HMvLL_synG1_hal_DEG_flt <- filter(HMvLL_synG1_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HMvLL_synG1_lyr_DEG_flt <- filter(HMvLL_synG1_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HMvLL_synG1_hal_DEG_flt <- filter(HMvLL_synG1_hal_DEG_flt, !(geneID %in% LL_lowC_genes_hal$V1))
HMvLL_synG1_lyr_DEG_flt <- filter(HMvLL_synG1_lyr_DEG_flt, !(geneID %in% LL_lowC_genes_lyr$V1))

HMvLL_synG4_hal_DEG_flt <- filter(HMvLL_synG4_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HMvLL_synG4_lyr_DEG_flt <- filter(HMvLL_synG4_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HMvLL_synG4_hal_DEG_flt <- filter(HMvLL_synG4_hal_DEG_flt, !(geneID %in% LL_lowC_genes_hal$V1))
HMvLL_synG4_lyr_DEG_flt <- filter(HMvLL_synG4_lyr_DEG_flt, !(geneID %in% LL_lowC_genes_lyr$V1))

### Plot overlaps

setwd("~/OneDrive/PhD/Project/Chapter_3/Pictures/Paper_picsV5")

# Prepare a palette of 4 colors with R colorbrewer:
myCol <- brewer.pal(4, "Paired")

# Function to output venn diagram for two sets
venn <- function(set1, set2, title, filename){
  venn.diagram(
    x = list(set1, set2),
    category.names = c("G1" , "G4"),
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

## HMvLL syn G1 vs syn G4 - halleri side

venn_HMvLL_syn1Vsyn4_hal_DEG <- venn(set1 = HMvLL_synG1_hal_DEG_flt$geneID,
                                    set2 = HMvLL_synG4_hal_DEG_flt$geneID,
                                    title = "HMvLL_hal",
                                    filename = NULL)

ggsave(venn_HMvLL_syn1Vsyn4_hal_DEG, 
       filename = "HMvLL_syn1Vsyn4_hal_DEG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HMvLL syn G1 vs syn G4 - halleri side

venn_HMvLL_syn1Vsyn4_lyr_DEG <- venn(set1 = HMvLL_synG1_lyr_DEG_flt$geneID,
                                     set2 = HMvLL_synG4_lyr_DEG_flt$geneID,
                                     title = "HMvLL_lyr",
                                     filename = NULL)

ggsave(venn_HMvLL_syn1Vsyn4_lyr_DEG, 
       filename = "HMvLL_syn1Vsyn4_lyr_DEG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)
