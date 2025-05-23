# This script checks for homologous pairs across DMGs

library(tidyverse)
library(data.table)

# Read reciprocal best-hits files

rbh_hal <- read.delim("~/path/to/rbh_Ahal_TAIR10.txt", header=FALSE)
rbh_lyr <- read.delim("~/path/to/rbh_Alyr_TAIR10.txt", header=FALSE)

# Read DMGs for PvS1 and PvS4 in both conditions

##MILD

MILD_CG_PvS1_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
MILD_CHG_PvS1_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
MILD_CHH_PvS1_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

MILD_CG_PvS1_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
MILD_CHG_PvS1_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
MILD_CHH_PvS1_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")
  
MILD_CG_PvS4_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
MILD_CHG_PvS4_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
MILD_CHH_PvS4_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

MILD_CG_PvS4_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
MILD_CHG_PvS4_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
MILD_CHH_PvS4_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

##STRESS

STRESS_CG_PvS1_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
STRESS_CHG_PvS1_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
STRESS_CHH_PvS1_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

STRESS_CG_PvS1_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn1_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
STRESS_CHG_PvS1_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn1_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
STRESS_CHH_PvS1_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn1_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

STRESS_CG_PvS4_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
STRESS_CHG_PvS4_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
STRESS_CHH_PvS4_hal <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

STRESS_CG_PvS4_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn4_v_pro1/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
STRESS_CHG_PvS4_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn4_v_pro1/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
STRESS_CHH_PvS4_lyr <- read.delim("~EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/STRESS_syn4_v_pro1/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

## Compute homologs per comparison

find_homologs_pairs <- function(hal_DMGs, lyr_DMGs) {
  hal_geneID <- select(hal_DMGs, "geneID")
  lyr_geneID <- select(lyr_DMGs, "geneID")
  homeolog_hal <- rbh_hal$V1[rbh_hal$V2 %in% hal_geneID$geneID]
  homeolog_lyr <- rbh_lyr$V1[rbh_lyr$V2 %in% lyr_geneID$geneID]
  return(intersect(homeolog_lyr, homeolog_hal))
  }

total_unique_genes <- function(hal_DMGs, lyr_DMGs) {
  hal_geneID <- select(hal_DMGs, "geneID")
  lyr_geneID <- select(lyr_DMGs, "geneID")
  homeolog_hal <- rbh_hal$V1[rbh_hal$V2 %in% hal_geneID$geneID]
  homeolog_lyr <- rbh_lyr$V1[rbh_lyr$V2 %in% lyr_geneID$geneID]
  return(union(homeolog_lyr, homeolog_hal))
}

MILD_CG_PvS1_homologs <- find_homologs_pairs(MILD_CG_PvS1_hal,
                                             MILD_CG_PvS1_lyr)
MILD_CG_PvS1_total <- total_unique_genes(MILD_CG_PvS1_hal,
                                         MILD_CG_PvS1_lyr)

MILD_CHG_PvS1_homologs <- find_homologs_pairs(MILD_CHG_PvS1_hal,
                                             MILD_CHG_PvS1_lyr)
MILD_CHG_PvS1_total <- total_unique_genes(MILD_CHG_PvS1_hal,
                                         MILD_CHG_PvS1_lyr)

MILD_CHH_PvS1_homologs <- find_homologs_pairs(MILD_CHH_PvS1_hal,
                                              MILD_CHH_PvS1_lyr)
MILD_CHH_PvS1_total <- total_unique_genes(MILD_CHH_PvS1_hal,
                                          MILD_CHH_PvS1_lyr)

MILD_CG_PvS4_homologs <- find_homologs_pairs(MILD_CG_PvS4_hal,
                                             MILD_CG_PvS4_lyr)
MILD_CG_PvS4_total <- total_unique_genes(MILD_CG_PvS4_hal,
                                         MILD_CG_PvS4_lyr)

MILD_CHG_PvS4_homologs <- find_homologs_pairs(MILD_CHG_PvS4_hal,
                                             MILD_CHG_PvS4_lyr)
MILD_CHG_PvS4_total <- total_unique_genes(MILD_CHG_PvS4_hal,
                                         MILD_CHG_PvS4_lyr)

MILD_CHH_PvS4_homologs <- find_homologs_pairs(MILD_CHH_PvS4_hal,
                                              MILD_CHH_PvS4_lyr)
MILD_CHH_PvS4_total <- total_unique_genes(MILD_CHH_PvS4_hal,
                                          MILD_CHH_PvS4_lyr)



STRESS_CG_PvS1_homologs <- find_homologs_pairs(STRESS_CG_PvS1_hal,
                                             STRESS_CG_PvS1_lyr)
STRESS_CG_PvS1_total <- total_unique_genes(STRESS_CG_PvS1_hal,
                                         STRESS_CG_PvS1_lyr)

STRESS_CHG_PvS1_homologs <- find_homologs_pairs(STRESS_CHG_PvS1_hal,
                                              STRESS_CHG_PvS1_lyr)
STRESS_CHG_PvS1_total <- total_unique_genes(STRESS_CHG_PvS1_hal,
                                          STRESS_CHG_PvS1_lyr)

STRESS_CHH_PvS1_homologs <- find_homologs_pairs(STRESS_CHH_PvS1_hal,
                                              STRESS_CHH_PvS1_lyr)
STRESS_CHH_PvS1_total <- total_unique_genes(STRESS_CHH_PvS1_hal,
                                          STRESS_CHH_PvS1_lyr)

STRESS_CG_PvS4_homologs <- find_homologs_pairs(STRESS_CG_PvS4_hal,
                                             STRESS_CG_PvS4_lyr)
STRESS_CG_PvS4_total <- total_unique_genes(STRESS_CG_PvS4_hal,
                                         STRESS_CG_PvS4_lyr)

STRESS_CHG_PvS4_homologs <- find_homologs_pairs(STRESS_CHG_PvS4_hal,
                                              STRESS_CHG_PvS4_lyr)
STRESS_CHG_PvS4_total <- total_unique_genes(STRESS_CHG_PvS4_hal,
                                          STRESS_CHG_PvS4_lyr)

STRESS_CHH_PvS4_homologs <- find_homologs_pairs(STRESS_CHH_PvS4_hal,
                                              STRESS_CHH_PvS4_lyr)
STRESS_CHH_PvS4_total <- total_unique_genes(STRESS_CHH_PvS4_hal,
                                          STRESS_CHH_PvS4_lyr)

fwrite(as.list(MILD_CG_PvS1_homologs), file = "output/path/MILD_CG_PvS1_homologs.txt", sep = "\n")
fwrite(as.list(MILD_CHG_PvS1_homologs), file = "output/path/MILD_CHG_PvS1_homologs.txt", sep = "\n")
fwrite(as.list(MILD_CHH_PvS1_homologs), file = "output/path/MILD_CHH_PvS1_homologs.txt", sep = "\n")
fwrite(as.list(MILD_CG_PvS4_homologs), file = "output/path/MILD_CG_PvS4_homologs.txt", sep = "\n")
fwrite(as.list(MILD_CHG_PvS4_homologs), file = "output/path/MILD_CHG_PvS4_homologs.txt", sep = "\n")
fwrite(as.list(MILD_CHH_PvS4_homologs), file = "output/path/MILD_CHH_PvS4_homologs.txt", sep = "\n")

fwrite(as.list(STRESS_CG_PvS1_homologs), file = "output/path/STRESS_CG_PvS1_homologs.txt", sep = "\n")
fwrite(as.list(STRESS_CHG_PvS1_homologs), file = "output/path/STRESS_CHG_PvS1_homologs.txt", sep = "\n")
fwrite(as.list(STRESS_CHH_PvS1_homologs), file = "output/path/STRESS_CHH_PvS1_homologs.txt", sep = "\n")
fwrite(as.list(STRESS_CG_PvS4_homologs), file = "output/path/STRESS_CG_PvS4_homologs.txt", sep = "\n")
fwrite(as.list(STRESS_CHG_PvS4_homologs), file = "output/path/STRESS_CHG_PvS4_homologs.txt", sep = "\n")
fwrite(as.list(STRESS_CHH_PvS4_homologs), file = "output/path/STRESS_CHH_PvS4_homologs.txt", sep = "\n")
