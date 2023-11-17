### We plot heatmaps for expression patterns in all
### of our samples (tpm normalized expression)

## Import libraries

library(tidyverse)
library(data.table)
library(DESeq2)
library(ComplexHeatmap)

## Paths to define

# Paths to genome annotation files
Ahal_anno_path <- c("path/to/anno_hal.gff")
Alyr_anno_path <- c("path/to/anno_lyr.gff")
# Path to file with low coverage genes
low_coverage_genes <- c("path/to/downstream_analyses/main_results/data/Fig5a/")
# Path to count tables (featureCounts output)
count_tables <- c("path/to/count_tables_folder/")

## In the first part of the script we import
## the annotations to exclude genes falling
## into low coverage scaffolds

# We import annotations to find the corresponding scaffold of DEGs

Ahal_v2_2 <- read.delim(Ahal_anno_path,
                        header=FALSE,
                        col.names = c("chr",
                                      "tool", 
                                      "context",
                                      "start", 
                                      "end", 
                                      "number",
                                      "strand", 
                                      "dot", 
                                      "extra"))
Alyr_v2_2 <- read.delim(Alyr_anno_path, 
                        header=FALSE,
                        col.names = c("chr",
                                      "tool", 
                                      "context",
                                      "start", 
                                      "end", 
                                      "number",
                                      "strand", 
                                      "dot", 
                                      "extra"))

Ahal_anno <- select(Ahal_v2_2, 
                    chr, 
                    context, 
                    start, 
                    end, 
                    strand, 
                    extra)
rm(Ahal_v2_2)


Alyr_anno <- select(Alyr_v2_2, 
                    chr, 
                    context, 
                    start, 
                    end, 
                    strand, 
                    extra)
rm(Alyr_v2_2)

# We select only genes

Ahal_anno <- filter(Ahal_anno, context == "gene")
Alyr_anno <- filter(Alyr_anno, context == "gene")

# Import low coverage scaffolds

HM_hal_lowC_genes  <- read.table(paste0(low_coverage_genes, "HM_lowC_genes_hal.txt"), quote="\"", comment.char="")
HM_lyr_lowC_genes  <- read.table(paste0(low_coverage_genes, "HM_lowC_genes_lyr.txt"), quote="\"", comment.char="")

LL_hal_lowC_genes  <- read.table(paste0(low_coverage_genes, "LL_lowC_genes_hal.txt"), quote="\"", comment.char="")
LL_lyr_lowC_genes  <- read.table(paste0(low_coverage_genes, "LL_lowC_genes_lyr.txt"), quote="\"", comment.char="")

# Extract geneID from extra column

Ahal_anno_geneID <- as.data.frame(str_split_fixed(Ahal_anno$extra, "=", 2)[,2])
Ahal_anno_geneID <- as.data.frame(str_split_fixed(Ahal_anno_geneID[,1], ";", 2)[,1])

colnames(Ahal_anno_geneID) <- "geneID"

Alyr_anno_geneID <- as.data.frame(str_split_fixed(Alyr_anno$extra, "=", 2)[,2])
Alyr_anno_geneID <- as.data.frame(str_split_fixed(Alyr_anno_geneID[,1], ";", 2)[,1])

colnames(Alyr_anno_geneID) <- "geneID"

# Add geneID to annotation dataframe

Ahal_anno <- Ahal_anno %>% mutate(geneID = Ahal_anno_geneID$geneID)
Alyr_anno <- Alyr_anno %>% mutate(geneID = Alyr_anno_geneID$geneID)

# Filter out genes falling in low coverage regions

filter_hal_HM <- which(Ahal_anno$geneID %in%
                         HM_hal_lowC_genes$V1)

filter_hal_LL <- which(Ahal_anno$geneID %in%
                         LL_hal_lowC_genes$V1)

filter_lyr_HM <- which(Alyr_anno$geneID %in%
                         HM_lyr_lowC_genes$V1)

filter_lyr_LL <- which(Alyr_anno$geneID %in%
                         LL_lyr_lowC_genes$V1)


# Processing count tables 

colnames <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Counts")


# Read files (raw counts)

HM_hal_G1_1 <- fread(paste0(count_tables, "HM_hal_G1_1_counts.txt"), 
                     col.names = colnames)
HM_hal_G1_2 <- fread(paste0(count_tables, "HM_hal_G1_2_counts.txt"), 
                     col.names = colnames)
HM_hal_G1_3 <- fread(paste0(count_tables, "HM_hal_G1_3_counts.txt"), 
                     col.names = colnames)
HM_hal_G1_4 <- fread(paste0(count_tables, "HM_hal_G1_4_counts.txt"), 
                     col.names = colnames)
HM_hal_G4_1 <- fread(paste0(count_tables, "HM_hal_G4_1_counts.txt"), 
                     col.names = colnames)
HM_hal_G4_2 <- fread(paste0(count_tables, "HM_hal_G4_2_counts.txt"), 
                     col.names = colnames)
HM_hal_G4_3 <- fread(paste0(count_tables, "HM_hal_G4_3_counts.txt"), 
                     col.names = colnames)
HM_lyr_G1_1 <- fread(paste0(count_tables, "HM_lyr_G1_1_counts.txt"), 
                     col.names = colnames)
HM_lyr_G1_2 <- fread(paste0(count_tables, "HM_lyr_G1_2_counts.txt"), 
                     col.names = colnames)
HM_lyr_G1_3 <- fread(paste0(count_tables, "HM_lyr_G1_3_counts.txt"), 
                     col.names = colnames)
HM_lyr_G1_4 <- fread(paste0(count_tables, "HM_lyr_G1_4_counts.txt"), 
                     col.names = colnames)
HM_lyr_G4_1 <- fread(paste0(count_tables, "HM_lyr_G4_1_counts.txt"), 
                     col.names = colnames)
HM_lyr_G4_2 <- fread(paste0(count_tables, "HM_lyr_G4_2_counts.txt"), 
                     col.names = colnames)
HM_lyr_G4_3 <- fread(paste0(count_tables, "HM_lyr_G4_3_counts.txt"), 
                     col.names = colnames)

HM_RS7_G1_1_hal <- fread(paste0(count_tables, "HM_RS7_G1_1_halcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_1_lyr <- fread(paste0(count_tables, "HM_RS7_G1_1_lyrcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_2_hal <- fread(paste0(count_tables, "HM_RS7_G1_2_halcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_2_lyr <- fread(paste0(count_tables, "HM_RS7_G1_2_lyrcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_3_hal <- fread(paste0(count_tables, "HM_RS7_G1_3_halcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_3_lyr <- fread(paste0(count_tables, "HM_RS7_G1_3_lyrcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_4_hal <- fread(paste0(count_tables, "HM_RS7_G1_4_halcounts.txt"), 
                         col.names = colnames)
HM_RS7_G1_4_lyr <- fread(paste0(count_tables, "HM_RS7_G1_4_lyrcounts.txt"), 
                         col.names = colnames)

HM_RS7_G4_1_hal <- fread(paste0(count_tables, "HM_RS7_G4_1_halcounts.txt"), 
                         col.names = colnames)
HM_RS7_G4_1_lyr <- fread(paste0(count_tables, "HM_RS7_G4_1_lyrcounts.txt"), 
                         col.names = colnames)
HM_RS7_G4_2_hal <- fread(paste0(count_tables, "HM_RS7_G4_2_halcounts.txt"), 
                         col.names = colnames)
HM_RS7_G4_2_lyr <- fread(paste0(count_tables, "HM_RS7_G4_2_lyrcounts.txt"), 
                         col.names = colnames)

HM_ALK_G1_1_hal <- fread(paste0(count_tables, "HM_ALK_G1_1_halcounts.txt"), 
                         col.names = colnames)
HM_ALK_G1_1_lyr <- fread(paste0(count_tables, "HM_ALK_G1_1_lyrcounts.txt"), 
                         col.names = colnames)
HM_ALK_G1_2_hal <- fread(paste0(count_tables, "HM_ALK_G1_2_halcounts.txt"), 
                         col.names = colnames)
HM_ALK_G1_2_lyr <- fread(paste0(count_tables, "HM_ALK_G1_2_lyrcounts.txt"), 
                         col.names = colnames)
HM_ALK_G1_3_hal <- fread(paste0(count_tables, "HM_ALK_G1_3_halcounts.txt"), 
                         col.names = colnames)
HM_ALK_G1_3_lyr <- fread(paste0(count_tables, "HM_ALK_G1_3_lyrcounts.txt"), 
                         col.names = colnames)

HM_ALK_G4_1_hal <- fread(paste0(count_tables, "HM_ALK_G4_1_halcounts.txt"), 
                         col.names = colnames)
HM_ALK_G4_1_lyr <- fread(paste0(count_tables, "HM_ALK_G4_1_lyrcounts.txt"), 
                         col.names = colnames)
HM_ALK_G4_2_hal <- fread(paste0(count_tables, "HM_ALK_G4_2_halcounts.txt"), 
                         col.names = colnames)
HM_ALK_G4_2_lyr <- fread(paste0(count_tables, "HM_ALK_G4_2_lyrcounts.txt"), 
                         col.names = colnames)
HM_ALK_G4_3_hal <- fread(paste0(count_tables, "HM_ALK_G4_3_halcounts.txt"), 
                         col.names = colnames)
HM_ALK_G4_3_lyr <- fread(paste0(count_tables, "HM_ALK_G4_3_lyrcounts.txt"), 
                         col.names = colnames)

HM_TKS_G1_1_hal <- fread(paste0(count_tables, "HM_TKS_G1_1_halcounts.txt"), 
                         col.names = colnames)
HM_TKS_G1_1_lyr <- fread(paste0(count_tables, "HM_TKS_G1_1_lyrcounts.txt"), 
                         col.names = colnames)
HM_TKS_G1_2_hal <- fread(paste0(count_tables, "HM_TKS_G1_2_halcounts.txt"), 
                         col.names = colnames)
HM_TKS_G1_2_lyr <- fread(paste0(count_tables, "HM_TKS_G1_2_lyrcounts.txt"), 
                         col.names = colnames)
HM_TKS_G1_3_hal <- fread(paste0(count_tables, "HM_TKS_G1_3_halcounts.txt"), 
                         col.names = colnames)
HM_TKS_G1_3_lyr <- fread(paste0(count_tables, "HM_TKS_G1_3_lyrcounts.txt"), 
                         col.names = colnames)

HM_TKS_G5_1_hal <- fread(paste0(count_tables, "HM_TKS_G5_1_halcounts.txt"), 
                         col.names = colnames)
HM_TKS_G5_1_lyr <- fread(paste0(count_tables, "HM_TKS_G5_1_lyrcounts.txt"), 
                         col.names = colnames)
HM_TKS_G5_2_hal <- fread(paste0(count_tables, "HM_TKS_G5_2_halcounts.txt"), 
                         col.names = colnames)
HM_TKS_G5_2_lyr <- fread(paste0(count_tables, "HM_TKS_G5_2_lyrcounts.txt"), 
                         col.names = colnames)
HM_TKS_G5_3_hal <- fread(paste0(count_tables, "HM_TKS_G5_3_halcounts.txt"), 
                         col.names = colnames)
HM_TKS_G5_3_lyr <- fread(paste0(count_tables, "HM_TKS_G5_3_lyrcounts.txt"), 
                         col.names = colnames)


LL_hal_G1_1 <- fread(paste0(count_tables, "LL_hal_G1_1_counts.txt"), 
                     col.names = colnames)
LL_hal_G1_2 <- fread(paste0(count_tables, "LL_hal_G1_2_counts.txt"), 
                     col.names = colnames)
LL_hal_G1_4 <- fread(paste0(count_tables, "LL_hal_G1_4_counts.txt"), 
                     col.names = colnames)
LL_hal_G1_5 <- fread(paste0(count_tables, "LL_hal_G1_5_counts.txt"), 
                     col.names = colnames)
LL_hal_G1_6 <- fread(paste0(count_tables, "LL_hal_G1_6_counts.txt"), 
                     col.names = colnames)
LL_hal_G4_1 <- fread(paste0(count_tables, "LL_hal_G4_1_counts.txt"), 
                     col.names = colnames)
LL_hal_G4_2 <- fread(paste0(count_tables, "LL_hal_G4_2_counts.txt"), 
                     col.names = colnames)

LL_lyr_G1_1 <- fread(paste0(count_tables, "LL_lyr_G1_1_counts.txt"), 
                     col.names = colnames)
LL_lyr_G1_2 <- fread(paste0(count_tables, "LL_lyr_G1_2_counts.txt"), 
                     col.names = colnames)
LL_lyr_G1_4 <- fread(paste0(count_tables, "LL_lyr_G1_4_counts.txt"), 
                     col.names = colnames)
LL_lyr_G1_5 <- fread(paste0(count_tables, "LL_lyr_G1_5_counts.txt"), 
                     col.names = colnames)
LL_lyr_G1_6 <- fread(paste0(count_tables, "LL_lyr_G1_6_counts.txt"), 
                     col.names = colnames)
LL_lyr_G4_1 <- fread(paste0(count_tables, "LL_lyr_G4_1_counts.txt"), 
                     col.names = colnames)
LL_lyr_G4_2 <- fread(paste0(count_tables, "LL_lyr_G4_2_counts.txt"), 
                     col.names = colnames)
LL_lyr_G4_3 <- fread(paste0(count_tables, "LL_lyr_G4_3_counts.txt"), 
                     col.names = colnames)
LL_lyr_G4_4 <- fread(paste0(count_tables, "LL_lyr_G4_4_counts.txt"), 
                     col.names = colnames)
LL_lyr_G4_5 <- fread(paste0(count_tables, "LL_lyr_G4_5_counts.txt"), 
                     col.names = colnames)

LL_RS7_G1_1_hal <- fread(paste0(count_tables, "LL_RS7_G1_1_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_1_lyr <- fread(paste0(count_tables, "LL_RS7_G1_1_lyrcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_2_hal <- fread(paste0(count_tables, "LL_RS7_G1_2_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_2_lyr <- fread(paste0(count_tables, "LL_RS7_G1_2_lyrcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_3_hal <- fread(paste0(count_tables, "LL_RS7_G1_3_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_3_lyr <- fread(paste0(count_tables, "LL_RS7_G1_3_lyrcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_4_hal <- fread(paste0(count_tables, "LL_RS7_G1_4_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G1_4_lyr <- fread(paste0(count_tables, "LL_RS7_G1_4_lyrcounts.txt"), 
                         col.names = colnames)

LL_RS7_G4_1_hal <- fread(paste0(count_tables, "LL_RS7_G4_1_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G4_1_lyr <- fread(paste0(count_tables, "LL_RS7_G4_1_lyrcounts.txt"), 
                         col.names = colnames)
LL_RS7_G4_2_hal <- fread(paste0(count_tables, "LL_RS7_G4_2_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G4_2_lyr <- fread(paste0(count_tables, "LL_RS7_G4_2_lyrcounts.txt"), 
                         col.names = colnames)
LL_RS7_G4_3_hal <- fread(paste0(count_tables, "LL_RS7_G4_3_halcounts.txt"), 
                         col.names = colnames)
LL_RS7_G4_3_lyr <- fread(paste0(count_tables, "LL_RS7_G4_3_lyrcounts.txt"), 
                         col.names = colnames)

LL_ALK_G1_1_hal <- fread(paste0(count_tables, "LL_ALK_G1_1_halcounts.txt"), 
                         col.names = colnames)
LL_ALK_G1_1_lyr <- fread(paste0(count_tables, "LL_ALK_G1_1_lyrcounts.txt"), 
                         col.names = colnames)
LL_ALK_G1_2_hal <- fread(paste0(count_tables, "LL_ALK_G1_2_halcounts.txt"), 
                         col.names = colnames)
LL_ALK_G1_2_lyr <- fread(paste0(count_tables, "LL_ALK_G1_2_lyrcounts.txt"), 
                         col.names = colnames)
LL_ALK_G1_3_hal <- fread(paste0(count_tables, "LL_ALK_G1_3_halcounts.txt"), 
                         col.names = colnames)
LL_ALK_G1_3_lyr <- fread(paste0(count_tables, "LL_ALK_G1_3_lyrcounts.txt"), 
                         col.names = colnames)

LL_ALK_G4_1_hal <- fread(paste0(count_tables, "LL_ALK_G4_1_halcounts.txt"), 
                         col.names = colnames)
LL_ALK_G4_1_lyr <- fread(paste0(count_tables, "LL_ALK_G4_1_lyrcounts.txt"), 
                         col.names = colnames)
LL_ALK_G4_2_hal <- fread(paste0(count_tables, "LL_ALK_G4_2_halcounts.txt"), 
                         col.names = colnames)
LL_ALK_G4_2_lyr <- fread(paste0(count_tables, "LL_ALK_G4_2_lyrcounts.txt"), 
                         col.names = colnames)
LL_ALK_G4_3_hal <- fread(paste0(count_tables, "LL_ALK_G4_3_halcounts.txt"), 
                         col.names = colnames)
LL_ALK_G4_3_lyr <- fread(paste0(count_tables, "LL_ALK_G4_3_lyrcounts.txt"), 
                         col.names = colnames)

LL_TKS_G1_1_hal <- fread(paste0(count_tables, "LL_TKS_G1_1_halcounts.txt"), 
                         col.names = colnames)
LL_TKS_G1_1_lyr <- fread(paste0(count_tables, "LL_TKS_G1_1_lyrcounts.txt"), 
                         col.names = colnames)
LL_TKS_G1_2_hal <- fread(paste0(count_tables, "LL_TKS_G1_2_halcounts.txt"), 
                         col.names = colnames)
LL_TKS_G1_2_lyr <- fread(paste0(count_tables, "LL_TKS_G1_2_lyrcounts.txt"), 
                         col.names = colnames)
LL_TKS_G1_3_hal <- fread(paste0(count_tables, "LL_TKS_G1_3_halcounts.txt"), 
                         col.names = colnames)
LL_TKS_G1_3_lyr <- fread(paste0(count_tables, "LL_TKS_G1_3_lyrcounts.txt"), 
                         col.names = colnames)

LL_TKS_G5_1_hal <- fread(paste0(count_tables, "LL_TKS_G5_1_halcounts.txt"), 
                         col.names = colnames)
LL_TKS_G5_1_lyr <- fread(paste0(count_tables, "LL_TKS_G5_1_lyrcounts.txt"), 
                         col.names = colnames)
LL_TKS_G5_2_hal <- fread(paste0(count_tables, "LL_TKS_G5_2_halcounts.txt"), 
                         col.names = colnames)
LL_TKS_G5_2_lyr <- fread(paste0(count_tables, "LL_TKS_G5_2_lyrcounts.txt"), 
                         col.names = colnames)
LL_TKS_G5_3_hal <- fread(paste0(count_tables, "LL_TKS_G5_3_halcounts.txt"), 
                         col.names = colnames)
LL_TKS_G5_3_lyr <- fread(paste0(count_tables, "LL_TKS_G5_3_lyrcounts.txt"), 
                         col.names = colnames)

## We create a data frame with all raw counts,
## filter them from low coverage genes
## and apply a variance stabilizing transformation
## to prevent genes with high counts driving all
## the variance

all_samples_HM_hal <- cbind(HM_hal_G1_1$Counts,
                            HM_hal_G1_2$Counts,
                            HM_hal_G1_3$Counts,
                            HM_hal_G1_4$Counts,
                            HM_hal_G4_1$Counts,
                            HM_hal_G4_2$Counts,
                            HM_hal_G4_3$Counts,
                            HM_RS7_G1_1_hal$Counts,
                            HM_RS7_G1_2_hal$Counts,
                            HM_RS7_G1_3_hal$Counts,
                            HM_RS7_G1_4_hal$Counts,
                            HM_RS7_G4_1_hal$Counts,
                            HM_RS7_G4_2_hal$Counts,
                            HM_ALK_G1_1_hal$Counts,
                            HM_ALK_G1_2_hal$Counts,
                            HM_ALK_G1_3_hal$Counts,
                            HM_ALK_G4_1_hal$Counts,
                            HM_ALK_G4_2_hal$Counts,
                            HM_ALK_G4_3_hal$Counts,
                            HM_TKS_G1_1_hal$Counts,
                            HM_TKS_G1_2_hal$Counts,
                            HM_TKS_G1_3_hal$Counts,
                            HM_TKS_G5_1_hal$Counts,
                            HM_TKS_G5_2_hal$Counts,
                            HM_TKS_G5_3_hal$Counts)

all_samples_HM_hal <- dplyr::slice(as.data.frame(all_samples_HM_hal),
                             -filter_hal_HM)

all_samples_HM_hal <- as.data.frame(vst(as.matrix(all_samples_HM_hal)))

names(all_samples_HM_hal) <- c("HM_hal_G1_1",
                               "HM_hal_G1_2",
                               "HM_hal_G1_3",
                               "HM_hal_G1_4",
                               "HM_hal_G4_1",
                               "HM_hal_G4_2",
                               "HM_hal_G4_3",
                               "HM_RS7_G1_1",
                               "HM_RS7_G1_2",
                               "HM_RS7_G1_3",
                               "HM_RS7_G1_4",
                               "HM_RS7_G4_1",
                               "HM_RS7_G4_2",
                               "HM_ALK_G1_1",
                               "HM_ALK_G1_2",
                               "HM_ALK_G1_3",
                               "HM_ALK_G4_1",
                               "HM_ALK_G4_2",
                               "HM_ALK_G4_3",
                               "HM_TKS_G1_1",
                               "HM_TKS_G1_2",
                               "HM_TKS_G1_3",
                               "HM_TKS_G5_1",
                               "HM_TKS_G5_2",
                               "HM_TKS_G5_3")

all_samples_LL_hal <- cbind(LL_hal_G1_1$Counts,
                            LL_hal_G1_2$Counts,
                            LL_hal_G1_4$Counts,
                            LL_hal_G1_5$Counts,
                            LL_hal_G1_6$Counts,
                            LL_hal_G4_1$Counts,
                            LL_hal_G4_2$Counts,
                            LL_RS7_G1_1_hal$Counts,
                            LL_RS7_G1_2_hal$Counts,
                            LL_RS7_G1_3_hal$Counts,
                            LL_RS7_G1_4_hal$Counts,
                            LL_RS7_G4_1_hal$Counts,
                            LL_RS7_G4_2_hal$Counts,
                            LL_RS7_G4_3_hal$Counts,
                            LL_ALK_G1_1_hal$Counts,
                            LL_ALK_G1_2_hal$Counts,
                            LL_ALK_G1_3_hal$Counts,
                            LL_ALK_G4_1_hal$Counts,
                            LL_ALK_G4_2_hal$Counts,
                            LL_ALK_G4_3_hal$Counts,
                            LL_TKS_G1_1_hal$Counts,
                            LL_TKS_G1_2_hal$Counts,
                            LL_TKS_G1_3_hal$Counts,
                            LL_TKS_G5_1_hal$Counts,
                            LL_TKS_G5_2_hal$Counts,
                            LL_TKS_G5_3_hal$Counts)

all_samples_LL_hal <- dplyr::slice(as.data.frame(all_samples_LL_hal),
                                   -filter_hal_LL)

all_samples_LL_hal <- as.data.frame(vst(as.matrix(all_samples_LL_hal)))

names(all_samples_LL_hal) <- c("LL_hal_G1_1",
                               "LL_hal_G1_2",
                               "LL_hal_G1_4",
                               "LL_hal_G1_5",
                               "LL_hal_G1_6",
                               "LL_hal_G4_1",
                               "LL_hal_G4_2",
                               "LL_RS7_G1_1",
                               "LL_RS7_G1_2",
                               "LL_RS7_G1_3",
                               "LL_RS7_G1_4",
                               "LL_RS7_G4_1",
                               "LL_RS7_G4_2",
                               "LL_RS7_G4_3",
                               "LL_ALK_G1_1",
                               "LL_ALK_G1_2",
                               "LL_ALK_G1_3",
                               "LL_ALK_G4_1",
                               "LL_ALK_G4_2",
                               "LL_ALK_G4_3",
                               "LL_TKS_G1_1",
                               "LL_TKS_G1_2",
                               "LL_TKS_G1_3",
                               "LL_TKS_G5_1",
                               "LL_TKS_G5_2",
                               "LL_TKS_G5_3")

all_samples_HM_lyr <- cbind(HM_lyr_G1_1$Counts,
                            HM_lyr_G1_2$Counts,
                            HM_lyr_G1_3$Counts,
                            HM_lyr_G1_4$Counts,
                            HM_lyr_G4_1$Counts,
                            HM_lyr_G4_2$Counts,
                            HM_lyr_G4_3$Counts,
                            HM_RS7_G1_1_lyr$Counts,
                            HM_RS7_G1_2_lyr$Counts,
                            HM_RS7_G1_3_lyr$Counts,
                            HM_RS7_G1_4_lyr$Counts,
                            HM_RS7_G4_1_lyr$Counts,
                            HM_RS7_G4_2_lyr$Counts,
                            HM_ALK_G1_1_lyr$Counts,
                            HM_ALK_G1_2_lyr$Counts,
                            HM_ALK_G1_3_lyr$Counts,
                            HM_ALK_G4_1_lyr$Counts,
                            HM_ALK_G4_2_lyr$Counts,
                            HM_ALK_G4_3_lyr$Counts,
                            HM_TKS_G1_1_lyr$Counts,
                            HM_TKS_G1_2_lyr$Counts,
                            HM_TKS_G1_3_lyr$Counts,
                            HM_TKS_G5_1_lyr$Counts,
                            HM_TKS_G5_2_lyr$Counts,
                            HM_TKS_G5_3_lyr$Counts)

all_samples_HM_lyr <- dplyr::slice(as.data.frame(all_samples_HM_lyr),
                                   -filter_lyr_HM)

all_samples_HM_lyr <- as.data.frame(vst(as.matrix(all_samples_HM_lyr)))

names(all_samples_HM_lyr) <- c("HM_lyr_G1_1",
                               "HM_lyr_G1_2",
                               "HM_lyr_G1_3",
                               "HM_lyr_G1_4",
                               "HM_lyr_G4_1",
                               "HM_lyr_G4_2",
                               "HM_lyr_G4_3",
                               "HM_RS7_G1_1",
                               "HM_RS7_G1_2",
                               "HM_RS7_G1_3",
                               "HM_RS7_G1_4",
                               "HM_RS7_G4_1",
                               "HM_RS7_G4_2",
                               "HM_ALK_G1_1",
                               "HM_ALK_G1_2",
                               "HM_ALK_G1_3",
                               "HM_ALK_G4_1",
                               "HM_ALK_G4_2",
                               "HM_ALK_G4_3",
                               "HM_TKS_G1_1",
                               "HM_TKS_G1_2",
                               "HM_TKS_G1_3",
                               "HM_TKS_G5_1",
                               "HM_TKS_G5_2",
                               "HM_TKS_G5_3")

all_samples_LL_lyr <- cbind(LL_lyr_G1_1$Counts,
                            LL_lyr_G1_2$Counts,
                            LL_lyr_G1_4$Counts,
                            LL_lyr_G1_5$Counts,
                            LL_lyr_G1_6$Counts,
                            LL_lyr_G4_1$Counts,
                            LL_lyr_G4_2$Counts,
                            LL_lyr_G4_3$Counts,
                            LL_lyr_G4_4$Counts,
                            LL_lyr_G4_5$Counts,
                            LL_RS7_G1_1_lyr$Counts,
                            LL_RS7_G1_2_lyr$Counts,
                            LL_RS7_G1_3_lyr$Counts,
                            LL_RS7_G1_4_lyr$Counts,
                            LL_RS7_G4_1_lyr$Counts,
                            LL_RS7_G4_2_lyr$Counts,
                            LL_RS7_G4_3_lyr$Counts,
                            LL_ALK_G1_1_lyr$Counts,
                            LL_ALK_G1_2_lyr$Counts,
                            LL_ALK_G1_3_lyr$Counts,
                            LL_ALK_G4_1_lyr$Counts,
                            LL_ALK_G4_2_lyr$Counts,
                            LL_ALK_G4_3_lyr$Counts,
                            LL_TKS_G1_1_lyr$Counts,
                            LL_TKS_G1_2_lyr$Counts,
                            LL_TKS_G1_3_lyr$Counts,
                            LL_TKS_G5_1_lyr$Counts,
                            LL_TKS_G5_2_lyr$Counts,
                            LL_TKS_G5_3_lyr$Counts)

all_samples_LL_lyr <- dplyr::slice(as.data.frame(all_samples_LL_lyr),
                                   -filter_lyr_LL)

all_samples_LL_lyr <- as.data.frame(vst(as.matrix(all_samples_LL_lyr)))

names(all_samples_LL_lyr) <- c("LL_lyr_G1_1",
                               "LL_lyr_G1_2",
                               "LL_lyr_G1_4",
                               "LL_lyr_G1_5",
                               "LL_lyr_G1_6",
                               "LL_lyr_G4_1",
                               "LL_lyr_G4_2",
                               "LL_lyr_G4_3",
                               "LL_lyr_G4_4",
                               "LL_lyr_G4_5",
                               "LL_RS7_G1_1",
                               "LL_RS7_G1_2",
                               "LL_RS7_G1_3",
                               "LL_RS7_G1_4",
                               "LL_RS7_G4_1",
                               "LL_RS7_G4_2",
                               "LL_RS7_G4_3",
                               "LL_ALK_G1_1",
                               "LL_ALK_G1_2",
                               "LL_ALK_G1_3",
                               "LL_ALK_G4_1",
                               "LL_ALK_G4_2",
                               "LL_ALK_G4_3",
                               "LL_TKS_G1_1",
                               "LL_TKS_G1_2",
                               "LL_TKS_G1_3",
                               "LL_TKS_G5_1",
                               "LL_TKS_G5_2",
                               "LL_TKS_G5_3")


## We select the top1000 most variable genes

top1000_HM_hal <- all_samples_HM_hal %>%
  rowwise() %>%
  mutate(variance = c_across(everything()) %>% var()) %>%
  ungroup() %>%
  slice_max(n = 1000, order_by = variance)
top1000_HM_hal <- top1000_HM_hal[,1:25]

top1000_LL_hal <- all_samples_LL_hal %>%
  rowwise() %>%
  mutate(variance = c_across(everything()) %>% var()) %>%
  ungroup() %>%
  slice_max(n = 1000, order_by = variance)
top1000_LL_hal <- top1000_LL_hal[,1:26]

top1000_HM_lyr <- all_samples_HM_lyr %>%
  rowwise() %>%
  mutate(variance = c_across(everything()) %>% var()) %>%
  ungroup() %>%
  slice_max(n = 1000, order_by = variance)
top1000_HM_lyr <- top1000_HM_lyr[,1:25]

top1000_LL_lyr <- all_samples_LL_lyr %>%
  rowwise() %>%
  mutate(variance = c_across(everything()) %>% var()) %>%
  ungroup() %>%
  slice_max(n = 1000, order_by = variance)
top1000_LL_lyr <- top1000_LL_lyr[,1:29]


## We plot the heatmap for each condition separately

Heatmap(top1000_HM_hal)
Heatmap(top1000_HM_lyr)
Heatmap(top1000_LL_hal)
Heatmap(top1000_LL_lyr)
