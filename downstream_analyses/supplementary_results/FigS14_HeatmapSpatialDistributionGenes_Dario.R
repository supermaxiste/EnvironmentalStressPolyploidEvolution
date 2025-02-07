### This script aims at plotting the location
### of all genes in lyrata and halleri along
### chromosomes

## Import libraries

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(karyoploteR)

## Import annotations

Ahal_anno <- read.delim("~/OneDrive/PhD/Project/bs_data/Ahal_genome/Ahal_v2_2.gff", header=FALSE)
names(Ahal_anno) <- c("scaffold", "source", "feature", 
                      "start", "end", "score", "strand",
                      "frame", "attribute")
Alyr_anno <- read.delim("~/OneDrive/PhD/Project/bs_data/Alyr_genome/Alyr_v2_2_renamed.gff", header=FALSE)
names(Alyr_anno) <- c("scaffold", "source", "feature", 
                      "start", "end", "score", "strand",
                      "frame", "attribute")

## Select genes only

Ahal_genes <- filter(Ahal_anno, feature == "gene")
Alyr_genes <- filter(Alyr_anno, feature == "gene")

## Remove full annotation

rm(Ahal_anno)
rm(Alyr_anno)

## Import agp files for coordinate conversion step

agp_hal <- read.delim("~/OneDrive/PhD/Project/Chapter2/AssemblyMapping/Ahal_mum_tiling.agp", header=FALSE, comment.char="#")

names(agp_hal) <- c("new_scaffold", "start", "end", "n", "WU", "old_scaffold", "type", "length", "orientiation")

agp_lyr <- read.delim("~/OneDrive/PhD/Project/Chapter2/AssemblyMapping/Alyr_mum_tiling.agp", header=FALSE, comment.char="#")

names(agp_lyr) <- c("new_scaffold", "start", "end", "n", "WU", "old_scaffold", "type", "length", "orientiation")

## Convert scaffold names in lyrata agp to matching names

agp_scaffold_renamer <- function(agp){
  
  # Get old scaffold numbers
  old_scaffolds <- agp$old_scaffold
  
  # Split the string to just get the number
  scaffold_number <- as.integer(str_split_fixed(old_scaffolds, "_", 2)[,2])
  
  # Add 2249 for the new numbering
  new_scaffold_number <- scaffold_number + 2239
  
  # Create new set of strings to put in the agp file
  new_scaffolds <- paste0("scaffold_", new_scaffold_number)
  
  # Replace old names with new names
  agp$old_scaffold <- new_scaffolds
  
  return(agp)
}
agp_lyr_renamed <- agp_scaffold_renamer(agp_lyr)

## Function to convert coordinates

coordinate_converter <- function(dmrseq_output, agp){
  
  ### Match old scaffolds to new chromosomes
  
  scaffold_match <- match(dmrseq_output$scaffold, agp$old_scaffold)
  
  ### Create vector with new scaffolds and new coordinates
  
  new_scaffolds <- agp$new_scaffold[scaffold_match]
  new_start <- agp$start[scaffold_match] + dmrseq_output$start
  new_end <- agp$start[scaffold_match] + dmrseq_output$end
  
  ### Modify dmrseq output with new coordinates
  
  new_dmrseq_output <- dmrseq_output
  new_dmrseq_output$seqnames <- new_scaffolds
  new_dmrseq_output$start <- new_start
  new_dmrseq_output$end <- new_end
  
  # Selecting only 8 chromosomes
  
  new_dmrseq_output <- new_dmrseq_output %>%
    dplyr::filter(seqnames=="scaffold_1_RagTag" | 
                    seqnames=="scaffold_2_RagTag" | 
                    seqnames=="scaffold_3_RagTag" | 
                    seqnames=="scaffold_4_RagTag" | 
                    seqnames=="scaffold_5_RagTag" | 
                    seqnames=="scaffold_6_RagTag" | 
                    seqnames=="scaffold_7_RagTag" | 
                    seqnames=="scaffold_8_RagTag")
  
  ### return file with new coordinates
  
  return(new_dmrseq_output)
}

Ahal_genes_chromCoord <- coordinate_converter(Ahal_genes, agp = agp_hal)
Alyr_genes_chromCoord <- coordinate_converter(Alyr_genes, agp = agp_lyr_renamed)

## Convert everything into GRanges

Ahal_genes_df <- data.frame(Chrom = as.character(Ahal_genes_chromCoord$seqnames), 
                                    Start = as.integer(Ahal_genes_chromCoord$start),
                                    End = as.integer(Ahal_genes_chromCoord$end),
                                    Name = as.character(Ahal_genes_chromCoord$seqnames),
                                    gieStain = "black")
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_1_RagTag"] <- 1
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_2_RagTag"] <- 2
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_3_RagTag"] <- 3
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_4_RagTag"] <- 4
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_5_RagTag"] <- 5
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_6_RagTag"] <- 6
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_7_RagTag"] <- 7
Ahal_genes_df$Chrom[Ahal_genes_df$Chrom == "scaffold_8_RagTag"] <- 8

Alyr_genes_df <- data.frame(Chrom = as.character(Alyr_genes_chromCoord$seqnames), 
                            Start = as.integer(Alyr_genes_chromCoord$start),
                            End = as.integer(Alyr_genes_chromCoord$end),
                            Name = as.character(Alyr_genes_chromCoord$seqnames),
                            gieStain = "black")

Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_1_RagTag"] <- 1
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_2_RagTag"] <- 2
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_3_RagTag"] <- 3
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_4_RagTag"] <- 4
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_5_RagTag"] <- 5
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_6_RagTag"] <- 6
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_7_RagTag"] <- 7
Alyr_genes_df$Chrom[Alyr_genes_df$Chrom == "scaffold_8_RagTag"] <- 8

# Define chromosome sizes halleri and chromosome names
# then turn into GRanges

chrom_sizes_hal <- data.frame(Chrom = c(1:8),
                          Start = c(rep(1, 8)),
                          End = c(29853912, 19430572, 22405798,
                                  22332524, 17489651, 21293240,
                                  25127945, 21784696),
                          Name = rep(NA, 8),
                          Colors = rep("lightgrey", 8))

chrom_names_hal <- chrom_sizes_hal$End
names(chrom_names_hal) <- chrom_sizes_hal$Chrom

chrom_sizes_hal_G <- GRanges(seqnames = chrom_sizes_hal$Chrom,
                         ranges = IRanges(start = chrom_sizes_hal$Start,
                                          end = chrom_sizes_hal$End),
                         seqlengths = chrom_names_hal)

# Define chromosome sizes lyrata and chromosome names
# then turn into GRanges

chrom_sizes_lyr <- data.frame(Chrom = c(1:8),
                          Start = c(rep(1, 8)),
                          End = c(26235489, 17728226, 22610115,
                                  18792104, 21074452, 18886653,
                                  22122096, 18984421),
                          Name = rep(NA, 8),
                          Colors = rep("lightgrey", 8))

chrom_names_lyr <- chrom_sizes_lyr$End
names(chrom_names_lyr) <- chrom_sizes_lyr$Chrom

chrom_sizes_lyr_G <- GRanges(seqnames = chrom_sizes_lyr$Chrom,
                             ranges = IRanges(start = chrom_sizes_lyr$Start,
                                              end = chrom_sizes_lyr$End),
                             seqlengths = chrom_names_lyr)

## Define windows to calculate average coverage

windows_hal <- tileGenome(chrom_names_hal,
                      tilewidth = 5e5, cut.last.tile.in.chrom = TRUE)

windows_lyr <- tileGenome(chrom_names_lyr,
                          tilewidth = 5e5, cut.last.tile.in.chrom = TRUE)

## Count how many genes are present in each window

geneDens_hal <- GenomicRanges::countOverlaps(windows_hal, toGRanges(Ahal_genes_df))

geneDens_lyr <- GenomicRanges::countOverlaps(windows_hal, toGRanges(Alyr_genes_df))


## Find which genes are part of which windows

geneWindow_hal <- GenomicRanges::findOverlaps(windows_hal, toGRanges(Ahal_genes_df))

geneWindow_lyr <- GenomicRanges::findOverlaps(windows_lyr, toGRanges(Alyr_genes_df))


# Path to count tables (featureCounts output)

folder_path <- "~/OneDrive/PhD/Project/Chapter2/RNAseq/04_count_tables/"
colnames <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "Counts")


# Read files (RPKM)

HM_hal_G1_1 <- fread(paste0(folder_path, "HM_hal_G1_1_counts_tpm.txt"), 
                     col.names = colnames)
HM_hal_G1_2 <- fread(paste0(folder_path, "HM_hal_G1_2_counts_tpm.txt"), 
                     col.names = colnames)
HM_lyr_G1_1 <- fread(paste0(folder_path, "HM_lyr_G1_1_counts_tpm.txt"), 
                     col.names = colnames)
HM_lyr_G1_2 <- fread(paste0(folder_path, "HM_lyr_G1_2_counts_tpm.txt"), 
                     col.names = colnames)

HM_RS7_G1_1_hal <- fread(paste0(folder_path, "HM_RS7_G1_1_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_1_lyr <- fread(paste0(folder_path, "HM_RS7_G1_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_2_hal <- fread(paste0(folder_path, "HM_RS7_G1_2_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_2_lyr <- fread(paste0(folder_path, "HM_RS7_G1_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_3_hal <- fread(paste0(folder_path, "HM_RS7_G1_3_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_3_lyr <- fread(paste0(folder_path, "HM_RS7_G1_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_4_hal <- fread(paste0(folder_path, "HM_RS7_G1_4_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G1_4_lyr <- fread(paste0(folder_path, "HM_RS7_G1_4_lyrcounts_tpm.txt"), 
                         col.names = colnames)

HM_RS7_G4_1_hal <- fread(paste0(folder_path, "HM_RS7_G4_1_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G4_1_lyr <- fread(paste0(folder_path, "HM_RS7_G4_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G4_2_hal <- fread(paste0(folder_path, "HM_RS7_G4_2_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_RS7_G4_2_lyr <- fread(paste0(folder_path, "HM_RS7_G4_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)

HM_ALK_G1_1_hal <- fread(paste0(folder_path, "HM_ALK_G1_1_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G1_1_lyr <- fread(paste0(folder_path, "HM_ALK_G1_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G1_2_hal <- fread(paste0(folder_path, "HM_ALK_G1_2_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G1_2_lyr <- fread(paste0(folder_path, "HM_ALK_G1_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G1_3_hal <- fread(paste0(folder_path, "HM_ALK_G1_3_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G1_3_lyr <- fread(paste0(folder_path, "HM_ALK_G1_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

HM_ALK_G4_1_hal <- fread(paste0(folder_path, "HM_ALK_G4_1_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G4_1_lyr <- fread(paste0(folder_path, "HM_ALK_G4_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G4_2_hal <- fread(paste0(folder_path, "HM_ALK_G4_2_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G4_2_lyr <- fread(paste0(folder_path, "HM_ALK_G4_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G4_3_hal <- fread(paste0(folder_path, "HM_ALK_G4_3_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_ALK_G4_3_lyr <- fread(paste0(folder_path, "HM_ALK_G4_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

HM_TKS_G1_1_hal <- fread(paste0(folder_path, "HM_TKS_G1_1_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G1_1_lyr <- fread(paste0(folder_path, "HM_TKS_G1_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G1_2_hal <- fread(paste0(folder_path, "HM_TKS_G1_2_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G1_2_lyr <- fread(paste0(folder_path, "HM_TKS_G1_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G1_3_hal <- fread(paste0(folder_path, "HM_TKS_G1_3_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G1_3_lyr <- fread(paste0(folder_path, "HM_TKS_G1_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

HM_TKS_G5_1_hal <- fread(paste0(folder_path, "HM_TKS_G5_1_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G5_1_lyr <- fread(paste0(folder_path, "HM_TKS_G5_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G5_2_hal <- fread(paste0(folder_path, "HM_TKS_G5_2_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G5_2_lyr <- fread(paste0(folder_path, "HM_TKS_G5_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G5_3_hal <- fread(paste0(folder_path, "HM_TKS_G5_3_halcounts_tpm.txt"), 
                         col.names = colnames)
HM_TKS_G5_3_lyr <- fread(paste0(folder_path, "HM_TKS_G5_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)


LL_hal_G1_1 <- fread(paste0(folder_path, "LL_hal_G1_1_counts_tpm.txt"), 
                     col.names = colnames)
LL_hal_G1_2 <- fread(paste0(folder_path, "LL_hal_G1_2_counts_tpm.txt"), 
                     col.names = colnames)
LL_hal_G4_1 <- fread(paste0(folder_path, "LL_hal_G4_1_counts_tpm.txt"), 
                     col.names = colnames)
LL_hal_G4_2 <- fread(paste0(folder_path, "LL_hal_G4_2_counts_tpm.txt"), 
                     col.names = colnames)

LL_lyr_G1_1 <- fread(paste0(folder_path, "LL_lyr_G1_1_counts_tpm.txt"), 
                     col.names = colnames)
LL_lyr_G1_2 <- fread(paste0(folder_path, "LL_lyr_G1_2_counts_tpm.txt"), 
                     col.names = colnames)
LL_lyr_G4_1 <- fread(paste0(folder_path, "LL_lyr_G4_1_counts_tpm.txt"), 
                     col.names = colnames)
LL_lyr_G4_2 <- fread(paste0(folder_path, "LL_lyr_G4_2_counts_tpm.txt"), 
                     col.names = colnames)
LL_lyr_G4_3 <- fread(paste0(folder_path, "LL_lyr_G4_3_counts_tpm.txt"), 
                     col.names = colnames)
LL_lyr_G4_4 <- fread(paste0(folder_path, "LL_lyr_G4_4_counts_tpm.txt"), 
                     col.names = colnames)

LL_RS7_G1_1_hal <- fread(paste0(folder_path, "LL_RS7_G1_1_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_1_lyr <- fread(paste0(folder_path, "LL_RS7_G1_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_2_hal <- fread(paste0(folder_path, "LL_RS7_G1_2_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_2_lyr <- fread(paste0(folder_path, "LL_RS7_G1_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_3_hal <- fread(paste0(folder_path, "LL_RS7_G1_3_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_3_lyr <- fread(paste0(folder_path, "LL_RS7_G1_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_4_hal <- fread(paste0(folder_path, "LL_RS7_G1_4_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G1_4_lyr <- fread(paste0(folder_path, "LL_RS7_G1_4_lyrcounts_tpm.txt"), 
                         col.names = colnames)

LL_RS7_G4_1_hal <- fread(paste0(folder_path, "LL_RS7_G4_1_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G4_1_lyr <- fread(paste0(folder_path, "LL_RS7_G4_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G4_2_hal <- fread(paste0(folder_path, "LL_RS7_G4_2_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G4_2_lyr <- fread(paste0(folder_path, "LL_RS7_G4_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G4_3_hal <- fread(paste0(folder_path, "LL_RS7_G4_3_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_RS7_G4_3_lyr <- fread(paste0(folder_path, "LL_RS7_G4_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

LL_ALK_G1_1_hal <- fread(paste0(folder_path, "LL_ALK_G1_1_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G1_1_lyr <- fread(paste0(folder_path, "LL_ALK_G1_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G1_2_hal <- fread(paste0(folder_path, "LL_ALK_G1_2_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G1_2_lyr <- fread(paste0(folder_path, "LL_ALK_G1_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G1_3_hal <- fread(paste0(folder_path, "LL_ALK_G1_3_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G1_3_lyr <- fread(paste0(folder_path, "LL_ALK_G1_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

LL_ALK_G4_1_hal <- fread(paste0(folder_path, "LL_ALK_G4_1_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G4_1_lyr <- fread(paste0(folder_path, "LL_ALK_G4_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G4_2_hal <- fread(paste0(folder_path, "LL_ALK_G4_2_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G4_2_lyr <- fread(paste0(folder_path, "LL_ALK_G4_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G4_3_hal <- fread(paste0(folder_path, "LL_ALK_G4_3_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_ALK_G4_3_lyr <- fread(paste0(folder_path, "LL_ALK_G4_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

LL_TKS_G1_1_hal <- fread(paste0(folder_path, "LL_TKS_G1_1_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G1_1_lyr <- fread(paste0(folder_path, "LL_TKS_G1_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G1_2_hal <- fread(paste0(folder_path, "LL_TKS_G1_2_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G1_2_lyr <- fread(paste0(folder_path, "LL_TKS_G1_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G1_3_hal <- fread(paste0(folder_path, "LL_TKS_G1_3_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G1_3_lyr <- fread(paste0(folder_path, "LL_TKS_G1_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

LL_TKS_G5_1_hal <- fread(paste0(folder_path, "LL_TKS_G5_1_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G5_1_lyr <- fread(paste0(folder_path, "LL_TKS_G5_1_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G5_2_hal <- fread(paste0(folder_path, "LL_TKS_G5_2_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G5_2_lyr <- fread(paste0(folder_path, "LL_TKS_G5_2_lyrcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G5_3_hal <- fread(paste0(folder_path, "LL_TKS_G5_3_halcounts_tpm.txt"), 
                         col.names = colnames)
LL_TKS_G5_3_lyr <- fread(paste0(folder_path, "LL_TKS_G5_3_lyrcounts_tpm.txt"), 
                         col.names = colnames)

# Compute average coverage for each window for each sample

computeCoverage <- function(sample, geneDens, geneWindows){
  coverage <- rep(0, length(geneDens))
  for (i in c(1:length(geneDens))){
    coverage[i] <-  sum(sample$Counts[geneWindows@to[geneWindows@from == i]]) / geneDens[i]
  }
  return(coverage)
}

HM_hal_G1_1_coverage <- computeCoverage(HM_hal_G1_1, geneDens_hal, geneWindow_hal)

HM_hal_G1_2_coverage <- computeCoverage(HM_hal_G1_2, geneDens_hal, geneWindow_hal)

HM_lyr_G1_1_coverage <- computeCoverage(HM_lyr_G1_1, geneDens_lyr, geneWindow_lyr)

HM_lyr_G1_2_coverage <- computeCoverage(HM_lyr_G1_2, geneDens_lyr, geneWindow_lyr)


HM_RS7_G1_1_hal_coverage <- computeCoverage(HM_RS7_G1_1_hal, geneDens_hal, geneWindow_hal)

HM_RS7_G1_1_lyr_coverage <- computeCoverage(HM_RS7_G1_1_lyr, geneDens_lyr, geneWindow_lyr)

HM_RS7_G1_2_hal_coverage <- computeCoverage(HM_RS7_G1_2_hal, geneDens_hal, geneWindow_hal)

HM_RS7_G1_2_lyr_coverage <- computeCoverage(HM_RS7_G1_2_lyr, geneDens_lyr, geneWindow_lyr)

HM_RS7_G1_3_hal_coverage <- computeCoverage(HM_RS7_G1_3_hal, geneDens_hal, geneWindow_hal)

HM_RS7_G1_3_lyr_coverage <- computeCoverage(HM_RS7_G1_3_lyr, geneDens_lyr, geneWindow_lyr)

HM_RS7_G1_4_hal_coverage <- computeCoverage(HM_RS7_G1_4_hal, geneDens_hal, geneWindow_hal)

HM_RS7_G1_4_lyr_coverage <- computeCoverage(HM_RS7_G1_4_lyr, geneDens_lyr, geneWindow_lyr)


HM_RS7_G4_1_hal_coverage <- computeCoverage(HM_RS7_G4_1_hal, geneDens_hal, geneWindow_hal)

HM_RS7_G4_1_lyr_coverage <- computeCoverage(HM_RS7_G4_1_lyr, geneDens_lyr, geneWindow_lyr)

HM_RS7_G4_2_hal_coverage <- computeCoverage(HM_RS7_G4_2_hal, geneDens_hal, geneWindow_hal)

HM_RS7_G4_2_lyr_coverage <- computeCoverage(HM_RS7_G4_2_lyr, geneDens_lyr, geneWindow_lyr)


HM_ALK_G1_1_hal_coverage <- computeCoverage(HM_ALK_G1_1_hal, geneDens_hal, geneWindow_hal)

HM_ALK_G1_1_lyr_coverage <- computeCoverage(HM_ALK_G1_1_lyr, geneDens_lyr, geneWindow_lyr)

HM_ALK_G1_2_hal_coverage <- computeCoverage(HM_ALK_G1_2_hal, geneDens_hal, geneWindow_hal)

HM_ALK_G1_2_lyr_coverage <- computeCoverage(HM_ALK_G1_2_lyr, geneDens_lyr, geneWindow_lyr)

HM_ALK_G1_3_hal_coverage <- computeCoverage(HM_ALK_G1_3_hal, geneDens_hal, geneWindow_hal)

HM_ALK_G1_3_lyr_coverage <- computeCoverage(HM_ALK_G1_3_lyr, geneDens_lyr, geneWindow_lyr)

HM_ALK_G4_1_hal_coverage <- computeCoverage(HM_ALK_G4_1_hal, geneDens_hal, geneWindow_hal)

HM_ALK_G4_1_lyr_coverage <- computeCoverage(HM_ALK_G4_1_lyr, geneDens_lyr, geneWindow_lyr)

HM_ALK_G4_2_hal_coverage <- computeCoverage(HM_ALK_G4_2_hal, geneDens_hal, geneWindow_hal)

HM_ALK_G4_2_lyr_coverage <- computeCoverage(HM_ALK_G4_2_lyr, geneDens_lyr, geneWindow_lyr)

HM_ALK_G4_3_hal_coverage <- computeCoverage(HM_ALK_G4_3_hal, geneDens_hal, geneWindow_hal)

HM_ALK_G4_3_lyr_coverage <- computeCoverage(HM_ALK_G4_3_lyr, geneDens_lyr, geneWindow_lyr)


HM_TKS_G1_1_hal_coverage <- computeCoverage(HM_TKS_G1_1_hal, geneDens_hal, geneWindow_hal)

HM_TKS_G1_1_lyr_coverage <- computeCoverage(HM_TKS_G1_1_lyr, geneDens_lyr, geneWindow_lyr)

HM_TKS_G1_2_hal_coverage <- computeCoverage(HM_TKS_G1_2_hal, geneDens_hal, geneWindow_hal)

HM_TKS_G1_2_lyr_coverage <- computeCoverage(HM_TKS_G1_2_lyr, geneDens_lyr, geneWindow_lyr)

HM_TKS_G1_3_hal_coverage <- computeCoverage(HM_TKS_G1_3_hal, geneDens_hal, geneWindow_hal)

HM_TKS_G1_3_lyr_coverage <- computeCoverage(HM_TKS_G1_3_lyr, geneDens_lyr, geneWindow_lyr)


HM_TKS_G5_1_hal_coverage <- computeCoverage(HM_TKS_G5_1_hal, geneDens_hal, geneWindow_hal)

HM_TKS_G5_1_lyr_coverage <- computeCoverage(HM_TKS_G5_1_lyr, geneDens_lyr, geneWindow_lyr)

HM_TKS_G5_2_hal_coverage <- computeCoverage(HM_TKS_G5_2_hal, geneDens_hal, geneWindow_hal)

HM_TKS_G5_2_lyr_coverage <- computeCoverage(HM_TKS_G5_2_lyr, geneDens_lyr, geneWindow_lyr)

HM_TKS_G5_3_hal_coverage <- computeCoverage(HM_TKS_G5_3_hal, geneDens_hal, geneWindow_hal)

HM_TKS_G5_3_lyr_coverage <- computeCoverage(HM_TKS_G5_3_lyr, geneDens_lyr, geneWindow_lyr)


LL_hal_G1_1_coverage <- computeCoverage(LL_hal_G1_1, geneDens_hal, geneWindow_hal)

LL_hal_G1_2_coverage <- computeCoverage(LL_hal_G1_2, geneDens_hal, geneWindow_hal)

LL_hal_G4_1_coverage <- computeCoverage(LL_hal_G4_1, geneDens_hal, geneWindow_hal)

LL_hal_G4_2_coverage <- computeCoverage(LL_hal_G4_2, geneDens_hal, geneWindow_hal)


LL_lyr_G1_1_coverage <- computeCoverage(LL_lyr_G1_1, geneDens_lyr, geneWindow_lyr)

LL_lyr_G1_2_coverage <- computeCoverage(LL_lyr_G1_2, geneDens_lyr, geneWindow_lyr)

LL_lyr_G4_1_coverage <- computeCoverage(LL_lyr_G4_1, geneDens_lyr, geneWindow_lyr)

LL_lyr_G4_2_coverage <- computeCoverage(LL_lyr_G4_2, geneDens_lyr, geneWindow_lyr)

LL_lyr_G4_3_coverage <- computeCoverage(LL_lyr_G4_3, geneDens_lyr, geneWindow_lyr)

LL_lyr_G4_4_coverage <- computeCoverage(LL_lyr_G4_4, geneDens_lyr, geneWindow_lyr)

LL_RS7_G1_1_hal_coverage <- computeCoverage(LL_RS7_G1_1_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G1_1_lyr_coverage <- computeCoverage(LL_RS7_G1_1_lyr, geneDens_lyr, geneWindow_lyr)

LL_RS7_G1_2_hal_coverage <- computeCoverage(LL_RS7_G1_2_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G1_2_lyr_coverage <- computeCoverage(LL_RS7_G1_2_lyr, geneDens_lyr, geneWindow_lyr)

LL_RS7_G1_3_hal_coverage <- computeCoverage(LL_RS7_G1_3_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G1_3_lyr_coverage <- computeCoverage(LL_RS7_G1_3_lyr, geneDens_lyr, geneWindow_lyr)

LL_RS7_G1_4_hal_coverage <- computeCoverage(LL_RS7_G1_4_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G1_4_lyr_coverage <- computeCoverage(LL_RS7_G1_4_lyr, geneDens_lyr, geneWindow_lyr)


LL_RS7_G4_1_hal_coverage <- computeCoverage(LL_RS7_G4_1_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G4_1_lyr_coverage <- computeCoverage(LL_RS7_G4_1_lyr, geneDens_lyr, geneWindow_lyr)

LL_RS7_G4_2_hal_coverage <- computeCoverage(LL_RS7_G4_2_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G4_2_lyr_coverage <- computeCoverage(LL_RS7_G4_2_lyr, geneDens_lyr, geneWindow_lyr)

LL_RS7_G4_3_hal_coverage <- computeCoverage(LL_RS7_G4_3_hal, geneDens_hal, geneWindow_hal)

LL_RS7_G4_3_lyr_coverage <- computeCoverage(LL_RS7_G4_3_lyr, geneDens_lyr, geneWindow_lyr)

LL_ALK_G1_1_hal_coverage <- computeCoverage(LL_ALK_G1_1_hal, geneDens_hal, geneWindow_hal)

LL_ALK_G1_1_lyr_coverage <- computeCoverage(LL_ALK_G1_1_lyr, geneDens_lyr, geneWindow_lyr)

LL_ALK_G1_2_hal_coverage <- computeCoverage(LL_ALK_G1_2_hal, geneDens_hal, geneWindow_hal)

LL_ALK_G1_2_lyr_coverage <- computeCoverage(LL_ALK_G1_2_lyr, geneDens_lyr, geneWindow_lyr)

LL_ALK_G1_3_hal_coverage <- computeCoverage(LL_ALK_G1_3_hal, geneDens_hal, geneWindow_hal)

LL_ALK_G1_3_lyr_coverage <- computeCoverage(LL_ALK_G1_3_lyr, geneDens_lyr, geneWindow_lyr)


LL_ALK_G4_1_hal_coverage <- computeCoverage(LL_ALK_G4_1_hal, geneDens_hal, geneWindow_hal)

LL_ALK_G4_1_lyr_coverage <- computeCoverage(LL_ALK_G4_1_lyr, geneDens_lyr, geneWindow_lyr)

LL_ALK_G4_2_hal_coverage <- computeCoverage(LL_ALK_G4_2_hal, geneDens_hal, geneWindow_hal)

LL_ALK_G4_2_lyr_coverage <- computeCoverage(LL_ALK_G4_2_lyr, geneDens_lyr, geneWindow_lyr)

LL_ALK_G4_3_hal_coverage <- computeCoverage(LL_ALK_G4_3_hal, geneDens_hal, geneWindow_hal)

LL_ALK_G4_3_lyr_coverage <- computeCoverage(LL_ALK_G4_3_lyr, geneDens_lyr, geneWindow_lyr)


LL_TKS_G1_1_hal_coverage <- computeCoverage(LL_TKS_G1_1_hal, geneDens_hal, geneWindow_hal)

LL_TKS_G1_1_lyr_coverage <- computeCoverage(LL_TKS_G1_1_lyr, geneDens_lyr, geneWindow_lyr)

LL_TKS_G1_2_hal_coverage <- computeCoverage(LL_TKS_G1_2_hal, geneDens_hal, geneWindow_hal)

LL_TKS_G1_2_lyr_coverage <- computeCoverage(LL_TKS_G1_2_lyr, geneDens_lyr, geneWindow_lyr)

LL_TKS_G1_3_hal_coverage <- computeCoverage(LL_TKS_G1_3_hal, geneDens_hal, geneWindow_hal)

LL_TKS_G1_3_lyr_coverage <- computeCoverage(LL_TKS_G1_3_lyr, geneDens_lyr, geneWindow_lyr)


LL_TKS_G5_1_hal_coverage <- computeCoverage(LL_TKS_G5_1_hal, geneDens_hal, geneWindow_hal)

LL_TKS_G5_1_lyr_coverage <- computeCoverage(LL_TKS_G5_1_lyr, geneDens_lyr, geneWindow_lyr)

LL_TKS_G5_2_hal_coverage <- computeCoverage(LL_TKS_G5_2_hal, geneDens_hal, geneWindow_hal)

LL_TKS_G5_2_lyr_coverage <- computeCoverage(LL_TKS_G5_2_lyr, geneDens_lyr, geneWindow_lyr)

LL_TKS_G5_3_hal_coverage <- computeCoverage(LL_TKS_G5_3_hal, geneDens_hal, geneWindow_hal)

LL_TKS_G5_3_lyr_coverage <- computeCoverage(LL_TKS_G5_3_lyr, geneDens_lyr, geneWindow_lyr)


## karyoploteR heatmap of all samples coverage for both sides

## create dataframe with all values

all_samples_coverage_hal <- c(HM_hal_G1_1_coverage,
                          HM_hal_G1_2_coverage,
                          HM_RS7_G1_1_hal_coverage,
                          HM_RS7_G1_2_hal_coverage,
                          HM_RS7_G1_3_hal_coverage,
                          HM_RS7_G1_4_hal_coverage,
                          HM_RS7_G4_1_hal_coverage,
                          HM_RS7_G4_2_hal_coverage,
                          HM_ALK_G1_1_hal_coverage,
                          HM_ALK_G1_2_hal_coverage,
                          HM_ALK_G1_3_hal_coverage,
                          HM_ALK_G4_1_hal_coverage,
                          HM_ALK_G4_2_hal_coverage,
                          HM_ALK_G4_3_hal_coverage,
                          HM_TKS_G1_1_hal_coverage,
                          HM_TKS_G1_2_hal_coverage,
                          HM_TKS_G1_3_hal_coverage,
                          HM_TKS_G5_1_hal_coverage,
                          HM_TKS_G5_2_hal_coverage,
                          HM_TKS_G5_3_hal_coverage,
                          LL_hal_G1_1_coverage,
                          LL_hal_G1_2_coverage,
                          LL_hal_G4_1_coverage,
                          LL_hal_G4_2_coverage,
                          LL_RS7_G1_1_hal_coverage,
                          LL_RS7_G1_2_hal_coverage,
                          LL_RS7_G1_3_hal_coverage,
                          LL_RS7_G1_4_hal_coverage,
                          LL_RS7_G4_1_hal_coverage,
                          LL_RS7_G4_2_hal_coverage,
                          LL_RS7_G4_3_hal_coverage,
                          LL_ALK_G1_1_hal_coverage,
                          LL_ALK_G1_2_hal_coverage,
                          LL_ALK_G1_3_hal_coverage,
                          LL_ALK_G4_1_hal_coverage,
                          LL_ALK_G4_2_hal_coverage,
                          LL_ALK_G4_3_hal_coverage,
                          LL_TKS_G1_1_hal_coverage,
                          LL_TKS_G1_2_hal_coverage,
                          LL_TKS_G1_3_hal_coverage,
                          LL_TKS_G5_1_hal_coverage,
                          LL_TKS_G5_2_hal_coverage,
                          LL_TKS_G5_3_hal_coverage)

all_samples_coverage_lyr <- c(HM_lyr_G1_1_coverage,
                              HM_lyr_G1_2_coverage,
                              HM_RS7_G1_1_lyr_coverage,
                              HM_RS7_G1_2_lyr_coverage,
                              HM_RS7_G1_3_lyr_coverage,
                              HM_RS7_G1_4_lyr_coverage,
                              HM_RS7_G4_1_lyr_coverage,
                              HM_RS7_G4_2_lyr_coverage,
                              HM_ALK_G1_1_lyr_coverage,
                              HM_ALK_G1_2_lyr_coverage,
                              HM_ALK_G1_3_lyr_coverage,
                              HM_ALK_G4_1_lyr_coverage,
                              HM_ALK_G4_2_lyr_coverage,
                              HM_ALK_G4_3_lyr_coverage,
                              HM_TKS_G1_1_lyr_coverage,
                              HM_TKS_G1_2_lyr_coverage,
                              HM_TKS_G1_3_lyr_coverage,
                              HM_TKS_G5_1_lyr_coverage,
                              HM_TKS_G5_2_lyr_coverage,
                              HM_TKS_G5_3_lyr_coverage,
                              LL_lyr_G1_1_coverage,
                              LL_lyr_G1_2_coverage,
                              LL_lyr_G4_1_coverage,
                              LL_lyr_G4_2_coverage,
                              LL_RS7_G1_1_lyr_coverage,
                              LL_RS7_G1_2_lyr_coverage,
                              LL_RS7_G1_3_lyr_coverage,
                              LL_RS7_G1_4_lyr_coverage,
                              LL_RS7_G4_1_lyr_coverage,
                              LL_RS7_G4_2_lyr_coverage,
                              LL_RS7_G4_3_lyr_coverage,
                              LL_ALK_G1_1_lyr_coverage,
                              LL_ALK_G1_2_lyr_coverage,
                              LL_ALK_G1_3_lyr_coverage,
                              LL_ALK_G4_1_lyr_coverage,
                              LL_ALK_G4_2_lyr_coverage,
                              LL_ALK_G4_3_lyr_coverage,
                              LL_TKS_G1_1_lyr_coverage,
                              LL_TKS_G1_2_lyr_coverage,
                              LL_TKS_G1_3_lyr_coverage,
                              LL_TKS_G5_1_lyr_coverage,
                              LL_TKS_G5_2_lyr_coverage,
                              LL_TKS_G5_3_lyr_coverage)


# replace NA with zeros

all_samples_coverage_hal[is.na(all_samples_coverage_hal)] <- 0
all_samples_coverage_lyr[is.na(all_samples_coverage_lyr)] <- 0

# define colors

col_samples_hal <- colByValue(all_samples_coverage_hal, colors = c("white", "black"), min = 0, max = 100)
col_samples_lyr <- colByValue(all_samples_coverage_lyr, colors = c("white", "black"), min = 0, max = 100)

# plot heatmap

nsamples_hal <- 45
nsamples_lyr <- 45

sampleNames_hal <- c("HM_hal_G1_1",
                     "HM_hal_G1_2",
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
                     "HM_TKS_G5_3",
                     "LL_hal_G1_1",
                     "LL_hal_G1_2",
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

sampleNames_lyr <- c("HM_lyr_G1_1",
                     "HM_lyr_G1_2",
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
                     "HM_TKS_G5_3",
                     "LL_lyr_G1_1",
                     "LL_lyr_G1_2",
                     "LL_lyr_G4_1",
                     "LL_lyr_G4_2",
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

pp <- getDefaultPlotParams(4)
pp$bottommargin <- 10
pp$leftmargin <- 0.02
pp$data1inmargin <- 3
pp$topmargin <- 15
pp$rightmargin <- 0.2

## plot halleri side 

kp <- plotKaryotype(genome = chrom_sizes_hal_G, labels.plotter = NULL, plot.type=4, main = "halleri-side (TPM)", plot.params = pp)
kpDataBackground(kp, col=colByChr(kp$genome, colors = c("#EEEEEE", "#999999")), data.panel = "ideogram")
kpText(kp, chr=kp$chromosomes, x=mid(kp$genome), y=0.5, labels = gsub("chr", "", kp$chromosomes), data.panel = "ideogram")
range <- 1
for(i in 1:nsamples_hal) {
  kpPlotRegions(kp, windows_hal, 
                col = col_samples_hal[range:(range+362)],
                r0=autotrack(i, nsamples_hal))
  kpAddLabels(kp, labels = sampleNames_hal[i], side = "right", r0=autotrack(i, nsamples_hal), cex=0.75, label.margin = 0.02)
  range <- range + 362
}

## plot lyrata side 

kp <- plotKaryotype(genome = chrom_sizes_lyr_G, labels.plotter = NULL, plot.type=4, main = "lyrata-side (TPM)", plot.params = pp)
kpDataBackground(kp, col=colByChr(kp$genome, colors = c("#EEEEEE", "#999999")), data.panel = "ideogram")
kpText(kp, chr=kp$chromosomes, x=mid(kp$genome), y=0.5, labels = gsub("chr", "", kp$chromosomes), data.panel = "ideogram")
range <- 1
for(i in 1:nsamples_lyr) {
  kpPlotRegions(kp, windows_lyr, 
                col = col_samples_lyr[range:(range+362)],
                r0=autotrack(i, nsamples_lyr))
  kpAddLabels(kp, labels = sampleNames_lyr[i], side = "right", r0=autotrack(i, nsamples_lyr), cex=0.75, label.margin = 0.02)
  range <- range + 362
}

