## This script will import methylation %
## of all samples for all overlapping Cs
## and output the average methylation 
## level in 100kb windows

### Import libraries 

library(data.table)
library(GenomicRanges)
library(stringr)
library(tidyverse)

### Command line arguments

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: cov gz file of interest
covgz_path <- comm_args[1]

# Third argument: progenitor side, hal or lyr
hal_or_lyr <- comm_args[2]

# Third argument: output name (with extension)
output_name <- comm_args[3]
  
### Import file to convert

covgz <- fread(covgz_path)

name_file <- function(covgz, hal_or_lyr){
  if (hal_or_lyr == "hal"){
    # Define chromosome sizes halleri
    names(covgz) <- c("coords_full", 
                      "hal_G1_1_perc", 
                      "hal_G1_1_totC",
                      "hal_G1_2_perc", 
                      "hal_G1_2_totC",
                      "hal_G1_3_perc", 
                      "hal_G1_3_totC",
                      "hal_G4_1_perc", 
                      "hal_G4_1_totC",
                      "hal_G4_2_perc", 
                      "hal_G4_2_totC",
                      "hal_G4_3_perc", 
                      "hal_G4_3_totC",
                      "RS7_G1_1_perc", 
                      "RS7_G1_1_totC",
                      "RS7_G1_2_perc", 
                      "RS7_G1_2_totC",
                      "RS7_G1_3_perc", 
                      "RS7_G1_3_totC",
                      "RS7_G4_1_perc", 
                      "RS7_G4_1_totC",
                      "RS7_G4_2_perc", 
                      "RS7_G4_2_totC",
                      "RS7_G4_3_perc", 
                      "RS7_G4_3_totC",
                      "ALK_G1_1_perc", 
                      "ALK_G1_1_totC",
                      "ALK_G1_2_perc", 
                      "ALK_G1_2_totC",
                      "ALK_G1_3_perc", 
                      "ALK_G1_3_totC",
                      "ALK_G4_1_perc", 
                      "ALK_G4_1_totC",
                      "ALK_G4_2_perc", 
                      "ALK_G4_2_totC",
                      "ALK_G4_3_perc", 
                      "ALK_G4_3_totC",
                      "TKS_G1_1_perc", 
                      "TKS_G1_1_totC",
                      "TKS_G1_2_perc", 
                      "TKS_G1_2_totC",
                      "TKS_G1_3_perc", 
                      "TKS_G1_3_totC",
                      "TKS_G5_1_perc", 
                      "TKS_G5_1_totC",
                      "TKS_G5_2_perc", 
                      "TKS_G5_2_totC",
                      "TKS_G5_3_perc", 
                      "TKS_G5_3_totC",
		      "context")
    return(covgz)
  }
  else {
    names(covgz) <- c("coords_full", 
                      "lyr_G1_1_perc", 
                      "lyr_G1_1_totC",
                      "lyr_G1_2_perc", 
                      "lyr_G1_2_totC",
                      "lyr_G1_3_perc", 
                      "lyr_G1_3_totC",
                      "lyr_G4_1_perc", 
                      "lyr_G4_1_totC",
                      "lyr_G4_2_perc", 
                      "lyr_G4_2_totC",
                      "lyr_G4_3_perc", 
                      "lyr_G4_3_totC",
                      "RS7_G1_1_perc", 
                      "RS7_G1_1_totC",
                      "RS7_G1_2_perc", 
                      "RS7_G1_2_totC",
                      "RS7_G1_3_perc", 
                      "RS7_G1_3_totC",
                      "RS7_G4_1_perc", 
                      "RS7_G4_1_totC",
                      "RS7_G4_2_perc", 
                      "RS7_G4_2_totC",
                      "RS7_G4_3_perc", 
                      "RS7_G4_3_totC",
                      "ALK_G1_1_perc", 
                      "ALK_G1_1_totC",
                      "ALK_G1_2_perc", 
                      "ALK_G1_2_totC",
                      "ALK_G1_3_perc", 
                      "ALK_G1_3_totC",
                      "ALK_G4_1_perc", 
                      "ALK_G4_1_totC",
                      "ALK_G4_2_perc", 
                      "ALK_G4_2_totC",
                      "ALK_G4_3_perc", 
                      "ALK_G4_3_totC",
                      "TKS_G1_1_perc", 
                      "TKS_G1_1_totC",
                      "TKS_G1_2_perc", 
                      "TKS_G1_2_totC",
                      "TKS_G1_3_perc", 
                      "TKS_G1_3_totC",
                      "TKS_G5_1_perc", 
                      "TKS_G5_1_totC",
                      "TKS_G5_2_perc", 
                      "TKS_G5_2_totC",
                      "TKS_G5_3_perc", 
                      "TKS_G5_3_totC",
		      "context")
    return(covgz)
  }
}

covgz <- name_file(covgz, hal_or_lyr)

## Split first coordinate column into three columns

split_coords <- str_split(covgz$coords_full, "_", n = 3, simplify = TRUE)

new_covgz <- mutate(covgz, scaffold = split_coords[,1],
                    start = as.numeric(split_coords[,2]),
                    end = as.numeric(split_coords[,3]))


# We select only the columns we will need for plotting to save memory

select_mem <- function(new_covgz, hal_or_lyr){
  if (hal_or_lyr == "hal"){
    covgz_final <- select(new_covgz, 
                          scaffold, 
                          start, 
                          end,
                          hal_G1_1_perc, 
                          hal_G1_2_perc, 
                          hal_G1_3_perc,
                          hal_G4_1_perc, 
                          hal_G4_2_perc, 
                          hal_G4_3_perc,
                          RS7_G1_1_perc, 
                          RS7_G1_2_perc, 
                          RS7_G1_3_perc,
                          RS7_G4_1_perc, 
                          RS7_G4_2_perc, 
                          RS7_G4_3_perc,
                          ALK_G1_1_perc, 
                          ALK_G1_2_perc, 
                          ALK_G1_3_perc,
                          ALK_G4_1_perc, 
                          ALK_G4_2_perc, 
                          ALK_G4_3_perc,
                          TKS_G1_1_perc, 
                          TKS_G1_2_perc, 
                          TKS_G1_3_perc,
                          TKS_G5_1_perc, 
                          TKS_G5_2_perc, 
                          TKS_G5_3_perc)
    return(covgz_final)
  }
  else{
    covgz_final <- select(new_covgz, 
                          scaffold, 
                          start, 
                          end,
                          lyr_G1_1_perc, 
                          lyr_G1_2_perc, 
                          lyr_G1_3_perc,
                          lyr_G4_1_perc, 
                          lyr_G4_2_perc, 
                          lyr_G4_3_perc,
                          RS7_G1_1_perc, 
                          RS7_G1_2_perc, 
                          RS7_G1_3_perc,
                          RS7_G4_1_perc, 
                          RS7_G4_2_perc, 
                          RS7_G4_3_perc,
                          ALK_G1_1_perc, 
                          ALK_G1_2_perc, 
                          ALK_G1_3_perc,
                          ALK_G4_1_perc, 
                          ALK_G4_2_perc, 
                          ALK_G4_3_perc,
                          TKS_G1_1_perc, 
                          TKS_G1_2_perc, 
                          TKS_G1_3_perc,
                          TKS_G5_1_perc, 
                          TKS_G5_2_perc, 
                          TKS_G5_3_perc)
    return(covgz_final)
    
  }
}

covgz_final <- select_mem(new_covgz, hal_or_lyr)

rm(covgz)
rm(new_covgz)

############## PART 2: AVERAGE METHYLATION ################ 

chrom_size <- function(hal_or_lyr){
  if (hal_or_lyr == "hal"){
    # Define chromosome sizes halleri
    chrom_sizes <- data.frame(Chrom = c("chr1",
                                        "chr2",
                                        "chr3",
                                        "chr4",
                                        "chr5",
                                        "chr6",
                                        "chr7",
                                        "chr8"),
                              Start = c(rep(1, 8)),
                              End = c(37868808, 19340014, 30847184,
                                      27260850, 25233939, 28041738,
                                      30089577, 26531562),
                              Name = c("chr1",
                                       "chr2",
                                       "chr3",
                                       "chr4",
                                       "chr5",
                                       "chr6",
                                       "chr7",
                                       "chr8"),
                              Colors = rep("lightgrey", 8))
    return(chrom_sizes)
  }
  else {
    # Define chromosome sizes lyrata
    chrom_sizes <- data.frame(Chrom = c("chr9",
                                        "chr10",
                                        "chr11",
                                        "chr12",
                                        "chr13",
                                        "chr14",
                                        "chr15",
                                        "chr16"),
                              Start = c(rep(1, 8)),
                              End = c(28718688, 21180362, 27642951,
                                      24725123, 23209782, 25845826,
                                      27879834, 20443964),
                              Name = c("chr9",
                                       "chr10",
                                       "chr11",
                                       "chr12",
                                       "chr13",
                                       "chr14",
                                       "chr15",
                                       "chr16"),
                              Colors = rep("lightgrey", 8))
    return(chrom_sizes)
  }
  
}

chrom_sizes <- chrom_size(hal_or_lyr = hal_or_lyr)

# Define vector with chromosome names

chrom_names <- chrom_sizes$End
names(chrom_names) <- chrom_sizes$Name

chrom_sizes_G <- GRanges(seqnames = chrom_sizes$Chrom,
                         ranges = IRanges(start = chrom_sizes$Start,
                                          end = chrom_sizes$End),
                         seqlengths = chrom_names)


windows <- tileGenome(chrom_names,
                      tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)

# Convert covgz final to GRanges to compute overlap

covgz_final_G <- GRanges(seqnames = covgz_final$scaffold,
                       ranges = IRanges(start = covgz_final$start,
                                        end = covgz_final$end)
                       )

dens <- GenomicRanges::countOverlaps(windows, covgz_final_G)

overlaps <- GenomicRanges::findOverlaps(windows, covgz_final_G)

window_scores_name <- function(hal_or_lyr){
  if (hal_or_lyr == "hal"){
    window_scores_df <- data.frame(hal_G1_1 = rep(0, length(dens)),
                                   hal_G1_2 = rep(0, length(dens)),
                                   hal_G1_3 = rep(0, length(dens)),
                                   hal_G4_1 = rep(0, length(dens)),
                                   hal_G4_2 = rep(0, length(dens)),
                                   hal_G4_3 = rep(0, length(dens)),
                                   RS7_G1_1 = rep(0, length(dens)),
                                   RS7_G1_2 = rep(0, length(dens)),
                                   RS7_G1_3 = rep(0, length(dens)),
                                   RS7_G4_1 = rep(0, length(dens)),
                                   RS7_G4_2 = rep(0, length(dens)),
                                   RS7_G4_3 = rep(0, length(dens)),
                                   ALK_G1_1 = rep(0, length(dens)),
                                   ALK_G1_2 = rep(0, length(dens)),
                                   ALK_G1_3 = rep(0, length(dens)),
                                   ALK_G4_1 = rep(0, length(dens)),
                                   ALK_G4_2 = rep(0, length(dens)),
                                   ALK_G4_3 = rep(0, length(dens)),
                                   TKS_G1_1 = rep(0, length(dens)),
                                   TKS_G1_2 = rep(0, length(dens)),
                                   TKS_G1_3 = rep(0, length(dens)),
                                   TKS_G5_1 = rep(0, length(dens)),
                                   TKS_G5_2 = rep(0, length(dens)),
                                   TKS_G5_3 = rep(0, length(dens)))
    return(window_scores_df)
  }
  
  else {
    window_scores_df <- data.frame(lyr_G1_1 = rep(0, length(dens)),
                                   lyr_G1_2 = rep(0, length(dens)),
                                   lyr_G1_3 = rep(0, length(dens)),
                                   lyr_G4_1 = rep(0, length(dens)),
                                   lyr_G4_2 = rep(0, length(dens)),
                                   lyr_G4_3 = rep(0, length(dens)),
                                   RS7_G1_1 = rep(0, length(dens)),
                                   RS7_G1_2 = rep(0, length(dens)),
                                   RS7_G1_3 = rep(0, length(dens)),
                                   RS7_G4_1 = rep(0, length(dens)),
                                   RS7_G4_2 = rep(0, length(dens)),
                                   RS7_G4_3 = rep(0, length(dens)),
                                   ALK_G1_1 = rep(0, length(dens)),
                                   ALK_G1_2 = rep(0, length(dens)),
                                   ALK_G1_3 = rep(0, length(dens)),
                                   ALK_G4_1 = rep(0, length(dens)),
                                   ALK_G4_2 = rep(0, length(dens)),
                                   ALK_G4_3 = rep(0, length(dens)),
                                   TKS_G1_1 = rep(0, length(dens)),
                                   TKS_G1_2 = rep(0, length(dens)),
                                   TKS_G1_3 = rep(0, length(dens)),
                                   TKS_G5_1 = rep(0, length(dens)),
                                   TKS_G5_2 = rep(0, length(dens)),
                                   TKS_G5_3 = rep(0, length(dens)))
    return(window_scores_df)
  }
}

window_scores_df <- window_scores_name(hal_or_lyr)

covgz_final <- as.data.frame(covgz_final)

for(i in 4:ncol(covgz_final)){
  print(paste0(round((i / ncol(covgz_final))*100, 1), "%"))
  for (j in 1:nrow(window_scores_df)){
    window_scores_df[j, i - 3] <- sum((covgz_final[overlaps[overlaps@from == j]@to, i])/dens[j])
  }
}

covgz_avg_fun <- function(hal_or_lyr){
  if (hal_or_lyr == "hal"){
    covgz_average <- data.frame(hal_G1 = apply(window_scores_df[,1:3], 1, mean),
                                hal_G4 = apply(window_scores_df[,4:6], 1, mean),
                                RS7_G1 = apply(window_scores_df[,7:9], 1, mean),
                                RS7_G4 = apply(window_scores_df[,10:12], 1, mean),
                                ALK_G1 = apply(window_scores_df[,13:15], 1, mean),
                                ALK_G4 = apply(window_scores_df[,16:18], 1, mean),
                                TKS_G1 = apply(window_scores_df[,19:21], 1, mean),
                                TKS_G5 = apply(window_scores_df[,22:24], 1, mean))
    return(covgz_average)
  }
  else{
    covgz_average <- data.frame(lyr_G1 = apply(window_scores_df[,1:3], 1, mean),
                                lyr_G4 = apply(window_scores_df[,4:6], 1, mean),
                                RS7_G1 = apply(window_scores_df[,7:9], 1, mean),
                                RS7_G4 = apply(window_scores_df[,10:12], 1, mean),
                                ALK_G1 = apply(window_scores_df[,13:15], 1, mean),
                                ALK_G4 = apply(window_scores_df[,16:18], 1, mean),
                                TKS_G1 = apply(window_scores_df[,19:21], 1, mean),
                                TKS_G5 = apply(window_scores_df[,22:24], 1, mean))
    return(covgz_average)
  }
}

covgz_average <- covgz_avg_fun(hal_or_lyr)

# write output 

fwrite(covgz_average, output_name, row.names = F,
            col.names = T, quote = F, sep = "\t")
