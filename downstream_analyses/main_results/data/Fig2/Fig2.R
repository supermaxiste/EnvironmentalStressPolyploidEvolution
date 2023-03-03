## This script will plot average methylation levels and coverage
## for our samples of interest

# Import libraries

library(circlize)
library(GenomicRanges)
library(ComplexHeatmap)

# Import data

data <- "/path/to/data"

HM_allhal_CG <- read.delim(paste0(data, "HM_allhal_avg_methylation_CG.txt"))
HM_allhal_CHG <- read.delim(paste0(data, "HM_allhal_avg_methylation_CHG.txt"))
HM_allhal_CHH <- read.delim(paste0(data, "HM_allhal_avg_methylation_CHH.txt"))

HM_alllyr_CG <- read.delim(paste0(data, "HM_alllyr_avg_methylation_CG.txt"))
HM_alllyr_CHG <- read.delim(paste0(data, "HM_alllyr_avg_methylation_CHG.txt"))
HM_alllyr_CHH <- read.delim(paste0(data, "HM_alllyr_avg_methylation_CHH.txt"))

LL_allhal_CG <- read.delim(paste0(data, "LL_allhal_avg_methylation_CG.txt"))
LL_allhal_CHG <- read.delim(paste0(data, "LL_allhal_avg_methylation_CHG.txt"))
LL_allhal_CHH <- read.delim(paste0(data, "LL_allhal_avg_methylation_CHH.txt"))

LL_alllyr_CG <- read.delim(paste0(data, "LL_alllyr_avg_methylation_CG.txt"))
LL_alllyr_CHG <- read.delim(paste0(data, "LL_alllyr_avg_methylation_CHG.txt"))
LL_alllyr_CHH <- read.delim(paste0(data, "LL_alllyr_avg_methylation_CHH.txt"))

# rename progenitor columns to be able to merge data

rename_cols <- function(df){
  names(df) <- c("pro_G1", "pro_G4",
                 "RS7_G1", "RS7_G4",
                 "ALK_G1", "ALK_G4",
                 "TKS_G1", "TKS_G5")
  return(df)
}

HM_allhal_CG <- rename_cols(HM_allhal_CG)
HM_allhal_CHG <- rename_cols(HM_allhal_CHG)
HM_allhal_CHH <- rename_cols(HM_allhal_CHH)

HM_alllyr_CG <- rename_cols(HM_alllyr_CG)
HM_alllyr_CHG <- rename_cols(HM_alllyr_CHG)
HM_alllyr_CHH <- rename_cols(HM_alllyr_CHH)

LL_allhal_CG <- rename_cols(LL_allhal_CG)
LL_allhal_CHG <- rename_cols(LL_allhal_CHG)
LL_allhal_CHH <- rename_cols(LL_allhal_CHH)

LL_alllyr_CG <- rename_cols(LL_alllyr_CG)
LL_alllyr_CHG <- rename_cols(LL_alllyr_CHG)
LL_alllyr_CHH <- rename_cols(LL_alllyr_CHH)

# Prepare data for lines in circlize

## halleri

chrom_names <- c(37868808, 19340014, 30847184,
                 27260850, 25233939, 28041738,
                 30089577, 26531562)

names(chrom_names) <- c(1:8)

chrom_sizes_G <- GRanges(seqnames = c(1:8),
                         ranges = IRanges(start = 1,
                                          end = chrom_names),
                         seqlengths = chrom_names)


windows <- tileGenome(chrom_names,
                      tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)

## lyrata

chrom_names_lyr <- c(28718688, 21180362, 27642951,
                     24725123, 23209782, 25845826,
                     27879834, 20443964)

names(chrom_names_lyr) <- c(1:8)

chrom_sizes_Glyr <- GRanges(seqnames = c(1:8),
                            ranges = IRanges(start = 1,
                                             end = chrom_names_lyr),
                            seqlengths = chrom_names_lyr)


windows_lyr <- tileGenome(chrom_names_lyr,
                          tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)

# We import coverage data to add coverage for 4th generation synthetics
# for both mild and stressful conditions

folder_path <- c("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/")

HM_RS7_G4_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_hal_dens.txt"))
HM_RS7_G4_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_hal_dens.txt"))
HM_RS7_G4_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_hal_dens.txt"))

HM_RS7_G4_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_lyr_dens.txt"))
HM_RS7_G4_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_lyr_dens.txt"))
HM_RS7_G4_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_lyr_dens.txt"))

LL_RS7K_G4_1_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_1_hal_dens.txt"))
LL_RS7K_G4_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_2_hal_dens.txt"))
LL_RS7K_G4_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_3_hal_dens.txt"))

LL_RS7K_G4_1_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_1_lyr_dens.txt"))
LL_RS7K_G4_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_2_lyr_dens.txt"))
LL_RS7K_G4_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_3_lyr_dens.txt"))

# total cytonsine density

hal_totC_dens <- read.delim(paste0(folder_path,"hal_totC_dens_Dario.txt"))
lyr_totC_dens <- read.delim(paste0(folder_path,"lyr_totC_dens_Dario.txt"))

### We first compute the real cytosine coverage per window for all samples

HM_RS7_G4_1_hal_real_coverage <- HM_RS7_G4_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_2_hal_real_coverage <- HM_RS7_G4_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_3_hal_real_coverage <- HM_RS7_G4_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7_G4_1_lyr_real_coverage <- HM_RS7_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_2_lyr_real_coverage <- HM_RS7_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_3_lyr_real_coverage <- HM_RS7_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7K_G4_1_hal_real_coverage <- LL_RS7K_G4_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_2_hal_real_coverage <- LL_RS7K_G4_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_3_hal_real_coverage <- LL_RS7K_G4_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7K_G4_1_lyr_real_coverage <- LL_RS7K_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_2_lyr_real_coverage <- LL_RS7K_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_3_lyr_real_coverage <- LL_RS7K_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

### Compute the average for both synthetics

HM_RS7_G4_hal_avg_coverage <- c((HM_RS7_G4_1_hal_real_coverage +
                                   HM_RS7_G4_2_hal_real_coverage +
                                   HM_RS7_G4_3_hal_real_coverage) / 3)

HM_RS7_G4_lyr_avg_coverage <- c((HM_RS7_G4_1_lyr_real_coverage +
                                   HM_RS7_G4_2_lyr_real_coverage +
                                   HM_RS7_G4_3_lyr_real_coverage) / 3)

LL_RS7K_G4_hal_avg_coverage <- c((LL_RS7K_G4_1_hal_real_coverage +
                                    LL_RS7K_G4_2_hal_real_coverage +
                                    LL_RS7K_G4_3_hal_real_coverage) / 3)

LL_RS7K_G4_lyr_avg_coverage <- c((LL_RS7K_G4_1_lyr_real_coverage +
                                    LL_RS7K_G4_2_lyr_real_coverage +
                                    LL_RS7K_G4_3_lyr_real_coverage) / 3)


# create groups for hal and lyr


groups_hal <- c(rep(1, windows@seqnames@lengths[1]),
                rep(2, windows@seqnames@lengths[2]),
                rep(3, windows@seqnames@lengths[3]),
                rep(4, windows@seqnames@lengths[4]),
                rep(5, windows@seqnames@lengths[5]),
                rep(6, windows@seqnames@lengths[6]),
                rep(7, windows@seqnames@lengths[7]),
                rep(8, windows@seqnames@lengths[8]))

groups_lyr <- c(rep(1, windows_lyr@seqnames@lengths[1]),
                rep(2, windows_lyr@seqnames@lengths[2]),
                rep(3, windows_lyr@seqnames@lengths[3]),
                rep(4, windows_lyr@seqnames@lengths[4]),
                rep(5, windows_lyr@seqnames@lengths[5]),
                rep(6, windows_lyr@seqnames@lengths[6]),
                rep(7, windows_lyr@seqnames@lengths[7]),
                rep(8, windows_lyr@seqnames@lengths[8]))

groups_all <-c(rep(1, windows@seqnames@lengths[1]),
               rep(2, windows@seqnames@lengths[2]),
               rep(3, windows@seqnames@lengths[3]),
               rep(4, windows@seqnames@lengths[4]),
               rep(5, windows@seqnames@lengths[5]),
               rep(6, windows@seqnames@lengths[6]),
               rep(7, windows@seqnames@lengths[7]),
               rep(8, windows@seqnames@lengths[8]),
               rep(9, windows_lyr@seqnames@lengths[1]),
               rep(10, windows_lyr@seqnames@lengths[2]),
               rep(11, windows_lyr@seqnames@lengths[3]),
               rep(12, windows_lyr@seqnames@lengths[4]),
               rep(13, windows_lyr@seqnames@lengths[5]),
               rep(14, windows_lyr@seqnames@lengths[6]),
               rep(15, windows_lyr@seqnames@lengths[7]),
               rep(16, windows_lyr@seqnames@lengths[8]))

groups_all <- factor(groups_all, levels = 
                       c(1:16))

### We build the circlize plot for mild conditions first

chromosomes <- as.character(c(1:16))
scaffold_length <- as.matrix(cbind(c(1:16), c(chrom_names, 
                                              chrom_names_lyr)))

# We split the plot in two parts

par(mar = c(2, 2, 2, 2), mfrow = c(3, 2))

# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")
# circos.track(chromosomes, ylim = c(-20, 20), bg.border = "black")

# Define loop to add heatmaps to each track and sector

HM_avg <- c(HM_RS7_G4_lyr_avg_coverage,
            HM_RS7_G4_hal_avg_coverage)

LL_avg <- c(LL_RS7K_G4_lyr_avg_coverage,
            LL_RS7K_G4_hal_avg_coverage)

HM_all_CG <- rbind(HM_alllyr_CG, HM_allhal_CG)
HM_all_CHG <- rbind(HM_alllyr_CHG, HM_allhal_CHG)
HM_all_CHH <- rbind(HM_alllyr_CHH, HM_allhal_CHH)

LL_all_CG <- rbind(LL_alllyr_CG, LL_allhal_CG)
LL_all_CHG <- rbind(LL_alllyr_CHG, LL_allhal_CHG)
LL_all_CHH <- rbind(LL_alllyr_CHH, LL_allhal_CHH)

col_meth <- colorRamp2(c(0, 50, 100), c("blue", "yellow", "red"))
col_meth2 <- colorRamp2(c(0, 5), c("white", "black"))

circos.clear()

# CG, cold conditions
# We thin the tracks

circos.par(track.height = 0.05, 
           gap.degree = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,20),
           start.degree = 87,
           cell.padding = c(0, 0, 0, 0))

# Initialize heatmap to show progenitors sides

circos.heatmap.initialize(HM_all_CG, split = as.factor(groups_all),
                          cluster = F)

# Add nine tracks

circos.track(chromosomes, ylim = c(-10,10), bg.border = "black")

for (i in c(1:8)){
  highlight.sector(as.factor(i), col = "#CD9600")
}

for (i in c(9:16)){
  highlight.sector(as.factor(i), col = "#F8766D")
}

circos.heatmap(mat = HM_avg, split = as.factor(groups_all),
               col = col_meth2, track.height = 0.05,
               show.sector.labels = F, cluster = F)

#lgd = Legend(title = "coverage", col_fun = col_meth2)
#grid.draw(lgd)

for (i in c(1:8)){
  mat <- as.matrix(HM_all_CG[,i])
    circos.heatmap(mat, split = as.factor(groups_all),
                   col = col_meth, track.height = 0.05,
                   bg.border = "black", show.sector.labels = F,
                   cluster = F)
  
}

circos.clear()

# CG, hot conditions

# We thin the tracks

circos.par(track.height = 0.05, 
           gap.degree = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,20),
           start.degree = 87,
           cell.padding = c(0, 0, 0, 0))

# Initialize heatmap to show progenitors sides

circos.heatmap.initialize(LL_all_CG, split = as.factor(groups_all), 
                          cluster = F)

# Add nine tracks

circos.track(chromosomes, ylim = c(-10,10), bg.border = "black")
for (i in c(1:8)){
  highlight.sector(as.factor(i), col = "#CD9600")
}

for (i in c(9:16)){
  highlight.sector(as.factor(i), col = "#F8766D")
}

circos.heatmap(mat = LL_avg, split = as.factor(groups_all),
               col = col_meth2, track.height = 0.05,
               show.sector.labels = F, cluster = F)

for (i in c(1:8)){
  mat <- as.matrix(LL_all_CG[,i])
  circos.heatmap(mat, split = as.factor(groups_all),
                 col = col_meth, track.height = 0.05,
                 bg.border = "black", show.sector.labels = F,
                 cluster = F)
  
}

circos.clear()

# CHG, cold conditions

# We thin the tracks

circos.par(track.height = 0.05, 
           gap.degree = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,20),
           start.degree = 87,
           cell.padding = c(0, 0, 0, 0))

# Initialize heatmap to show progenitors sides

circos.heatmap.initialize(HM_all_CHG, split = as.factor(groups_all),
                          cluster = F)

# Add nine tracks

circos.track(chromosomes, ylim = c(-10,10), bg.border = "black")
for (i in c(1:8)){
  highlight.sector(as.factor(i), col = "#CD9600")
}

for (i in c(9:16)){
  highlight.sector(as.factor(i), col = "#F8766D")
}

circos.heatmap(mat = HM_avg, split = as.factor(groups_all),
               col = col_meth2, track.height = 0.05,
               show.sector.labels = F, cluster = F)

for (i in c(1:8)){
  mat <- as.matrix(HM_all_CHG[,i])
  circos.heatmap(mat, split = as.factor(groups_all),
                 col = col_meth, track.height = 0.05,
                 bg.border = "black", show.sector.labels = F,
                 cluster = F)
  
}

circos.clear()

# CHG, hot conditions

# We thin the tracks

circos.par(track.height = 0.05, 
           gap.degree = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,20),
           start.degree = 87,
           cell.padding = c(0, 0, 0, 0))

# Initialize heatmap to show progenitors sides

circos.heatmap.initialize(LL_all_CHG, split = as.factor(groups_all),
                          cluster = F)

# Add nine tracks

circos.track(chromosomes, ylim = c(-10,10), bg.border = "black")
for (i in c(1:8)){
  highlight.sector(as.factor(i), col = "#CD9600")
}

for (i in c(9:16)){
  highlight.sector(as.factor(i), col = "#F8766D")
}

circos.heatmap(mat = LL_avg, split = as.factor(groups_all),
               col = col_meth2, track.height = 0.05,
               show.sector.labels = F,
               cluster = F)

for (i in c(1:8)){
  mat <- as.matrix(LL_all_CHG[,i])
  circos.heatmap(mat, split = as.factor(groups_all),
                 col = col_meth, track.height = 0.05,
                 bg.border = "black", show.sector.labels = F)
  
}

circos.clear()

# CHH, cold conditions
# We thin the tracks

circos.par(track.height = 0.05, 
           gap.degree = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,20),
           start.degree = 87,
           cell.padding = c(0, 0, 0, 0))

# Initialize heatmap to show progenitors sides

circos.heatmap.initialize(HM_all_CHH, split = as.factor(groups_all),
                          cluster = F)

# Add nine tracks

circos.track(chromosomes, ylim = c(-10,10), bg.border = "black")
for (i in c(1:8)){
  highlight.sector(as.factor(i), col = "#CD9600")
}

for (i in c(9:16)){
  highlight.sector(as.factor(i), col = "#F8766D")
}

circos.heatmap(mat = HM_avg, split = as.factor(groups_all),
               col = col_meth2, track.height = 0.05,
               show.sector.labels = F, cluster = F)

for (i in c(1:8)){
  mat <- as.matrix(HM_all_CHH[,i])
  circos.heatmap(mat, split = as.factor(groups_all),
                 col = col_meth, track.height = 0.05,
                 bg.border = "black", show.sector.labels = F,
                 cluster = F)
  
}

# CHH, hot conditions
# We thin the tracks

circos.par(track.height = 0.05, 
           gap.degree = c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,20),
           start.degree = 87,
           cell.padding = c(0, 0, 0, 0))

# Initialize heatmap to show progenitors sides

circos.heatmap.initialize(LL_all_CHH, split = as.factor(groups_all),
                          cluster = F)

# Add nine tracks

circos.track(chromosomes, ylim = c(-10,10), bg.border = "black")
for (i in c(1:8)){
  highlight.sector(as.factor(i), col = "#CD9600")
}

for (i in c(9:16)){
  highlight.sector(as.factor(i), col = "#F8766D")
}

circos.heatmap(mat = LL_avg, split = as.factor(groups_all),
               col = col_meth2, track.height = 0.05,
               show.sector.labels = F, cluster = F)

for (i in c(1:8)){
  mat <- as.matrix(LL_all_CHH[,i])
  circos.heatmap(mat, split = as.factor(groups_all),
                 col = col_meth, track.height = 0.05,
                 bg.border = "black", show.sector.labels = F,
                 cluster = F)
  
}

