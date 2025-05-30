### DNA spatial distribution

### Import libraries

library(tidyverse)
library(karyoploteR)
library(data.table)

## Path to density data folder

folder_path <- c("~/OneDrive/PhD/Project/Chapter_3/new_coverage_analysis_v2/00_data/")

## Import all densities

HM_hal_G1_1_dens <- read.csv(paste0(folder_path, "HM_hal_G1_1_dens.txt"))

HM_lyr_G1_1_dens <- read.csv(paste0(folder_path, "HM_lyr_G1_1_dens.txt"))

## polyploid data

HM_RS7K_G1_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_1_hal_dens.txt"))
HM_RS7K_G1_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_2_hal_dens.txt"))
HM_RS7K_G1_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_3_hal_dens.txt"))

HM_RS7K_G1_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_1_lyr_dens.txt"))
HM_RS7K_G1_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_2_lyr_dens.txt"))
HM_RS7K_G1_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_3_lyr_dens.txt"))

HM_RS7_G4_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_hal_dens.txt"))
HM_RS7_G4_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_hal_dens.txt"))
HM_RS7_G4_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_hal_dens.txt"))

HM_RS7_G4_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_lyr_dens.txt"))
HM_RS7_G4_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_lyr_dens.txt"))
HM_RS7_G4_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_lyr_dens.txt"))

LL_RS7_G1_1_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_1_hal_dens.txt"))
LL_RS7_G1_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_2_hal_dens.txt"))
LL_RS7_G1_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_3_hal_dens.txt"))

LL_RS7_G1_1_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_1_lyr_dens.txt"))
LL_RS7_G1_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_2_lyr_dens.txt"))
LL_RS7_G1_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_3_lyr_dens.txt"))

LL_RS7K_G4_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_2_hal_dens.txt"))
LL_RS7K_G4_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_3_hal_dens.txt"))

LL_RS7K_G4_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_2_lyr_dens.txt"))
LL_RS7K_G4_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G4_3_lyr_dens.txt"))

## total density

hal_totC_dens <- read.delim("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/hal_totC_dens_Dario.txt")
lyr_totC_dens <- read.delim("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/lyr_totC_dens_Dario.txt")

### We first compute the real cytosine coverage per window for all samples

## polyploid data

HM_RS7K_G1_1_hal_real_coverage <- HM_RS7K_G1_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7K_G1_2_hal_real_coverage <- HM_RS7K_G1_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7K_G1_3_hal_real_coverage <- HM_RS7K_G1_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7K_G1_1_lyr_real_coverage <- HM_RS7K_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7K_G1_2_lyr_real_coverage <- HM_RS7K_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7K_G1_3_lyr_real_coverage <- HM_RS7K_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_RS7_G4_1_hal_real_coverage <- HM_RS7_G4_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_2_hal_real_coverage <- HM_RS7_G4_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_3_hal_real_coverage <- HM_RS7_G4_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7_G4_1_lyr_real_coverage <- HM_RS7_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_2_lyr_real_coverage <- HM_RS7_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_3_lyr_real_coverage <- HM_RS7_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7_G1_1_hal_real_coverage <- LL_RS7_G1_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G1_2_hal_real_coverage <- LL_RS7_G1_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G1_3_hal_real_coverage <- LL_RS7_G1_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7_G1_1_lyr_real_coverage <- LL_RS7_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G1_2_lyr_real_coverage <- LL_RS7_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G1_3_lyr_real_coverage <- LL_RS7_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7K_G4_2_hal_real_coverage <- LL_RS7K_G4_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_3_hal_real_coverage <- LL_RS7K_G4_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7K_G4_2_lyr_real_coverage <- LL_RS7K_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_3_lyr_real_coverage <- LL_RS7K_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

# Compute average for each sample

HM_RS7K_G1_hal_avg_coverage <- c((HM_RS7K_G1_1_hal_real_coverage +
                                    HM_RS7K_G1_2_hal_real_coverage +
                                    HM_RS7K_G1_3_hal_real_coverage) / 3)

HM_RS7K_G1_lyr_avg_coverage <- c((HM_RS7K_G1_1_lyr_real_coverage +
                                    HM_RS7K_G1_2_lyr_real_coverage +
                                    HM_RS7K_G1_3_lyr_real_coverage) / 3)

HM_RS7_G4_hal_avg_coverage <- c((HM_RS7_G4_1_hal_real_coverage +
                                   HM_RS7_G4_2_hal_real_coverage +
                                   HM_RS7_G4_3_hal_real_coverage) / 3)

HM_RS7_G4_lyr_avg_coverage <- c((HM_RS7_G4_1_lyr_real_coverage +
                                   HM_RS7_G4_2_lyr_real_coverage +
                                   HM_RS7_G4_3_lyr_real_coverage) / 3)

LL_RS7_G1_hal_avg_coverage <- c((LL_RS7_G1_1_hal_real_coverage +
                                   LL_RS7_G1_2_hal_real_coverage +
                                   LL_RS7_G1_3_hal_real_coverage) / 3)

LL_RS7_G1_lyr_avg_coverage <- c((LL_RS7_G1_1_lyr_real_coverage +
                                   LL_RS7_G1_2_lyr_real_coverage +
                                   LL_RS7_G1_3_lyr_real_coverage) / 3)

LL_RS7K_G4_hal_avg_coverage <- c((LL_RS7K_G4_2_hal_real_coverage +
                                    LL_RS7K_G4_3_hal_real_coverage) / 2)

LL_RS7K_G4_lyr_avg_coverage <- c((LL_RS7K_G4_2_lyr_real_coverage +
                                    LL_RS7K_G4_3_lyr_real_coverage) / 2)


### We define a function to plot DEG 
### on chromosomes with karyoploteR


DEGs_karyoploter <- function(deg_output, 
                             side, 
                             coverage, 
                             chrom_length,
                             homeo_regions,
                             nonhomeo_regions){

  direction <- ifelse(deg_output$dir>0, "increase", "decrease")

  deg_output <- mutate(deg_output, direction = direction)
  
  # We import the annotations to add genomic coordinates
  # for each DEG based on their geneID
  
  annotation_import <- function(side){
    if(side=="hal"){
      anno <- read.delim("~/OneDrive/PhD/Project/bs_data/Ahal_genome_Dario/hal_MAKER_PASA_46264_noTEs.gff",
                         header=FALSE,
                         col.names = c("chromosome", 
                                       "tool", 
                                       "context", 
                                       "start", 
                                       "end", 
                                       "number",
                                       "strand", 
                                       "dot", 
                                       "extra"))
    }
    else {
      anno <- read.delim("~/OneDrive/PhD/Project/bs_data/Alyr_genome_Dario/lyr_r2_MAKER_28737_8chr_renamed.gff",
                         header=FALSE,
                         col.names = c("chromosome", 
                                       "tool", 
                                       "context", 
                                       "start", 
                                       "end", 
                                       "number",
                                       "strand", 
                                       "dot", 
                                       "extra"))
    }
    return(anno)
  }
  
  anno <- annotation_import(side)
  
  anno <- filter(anno, context=="gene")
  
  # we get the geneID from the extra column and add it to the annotation
  
  geneID_anno <- str_split_fixed( # get all the characters after the =
    as.vector(
    # get all the characters before the ;
    str_split_fixed(anno$extra, ";", 2)[,1]
    ), "=", 2
  )[,2]
  
  anno <- anno %>%
    mutate(geneID = geneID_anno)
  
  deg_coordinates <- anno %>%
    filter(geneID %in% deg_output$geneID)
  
  new_deg_output <- deg_output %>%
    filter(geneID %in% deg_coordinates$geneID) %>%
    mutate(chromosome = deg_coordinates$chromosome,
           start = deg_coordinates$start,
           end = deg_coordinates$end)


  # Define chromosome sizes
  
  import_chromSizes <- function(side){
    if(side=="hal"){
      chrom_sizes <- data.frame(Chrom = c(1:8),
                                Start = c(rep(1, 8)),
                                End = c(37868808, 
                                        19340014, 
                                        30847184,
                                        27260850, 
                                        25233939, 
                                        28041738,
                                        30089577, 
                                        26531562),
                                Name = rep(NA, 8),
                                Colors = rep("lightgrey", 
                                             8))
      
    }
    
    else{
      # Define chromosome sizes lyrata
      chrom_sizes <- data.frame(Chrom = c(9:16),
                                Start = c(rep(1, 8)),
                                End = c(28718688, 
                                        21180362, 
                                        27642951,
                                        24725123, 
                                        23209782, 
                                        25845826,
                                        27879834, 
                                        20443964),
                                Name = rep(NA, 8),
                                Colors = rep("lightgrey", 
                                             8))
    }
    return(chrom_sizes)
  }
  
  # Import chromosome sizes
  
  chrom_sizes <- import_chromSizes(side)

  # Define vector with chromosome names
  chrom_names <- chrom_sizes$End
  names(chrom_names) <- chrom_sizes$Chrom

  chrom_sizes_G <- GRanges(seqnames = chrom_sizes$Chrom,
                           ranges = IRanges(start = chrom_sizes$Start,
                                            end = chrom_sizes$End),
                           seqlengths = chrom_names)

  # Turn chromosome scaffolds into simple numbers
  
  new_deg_output$chromosome <- as.numeric(
    str_replace(new_deg_output$chromosome, "chr", ""))

  # Turn deg output into GRanges object
  
  new_deg_output_G <- GRanges(seqnames = new_deg_output$chromosome,
                              ranges = IRanges(start = new_deg_output$start,
                                               end = new_deg_output$end))
  new_deg_output_G <- GenomicRanges::reduce(new_deg_output_G)
  
  get_heatmap_data <- function(side){
    if(side=="hal"){
      seqnames_dens <- as.numeric(str_replace_all(HM_hal_G1_1_dens$seqnames, "chr", ""))
      heatmap_data <- GRanges(seqnames = seqnames_dens,
                              IRanges(start = HM_hal_G1_1_dens$start,
                                      end = HM_hal_G1_1_dens$end))
      return(heatmap_data)
    }
    else{
      seqnames_dens <- as.numeric(str_replace_all(HM_lyr_G1_1_dens$seqnames, "chr", ""))
      heatmap_data <- GRanges(seqnames = seqnames_dens,
                              IRanges(start = HM_lyr_G1_1_dens$start,
                                      end = HM_lyr_G1_1_dens$end))
      return(heatmap_data)
    }
  }
  
  heatmap_data <- get_heatmap_data(side = side)
  
  adjust_chr <- function(chr_numbers){
    if(side=="lyr"){
      chr_numbers = chr_numbers + 8
      return(chr_numbers)
    }
    return(chr_numbers)
  }
  
  homeo_regions_G <- GRanges(seqnames = adjust_chr(as.numeric(
    str_replace(homeo_regions$seqnames, "chr", ""))),
                             IRanges(start = homeo_regions$start,
                                     end = homeo_regions$end))
  
  nonhomeo_regions_G <- GRanges(seqnames = adjust_chr(as.numeric(
    str_replace(nonhomeo_regions$start_scaffold, "chr", ""))),
    IRanges(start = nonhomeo_regions$start_coordinate,
            end = nonhomeo_regions$end_coordinate))
  
  coverage[coverage > 40] <- 40
  
  # Plot with karyoploter

  # First set parameters

  plot.params <- getDefaultPlotParams(plot.type=2)
  plot.params$data1height <- 50
  plot.params$data2height <- 50
  plot.params$ideogramheight <- 50

  kp <- plotKaryotype(genome = chrom_sizes_G, plot.type = 2, plot.params = plot.params, chromosomes = "all")
  kpDataBackground(kp, data.panel = 1, color = "grey")
  kpDataBackground(kp, data.panel = 2, color = "white")
  kpDataBackground(kp, data.panel = "ideogram", color = "white")
  
  # kpLines(kp, data.panel = 1, chr = c(1:8), 
  #         x = average_chromosomes, 
  #         y = coverage, 
  #         col = "black", 
  #         ymax = 20)
  # kpAxis(kp, data.panel = 1, 
  #        side = 2, 
  #        numticks = 2, 
  #        tick.pos = c(0, 1), 
  #        labels = c(0, 20), 
  #        cex = 0.5)
  # 
  kpHeatmap(kp, data = heatmap_data, y = coverage, colors = c("white", "black"))
  kpPlotRegions(kp, data = new_deg_output_G, 
                col = c("blue"), border = NA, 
                data.panel = "ideogram")
  kpPlotRegions(kp, data = homeo_regions_G,
                col = c("#CC79A7"), border = NA, 
                data.panel = 2)
  kpPlotRegions(kp, data = nonhomeo_regions_G, 
                col = c("#009E73"), border = NA, 
                data.panel = 2)
  
}

homeo_region_HM <- read.delim("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/homeo_HM_overlap_continuous.txt", sep = " ")
homeo_region_LL <- read.delim("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/homeo_LL_overlap_continuous.txt", sep = " ")
nonhomeo_region_HM <- read.delim("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/nonhomeo_HM_overlap_continuous.txt", sep = " ")
nonhomeo_region_LL <- read.delim("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis_v2/00_data/nonhomeo_LL_overlap_continuous.txt", sep = " ")


homeo_region_empty <- data.frame(seqnames = NULL,
                                 start = NULL,
                                 end = NULL)
nonhomeo_region_empty <- data.frame(start_scaffold = NULL,
                                    start_coordinate = NULL,
                                    end_coordinate = NULL)

HM_halG1_v_RS7G1_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/HM_halG1_v_RS7G1_DEG.txt")
HM_halG1_v_RS7G4_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/HM_halG1_v_RS7G4_DEG.txt")
HM_lyrG1_v_RS7G1_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyrG1_v_RS7G4_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/HM_lyrG1_v_RS7G4_DEG.txt")

DEGs_karyoploter(HM_halG1_v_RS7G1_DEG, 
                 "hal", 
                 HM_RS7K_G1_hal_avg_coverage,
                 chrom_length_hal,
                 homeo_region_empty,
                 nonhomeo_region_empty)
DEGs_karyoploter(HM_halG1_v_RS7G4_DEG, 
                 "hal",
                 HM_RS7_G4_hal_avg_coverage,
                 chrom_length_hal,
                 homeo_region_HM,
                 nonhomeo_region_HM)
DEGs_karyoploter(HM_lyrG1_v_RS7G1_DEG, 
                 "lyr",
                 HM_RS7K_G1_lyr_avg_coverage,
                 chrom_length_lyr,
                 homeo_region_empty,
                 nonhomeo_region_empty)
DEGs_karyoploter(HM_lyrG1_v_RS7G4_DEG, 
                 "lyr",
                 HM_RS7_G4_lyr_avg_coverage,
                 chrom_length_lyr,
                 homeo_region_HM,
                 nonhomeo_region_HM)

LL_halG1_v_RS7G1_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/LL_halG1_v_RS7G1_DEG.txt")
LL_halG1_v_RS7G4_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/LL_halG1_v_RS7G4_DEG.txt")
LL_lyrG1_v_RS7G1_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyrG1_v_RS7G4_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG_v2/LL_lyrG1_v_RS7G4_DEG.txt")

DEGs_karyoploter(LL_halG1_v_RS7G1_DEG, 
                 "hal", 
                 LL_RS7_G1_hal_avg_coverage,
                 chrom_length_hal,
                 homeo_region_empty,
                 nonhomeo_region_empty)
DEGs_karyoploter(LL_halG1_v_RS7G4_DEG, 
                 "hal", 
                 LL_RS7K_G4_hal_avg_coverage,
                 chrom_length_hal,
                 homeo_region_LL,
                 nonhomeo_region_LL)
DEGs_karyoploter(LL_lyrG1_v_RS7G1_DEG, 
                 "lyr",
                 LL_RS7_G1_1_lyr_real_coverage,
                 chrom_length_lyr,
                 homeo_region_empty,
                 nonhomeo_region_empty)
DEGs_karyoploter(LL_lyrG1_v_RS7G4_DEG, 
                 "lyr",
                 LL_RS7K_G4_lyr_avg_coverage,
                 chrom_length_lyr,
                 homeo_region_LL,
                 nonhomeo_region_LL)

HM_RS7G1_v_G4_hal_DEG <- read.delim("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG/HM_RS7G1_v_G4_hal_DEG.txt")

DEGs_karyoploter(HM_RS7G1_v_G4_hal_DEG, 
                 "hal", 
                 HM_RS7_G4_1_hal_real_coverage,
                 chrom_length_hal,
                 average_chromosomes_hal)

## We also plot the DMR density but we change
## coverage thresholds

# Import DMR density files

MILD_PvsS1_hal <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/HM_SYN1vPRO1_hal_dens.txt")
MILD_PvsS4_hal <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/HM_SYN4vPRO1_hal_dens.txt")
MILD_PvsS1_lyr <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/HM_SYN1vPRO1_lyr_dens.txt")
MILD_PvsS4_lyr <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/HM_SYN4vPRO1_lyr_dens.txt")

STRESS_PvsS1_hal <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/LL_SYN1vPRO1_hal_dens.txt")
STRESS_PvsS4_hal <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/LL_SYN4vPRO1_hal_dens.txt")
STRESS_PvsS1_lyr <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/LL_SYN1vPRO1_lyr_dens.txt")
STRESS_PvsS4_lyr <- fread("~/Dropbox/2022_Akam_Epigenetic_Stefan/Manuscript/Supplementary/Figures/newFigS2/LL_SYN4vPRO1_lyr_dens.txt")

# Redefine the function to threshold set at X

### We define a function to plot DMRs 
### on chromosomes with karyoploteR


DMRs_karyoploter <- function(deg_output, 
                             side, 
                             coverage, 
                             chrom_length,
                             homeo_regions,
                             nonhomeo_regions){
  
  coverage <- coverage$dens
  
  direction <- ifelse(deg_output$dir>0, "increase", "decrease")
  
  deg_output <- mutate(deg_output, direction = direction)
  
  # We import the annotations to add genomic coordinates
  # for each DEG based on their geneID
  
  annotation_import <- function(side){
    if(side=="hal"){
      anno <- read.delim("~/OneDrive/PhD/Project/bs_data/Ahal_genome_Dario/hal_MAKER_PASA_46264_noTEs.gff",
                         header=FALSE,
                         col.names = c("chromosome", 
                                       "tool", 
                                       "context", 
                                       "start", 
                                       "end", 
                                       "number",
                                       "strand", 
                                       "dot", 
                                       "extra"))
    }
    else {
      anno <- read.delim("~/OneDrive/PhD/Project/bs_data/Alyr_genome_Dario/lyr_r2_MAKER_28737_8chr_renamed.gff",
                         header=FALSE,
                         col.names = c("chromosome", 
                                       "tool", 
                                       "context", 
                                       "start", 
                                       "end", 
                                       "number",
                                       "strand", 
                                       "dot", 
                                       "extra"))
    }
    return(anno)
  }
  
  anno <- annotation_import(side)
  
  anno <- filter(anno, context=="gene")
  
  # we get the geneID from the extra column and add it to the annotation
  
  geneID_anno <- str_split_fixed( # get all the characters after the =
    as.vector(
      # get all the characters before the ;
      str_split_fixed(anno$extra, ";", 2)[,1]
    ), "=", 2
  )[,2]
  
  anno <- anno %>%
    mutate(geneID = geneID_anno)
  
  deg_coordinates <- anno %>%
    filter(geneID %in% deg_output$geneID)
  
  new_deg_output <- deg_output %>%
    filter(geneID %in% deg_coordinates$geneID) %>%
    mutate(chromosome = deg_coordinates$chromosome,
           start = deg_coordinates$start,
           end = deg_coordinates$end)
  
  
  # Define chromosome sizes
  
  import_chromSizes <- function(side){
    if(side=="hal"){
      chrom_sizes <- data.frame(Chrom = c(1:8),
                                Start = c(rep(1, 8)),
                                End = c(37868808, 
                                        19340014, 
                                        30847184,
                                        27260850, 
                                        25233939, 
                                        28041738,
                                        30089577, 
                                        26531562),
                                Name = rep(NA, 8),
                                Colors = rep("lightgrey", 
                                             8))
      
    }
    
    else{
      # Define chromosome sizes lyrata
      chrom_sizes <- data.frame(Chrom = c(9:16),
                                Start = c(rep(1, 8)),
                                End = c(28718688, 
                                        21180362, 
                                        27642951,
                                        24725123, 
                                        23209782, 
                                        25845826,
                                        27879834, 
                                        20443964),
                                Name = rep(NA, 8),
                                Colors = rep("lightgrey", 
                                             8))
    }
    return(chrom_sizes)
  }
  
  # Import chromosome sizes
  
  chrom_sizes <- import_chromSizes(side)
  
  # Define vector with chromosome names
  chrom_names <- chrom_sizes$End
  names(chrom_names) <- chrom_sizes$Chrom
  
  chrom_sizes_G <- GRanges(seqnames = chrom_sizes$Chrom,
                           ranges = IRanges(start = chrom_sizes$Start,
                                            end = chrom_sizes$End),
                           seqlengths = chrom_names)
  
  # Turn chromosome scaffolds into simple numbers
  
  new_deg_output$chromosome <- as.numeric(
    str_replace(new_deg_output$chromosome, "chr", ""))
  
  # Turn deg output into GRanges object
  
  new_deg_output_G <- GRanges(seqnames = new_deg_output$chromosome,
                              ranges = IRanges(start = new_deg_output$start,
                                               end = new_deg_output$end))
  new_deg_output_G <- GenomicRanges::reduce(new_deg_output_G)
  
  get_heatmap_data <- function(side){
    if(side=="hal"){
      seqnames_dens <- as.numeric(str_replace_all(HM_hal_G1_1_dens$seqnames, "chr", ""))
      heatmap_data <- GRanges(seqnames = seqnames_dens,
                              IRanges(start = HM_hal_G1_1_dens$start,
                                      end = HM_hal_G1_1_dens$end))
      return(heatmap_data)
    }
    else{
      seqnames_dens <- as.numeric(str_replace_all(HM_lyr_G1_1_dens$seqnames, "chr", ""))
      heatmap_data <- GRanges(seqnames = seqnames_dens,
                              IRanges(start = HM_lyr_G1_1_dens$start,
                                      end = HM_lyr_G1_1_dens$end))
      return(heatmap_data)
    }
  }
  
  heatmap_data <- get_heatmap_data(side = side)
  
  adjust_chr <- function(chr_numbers){
    if(side=="lyr"){
      chr_numbers = chr_numbers + 8
      return(chr_numbers)
    }
    return(chr_numbers)
  }
  
  homeo_regions_G <- GRanges(seqnames = adjust_chr(as.numeric(
    str_replace(homeo_regions$seqnames, "chr", ""))),
    IRanges(start = homeo_regions$start,
            end = homeo_regions$end))
  
  nonhomeo_regions_G <- GRanges(seqnames = adjust_chr(as.numeric(
    str_replace(nonhomeo_regions$start_scaffold, "chr", ""))),
    IRanges(start = nonhomeo_regions$start_coordinate,
            end = nonhomeo_regions$end_coordinate))
  
  coverage[coverage > 40] <- 40
  
  # Plot with karyoploter
  
  # First set parameters
  
  plot.params <- getDefaultPlotParams(plot.type=2)
  plot.params$data1height <- 50
  plot.params$data2height <- 50
  plot.params$ideogramheight <- 50
  
  kp <- plotKaryotype(genome = chrom_sizes_G, plot.type = 2, plot.params = plot.params)
  kpDataBackground(kp, data.panel = 1, color = "grey")
  kpDataBackground(kp, data.panel = 2, color = "white")
  kpDataBackground(kp, data.panel = "ideogram", color = "white")
  
  # kpLines(kp, data.panel = 1, chr = c(1:8), 
  #         x = average_chromosomes, 
  #         y = coverage, 
  #         col = "black", 
  #         ymax = 20)
  # kpAxis(kp, data.panel = 1, 
  #        side = 2, 
  #        numticks = 2, 
  #        tick.pos = c(0, 1), 
  #        labels = c(0, 20), 
  #        cex = 0.5)
  # 
  kpHeatmap(kp, data = heatmap_data, y = coverage, colors = c("white", "grey", "black"))
  kpPlotRegions(kp, data = new_deg_output_G, 
                col = c("blue"), border = NA, 
                data.panel = "ideogram")
  kpPlotRegions(kp, data = homeo_regions_G,
                col = c("#CC79A7"), border = NA, 
                data.panel = 2)
  kpPlotRegions(kp, data = nonhomeo_regions_G, 
                col = c("#009E73"), border = NA, 
                data.panel = 2)
  
}

# halleri side plots
DMRs_karyoploter(HM_halG1_v_RS7G1_DEG, 
                 "hal", 
                 MILD_PvsS1_hal,
                 chrom_length_hal,
                 homeo_region_empty,
                 nonhomeo_region_empty)

DMRs_karyoploter(HM_halG1_v_RS7G1_DEG, 
                 "hal", 
                 MILD_PvsS4_hal,
                 chrom_length_hal,
                 homeo_region_empty,
                 nonhomeo_region_empty)

DMRs_karyoploter(HM_halG1_v_RS7G1_DEG, 
                 "hal", 
                 STRESS_PvsS1_hal,
                 chrom_length_hal,
                 homeo_region_empty,
                 nonhomeo_region_empty)

DMRs_karyoploter(HM_halG1_v_RS7G1_DEG, 
                 "hal", 
                 STRESS_PvsS4_hal,
                 chrom_length_hal,
                 homeo_region_empty,
                 nonhomeo_region_empty)

#lyrata side plots

DMRs_karyoploter(HM_lyrG1_v_RS7G1_DEG, 
                 "lyr", 
                 MILD_PvsS1_lyr,
                 chrom_length_lyr,
                 homeo_region_empty,
                 nonhomeo_region_empty)

DMRs_karyoploter(HM_lyrG1_v_RS7G1_DEG, 
                 "lyr", 
                 MILD_PvsS4_lyr,
                 chrom_length_lyr,
                 homeo_region_empty,
                 nonhomeo_region_empty)

DMRs_karyoploter(HM_lyrG1_v_RS7G1_DEG, 
                 "lyr", 
                 STRESS_PvsS1_lyr,
                 chrom_length_lyr,
                 homeo_region_empty,
                 nonhomeo_region_empty)

DMRs_karyoploter(HM_lyrG1_v_RS7G1_DEG, 
                 "lyr", 
                 STRESS_PvsS4_lyr,
                 chrom_length_lyr,
                 homeo_region_empty,
                 nonhomeo_region_empty)
