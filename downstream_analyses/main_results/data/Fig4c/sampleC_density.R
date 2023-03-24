### This script:
### 1) Converts genomic coordinates from Bismark cov.gz
### output to the genomic coordinates of our remapped genome
### assemblies *
### 2) Computes the cytosine coverage in 100kb windows
### 3) Outputs the coverage for each genomic bin


### * to do that, we use the agp file resulting
### from RagTag, describing the new position of all scaffolds
### from the old assembly.

### Import libraries 

library(data.table)
library(GenomicRanges)


### Command line arguments

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: cov gz file of interest
covgz_path <- comm_args[1]

# Second argument: progenitor side, hal or lyr
hal_or_lyr <- comm_args[2]

# Third argument: output name (with extension)
output_name <- comm_args[3]

### Import file to convert

covgz <- fread(covgz_path)
names(covgz) <- c("seqnames", "start", "end", "methylated", "unmethylated", "total")

print("All files imported, starting cytosine coverage computation...")

############################## PART 2: CYTOSINE COVERAGE ############################## 

# Define chromosome sizes for each progenitor

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

# Convert new covgz to GRanges to compute overlap

new_covgz_G <- GRanges(seqnames = covgz$seqnames,
                       ranges = IRanges(start = covgz$start,
                                        end = covgz$end),
                       score = covgz$total
                       )

dens <- GenomicRanges::countOverlaps(windows, new_covgz_G)

overlaps <- GenomicRanges::findOverlaps(windows, new_covgz_G)

window_scores <- rep(0, length(dens))

for (i in c(1:length(windows))){
  print(paste0(round((i / length(windows))*100, 1), "%"))
  window_scores[i] <- sum(covgz[overlaps[overlaps@from == i]@to,]$total)
}

output <- cbind(as.data.frame(windows), 
                dens = dens, 
                window_scores = window_scores,
                coverage = window_scores / dens)

### Write to output

# Filename

fwrite(output, file=output_name, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

print("Done")
