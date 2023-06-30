### This script will use the following inputs:
### 1) CX_report file
### 2) low coverage regions file
### 3) output file name

### and will produce one CX file with
### cytosines from regions with coverage
### over 2X

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: CX_report file
CX_file <- comm_args[1]

# Second argument: low cov scaffolds file
low_file <- comm_args[2]

# Fourth argument: output name (without extension)
output <- comm_args[3]

### Import libraries

library(data.table)
library(GenomicRanges)

### Read files

CX <- fread(CX_file)

low <- fread(low_file, header = F)

### Turn both files into GRanges
### to find overlaps

cytosines <- GRanges(seqnames = CX$V1, 
                     ranges = IRanges(start = CX$V2,
                                      end = CX$V2))

lowCov_regions <- GRanges(seqnames = low$V1,
                          ranges = IRanges(start = low$V2,
                                           end = low$V3))

overlaps <- findOverlaps(cytosines, lowCov_regions)@from

### Remove overlaps which represent
### cytosines in low coverage regions

CX_new <- CX[-overlaps, ]

### Write to output

fwrite(CX_new,
       file = paste0(output, "_noLow.txt"), 
       col.names = F, 
       row.names = F,
       quote = F, 
       sep = "\t")

