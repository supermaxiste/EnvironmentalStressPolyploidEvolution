### This script provides a list of Arabidopsis thaliana
### gene IDs with a list of differentially methylated genes
### as input (output from ARPEGGIO) together with a two-way
### blast hit list between genes of a species against 
### A. thaliana genes

### import libraries

library(stringr)

### Set up command line arguments

args <- commandArgs(trailingOnly = TRUE)
dmr_path <- normalizePath(args[1])
blast_path <- normalizePath(args[2])
output_path <- normalizePath(args[3])
output_name <- args[4]

### Read files

dmr <- read.csv(dmr_path, sep="")
blast <- read.delim(blast_path)
names(blast) <- c("thalianaID", "geneID")

### match and pick corresponding thaliana ID

matches <- match(dmr$geneID, blast$geneID)
no_matches <- sum(is.na(matches))

thalianaID_output <- blast$thalianaID[na.omit(matches)]

print(paste0("Total genes: ", length(unique(dmr$geneID))))
print(paste0("Total matches: ", length(unique(thalianaID_output)), " (", round(((length(unique(thalianaID_output))) / length(unique(dmr$geneID)))*100, 1), "%)" ))

### Save thaliana gene IDs

setwd(output_path)

write.table(unique(thalianaID_output), file = output_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
