# load GENESPACE library
library(GENESPACE) #GENESPACE v1.3.1: synteny and orthology constrained comparative genomics
library(svglite)
library(ggplot2)

# Get the working directory from which the script is called:
script_wd <- getwd()

# Parse command line arguments:
args <- commandArgs(trailingOnly = TRUE)
wd <- paste0(script_wd,args[1])
path2mcscanx <- paste0(script_wd,args[2])

# Init genespace: 
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  ploidy = c(1,1,1,1),
  genomeIDs = c("A.halleri","A.lyrata","A.lyrata_jgi","A.thaliana")
  )
# Run genespace: 
gpar <- run_genespace(gsParam = gpar) 

# NOTE: 8. Constructing syntenic pan-gene set failed in our run due to "Error: vector memory exhausted (limit reached?)"



