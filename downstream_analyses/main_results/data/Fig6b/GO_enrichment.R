#############
# Trying to automate GO enrichment on gene lists
# First argument of the script is a folder
###
commandArgs()
refList <- commandArgs()[6]
inputdir <- commandArgs()[7]
outdir <- commandArgs()[8]


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("topGO", "biomaRt", "Rgraphviz")){

  if (pkg %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(pkg, dependencies = T)
}
}

library(topGO)
library(devtools)
if (packageVersion("dbplyr") != "2.3.4") {
  devtools::install_version("dbplyr", version = "2.3.4", repos = "https://stat.ethz.ch/CRAN/")
}
library(biomaRt)
library(data.table)
library(usethis)


## Make the database of TAIR genes from ensembl

#collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'https://plants.ensembl.org')
# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c('ensembl_gene_id', 
                                       'go_id'), 
                        mart = mart)

#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]

# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))

## Run the analysis using "interesting" genes
all.genes <- names(geneID2GO)

# Reduce gene universe to include only genes that were expressed
avail.genes <- scan(refList, what = "character")
gene.universe <- geneID2GO[avail.genes]

# we create a function for writing gene sets for GO categories

writeGeneSets <- function(genesets, filename) {
  if (missing(filename)) {
    conn <- stdout()
  } else {
    conn <- file(filename, open="w")
    on.exit( close(conn) )
  }
  for (name in names(genesets)) {
    writeLines(paste(">", name), conn)
    writeLines(paste(genesets[[name]], collapse="\n"), conn)
  }
}

## Now let's make it a function
get_GO_enrichment <- function(geneFile){
  geneFile_path <- paste(inputdir,
                         geneFile,
                         sep = "/")
  int.genes_names <- scan(geneFile_path, what = "character")
  int.genes <- factor(as.integer(all.genes %in% int.genes_names), levels = c(0,1))

  names(int.genes) <- all.genes
  topGo.obj <- paste("go.obj", "BP", sep = ".")

 
  try({topGo.obj <- new("topGOdata", ontology="BP", 
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = gene.universe,
                   nodeSize = 10)
  
    show(topGo.obj)
    results.elim <- runTest(topGo.obj, 
                            algorithm = "elim", 
                            statistic = "fisher")
    pValue.elim <- score(results.elim)
    sigGOterms <- names(which(pValue.elim < 0.05))
    results.tab <- GenTable(object = topGo.obj, 
                            elimFisher = results.elim, 
                            orderBy = "elimFisher", 
                            ranksOf = "elimFisher", 
                            topNodes = length(sigGOterms),
                            numChar = 200)
    outTable <- paste(outdir, 
                      paste(tools::file_path_sans_ext(basename(geneFile)), 
                            "BP", 
                            "resultsTopGO.txt", 
                            sep = "_"), 
                      sep = "/")
    outGenes <- paste(outdir, 
                      paste(tools::file_path_sans_ext(basename(geneFile)), 
                            "BP", 
                            "resultsTopGO_genes.txt",
                            sep = "_"), 
                      sep = "/")
    allGO <- genesInTerm(topGo.obj)
    genes_per_category <- lapply(allGO,function(x) x[x %in% int.genes_names])
    genes_category_list <- sapply(sigGOterms, function(x) genes_per_category[[x]])
    writeGeneSets(genes_category_list, outGenes)
    write.table(results.tab, 
                file=outTable, 
                quote=F, 
                sep='\t', 
                row.names = F, 
                col.names = T)
    rm(results.tab)
    rm(topGo.obj)
    })
}

setwd(".")
fileList <- dir(path = inputdir, pattern = "*.txt")
lapply(fileList, get_GO_enrichment)

rm(list = ls(all.names = TRUE))
