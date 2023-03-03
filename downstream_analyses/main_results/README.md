# Main results

To find which script is related to which result and the data required, use the index below to navigate the sections.

 - [Figure 2](#figure-2): circlize plots of genome-wide methylation levels
 - [Figure 3](#figure-3): global methylation levels, PCA plot and average methylation levels around gene bodies
 - [Figure 4](#figure-4): differentially methylated regions (DMRs) between synthetic and progenitors, genomic context of the DMRs and DMRs between natural species and progenitors/synthetics
 - [Figure 5](#figure-5): heatmap of highly variable genes and differentially expressed genes across several comparisons
 - [Figure 6](#figure-6): Venn diagram of differentially expressed and differentially methylated genes, and visualization of gene overrepresentation analyses


## Figure 2

Package requirements: `circlize`, `GenomicRanges` and `ComplexHeatmap`
Script: `Fig2_CirclizeHeatmap.R`
Data: `data/Fig2`
Details: the input data used for this figure is from data for Fig.3a ( link)[https://figshare.com/projects/Data_for_MDS_analyses/134765] and Bismark coverage files that can be found (here)[https://doi.org/10.5281/zenodo.7323783]. The pipeline to generate this figure is explained in detail in the `data/Fig2` folder.

## Figure 3

### Figure 3a

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`
Script: `Fig3a_MDSplot.R`
Data: on Figshare following (this link)[https://figshare.com/projects/Data_for_MDS_analyses/134765].
Details: To find out how to obtain the data used, refer to the README in `data/Fig3a`

### Figure 3b

Package requirements: `tidyverse`, `data.table` and `patchwork`
Script: `Fig3b_DotplotGlobalMethylationLevel.R`
Data: `data/Fig3b`
Details: To find out how to obtain the data used, refer to the `METHImpute_analyses` folder from the main repository ((direct link)[https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/METHImpute_analyses]).

## Figure 4

TBD

## Figure 5

TBD

## Figure 6

TBD
