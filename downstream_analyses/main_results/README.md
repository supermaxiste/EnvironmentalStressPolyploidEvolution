# Main results

To find which script is related to which result and the data required, use the index below to navigate the sections.

 - [Figure 1d](#figure-1d): circlize plots of genome-wide methylation levels
 - [Figure 2](#figure-2): global methylation levels, PCA plot and average methylation levels around gene bodies
 - [Figure 4](#figure-4): differentially methylated regions (DMRs) between synthetic and progenitors, genomic context of the DMRs and DMRs between natural species and progenitors/synthetics
 - [Figure 5](#figure-5): heatmap of highly variable genes and differentially expressed genes across several comparisons
 - [Figure 6](#figure-6): Venn diagram of differentially expressed and differentially methylated genes, and visualization of gene overrepresentation analyses


## Figure 1d

Package requirements: `circlize`, `GenomicRanges` and `ComplexHeatmap`  
Script: `Fig1d_CirclizeHeatmap.R`  
Data: `data/Fig1d`  
Details: the input data used for this figure is from data for Fig.3a [link](https://figshare.com/projects/Data_for_MDS_analyses/134765) and Bismark coverage files that can be found [here](https://doi.org/10.5281/zenodo.7323783). The pipeline to generate this figure is explained in detail in the `data/Fig1d` folder.

## Figure 2

### Figure 2a

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`  
Script: `Fig2a_MDSplot.R`  
Data: on Figshare following [this link](https://figshare.com/projects/Data_for_MDS_analyses/134765).   
Details: To find out how to obtain the data used, refer to the README in `data/Fig2a`

### Figure 2b

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig2b_DotplotGlobalMethylationLevel.R`  
Data: `data/Fig2b`  
Details: To find out how to obtain the data used, refer to the `METHImpute_analyses` folder from the main repository ([direct link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/METHImpute_analyses)).

### Figure 2c

Package requirements: `tidyverse`, `data.table`, `scales` and `patchwork`  
Script: `Fig2c_enrichmentPlot.R`  
Data: `data/Fig2c`  
Details: To find out how to obtain the data used, refer to the `METHImpute_analyses` folder from the main repository ([direct link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/METHImpute_analyses)).

## Figure 3

### Figure 3a

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig3a_DMRBarplots_conditions.R`.  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically synthetics generations  vs generation 4 (both conditions):
  - [MILD_syn1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_syn4)
  - [STRESS_syn1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_syn4)

### Figure 3b

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig3b_DMRBarplot.R`.  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically synthetics generations 1 & 4 vs progenitors from generation 1 (both conditions):
  - [MILD_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_pro1)
  - [MILD_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn4_v_pro1)
  - [STRESS_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_pro1)
  - [STRESS_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn4_v_pro1).  
Details: the script changes working directory 4 times, one for each of the outputs outlined above. Make sure to set the correct paths for lines 13, 62, 155 and 202

### Figure 3c

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig3c_DMRLineplot.R`.  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically alk & tks generation 1 vs progenitors generation 1 and synthetics generation 1 & 4 (both conditions):
  - [MILD_alk1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_alk1_v_pro1)
  - [MILD_alk1_v_syn1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_alk1_v_syn1)
  - [MILD_alk1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_alk1_v_syn4)
  - [MILD_tks1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_tks1_v_pro1)
  - [MILD_tks1_v_syn1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_tks1_v_syn1)
  - [MILD_tks1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_tks1_v_syn4)
  - [STRESS_alk1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_alk1_v_pro1)
  - [STRESS_alk1_v_syn1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_alk1_v_syn1)
  - [STRESS_alk1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_alk1_v_syn4)
  - [STRESS_tks1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_tks1_v_pro1)
  - [STRESS_tks1_v_syn1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_tks1_v_syn1)
  - [STRESS_tks1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_tks1_v_syn4)   
Details: the script changes working directory 12 times, one for each of the outputs outlined above and requires 4 paths to the data in `data/Fig4c`. Make sure to set the correct paths for lines 13, 33, 35, 99, 189, 329, 408, 497, 638, 658, 660, 725, 814, 954, 1033 and 1122.

## Figure 5

Package requirements: `grid`, `gridExtra`, `tidyverse`, `data.table`, `GenomicRanges` and `patchwork`  
Script: `Fig5_DonutPlot.R`.  
Data: 500bp promoter annotation files in `data/Fig5`, classical genome annotation files from TBD and the same `drmseq` output data for `Fig3b`:  
  - [MILD_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_pro1)
  - [MILD_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn4_v_pro1)
  - [STRESS_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_pro1)
  - [STRESS_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn4_v_pro1)   
Details: the script changes working directory 4 times, one for each of the outputs outlined above and requires 4 paths: 2 for annotations in `data/Fig3c` and 2 for the classical annotations. Make sure to set the correct paths for lines 24, 30, 208, 220, 723, 743, 777, 797.


### Figure 5bcd

Package requirements: `tidyverse` and `patchwork` 
Script: `Fig5bcd_DEGBarplots.R`  
Data: `data/Fig5a` and files with differentially expressed genes from TBD   
Details: the script requires to define paths in L11-14 for the folder with DEG files and data from `data/Fig5a`.

## Figure 6

### Figure 6a

Package requirements: `tidyverse`, `patchwork` , `VennDiagram`, `ggVennDiagram`, `RColorBrewer` and `svglite`
Script: `Fig8a_OverlapDEGDMG.R`  
Data: `data/Fig8a`
Details: L19 set path to `data/Fig8a`, L184 and L457 requires low coverage genes and expression data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`. L802 L1554 require setting the output folder for the output figures.

TBD

## Figure 7

Package requirements: `tidyverse`, `data.table`, `DESeq2` and `ComplexHeatmap`  
Script: `Fig7_ExpressionHeatmaps.R`  
Data: `data/Fig7`, classical genome annotation files from TBD and count tables from TBD   
Details: the script requires to define paths in L13-19 for the annotations, count tables and data from `data/Fig7`.

## Figure 8

### Figure 8a

Package requirements: `tidyverse`, `patchwork` , `VennDiagram`, `ggVennDiagram`, `RColorBrewer` and `svglite`  
Script: `Fig8a_OverlapDEGDMG.R`  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically synthetics generations 1 & 4 vs progenitors from generation 1 (both conditions):
  - [MILD_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_pro1)
  - [MILD_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn4_v_pro1)
  - [STRESS_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_pro1)
  - [STRESS_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn4_v_pro1).  
  Details: L19 set path to `data/Fig8a`, L184 and L457 requires low coverage genes and expression data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`. L802 L1554 require setting the output folder for the output figures. The whole script includes more comparisons than needed, you can exclude steps using the naturals `ALK` and `TKS`.

### Figure 8b

Package requirements: `tidyverse`, `patchwork` , `VennDiagram`, `ggVennDiagram`, `RColorBrewer` and `svglite`  
Script: `Fig8b_overlapDEGConditions.R`  
Data: expression data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`  
Details: L16 requires to set a path to DEGs and L18 a path to low coverage genes, all obtained from `RNAseq_edgeR_fullReport.Rmd`. L48 requires to set an output folder for the figures.  
