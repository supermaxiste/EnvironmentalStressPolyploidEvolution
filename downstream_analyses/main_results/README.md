# Main results

To find which script is related to which result and the data required, use the index below to navigate the sections.

 - [Figure 2](#figure-1d): circlize plots of genome-wide methylation levels
 - [Figure 3](#figure-2): global methylation levels, PCA plot and average methylation levels around gene bodies
 - [Figure 4](#figure-4): differentially methylated regions (DMRs) between synthetic and progenitors, genomic context of the DMRs and DMRs between natural species and progenitors/synthetics
 - [Figure 5](#figure-5): heatmap of highly variable genes and differentially expressed genes across several comparisons
 - [Figure 6](#figure-6): Venn diagram of differentially expressed and differentially methylated genes, and visualization of gene overrepresentation analyses


## Figure 2

Package requirements: `circlize`, `GenomicRanges` and `ComplexHeatmap`  
Script: `Fig2_CirclizeHeatmap.R`  
Data: `data/Fig1d`  
Details: the input data used for this figure is from data for Fig.3a [link](https://figshare.com/projects/Data_for_MDS_analyses/134765) and Bismark coverage files that can be found [here](https://doi.org/10.5281/zenodo.7323783). The pipeline to generate this figure is explained in detail in the `data/Fig1d` folder.

## Figure 3

### Figure 3a

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`  
Script: `Fig3a_MDSplot.R`  
Data: on Figshare following [this link](https://figshare.com/projects/Data_for_MDS_analyses/134765).   
Details: To find out how to obtain the data used, refer to the README in `data/Fig2a`

### Figure 3b

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig3b_DotplotGlobalMethylationLevel.R`  
Data: `data/Fig3b`  
Details: To find out how to obtain the data used, refer to the `METHImpute_analyses` folder from the main repository ([direct link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/METHImpute_analyses)).

### Figure 3c

Package requirements: `tidyverse`, `data.table`, `scales` and `patchwork`  
Script: `Fig3c_enrichmentPlot.R`  
Data: `data/Fig3c`  
Details: To find out how to obtain the data used, refer to the `METHImpute_analyses` folder from the main repository ([direct link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/METHImpute_analyses)).

## Figure 4

### Figure 4a

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig4a_DMRBarplots_conditions.R`.  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically synthetics generations  vs generation 4 (both conditions):
  - [MILD_syn1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_syn4)
  - [STRESS_syn1_v_syn4](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_syn4)

### Figure 4b

Package requirements: `tidyverse`, `data.table`, `GenomicRanges` and `patchwork`  
Script: `Fig4b_DMRs_LinePlot_progenitorsVsyn14nat1.R`.  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically progenitors generation 1 vs synthetics generation 1 & 4 and alk & tks generation 1 (both conditions):
  - [MILD_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_pro1)
  - [MILD_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn4_v_pro1)
  - [MILD_alk1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_alk1_v_pro1)
  - [MILD_tks1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_tks1_v_pro1)
  - [STRESS_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_pro1)
  - [STRESS_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn4_v_pro1)
  - [STRESS_alk1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_alk1_v_pro1)
  - [STRESS_tks1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_tks1_v_pro1)  
Details: the script changes working directory 12 times, one for each of the outputs outlined above and requires 4 paths to the data in `data/Fig4c`. Make sure to set the correct paths for lines 19, 39, 40, 130, 234, 389, 409, 410, 500, 604, 761, 872, 976, 1131, 1242 and 1346.


### Figure 4c

Package requirements: `tidyverse`, `data.table` and `patchwork`  
Script: `Fig4c_DMRLineplot.R`.  
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

### Figure 5a

Package requirements: `edgeR`, `data.table` and `ggplot2`  
Script: `RNAseq_edgeR_fullReport.Rmd`  
Data: output counts from the [RNAseq analyses](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/RNAseq)  
Details: the easiest way to run this R Markdown file is to do so interactively. L17, L41-45 and L53 require to set paths to the required files.  

### Figure 5bcde

Package requirements: `tidyverse` and `patchwork` 
Script: `Fig5bcde_DEGBarplots.R`  
Data: expression and low coverage genes data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`    
Details: the script requires to define paths in L11-14 for low coverage genes and expression data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`.  


## Figure 6

Package requirements: `tidyverse`, `data.table`, `DESeq2` and `ComplexHeatmap`  
Script: `Fig6_ExpressionHeatmaps.R`  
Data: `data/Fig6`, including classical genome annotation files and count tables  
Details: the script requires to define paths in L13-19 for the annotations, count tables and data from `data/Fig6`.

## Figure 7

### Figure 7a

Package requirements: `tidyverse`, `patchwork` , `VennDiagram`, `ggVennDiagram`, `RColorBrewer` and `svglite`  
Script: `Fig7a_OverlapDEGDMG.R`  
Data: `dmrseq` output from the ARPEGGIO pipeline, specifically synthetics generations 1 & 4 vs progenitors from generation 1 (both conditions):
  - [MILD_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn1_v_pro1)
  - [MILD_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/MILD_syn4_v_pro1)
  - [STRESS_syn1_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn1_v_pro1)
  - [STRESS_syn4_v_pro1](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses/STRESS_syn4_v_pro1).  
  Details: L19 set path to `data/Fig7a`, L184 and L457 requires low coverage genes and expression data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`. L802 L1554 require setting the output folder for the output figures. The whole script includes more comparisons than needed, you can exclude steps using the naturals `ALK` and `TKS`.

### Figure7b

Package requirements: `tidyverse`, `patchwork` , `VennDiagram`, `ggVennDiagram`, `RColorBrewer` and `svglite`  
Script: `Fig7b_overlapDEGConditions.R`  
Data: expression data obtained from `downstream_analyses/supplementary_results/RNAseq_edgeR_fullReport.Rmd`  
Details: L16 requires to set a path to DEGs and L18 a path to low coverage genes, all obtained from `RNAseq_edgeR_fullReport.Rmd`. L48 requires to set an output folder for the figures.  
