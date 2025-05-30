# Supplementary results

To find which script is related to which supplementary result and the data required, use the index below to navigate the sections.

 - [Figure S2a](#figure-s2a): BS-seq coverage per chromosome compared to DMR density, DEG density and HE,non-HE regions
 - [Figure S2b](#figure-s2b): DNA-seq coverage per chromosome compared to DMR density and HE,non-HE regions
 - [Figure S3](#figure-s3): BS-seq coverage of natural ALK & TKS samples
 - [Figure S4](#figure-s4): MDS plots of the methylation state in all samples for different cytosine thresholds.
 - [Figure S5abc](#figure-s5abc): Differentially mehylated regions for:
     - P G1 vs S G1-G4
     - ALK G1 vs S G1-G4  
     - TKS G1 vs S G1-G4  
     - P G1 vs ALK G1-G4  
     - P G1 vs TKS G1-G5
 - [Figure S5d](#figure-s5d): Differentially methylated regions for natural ALK G1 vs TKS G1
 - [Figure S6](#figure-s6): Differentially methylated regions for G1 vs G4-5 for progenitors, natural ALK & TKS
 - [Figure S7](#figure-s7): Stacked barplot with proportion of DMRs contexts for P G1 vs S G1-G4
 - [Figure S8](#figure-s8): Density distribution of DMRs across chromosomes (synthetic G4 vs natural ALK G4)
 - [Figure S9](#figure-s9): Differentially expressed genes for P G1 vs S G1-G4
 - [Figure S10](#figure-s10): Chisqure tests for overlaps between DMGs & DEGs
 - [Figure S11](#figure-s11): Relationship and correlation between expression and methylation changes P G1 vs S G1
 - [Figure S12](#figure-s14): Overlap across comparisons for tetraploid-specific genes




## Figure S2a

Package requirements: `karyoploteR`, `tidyverse` and `data.table`  
Script: `FigS2a_chromosomeDEG.R`  
Data: density files for Fig.4c [link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results/data/Fig4c/density_files), files from `data/FigS2b` and results from `main_results/RNAseq_edgeR_fullReport.Rmd` (without path).  
Details: to run the script set up the path to the density files on `L10` and the mentioned paths in `L54` `L55` `L148` `L161` `L344-347`, `L357-360`, `L387-390`, `L417`, `L430-438`, `L464`, `L477`.

## Figure S2b

Package requirements: `karyoploteR`, `tidyverse` and `data.table`  
Script: `FigS2a_chromosomeDNACov.R`  
Data: density files for Fig.4c [link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results/data/Fig4c/density_files), files from `data/FigS2b` and results from `main_results/RNAseq_edgeR_fullReport.Rmd` (without path).  
Details: to run the script set up the path to the density files on `L11` and the mentioned paths in `L53` `L54` `L143` `L156` `L339-342`, `L352-355`, `L382-385`, `L412`, `L425-433`, `L459`, `L472`.

## Figure S3

Package requirements: `karyoploteR`, `tidyverse`, `data.table` and `ggplot2`
Script: `FigS3_BSseq_coverage_ALK_TKS.R`  
Data: density files for Fig.4c [link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results/data/Fig4c/density_files), files from `data/FigS2b` and results from `main_results/RNAseq_edgeR_fullReport.Rmd` (without path).  
Details: to run the script set up the path to the density files on `L15`.

## Figure S4

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`  
Script: [`Fig3a_MDSplot.R`](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results)  
Data: on Figshare following [this link](https://figshare.com/projects/Data_for_MDS_analyses/134765).    
Details: on L359, L390, L422, L453 change the `top` parameter to the number of cytosines you'd like to include for the MDS plot.

## Figure S5

### FigS5a

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`  
Script: `FigS3a_DMRBarplot.R`  
Data: `data/FigS5a`   
Details: To run the script set the right data inputs on `L13`, `L62`, `L155` and `L202`

### FigS5c

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`  
Script: `FigS5c_DMRBarplots_ALKvTKS.R`  
Data: `data/FigS5c`   
Details: To run the script set the right data inputs on `L12` and `L96`

## Figure S6

Package requirements: `tidyverse`, `data.table`, `patchwork` and `GenomicRanges`
Script:  `FigS6_DMRs_g4v1.R`
Data: from [ARPEGGIO analyses](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses), specifically:  
  - `MILD_alk1_v_alk4` & `STRESS_alk1_v_alk4`
  - `MILD_tks1_v_tks5` & `STRESS_tks1_v_tks5`
  - `MILD_pro1_v_pro4` & `STRESS_pro1_v_pro4`
  
Details: `L15`, `L45-46`, `L159`, `L189-90`, `L303`, `L443`, `L582`, `L586`, `L711`, `L715` need to be changed to the path of the corresponding datasets.

### FigS7

Figure done with `Excel`. 
Data: `data/FigS5a`   

## Figure S8

Package requirements: `tidyverse` and `patchwork`  
Script: `FigS8_ChromosomeSpatialDensityDMRs.R`  
Data: `data/FigS8`
Details: Change `L10` and `40-41` to the specified path.


## Figure S9

Package requirements: edgeR, data.table, ggplot2
Script: `main_results/RNAseq_edgeR_fullReport.Rmd`   
Data: expression data from `RNASeq`

## Figure S10

Package requirements: `tidyverse`, `patchwork`, `VennDiagram`, `ggVennDiagram`, `RColorBrewer`
Script: `FigS10_chisquare.R`  
Data: from [ARPEGGIO analyses](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses) and [low coverage regions](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results/data/Fig4c).  
Details: L18 has to be set for the `ARPEGGIO_analyses` folder. L184-185 and L187-188 have paths to be changed for files about low coverage regions `main_results/data/Fig4c`.  

## Figure S11

Package requirements: `data.table`, `UpSetR`, `ggplot2`, `patchwork`, `tidyverse`
Script: `FigS11_RawDE_vs_RawDM_correlation.R`  
Data: `data/FigS5a`, [low coverage regions](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results/data/Fig4c) and results from `main_results/RNAseq_edgeR_fullReport.Rmd`.  
Details: L13 has to be set for the `data/FigS5a`. L359 should point to the results from `main_results/RNAseq_edgeR_fullReport.Rmd`.