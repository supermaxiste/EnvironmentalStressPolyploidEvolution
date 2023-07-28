# Supplementary results

To find which script is related to which supplementary result and the data required, use the index below to navigate the sections.

 - [Figure S1](#figure-s1): Genome-wide coverage on both conditions and both progenitor's sides highlighting homoeologous exchange (HE) events.
 - [Figure S2-3](#figure-s2-3): MDS plots of the methylation state in all samples for different cytosine thresholds.
 - [Figure S4](#figure-s4): Barplots of the DMRs obtained from comparing synthetic *A. kamchatica* from the same generation across conditions. 
 - [Figure S5](#figure-s5): Barplots of the DMRs obtained from comparing natural polyploids and diploids across generations. First generations were used as a reference against fourth/fifth generations.
 - [Figure S6-7](#figure-s6-7): TBD
 - [Figure S8](#figure-s8): TBD
 - [Figure S9](#figure-s9): TBD
 - [Figure S10](#figure-s10): TBD
 - [Figure S11](#figure-s11): TBD
 - [Figure S12-13](#figure-s12-13): TBD
 - [Figure S14](#figure-s14): TBD
 - [Figure S15](#figure-s15): TBD
 - [Figure S16-22](#figure-s16-22): TBD


## Figure S1

Package requirements: `karyoploteR`, `tidyverse` and `patchwork`  
Script: `FigS1_HECoverage.R`  
Data: density files for Fig.4c [link](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results/data/Fig4c/density_files).  
Details: to run the script set up the path to the density files on `L14`.

## Figure S2-3

Package requirements: `tidyverse`, `data.table`, `patchwork` and `limma`  
Script: [`Fig3a_MDSplot.R`](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/downstream_analyses/main_results)  
Data: on Figshare following [this link](https://figshare.com/projects/Data_for_MDS_analyses/134765).    
Details: on L359, L390, L422, L453 change the `top` parameter to the number of cytosines you'd like to include for the MDS plot.

## Figure S4

Package requirements: `tidyverse`, `data.table` and `patchwork`
Script: `FigS4_DMRBarplots_conditions.R`   
Data: ARPEGGIO results for HM vs LL synthetics G1 and G4 (link TBD)  
Details: L12 and L120 need to be changed to the path of the corresponding datasets.

## Figure S5

Package requirements: `tidyverse`, `data.table`, `patchwork` and `GenomicRanges`
Script:  `FigS5_DMRs_g4v1.R`
Data: from [ARPEGGIO analyses](https://github.com/supermaxiste/EnvironmentalStressPolyploidEvolution/tree/main/ARPEGGIO_analyses), specifically:  
  - `MILD_alk1_v_alk4` & `STRESS_alk1_v_alk4`
  - `MILD_tks1_v_tks5` & `STRESS_tks1_v_tks5`
  - `MILD_pro1_v_pro4` & `STRESS_pro1_v_pro4`
  
Details: L15, L45-46, L159, L189-90, L303, L443, L582, L586, L711, L715 need to be changed to the path of the corresponding datasets.

## Figure S6-7

Package requirements: `edgeR`, `data.table` and `ggplot2`  
Script: `RNAseq_edgeR_fullReport.Rmd` 
Data: TBD  
Details: L41-42, L44-45 and L53 need to be changed to the path of the low coverage genes files (TBD) and count files (TBD) respectively.

## Figure S8

Package requirements:  
Script:  
Data:   
Details: 

## Figure S9

Package requirements:  
Script:  
Data:   
Details: 

## Figure S10

Package requirements:  
Script:  
Data:   
Details: 

## Figure S11

Package requirements:  
Script:  
Data:   
Details: 

## Figure S12-13

Package requirements:  
Script:  
Data:   
Details: 

## Figure S14

Package requirements:  
Script:  
Data:   
Details: 

## Figure S15

Package requirements:  
Script: 
Data:   
Details: 

## Figure S16-22

Package requirements: 
Script:   
Data:   
Details: 