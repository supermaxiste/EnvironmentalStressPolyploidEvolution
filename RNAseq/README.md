# RNAseq analyses

In this folder you will find coordinates about the RNAseq data used (below), a `rna_env.yaml` file including all of the tools needed to run the full workflow together with shell scripts performing steps in the workflow.

### Setup and run the workflow

To modify scripts as little as possible and execute the analyses in an orderly fashion, we recommend creating the folders below and renaming files according to the names in the `Data coordinates` section below.

  - `00_data` where all the raw data (or symlinks) should be
  - `01_qualityCheck`
  - `02_alignment`
  - `03_read_sorting`
  - `04_count_tables`
  
The naming of the folders reflects the execution order of the scripts:

  1) `qualityCheck.sh` that quality checks raw data with `FastQC` (v0.11.8)
  2) `alignment_noclip.sh` aligns reads to the assemblies with `STAR` (v.2.7.3a)
  3) `read_sorting.sh` classifies polyploid reads with `EAGLE` (v1.1.3)
  4) `featureCounts.sh` creates count tables with featureCounts (part of `subread` v2.0.1)

Please check and modify the paths at the beginning of each script to make sure that input files/folders, output folders and resources used are set appropriately.

Before executing the scripts make sure to activate the `conda` environment created using `rna_env.yaml`. For `read_sorting.sh` script you will need to install [EAGLE](https://github.com/tony-kuo/eagle) (v1.1.3) and set a path to the `eagle` binary (the `conda` environment is not needed). Once all is set, run scripts as follows:

```
sh script.sh
```

### Data coordinates

All raw RNAseq data is bundled in [BioProject PRJDB12582](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJDB12582). Links and DRA accessions for individual files can be found below.

#### Cold conditions

_A. halleri_ samples:

 - `HM_hal_G1_1`: `DRX320737` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320737)
 - `HM_hal_G1_2`: `DRX320738` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320738)
 - `HM_hal_G1_3`: `DRX320739` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320739)
 - `HM_hal_G1_4`: `DRX320740` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320740)
 - `HM_hal_G4_1`: `DRX320741` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320741)
 - `HM_hal_G4_2`: `DRX320742` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320742)
 - `HM_hal_G4_3`: `DRX320743` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320743)
 
_A. lyrata_ samples:
 
 - `HM_lyr_G1_1`: `DRX320752` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320752)
 - `HM_lyr_G1_2`: `DRX320753` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320753)
 - `HM_lyr_G1_3`: `DRX320754` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320754)
 - `HM_lyr_G1_4`: `DRX320755` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320755)
 - `HM_lyr_G4_1`: `DRX320756` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320756)
 - `HM_lyr_G4_2`: `DRX320757` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320757)
 - `HM_lyr_G4_3`: `DRX320758` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320758)

_A. kamchatica_ synthetics samples:

 - `HM_SYN_G1_1`: `DRX320770` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320770)
 - `HM_SYN_G1_2`: `DRX320771` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320771)
 - `HM_SYN_G1_3`: `DRX320772` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320772)
 - `HM_SYN_G4_1`: `DRX320773` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320773)
 - `HM_SYN_G4_2`: `DRX320774` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320774)
 - `HM_SYN_G4_3`: `DRX320775` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320775)

_A. kamchatica_ Alaska line samples:

 - `HM_ALK_G1_1`: `DRX320782` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320782)
 - `HM_ALK_G1_2`: `DRX320783` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320783)
 - `HM_ALK_G1_3`: `DRX320784` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320784)
 - `HM_ALK_G4_1`: `DRX320785` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320785)
 - `HM_ALK_G4_2`: `DRX320786` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320786)
 - `HM_ALK_G4_3`: `DRX320787` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320787)

_A. kamchatica_ Takashima line samples:

 - `HM_TKS_G1_1`: `DRX320794` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320794)
 - `HM_TKS_G1_2`: `DRX320795` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320795)
 - `HM_TKS_G1_3`: `DRX320796` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320796)
 - `HM_TKS_G5_1`: `DRX320797` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320797)
 - `HM_TKS_G5_2`: `DRX320798` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320798)
 - `HM_TKS_G5_3`: `DRX320799` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320799)


#### Hot conditions

_A. halleri_ samples:

 - `LL_hal_G1_1`: `DRX320744` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320744)
 - `LL_hal_G1_2`: `DRX320745` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320745)
 - `LL_hal_G1_4`: `DRX320747` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320747) (this sample corresponds to `LL_hal_G1_1` in methylation data)
 - `LL_hal_G1_5`: `DRX320748` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320748) (this sample corresponds to `LL_hal_G1_2` in methylation data)
 - `LL_hal_G1_6`: `DRX320749` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320749) (this sample corresponds to `LL_hal_G1_3` in methylation data)
 - `LL_hal_G4_1`: `DRX320750` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320750)
 - `LL_hal_G4_2`: `DRX320751` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320751)
 
_A. lyrata_ samples:
 
 - `LL_lyr_G1_1`: `DRX320759` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320759)
  - `LL_lyr_G1_2`: `DRX320760` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320760)
  - `LL_lyr_G1_4`: `DRX320762` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320762) (this sample corresponds to `LL_lyr_G1_1` in methylation data)
  - `LL_lyr_G1_5`: `DRX320763` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320763) (this sample corresponds to `LL_lyr_G1_2` in methylation data)
  - `LL_lyr_G1_6`: `DRX320764` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320764) (this sample corresponds to `LL_lyr_G1_3` in methylation data)
  - `LL_lyr_G4_1`: `DRX320765` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320765) (this sample corresponds to `LL_lyr_G4_1` in methylation data)
  - `LL_lyr_G4_2`: `DRX320766` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320766) (this sample corresponds to `LL_lyr_G4_2` in methylation data)
  - `LL_lyr_G4_3`: `DRX320767` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320767)
  - `LL_lyr_G4_4`: `DRX320768` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320768)
  - `LL_lyr_G4_5`: `DRX320769` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320769) (this sample corresponds to `LL_lyr_G4_3` in methylation data)

_A. kamchatica_ synthetics samples:

 - `LL_SYN_G1_1`: `DRX320776` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320776)
 - `LL_SYN_G1_2`: `DRX320777` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320777)
 - `LL_SYN_G1_3`: `DRX320778` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320778)
 - `LL_SYN_G4_1`: `DRX320779` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320779)
 - `LL_SYN_G4_2`: `DRX320780` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320780)
 - `LL_SYN_G4_3`: `DRX320781` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320781)

_A. kamchatica_ Alaska line samples:

 - `LL_ALK_G1_1`: `DRX320788` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320788)
 - `LL_ALK_G1_2`: `DRX320789` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320789)
 - `LL_ALK_G1_3`: `DRX320790` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320790)
 - `LL_ALK_G4_1`: `DRX320791` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320791)
 - `LL_ALK_G4_2`: `DRX320792` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320792)
 - `LL_ALK_G4_3`: `DRX320793` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320793)

_A. kamchatica_ Takashima line samples:

 - `LL_TKS_G1_1`: `DRX320800` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320800)
 - `LL_TKS_G1_2`: `DRX320801` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320801)
 - `LL_TKS_G1_3`: `DRX320802` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320802)
 - `LL_TKS_G5_1`: `DRX320803` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320803)
 - `LL_TKS_G5_2`: `DRX320804` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320804)
 - `LL_TKS_G5_3`: `DRX320805` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320805)
