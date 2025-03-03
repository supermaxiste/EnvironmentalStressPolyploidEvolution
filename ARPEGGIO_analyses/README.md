# ARPEGGIO analyses

Here you will find 24 folders for all 24 comparisons ran with ARPEGGIO to compare DNA methylation patterns together with a summary on data quality. To navigate the folders use the following legend:

 - `MILD` and `STRESS` refer to the two conditions in which plants were grown

 - `1`, `4` and `5` represent the generation of a given group of plants

 - `pro` are the progenitor species _A. halleri_ and _A. lyrata_

 - `syn` are the synthetic _A. kamchatica_ created after crossing _A. halleri_ and _A. lyrata_

 - `alk` and `tks` are two natural lines of _A. kamchatica_. `alk` is short for Alaska and `tks` is short for Takashima. The Alaska line is useed to colder conditions while Takashima is used to warmer conditions.
 
In each comparison folder you will find the specific `config.yaml` and `metadata.txt` files used together with a `README` including links to the specific raw data used. All of the links can also be found below. All data is bundled in [BioProject PRJDB12567](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJDB12567).

## About `config.yaml` files

In short, parameters with paths require changes to the local system and for details about assemblies and annotations used, see the following paragraph.  

Conversion checks were carried out using the *Arabidopsis halleri* chloroplast genome (GenBank: KU764767.1, [NCBI link](https://www.ncbi.nlm.nih.gov/nuccore/KU764767.1/)) from [Assaf *et al.* (2017)](https://doi.org/10.1038/s41598-017-07891-5).  
Alignment was done used pre-processed assemblies:

  - For *Arabidopsis halleri* only the 8 main chromosomes were kept and renamed `chr1`-`chr8`
  - For *Arabidopsis lyrata* only the 8 main chromosomes were kept and renamed `chr9`-`chr16` (ARPEGGIO requirement to not have overlapping chromosome names)

Concerning the annotations of the above mentioned assemblies during the analyses for the paper:

  - For *Arabidopsis halleri* the annotation chromosomes were renamed `chr1`-`chr8`. Since the annotation was very recently made at the time of the analyses, the defult gene names were pretty long and not matching the official "released" version. The analyses are still reproducible but there might be a name mismatch between results presented in this repository and reproduced ones.
  - For *Arabidopsis lyrata* the annotation chromosomes were renamed `chr9`-`chr16`
 
### Cold conditions

_A. halleri_ samples:

 - `HM_hal_G1_1`: `DRX320070` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320070)
 - `HM_hal_G1_2`: `DRX320071` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320071)
 - `HM_hal_G1_3`: `DRX320072` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320072)
 - `HM_hal_G4_1`: `DRX320073` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320073)
 - `HM_hal_G4_2`: `DRX320074` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320074)
 - `HM_hal_G4_3`: `DRX320075` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320075)
 
_A. lyrata_ samples:
 
 - `HM_lyr_G1_1`: `DRX320081` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320081)
 - `HM_lyr_G1_2`: `DRX320082` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320082)
 - `HM_lyr_G1_3`: `DRX320083` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320083)
 - `HM_lyr_G4_1`: `DRX320084` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320084)
 - `HM_lyr_G4_2`: `DRX320085` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320085)
 - `HM_lyr_G4_3`: `DRX320086` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320086)

_A. kamchatica_ synthetics samples:

 - `HM_SYN_G1_1`: `DRX320093` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320093)
 - `HM_SYN_G1_2`: `DRX320094` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320094)
 - `HM_SYN_G1_3`: `DRX320095` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320095)
 - `HM_SYN_G4_1`: `DRX320096` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320096)
 - `HM_SYN_G4_2`: `DRX320097` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320097)
 - `HM_SYN_G4_3`: `DRX320098` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320098)

_A. kamchatica_ Alaska line samples:

 - `HM_ALK_G1_1`: `DRX320105` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320105)
 - `HM_ALK_G1_2`: `DRX320106` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320106)
 - `HM_ALK_G1_3`: `DRX320107` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320107)
 - `HM_ALK_G4_1`: `DRX320108` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320108)
 - `HM_ALK_G4_2`: `DRX320109` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320109)
 - `HM_ALK_G4_3`: `DRX320110` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320110)

_A. kamchatica_ Takashima line samples:

 - `HM_TKS_G1_1`: `DRX320117` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320117)
 - `HM_TKS_G1_2`: `DRX320118` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320118)
 - `HM_TKS_G1_3`: `DRX320119` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320119)
 - `HM_TKS_G5_1`: `DRX320120` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320120)
 - `HM_TKS_G5_2`: `DRX320121` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320121)
 - `HM_TKS_G5_3`: `DRX320122` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320122)


### Cold conditions

_A. halleri_ samples:

 - `LL_hal_G1_1`: `DRX320076`  [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320076)
 - `LL_hal_G1_2`: `DRX320077` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320077)
 - `LL_hal_G1_3`: `DRX320078` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320078)
 - `LL_hal_G4_1`: `DRX320079` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320079)
 - `LL_hal_G4_2`: `DRX320080` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320080)
 
_A. lyrata_ samples:
 
 - `LL_lyr_G1_1`: `DRX320087` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320087)
 - `LL_lyr_G1_2`: `DRX320088` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320088)
 - `LL_lyr_G1_3`: `DRX320089` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320089)
 - `LL_lyr_G4_1`: `DRX320090` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320090)
 - `LL_lyr_G4_2`: `DRX320091` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320091)
 - `LL_lyr_G4_3`: `DRX320092` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320092)

_A. kamchatica_ synthetics samples:

 - `LL_SYN_G1_1`: `DRX320099` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320099)
 - `LL_SYN_G1_2`: `DRX320100` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320100)
 - `LL_SYN_G1_3`: `DRX320101` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320101)
 - `LL_SYN_G4_1`: `DRX320102` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320102)
 - `LL_SYN_G4_2`: `DRX320103` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320103)
 - `LL_SYN_G4_3`: `DRX320104` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320104)

_A. kamchatica_ Alaska line samples:

 - `LL_ALK_G1_1`: `DRX320111` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320111)
 - `LL_ALK_G1_2`: `DRX320112` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320112)
 - `LL_ALK_G1_3`: `DRX320113` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320113)
 - `LL_ALK_G4_1`: `DRX320114` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320114)
 - `LL_ALK_G4_2`: `DRX320115` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320115)
 - `LL_ALK_G4_3`: `DRX320116` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320116)

_A. kamchatica_ Takashima line samples:

 - `LL_TKS_G1_1`: `DRX320123` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320123)
 - `LL_TKS_G1_2`: `DRX320124` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320124)
 - `LL_TKS_G1_3`: `DRX320125` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320125)
 - `LL_TKS_G5_1`: `DRX320126` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320126)
 - `LL_TKS_G5_2`: `DRX320127` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320127)
 - `LL_TKS_G5_3`: `DRX320128` [Link](https://www.ncbi.nlm.nih.gov/sra/?term=DRX320128)
