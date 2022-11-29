## MDS analyses

The starting point of MDS analyses are `cov` files from Bismark which can be obtained from ARPEGGIO runs. If these files are not available, we provide the final files used for MDS plotting (after performing all of the steps outlined below). To perform the analyses for all progenitors' sides and conditions, there are a total of 95 files. Cold conditions has 48 files and includes the following samples:

 - _A. halleri_ generation 1 (3 replicates)
 - _A. halleri_ generation 4 (3 replicates)
 - _A. lyrata_ generation 1 (3 replicates)
 - _A. lyrata_ generation 4 (3 replicates)
 - _A. kamchatica_ synthetic generation 1 (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ synthetic generation 4 (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 1, Alaska line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 4, Alaska line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 1, Takashima line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 4, Takashima line (3 replicates for each side: 6 replicates total)
 
 Hot conditions has 47 files and includes the following samples:

 - _A. halleri_ generation 1 (3 replicates)
 - _A. halleri_ generation 4 (2 replicates)
 - _A. lyrata_ generation 1 (3 replicates)
 - _A. lyrata_ generation 4 (3 replicates)
 - _A. kamchatica_ synthetic generation 1 (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ synthetic generation 4 (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 1, Alaska line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 4, Alaska line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 1, Takashima line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 4, Takashima line (3 replicates for each side: 6 replicates total)

 The following steps provide a guide to analyse the data:

  1) After gathering all `cov` files, in order to merge them, all need to be sorted. Assuming that we're working with a file called `filename_bismark.cov.gz`, the command to sort it is:

  ```
  zcat filename_bismark.cov.gz | sort -k1,1 -k2,2n | pigz > filename_bismark.cov.gz
  ```

  This command should be run for all samples and for _A. kamchatica_ samples make sure to name the files differently for each progenitors' side, for example `filename_hal_bismark.cov.gz` and `filename_lyr_bismark.cov.gz`

  2) Next we reformat the columns from the `cov` files by merging the first there columns with the spatial coordinates (providing unique names to each position) and summing the last two columns to obtain the total count for both methylated and unmethylated cytosines (coverage). We do this as follows:

  ```
  data=path/to/*cov.gz
  output=output/path/

  for filename in $data; do
          base=$(basename $filename _bismark.cov.gz)
          zcat $filename | awk '{print $1"_"$2"_"$3"\t"$4"\t"$5+$6}' | sort > output/path/${base}_minimal.cov
          pigz output/path/${base}_minimal.cov
  done
  ```
  Make sure to modify `data` and `output` paths.

  3) Files can now be merged based on their overlapping cytosines. For each sample, two values will be kept: the methylation proportion and the total number of cytosines. Assuming we have files from a samples named `filename_1_lyr_minimal.cov`, `filename_2_lyr_minimal.cov` and `filename_3_lyr_minimal.cov`, we merge them as follows:

  ```
  join -j 1 -o 1.1,1.2,1.3,2.2,2.3 filename_1_lyr_minimal.cov filename_2_lyr_minimal.cov > filename_12_lyr_minimal.cov
  join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 filename_12_lyr_minimal.cov filename_3_lyr_minimal.cov > filename_123_lyr_minimal.cov
  ```
  The structure of the output file `filename_123_lyr_minimal.cov` will be as follows:

   - Column 1: cytosine coordinates
   - Column 2-3: methylation proportion and total cytosines for replicate 1
   - Column 4-5: methylation proportion and total cytosines for replicate 2
   - Column 6-7: methylation proportion and total cytosines for replicate 3


   4) Merged files can now be merged further in four final files (2 conditions x 2 progenitors' sides). We provide our own scripts below for both progenitors' sides and both conditions (HM = cold, LL = hot) since the number of replicates is different:

```
data=../path/to/data/

# HM synthetics

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_RS7_G1_1hal_minimal.cov ${data}HM_RS7_G1_2hal_minimal.cov > HM_RS7_G1_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_RS7_G1_12hal_minimal.cov ${data}HM_RS7_G1_3hal_minimal.cov > HM_RS7_G1_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_RS7_G1_1lyr_minimal.cov ${data}HM_RS7_G1_2lyr_minimal.cov > HM_RS7_G1_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_RS7_G1_12lyr_minimal.cov ${data}HM_RS7_G1_3lyr_minimal.cov > HM_RS7_G1_123lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_RS7_G4_1hal_minimal.cov ${data}HM_RS7_G4_2hal_minimal.cov > HM_RS7_G4_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_RS7_G4_12hal_minimal.cov ${data}HM_RS7_G4_3hal_minimal.cov > HM_RS7_G4_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_RS7_G4_1lyr_minimal.cov ${data}HM_RS7_G4_2lyr_minimal.cov > HM_RS7_G4_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_RS7_G4_12lyr_minimal.cov ${data}HM_RS7_G4_3lyr_minimal.cov > HM_RS7_G4_123lyr_minimal.cov

# HM natural ALK

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_ALK_G1_1hal_minimal.cov ${data}HM_ALK_G1_2hal_minimal.cov > HM_ALK_G1_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_ALK_G1_12hal_minimal.cov ${data}HM_ALK_G1_3hal_minimal.cov > HM_ALK_G1_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_ALK_G1_1lyr_minimal.cov ${data}HM_ALK_G1_2lyr_minimal.cov > HM_ALK_G1_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_ALK_G1_12lyr_minimal.cov ${data}HM_ALK_G1_3lyr_minimal.cov > HM_ALK_G1_123lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_ALK_G4_1hal_minimal.cov ${data}HM_ALK_G4_2hal_minimal.cov > HM_ALK_G4_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_ALK_G4_12hal_minimal.cov ${data}HM_ALK_G4_3hal_minimal.cov > HM_ALK_G4_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_ALK_G4_1lyr_minimal.cov ${data}HM_ALK_G4_2lyr_minimal.cov > HM_ALK_G4_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_ALK_G4_12lyr_minimal.cov ${data}HM_ALK_G4_3lyr_minimal.cov > HM_ALK_G4_123lyr_minimal.cov

# HM natural TKS

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_TKS_G1_1hal_minimal.cov ${data}HM_TKS_G1_2hal_minimal.cov > HM_TKS_G1_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_TKS_G1_12hal_minimal.cov ${data}HM_TKS_G1_3hal_minimal.cov > HM_TKS_G1_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_TKS_G1_1lyr_minimal.cov ${data}HM_TKS_G1_2lyr_minimal.cov > HM_TKS_G1_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_TKS_G1_12lyr_minimal.cov ${data}HM_TKS_G1_3lyr_minimal.cov > HM_TKS_G1_123lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_TKS_G5_1hal_minimal.cov ${data}HM_TKS_G5_2hal_minimal.cov > HM_TKS_G5_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_TKS_G5_12hal_minimal.cov ${data}HM_TKS_G5_3hal_minimal.cov > HM_TKS_G5_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_TKS_G5_1lyr_minimal.cov ${data}HM_TKS_G5_2lyr_minimal.cov > HM_TKS_G5_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_TKS_G5_12lyr_minimal.cov ${data}HM_TKS_G5_3lyr_minimal.cov > HM_TKS_G5_123lyr_minimal.cov

# HM progenitors

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_hal_G1_1_minimal.cov ${data}HM_hal_G1_2_minimal.cov > HM_hal_G1_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_hal_G1_12_minimal.cov ${data}HM_hal_G1_3_minimal.cov > HM_hal_G1_123_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_lyr_G1_1_minimal.cov ${data}HM_lyr_G1_2_minimal.cov > HM_lyr_G1_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_lyr_G1_12_minimal.cov ${data}HM_lyr_G1_3_minimal.cov > HM_lyr_G1_123_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_hal_G4_1_minimal.cov ${data}HM_hal_G4_2_minimal.cov > HM_hal_G4_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_hal_G4_12_minimal.cov ${data}HM_hal_G4_3_minimal.cov > HM_hal_G4_123_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}HM_lyr_G4_1_minimal.cov ${data}HM_lyr_G4_2_minimal.cov > HM_lyr_G4_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 HM_lyr_G4_12_minimal.cov ${data}HM_lyr_G4_3_minimal.cov > HM_lyr_G4_123_minimal.cov

# LL synthetics

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_RS7_G1_1hal_minimal.cov ${data}LL_RS7_G1_2hal_minimal.cov > LL_RS7_G1_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_RS7_G1_12hal_minimal.cov ${data}LL_RS7_G1_3hal_minimal.cov > LL_RS7_G1_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_RS7_G1_1lyr_minimal.cov ${data}LL_RS7_G1_2lyr_minimal.cov > LL_RS7_G1_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_RS7_G1_12lyr_minimal.cov ${data}LL_RS7_G1_3lyr_minimal.cov > LL_RS7_G1_123lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_RS7_G4_1hal_minimal.cov ${data}LL_RS7_G4_2hal_minimal.cov > LL_RS7_G4_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_RS7_G4_12hal_minimal.cov ${data}LL_RS7_G4_3hal_minimal.cov > LL_RS7_G4_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_RS7_G4_1lyr_minimal.cov ${data}LL_RS7_G4_2lyr_minimal.cov > LL_RS7_G4_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_RS7_G4_12lyr_minimal.cov ${data}LL_RS7_G4_3lyr_minimal.cov > LL_RS7_G4_123lyr_minimal.cov

# LL natural ALK

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_ALK_G1_1hal_minimal.cov ${data}LL_ALK_G1_2hal_minimal.cov > LL_ALK_G1_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_ALK_G1_12hal_minimal.cov ${data}LL_ALK_G1_3hal_minimal.cov > LL_ALK_G1_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_ALK_G1_1lyr_minimal.cov ${data}LL_ALK_G1_2lyr_minimal.cov > LL_ALK_G1_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_ALK_G1_12lyr_minimal.cov ${data}LL_ALK_G1_3lyr_minimal.cov > LL_ALK_G1_123lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_ALK_G4_1hal_minimal.cov ${data}LL_ALK_G4_2hal_minimal.cov > LL_ALK_G4_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_ALK_G4_12hal_minimal.cov ${data}LL_ALK_G4_3hal_minimal.cov > LL_ALK_G4_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_ALK_G4_1lyr_minimal.cov ${data}LL_ALK_G4_2lyr_minimal.cov > LL_ALK_G4_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_ALK_G4_12lyr_minimal.cov ${data}LL_ALK_G4_3lyr_minimal.cov > LL_ALK_G4_123lyr_minimal.cov

# LL natural TKS

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_TKS_G1_1hal_minimal.cov ${data}LL_TKS_G1_2hal_minimal.cov > LL_TKS_G1_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_TKS_G1_12hal_minimal.cov ${data}LL_TKS_G1_3hal_minimal.cov > LL_TKS_G1_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_TKS_G1_1lyr_minimal.cov ${data}LL_TKS_G1_2lyr_minimal.cov > LL_TKS_G1_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_TKS_G1_12lyr_minimal.cov ${data}LL_TKS_G1_3lyr_minimal.cov > LL_TKS_G1_123lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_TKS_G5_1hal_minimal.cov ${data}LL_TKS_G5_2hal_minimal.cov > LL_TKS_G5_12hal_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_TKS_G5_12hal_minimal.cov ${data}LL_TKS_G5_3hal_minimal.cov > LL_TKS_G5_123hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_TKS_G5_1lyr_minimal.cov ${data}LL_TKS_G5_2lyr_minimal.cov > LL_TKS_G5_12lyr_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_TKS_G5_12lyr_minimal.cov ${data}LL_TKS_G5_3lyr_minimal.cov > LL_TKS_G5_123lyr_minimal.cov

# LL progenitors

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_hal_G1_1_minimal.cov ${data}LL_hal_G1_2_minimal.cov > LL_hal_G1_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_hal_G1_12_minimal.cov ${data}LL_hal_G1_3_minimal.cov > LL_hal_G1_123_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_lyr_G1_1_minimal.cov ${data}LL_lyr_G1_2_minimal.cov > LL_lyr_G1_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_lyr_G1_12_minimal.cov ${data}LL_lyr_G1_3_minimal.cov > LL_lyr_G1_123_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_hal_G4_1_minimal.cov ${data}LL_hal_G4_2_minimal.cov > LL_hal_G4_12_minimal.cov

join -j 1 -o 1.1,1.2,1.3,2.2,2.3 ${data}LL_lyr_G4_1_minimal.cov ${data}LL_lyr_G4_2_minimal.cov > LL_lyr_G4_12_minimal.cov
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 LL_lyr_G4_12_minimal.cov ${data}LL_lyr_G4_3_minimal.cov > LL_lyr_G4_123_minimal.cov

### Final merging

# halleri side, HM conditions

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_RS7_G1_123hal_minimal.cov ${data}HM_RS7_G4_123hal_minimal.cov > ${data}HM_synthe_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_ALK_G1_123hal_minimal.cov ${data}HM_ALK_G4_123hal_minimal.cov > ${data}HM_ALK_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_TKS_G1_123hal_minimal.cov ${data}HM_TKS_G5_123hal_minimal.cov > ${data}HM_TKS_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_hal_G1_123_minimal.cov ${data}HM_hal_G4_123_minimal.cov > ${data}HM_pro_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}HM_pro_allgen_hal_minimal.cov ${data}HM_synthe_allgen_hal_minimal.cov > ${data}HM_prosyn_hal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}HM_ALK_allgen_hal_minimal.cov ${data}HM_TKS_allgen_hal_minimal.cov > ${data}HM_nat_hal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25 ${data}HM_prosyn_hal.cov ${data}HM_nat_hal.cov > ${data}HM_allhal.cov

#Final ordering in the file:
#halleri_G1+G4	HM_RS7K_G1+G4	HM_ALK_G1+G4	HM_TKS_G1+5

# halleri side, LL conditions

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_RS7_G1_123hal_minimal.cov ${data}LL_RS7_G4_123hal_minimal.cov > ${data}LL_synthe_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_ALK_G1_123hal_minimal.cov ${data}LL_ALK_G4_123hal_minimal.cov > ${data}LL_ALK_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_TKS_G1_123hal_minimal.cov ${data}LL_TKS_G5_123hal_minimal.cov > ${data}LL_TKS_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5 LL_hal_G1_123_minimal.cov ${data}LL_hal_G4_12_minimal.cov > ${data}LL_pro_allgen_hal_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}LL_pro_allgen_hal_minimal.cov ${data}LL_synthe_allgen_hal_minimal.cov > ${data}LL_prosyn_hal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}LL_ALK_allgen_hal_minimal.cov ${data}LL_TKS_allgen_hal_minimal.cov > ${data}LL_nat_hal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25 ${data}LL_prosyn_hal.cov ${data}LL_nat_hal.cov > ${data}LL_allhal.cov

#Final ordering in the file:
#LL_halleri_G1+G4  LL_RS7K_G1+G4   LL_ALK_G1+G4    LL_TKS_G1+5

# lyrata side, HM conditions

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_RS7_G1_123lyr_minimal.cov ${data}HM_RS7_G4_123lyr_minimal.cov > ${data}HM_synthe_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_ALK_G1_123lyr_minimal.cov ${data}HM_ALK_G4_123lyr_minimal.cov > ${data}HM_ALK_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_TKS_G1_123lyr_minimal.cov ${data}HM_TKS_G5_123lyr_minimal.cov > ${data}HM_TKS_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}HM_lyr_G1_123_minimal.cov ${data}HM_lyr_G4_123_minimal.cov > ${data}HM_pro_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}HM_pro_allgen_lyr_minimal.cov ${data}HM_synthe_allgen_lyr_minimal.cov > ${data}HM_prosyn_lyr.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}HM_ALK_allgen_lyr_minimal.cov ${data}HM_TKS_allgen_lyr_minimal.cov > ${data}HM_nat_lyr.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25 ${data}HM_prosyn_lyr.cov ${data}HM_nat_lyr.cov > ${data}HM_alllyr.cov

#Final ordering in the file:
#lyrata_G1+G4  HM_RS7K_G1+G4   HM_ALK_G1+G4    HM_TKS_G1+5

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_RS7_G1_123lyr_minimal.cov ${data}LL_RS7_G4_123lyr_minimal.cov > ${data}LL_synthe_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_ALK_G1_123lyr_minimal.cov ${data}LL_ALK_G4_123lyr_minimal.cov > ${data}LL_ALK_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_TKS_G1_123lyr_minimal.cov ${data}LL_TKS_G5_123lyr_minimal.cov > ${data}LL_TKS_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 ${data}LL_lyr_G1_123_minimal.cov ${data}LL_lyr_G4_123_minimal.cov > ${data}LL_pro_allgen_lyr_minimal.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}LL_pro_allgen_lyr_minimal.cov ${data}LL_synthe_allgen_lyr_minimal.cov > ${data}LL_prosyn_lyr.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 ${data}LL_ALK_allgen_lyr_minimal.cov ${data}LL_TKS_allgen_lyr_minimal.cov > ${data}LL_nat_lyr.cov

join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25 ${data}LL_prosyn_lyr.cov ${data}LL_nat_lyr.cov > ${data}LL_alllyr.cov

#Final ordering in the file
#lyrata_G1+G4	HM_RS7_G1+G4	HM_ALK_G1+G4	HM_TKS_G1+5
```

   It's important to keep in mind the way files are merged to know which column represents which sample. In all of our files the order is the following. For cold conditions (both hal and lyr sides) and hot conditions (only lyr side):

  - Column 1: cytosine coordinates
  - Columns 2-7: Progenitor generation 1 (3 replicates)
  - Columns 8-13: Progenitor generation 4 (3 replicates)
  - Columns 14-19: _A. kamchatica_ synthetic generation 1 (3 replicates)
  - Columns 20-25: _A. kamchatica_ synthetic generation 4 (3 replicates)
  - Columns 26-31: _A. kamchatica_ natural generation 1, Alaska line (3 replicates)
  - Columns 32-37: _A. kamchatica_ natural generation 4, Alaska line (3 replicates)
  - Columns 38-43: _A. kamchatica_ natural generation 1, Takashima line (3 replicates)
  - Columns 44-49: _A. kamchatica_ natural generation 4, Takashima line (3 replicates)

For hot conditions (hal side):

   - Column 1: cytosine coordinates
   - Columns 2-7: Progenitor generation 1 (3 replicates)
   - Columns 8-11: Progenitor generation 4 (2 replicates)
   - Columns 12-17: _A. kamchatica_ synthetic generation 1 (3 replicates)
   - Columns 18-23: _A. kamchatica_ synthetic generation 4 (3 replicates)
   - Columns 24-29: _A. kamchatica_ natural generation 1, Alaska line (3 replicates)
   - Columns 30-35: _A. kamchatica_ natural generation 4, Alaska line (3 replicates)
   - Columns 36-41: _A. kamchatica_ natural generation 1, Takashima line (3 replicates)
   - Columns 42-47: _A. kamchatica_ natural generation 4, Takashima line (3 replicates)


  5) The last step is to filter all cytosines that have >=3 coverage. We use the same output files from the previous step:

```
data=../path/to/data/

#HM

awk -F " " '{ if(($3 >= 3) && ($5 >= 3) && ($7 >= 3) && ($9 >= 3) && ($11 >= 3) && ($13 >= 3) && ($15 >= 3) && ($17 >= 3) && ($19 >= 3) && ($21 >= 3) && ($23 >= 3) && ($25 >= 3) && ($27 >= 3) && ($29 >= 3) && ($31 >= 3) && ($33 >= 3) && ($35 >= 3) && ($37 >= 3) && ($39 >= 3) && ($41 >= 3) && ($43 >= 3) && ($45 >= 3) && ($47 >= 3) && ($49 >= 3)) { print } }' ${data}HM_allhal.cov > HM_allhal_filtered.cov

echo "halleri done (HM)"

awk -F " " '{ if(($3 >= 3) && ($5 >= 3) && ($7 >= 3) && ($9 >= 3) && ($11 >= 3) && ($13 >= 3) && ($15 >= 3) && ($17 >= 3) && ($19 >= 3) && ($21 >= 3) && ($23 >= 3) && ($25 >= 3) && ($27 >= 3) && ($29 >= 3) && ($31 >= 3) && ($33 >= 3) && ($35 >= 3) && ($37 >= 3) && ($39 >= 3) && ($41 >= 3) && ($43 >= 3) && ($45 >= 3) && ($47 >= 3) && ($49 >= 3)) { print } }' ${data}HM_alllyr.cov > HM_alllyr_filtered.cov

echo "lyrata done (HM)"

#LL

awk -F " " '{ if(($3 >= 3) && ($5 >= 3) && ($7 >= 3) && ($9 >= 3) && ($11 >= 3) && ($13 >= 3) && ($15 >= 3) && ($17 >= 3) && ($19 >= 3) && ($21 >= 3) && ($23 >= 3) && ($25 >= 3) && ($27 >= 3) && ($29 >= 3) && ($31 >= 3) && ($33 >= 3) && ($35 >= 3) && ($37 >= 3) && ($39 >= 3) && ($41 >= 3) && ($43 >= 3) && ($45 >= 3) && ($47 >= 3))  { print } }' ${data}LL_allhal.cov > LL_allhal_filtered.cov

echo "halleri done (LL)"

awk -F " " '{ if(($3 >= 3) && ($5 >= 3) && ($7 >= 3) && ($9 >= 3) && ($11 >= 3) && ($13 >= 3) && ($15 >= 3) && ($17 >= 3) && ($19 >= 3) && ($21 >= 3) && ($23 >= 3) && ($25 >= 3) && ($27 >= 3) && ($29 >= 3) && ($31 >= 3) && ($33 >= 3) && ($35 >= 3) && ($37 >= 3) && ($39 >= 3) && ($41 >= 3) && ($43 >= 3) && ($45 >= 3) && ($47 >= 3) && ($49 >= 3)) { print } }' ${data}LL_alllyr.cov > LL_alllyr_filtered.cov

echo "lyrata done (LL)"
```

The outputs from this step can be found at this link: (https://figshare.com/projects/Data_for_MDS_analyses/134765)[https://figshare.com/projects/Data_for_MDS_analyses/134765]
