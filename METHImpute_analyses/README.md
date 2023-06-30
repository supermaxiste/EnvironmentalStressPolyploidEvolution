# METHImpute analyses

All the analyses with [METHImpute](https://bioconductor.org/packages/3.14/bioc/html/methimpute.html) were done with the command line R script `METHImpute_analysisNoCategories.R`. To run this script you will need the following:

 1) Bismark outputs with information about all cytosines in the genome. These outputs are named `sample.CX_report.txt` and can be found at LINK TBD

 2) Annotation files in `gff` format for the species of interest. These can be downloaded here: LINK TBD
 
 3) The R package `methimpute` installed. To install it, find instructions [here](https://bioconductor.org/packages/release/bioc/html/methimpute.html) or create a `conda` environment using the `methimpute.yaml` file.

 After downloading all the necessary files, the R script can be run as follows:

 ```
 Rscript METHImpute_analysisNoCategories.R arg1 arg2 arg3 arg4
 ```

  - `arg1`: path to Bismark cytosine report for a given sample (e.g. `sample.CX_report.txt`), see point 1) above 
 
  - `arg2`: path to annotation file, see point 2) above 
 
  - `arg3`: prefix for output files 
 
  - `arg4`: path to output folder
  
**IMPORTANT**: line 59 is hardcoded to use 8 threads for the computation. If you're running this script on a laptop or a machine with fewer threads available, please modify the `num.threads` argument.

 ## Example use case

 Let's assume that we would like to analyze the sample `Arabi1` with `METHImpute`. We have a Bismark report `arabi1.CX_report.txt` and the annotation file `arabi1_anno.gff` in the same folder. We would like to keep the `arabi1` prefix and have all the output in the folder `arabi1_output`. Together with `arabi1.CX_report.txt` and `arabi1_anno.gff`, the `METHImpute_analysisNoCategories.R` script and `arabi1_output` are in the same folder. To run the analysis we would use the following command:

 ```
 Rscript METHImpute_analysisNoCategories.R arabi1.CX_report.txt arabi1_anno.gff arabi1 arabi1_output
 ```

 The output for this command would be the following:

  - `arabi1_correlationPlot.png`: a plot showing the correlation of methylation levels between neighbouring cytosines. It includes 6 sub-plots for each pairwise comparison between CG, CHG and CHH context.
  - `arabi1_conversionRates.txt`: imputed conversion rates for each context
  - `arabi1_MU.txt`: summary of the (imputed) total methylated and total unmethylated cytosines in the sample
  - `arabi1_high_conf_MU.txt`: summary of the (imputed) total methylated and total unmethylated cytosines in the sample with high confidence. Those are defined as cytosine with a posterior probability associated to imputed state equal or greater than `0.9`.
  - `arabi1_status.txt`: summary of the state of each cytosine in the genome. Each row in this file represents a cytosine. There are two columns: the first specifies the imputed methylation state (1 = unmethylated, 2 = methylated), while the second is the posterior probability (from 0 to 1) from the Hidden Markov Model associated to the given methylation state.
  - `arabi1_dataframe_meth.txt`: dataframe used to plot the average methylation level within and around gene bodies (500bp). This dataframe should be used in combination with the script `XX_v_XX_EnrichmentPlot.R` to produce the plots.
  - `arabi1_dataframe_rcmeth.txt`: dataframe with information about  recalibrated methylation level within and around gene bodies (500bp). This data is a predicted average methylation based on imputed methylation states.
