# snake-GENESPACE run

## INTRODUCTION

This directory contains code for the generation of synteny plots. 
Synteny analysis was performed using snake-GENESPACE (Yip Tong et al., in prep), a snakemake workflow around the tool GENESPACE (Lovell et al., 2022). 

snake-GENESPACE offers a reproducibility report which includes md5 sums of input and output files as well as a record of the version of each tool used by snake-GENESPACE. 

## DATA

The A. lyrata NT1 assembly and annotation were downloaded manually at:
https://figshare.com/projects/Arabidopsis_lyrata_genome_assemblies/162343

A.lyrata JGI and A.thaliana data was obtained from the Joint Genome Institut (JGI) database.
The downloaded data are fasta files with all peptide sequences and gff3 annotations of genes.
For A.lyrata: https://phytozome-next.jgi.doe.gov/info/Alyrata_v2_1
For A. thaliana: https://phytozome-next.jgi.doe.gov/info/Athaliana_TAIR10
Files were downloaded using this command:
```console
curl --cookie jgi_session=/api/sessions/068cf3e6e17b1df79495186cff4312a1 --output download.20240315.170054.zip -d "{\"ids\":{\"Phytozome-167\":[\"52b9c702166e730e43a34e56\",\"52b9c700166e730e43a34e53\"],\"Phytozome-384\":[\"585486937ded5e78cff8c522\",\"585486957ded5e78cff8c529\"]}}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/

```
Files were then processed manually to achieve the GENESPACE desired format. This was done prior to the creation of snake-GENESPACE, which is why we did not use snake-GENESPACE in that case. This already formated data can be directly fed into snake-GENESPACE, along with unformated gff annotations and fasta assemblies.
The peptide fasta and primary transcript bed files were processed (renaming sequences) using the following code:
```console
###################
## A. lyrata JGI ##
###################
## Peptide
mv genomes/lyr_jgi/Alyrata_384_v2.1.protein_primaryTranscriptOnly.fa genomes/lyr_jgi/A.lyrata_jgi_peptide.fa
awk -F '[.]' '{print $1}' genomes/lyr_jgi/A.lyrata_jgi_peptide.fa > genespace/peptide/A.lyrata_jgi.fa # A.lyrata jgi peptide fasta is ready. 
md5sum genespace/peptide/A.lyrata_jgi.fa # 94f0ffd10802ff01f468014ce17c0518
## Bed
mv genomes/lyr_jgi/Alyrata_384_v2.1.gene.gff3 genomes/lyr_jgi/A.lyrata_jgi_gene.gff3
# Converting gff to bed (Keep only scaffolds 1 to 8 as these match chromosomes 1-8 of our assembly):
awk -F '\t|;' '$1 ~ /scaffold_[1-8]$/ && $3 ~ /gene/ { print $1"\t"($4-1)"\t"$5"\t"$10 }' genomes/lyr_jgi/A.lyrata_jgi_gene.gff3 > genomes/lyr_jgi/A.lyrata_jgi_gene.bed
sed 's/Name=//g' genomes/lyr_jgi/A.lyrata_jgi_gene.bed > genespace/bed/A.lyrata_jgi.bed # A.lyrata jgi bed is ready.
md5sum genespace/bed/A.lyrata_jgi.bed # f4402b9b8297c79b5e42f3eb3f323bf5

#################
## A. thaliana ##
#################
## Peptide 
mv genomes/thaliana/Athaliana_167_protein_primaryTranscriptOnly.fa genomes/thaliana/A.thaliana_peptide.fa
awk -F '[.]' '{print $1}' genomes/thaliana/A.thaliana_peptide.fa > genespace/peptide/A.thaliana.fa # A.thaliana peptide fasta is ready.
md5sum genespace/peptide/A.thaliana.fa # a8b97ad8800141b2488d26d691bbcf6b
## Bed
mv genomes/thaliana/Athaliana_167_gene.gff3 genomes/thaliana/A.thaliana_gene.gff3
md5sum genomes/thaliana/A.thaliana_gene.gff3 # 6a96ce30bca76c16fd509decbf68ed4e
# Converting gff to bed:
awk -F '\t|;' '$3 ~ /gene/ { print $1"\t"($4-1)"\t"$5"\t"$10 }' genomes/thaliana/A.thaliana_gene.gff3 > genomes/thaliana/A.thaliana_gene.bed
sed 's/Name=//g' genomes/thaliana/A.thaliana_gene.bed > genespace/bed/A.thaliana.bed # A.thaliana bed is ready. 
md5sum genespace/bed/A.thaliana.bed # 978a349da8ea430d2005f6476b22a87a
```

The input directory for snake-GENESPACE was then created with a bed and peptide directory containing the A.thaliana and A.lyrata JGI data. It also had one directory each for A.halleri, A.lyrata (new assemblies and annotations) and A.lyrata_NT1 with fasta assembly and gff annotation. 

Some assemblies were also subsetted to keep only scaffolds identified as chromosomes:
```console

#############
# A.halleri #
#############
full_hal="./full_assemblies/A.halleri_genome.fa"
chr_only_hal="./input_dir/A.halleri/A.halleri_chr_only.fa"

md5sum $full_hal #f4127e868e9f03b382f1bf997457abd9

# Process the input file to keep only sequences starting with "chr"
grep -A 1 "^>chr" "$full_hal" > "$chr_only_hal"

md5sum $chr_only_hal #9df416da2ece1841934c2a180f57a6d5

############
# A.lyrata #
############
full_lyr="./full_assemblies/A.lyrata_genome.fa"
chr_only_lyr="./input_dir/A.lyrata/A.lyrata_chr_only.fa"

md5sum $full_lyr #5c6f86a940724a21075f4bcc98c2b5de

# Process the input file to keep only sequences starting with "chr"
grep -A 1 "^>chr" "$full_lyr" > "$chr_only_lyr"

md5sum $chr_only_lyr #2196214f04b5f474b4eec36eaf6a96d7

##############
# A.thaliana #
##############
full_bed_thal="./full_assemblies/A.thaliana_all_chromo.bed"
auto_chr_only_that="./input_dir/bed/A.thaliana.bed"

md5sum $full_bed_thal #978a349da8ea430d2005f6476b22a87a

# Process the input file to keep only sequences starting with "chr"
grep -E '^Chr[1-5]\b' $full_bed_thal > $auto_chr_only_that

md5sum $auto_chr_only_that #5785307472e7dd885cfc9de204799a72
```


## ANALYSIS 

snake-GENESPACE was downloaded from github:
```console
git clone https://github.com/kenji-yt/snake-GENESPACE.git"
```
The run was then performed on the unziped input directory:
```console
snakemake --use-conda --config INPUT=../input_dir -c20 
```
Additional scripts can be found in this directory:

- edit_figs.R is an R script used to create modified versions of the standard GENESPACE output plots. 

- modified_dotplot.R is an R script adapted for GENESPACE to allow white backgrounds. It is used in "edit_figs.R".     

## OUTPUTS: 

You can find the snake-GENESPACE reproducibility report and run logs in the results directory. Other snake-GENESPACE outputs were not included due to their large size.

