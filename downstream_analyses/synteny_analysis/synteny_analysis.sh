##################################
### THIS PART WAS RUN ON LINUX ###
##################################
# Assumes current wd has a genomes directory and gffread directory. 
# The genomes directory has 4 subdirectories. One for each genome: halleri/, lyrata/, lyr_jgi/, thaliana/. 
# gffread can be downloaded here: https://github.com/gpertea/gffread

# Mainly what we do here is get the data in the right place and format for GENESPACE to run. 
# We first make a genespace directory with sub directories for the bed files and peptide fasta files: 
mkdir -p genespace/bed genespace/peptide

# Expand all files:
find genomes -type f | xargs -I {} gunzip {}

#--------------------------------------#
# Phytozome assemblies and annotations #
#--------------------------------------#
# lyr_jgi/ and thaliana/ have data coming from the Joint Genome Institut (JGI) database.
# The directories contain fasta files with all peptide sequences andc gff3 annotations of gene position information.
# For A.lyrata: https://phytozome-next.jgi.doe.gov/info/Alyrata_v2_1
# For A. thaliana: https://phytozome-next.jgi.doe.gov/info/Athaliana_TAIR10
# Files were downloaded using this command:
# curl --cookie jgi_session=/api/sessions/068cf3e6e17b1df79495186cff4312a1 --output download.20240315.170054.zip -d "{\"ids\":{\"Phytozome-167\":[\"52b9c702166e730e43a34e56\",\"52b9c700166e730e43a34e53\"],\"Phytozome-384\":[\"585486937ded5e78cff8c522\",\"585486957ded5e78cff8c529\"]}}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/
# Files were downloaded and moved into the respective directories. Original name and md5sums can be found in the JGI_download_manifests/ directory. 
# Here we expand, rename and get the md5sum of each files. 

###############
## A. lyrata ##
###############
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

#----------------------------------#
# Novel assemblies and annotations #
#----------------------------------#
# Here we start with our newly generated assemblies and annotations in their respective directories.
# We extract peptide sequences and reformat for genespace.

################
## A. halleri ##
################

## Bed
md5sum genomes/halleri/A.halleri_genome.fa # f4127e868e9f03b382f1bf997457abd9
md5sum genomes/halleri/A.halleri.gff #8abc2179ea061fe36638b25311540478

# Converting gff to gene bed file only keeping chromosomes:
awk -F '\t|;|=' '$1 ~ /chr[0-8]/ && $3 ~ /gene/ { print $1"\t"($4-1)"\t"$5"\t"$10 }' genomes/halleri/A.halleri.gff > genespace/bed/A.halleri.bed
# A.halleri bed is ready. 
md5sum genespace/bed/A.halleri.bed # 7df46f4f3f344666fc0bc161f5daf78b

## Peptide
# Make file with transcript name and corresponding gene name (last two blocks removes splicing variants):
awk -F '\t|;' '$1 ~ /chr[0-8]/ && $3 ~ /mRNA/ { print $9"\t"$10 }' genomes/halleri/A.halleri.gff \
| sed 's/ID=//g' \
| sed 's/Parent=//'g \
| sed 's/mRNA-1/mRNA-1\t/g' \
|  awk -F '\t' '!seen[$1]++ {print $1"\t"$NF}' > genomes/halleri/transcript2gene.key

# Extract peptide sequence using gffread:
./gffread/gffread genomes/halleri/A.halleri.gff -g genomes/halleri/A.halleri_genome.fa -J -E -y genomes/halleri/A.halleri_peptide.fa
md5sum genomes/halleri/A.halleri_peptide.fa # 4f36193d33c37d2cc1d7f96a89106372
# Rename peptide sequences and keep only primary gene model 
# We do this with a custom script which should be in the same directory
# (This takes a bit of running time as the script is not very efficient :D)  
bash format_genespace.sh genomes/halleri/transcript2gene.key genomes/halleri/A.halleri_peptide.fa genespace/peptide/A.halleri.fa # A.halleri peptide fasta is ready 
md5sum genespace/peptide/A.halleri.fa # 5c65f39de961d84c5be3c407e4745957

###############
## A. lyrata ##
###############

## Bed
md5sum genomes/lyrata/A.lyrata_genome.fa # 5c6f86a940724a21075f4bcc98c2b5de
md5sum genomes/lyrata/A.lyrata.gff # 3bd6e399c26c28cd4a70c3b93c5e7698

# Converting gff to gene bed file only keeping chromosomes:
awk -F '\t|;|=' '$1 ~ /chr[0-8]/ && $3 ~ /gene/ { print $1"\t"($4-1)"\t"$5"\t"$10 }' genomes/lyrata/A.lyrata.gff > genespace/bed/A.lyrata.bed
# A.lyrata bed is ready. 
md5sum genespace/bed/A.lyrata.bed # 46925600b394c434ca8075b4a76395b5 


## Peptide 
# Make file with transcript name and corresponding gene name (last block removes splicing variants):
awk -F '\t|;' '$1 ~ /chr[0-8]/ && $3 ~ /mRNA/ { print $9"\t"$10 }' genomes/lyrata/A.lyrata.gff \
| sed 's/ID=//g' \
| sed 's/Parent=//'g \
|  awk -F '\t|[.]' '!seen[$1]++ {print $1"\t"$NF}' > genomes/lyrata/transcript2gene.key

# Extract peptide sequence using gffread:
./gffread/gffread genomes/lyrata/A.lyrata.gff -g genomes/lyrata/A.lyrata_genome.fa -J -E -y genomes/lyrata/A.lyrata_peptide.fa
md5sum genomes/lyrata/A.lyrata_peptide.fa # d85fc9b534af1a5e2d808c8aa9784cb5
# Rename peptide sequences and keep only primary gene modeel using the same script as for halleri:
bash format_genespace.sh genomes/lyrata/transcript2gene.key genomes/lyrata/A.lyrata_peptide.fa genespace/peptide/A.lyrata.fa # A.lyrata peptide fasta is ready 
md5sum genespace/peptide/A.lyrata.fa # 4f2e6c6812ae240427d9f15db20816c7 


####################################
### THIS PART WAS RUN ON MAC OS  ###
####################################

## Running the following part requires:
# conda environment called "genespace_env" with Orthofinder installed. 
# Orthofinder dependencies installed: mcl, fastme & diamond. 
# Rstudio. 
# GENESPACE installed in R. 
# Compiled MCScanX scripts in the current directory (in an MCScanX directory).
# MCScanX can be downloaded here: https://github.com/wyp1125/MCScanX 

# Activate genespace environment (Orthofinder is installed on it): 
conda activate genespace_env

# Make a directory to hold the final publication figures:
mkdir figures/

# Run genespace script with all trailing arguments:
Rscript genespace_commands.R /genespace /MCScanX 

#Note: Step 8 failed due to "out of memory" error but the output was not relevant to our study. 

# Run a command to edit final figures for publication:
Rscript edit_fig.R