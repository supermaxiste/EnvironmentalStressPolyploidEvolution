####### GENERAL PARAMETERS #########

DATA=00_data/
CORES=8
OUTPUT=02_alignment_new/
FOLDER_HAL=Ahal_genome/
FOLDER_LYR=Alyr_genome/
GENOME_HAL=Ahal_genome/Ahal_genome.fa
GENOME_LYR=Alyr_genome/Alyr_genome_renamed.fa

####### GENOME PREPARATION (INDEXING) ########

#STAR --runThreadN ${CORES} --runMode genomeGenerate --genomeDir ${FOLDER_HAL} --genomeFastaFiles ${GENOME_HAL}
#STAR --runThreadN ${CORES} --runMode genomeGenerate --genomeDir ${FOLDER_LYR} --genomeFastaFiles ${GENOME_LYR}


###### ALIGNMENT #######


for filename in ${DATA}/*SYN*R1.fastq.gz; do
	base=$(basename $filename _R1.fastq.gz)
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_HAL} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_hal --readFilesCommand zcat 
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_LYR} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_lyr --readFilesCommand zcat 
done

for filename in ${DATA}/*TKS*R1.fastq.gz; do
	base=$(basename $filename _R1.fastq.gz)
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_HAL} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_hal --readFilesCommand zcat
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_LYR} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_lyr --readFilesCommand zcat
done

for filename in ${DATA}/*ALK*R1.fastq.gz; do
	base=$(basename $filename _R1.fastq.gz)
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_HAL} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_hal --readFilesCommand zcat
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_LYR} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_lyr --readFilesCommand zcat
done


for filename in ${DATA}/*hal*R1.fastq.gz; do
	base=$(basename $filename _R1.fastq.gz)
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_HAL} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_hal --readFilesCommand zcat 
done

for filename in ${DATA}/*lyr*R1.fastq.gz; do
	base=$(basename $filename _R1.fastq.gz)
	STAR --runThreadN ${CORES} --genomeDir ${FOLDER_LYR} --readFilesIn ${DATA}${base}_R1.fastq.gz ${DATA}${base}_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT}${base}_lyr --readFilesCommand zcat 
done
	
