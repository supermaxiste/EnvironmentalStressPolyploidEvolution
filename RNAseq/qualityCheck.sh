##### GENERAL PARAMETERS #####

FAST_OUTPUT=01_qualityCheck/
CORES=4
DATA=00_data

##### QUALITY CONTROL #######

for filename in ${DATA}/*fastq.gz; do
	fastqc -o ${FAST_OUTPUT} -t ${CORES} $filename
done
