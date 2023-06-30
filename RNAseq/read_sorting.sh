#!/usr/bin/env bash


# Paths to alignment files, assemblies, EAGLE and output folder

bam=02_alignment/
hal=Ahal_genome/Ahal_genome.fa
lyr=Alyr_genome/Alyr_gemome_renamed.fa
eagle=/path/to/eagle/eagle
output=03_read_sorting/


# Commands

for filename in ${bam}*SYN*_halAligned.sortedByCoord.out.bam; do
	base=$(basename $filename _halAligned.sortedByCoord.out.bam)
	$eagle-rc --ngi --paired --ref1=$hal --bam1=${bam}${base}_halAligned.sortedByCoord.out.bam --ref2=$lyr --bam2=${bam}${base}_lyrAligned.sortedByCoord.out.bam -o ${output}${base}_classified > ${output}${base}_classified_reads.list
done

for filename in ${bam}*TKS*_halAligned.sortedByCoord.out.bam; do
	base=$(basename $filename _halAligned.sortedByCoord.out.bam)
	$eagle-rc --ngi --paired --ref1=$hal --bam1=${bam}${base}_halAligned.sortedByCoord.out.bam --ref2=$lyr --bam2=${bam}${base}_lyrAligned.sortedByCoord.out.bam -o ${output}${base}_classified > ${output}${base}_classified_reads.list
done

for filename in ${bam}*ALK*_halAligned.sortedByCoord.out.bam; do
	base=$(basename $filename _halAligned.sortedByCoord.out.bam)
	$eagle-rc --ngi --paired --ref1=$hal --bam1=${bam}${base}_halAligned.sortedByCoord.out.bam --ref2=$lyr --bam2=${bam}${base}_lyrAligned.sortedByCoord.out.bam -o ${output}${base}_classified > ${output}${base}_classified_reads.list
done

