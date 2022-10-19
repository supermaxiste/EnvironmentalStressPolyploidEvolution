#!/usr/bin/env bash 

######### GENERAL PARAMETERS ###########

bam_diplo=02_alignment/
bam_allo=03_read_sorting/
anno_hal=Ahal_genome/Ahal_genome_anno.gtf
anno_lyr=Alyr_genome/Alyr_genome_anno_renamed.gtf
output=04_count_tables/
cores=2

######### COMMANDS ##########

for filename in ${bam_diplo}HM_hal*.bam; do
	base=$(basename $filename _halAligned.sortedByCoord.out.bam)
	featureCounts -a $anno_hal -o ${output}${base}_counts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_diplo}HM_lyr*.bam; do
	base=$(basename $filename _lyrAligned.sortedByCoord.out.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_counts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_diplo}LL_hal*.bam; do
	base=$(basename $filename _halAligned.sortedByCoord.out.bam)
	featureCounts -a $anno_hal -o ${output}${base}_counts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_diplo}LL_lyr*.bam; do
	base=$(basename $filename _lyrAligned.sortedByCoord.out.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_counts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}HM_SYN*_classified1.ref.bam; do
	base=$(basename $filename _classified1.ref.bam)
	featureCounts -a $anno_hal -o ${output}${base}_halcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}HM_SYN*_classified2.ref.bam; do
	base=$(basename $filename _classified2.ref.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_lyrcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}LL_SYN*_classified1.ref.bam; do
	base=$(basename $filename _classified1.ref.bam)
	featureCounts -a $anno_hal -o ${output}${base}_halcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}LL_SYN*_classified2.ref.bam; do
	base=$(basename $filename _classified2.ref.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_lyrcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}HM_ALK*_classified1.ref.bam; do
	base=$(basename $filename _classified1.ref.bam)
	featureCounts -a $anno_hal -o ${output}${base}_halcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}HM_ALK*_classified2.ref.bam; do
	base=$(basename $filename _classified2.ref.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_lyrcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}LL_ALK*_classified1.ref.bam; do
	base=$(basename $filename _classified1.ref.bam)
	featureCounts -a $anno_hal -o ${output}${base}_halcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}LL_ALK*_classified2.ref.bam; do
	base=$(basename $filename _classified2.ref.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_lyrcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}HM_TKS*_classified1.ref.bam; do
	base=$(basename $filename _classified1.ref.bam)
	featureCounts -a $anno_hal -o ${output}${base}_halcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}HM_TKS*_classified2.ref.bam; do
	base=$(basename $filename _classified2.ref.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_lyrcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}LL_TKS*_classified1.ref.bam; do
	base=$(basename $filename _classified1.ref.bam)
	featureCounts -a $anno_hal -o ${output}${base}_halcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

for filename in ${bam_allo}LL_TKS*_classified2.ref.bam; do
	base=$(basename $filename _classified2.ref.bam)
	featureCounts -a $anno_lyr -o ${output}${base}_lyrcounts.txt --primary -p -T $cores --tmpDir $output $filename
done

