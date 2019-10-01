#!/bin/bash

infiles=(rsem/*.genes.results)
#countfile=csf_smartseq_counts.txt
tpmfile=csf_smartseq_tpm.txt

# first column gene ids
#cut -f1 ${infiles[1]} > $countfile
cut -f1 ${infiles[1]} > $tpmfile

for i in ${infiles[@]}; do
	sample=$(echo $i | sed 's/rsem\/\(.*\)\.genes\.results/\1/')
	echo $sample;
	#paste $countfile <(awk -v sample="$sample" 'NR==1 {print sample}; NR!=1{print $5}' $i) > tmp
	#mv tmp $countfile	

	paste $tpmfile <(awk -v sample="$sample" 'NR==1 {print sample}; NR!=1{print $6}' $i) > tmp
	mv tmp $tpmfile	
done

rm tmp
