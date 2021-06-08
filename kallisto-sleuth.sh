#!/bin/bash

threads=10
path=$1
ref_fasta=$2
gtf=$3

echo "building kallisto index"

kallisto index \
	-i kallisto_index \
	$ref_fasta 

echo "running kallisto"
mkdir kallisto_results
for i in `ls -d $path/*.gz | cut -d _ -f1,2,3,4 | sort | uniq`
	do echo "aligning $i" 
	rone="$i""_1.fastq.gz" 
	rtwo="$i""_2.fastq.gz"
	echo $rone $rtwo $i
	outname=`basename $i`
	kallisto quant -i kallisto_index \
				-b 100 \
				-o kallisto_results/$outname \
				-t $threads \
				-g $gtf \
				$rone \
				$rtwo 
	done

echo "running Sleuth"
rscript --vanilla run_sleuth.r .