#!/bin/bash

## name sort bams from STAR output for counting in featurecounts


for bc in `ls -1 | grep .*.bam`;

do

v=$(echo $bc | sed 's/\(.*\).bam/\1/') ;


samtools sort -n -o "$v"_sorted.bam -@ 8 $bc ;

done


## count feature expression from bam files using featurecount from subread package

#run with 8 cores, rev stranded, paired end

list=$(echo *_sorted.bam)

featureCounts -T 8 -s 2 -p -a "./sars2_jq1_counts.txt" \
 	$list


## clean count matrix for use in analysis with DESEQ2

#remove columns specifying gene position, chromosome etc.
cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21 sars2_jq1_counts.txt > sars2_jq1_counts.Rmatrix.txt

#remove first line showing command used to generate file
tail -n +2 sars2_jq1_counts.Rmatrix.txt > FILE.tmp && mv FILE.tmp sars2_jq1_counts.Rmatrix.txt

#remove .bam extension from sample name
sed 's/_sorted.bam//g' sars2_jq1_counts.Rmatrix.txt > FILE.tmp && mv FILE.tmp sars2_jq1_counts.Rmatrix.txt

