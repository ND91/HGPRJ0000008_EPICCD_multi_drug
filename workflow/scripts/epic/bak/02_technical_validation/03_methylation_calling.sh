#!/usr/bin/env bash
# Align Miseq-sequenced bisulfite-converted reads to GRCh37

workDir="/media/hdd1/ayliyim/bsseq"
samDir="${workDir}/sam"

outDir="${workDir}/methylation"
mkdir "${outDir}"

for sampleDir in ${samDir}/*
do
	sampleID=`basename ${sampleDir}`
	echo "${sampleID}"
	mkdir -p "${outDir}/${sampleID}"
    
	bamfile=${sampleDir}/*.bam
    
	#echo ${bamfile}
	
	bismark_methylation_extractor --gzip --bedGraph ${bamfile}
done