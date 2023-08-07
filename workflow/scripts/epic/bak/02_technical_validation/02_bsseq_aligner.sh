#!/usr/bin/env bash
# Align Miseq-sequenced bisulfite-converted reads to GRCh37

# Move to scratch
workDir="/media/hdd1/ayliyim/bsseq"
genomeDir="${workDir}/fasta"
outDir="${workDir}/aligned"
mkdir "${outDir}"

fastqDir="${HOME}/archive/raw_data/fastq/bsseq/GENDX_R210510"

cd "${workDir}"

for sampleDir in ${fastqDir}/*
do
	sampleID=`basename ${sampleDir}`
	echo "${sampleID}"
	mkdir -p "${workDir}/${sampleID}"
    
	reads1=${sampleDir}/*_R1_001.fastq.gz
    reads2=${sampleDir}/*_R2_001.fastq.gz
    
	echo ${reads1}
	echo ${reads2}
	
	bismark --bowtie2 -p 8 --temp_dir "${outDir}/${sampleID}" --output_dir "${outDir}/${sampleID}" ${genomeDir} -1 ${reads1} -2 ${reads2}
done