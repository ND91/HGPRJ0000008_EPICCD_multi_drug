#!/usr/bin/env bash
# Run FastQC on the .fastq files to ensure the files are of good quality

RUN_PATH="${HOME}/archive/projects/PRJ0000008_CDMEDRESP/anti_a4b7/data/bsseq/data"

for file in $(ls ${RUN_PATH}); do
	SAMPLE=`basename ${file}`
	OUTPUT=${RUN_PATH}/../fastqc/${SAMPLE}
	mkdir -p ${OUTPUT}
	
	fastqc -t 8 -o "${OUTPUT}" ${RUN_PATH}/${file} 
done
