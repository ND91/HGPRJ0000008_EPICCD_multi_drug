#!/usr/bin/env bash
# Create a bowtie genome index for GRCh37

# Move to scratch
workDir="/media/hdd1/ayliyim/bsseq"
mkdir -p "${workDir}"
cd "${workDir}"

if [ ! -d "${workDir}/fasta" ]; then
	echo "Copying human genome GRCh37 reference data to ${workDir}"
	cp -r "${HOME}/archive/common_data/genome/homo_sapiens/GRCh37/fasta" "${workDir}/fasta"
fi

bismark_genome_preparation --bowtie2 --verbose "${workDir}/fasta/"
#bismark_genome_preparation --path_to_aligner "/mnt/nfs_kgres/home/ayliyim/local/miniconda3/envs/bismark/bin/bowtie2/" --verbose "/media/hdd1/ayliyim/bsseq/fasta/"