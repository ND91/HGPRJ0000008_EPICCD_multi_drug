# Preparation

rule rnaseq_copy_fastqs:
  input:
    fastq=lambda w: (filedata[filedata.new_path == w.readfile].old_path + '.fastq.gz').tolist()
  output:
    "resources/rnaseq/fastq/{readfile}.fastq.gz"
  conda:
    "envs/rsync.yaml"
  message:
    "--- Copying {input.fastq} to scratch ---"
  shell:
    """
    rsync -a {input.fastq} {output}
    """

rule copy_star_reference:
  input:
    config['reference_dir']
  output:
    directory("resources/rnaseq/reference")
  conda:
    "envs/rsync.yaml"
  message:
    "--- Copying reference index genome to scratch ---"
  shell:
    """
    rsync -ar {input}/ {output}
    """

# QC

rule rnaseq_fastqc:
  input:
    "resources/rnaseq/fastq/{readfile}.fastq.gz"
  output:
    fastqcdir=directory("output/rnaseq/fastqc/{readfile}"),
    fastqchtml="output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.html",
    fastqczip="output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.zip",
  conda:
    "../envs/fastqc.yaml"
  log:
    "output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.log",
  message:
    "--- FastQC {wildcards.readfile} ---"
  threads: 
    8
  benchmark:
    "output/rnaseq/fastqc/{readfile}/{readfile}.txt",
  shell:
    """
    fastqc --outdir={output.fastqcdir} \
           --threads {threads} \
           {input} \
           &> {log}
    """
    
rule rnaseq_multiqc:
  input:
    fastqc=expand("output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.html", readfile=readfiles),
    featurecounts="output/rnaseq/counts/counts.txt",
  output:
    multiqcdir=directory('output/rnaseq/multiqc'),
  conda:
    "../envs/multiqc.yaml"
  log:
    "output/rnaseq/multiqc/multiqc.log"
  message:
    "--- MultiQC: Summarizing statistics ---"
  threads: 
    8
  shell:
    """    
    multiqc . -o {output.multiqcdir} &> {log}
    """

# Alignment and mapping

rule star_genome_load:
  input:
    reference="resources/rnaseq/reference"
  output:
    touch("output/rnaseq/genome_loaded.done")
  conda:
    "../envs/star.yaml"
  message:
    "--- STAR: load genome ---"
  threads: 
    8
  shell:
    """
    STAR --genomeDir {input.reference} \
         --genomeLoad LoadAndExit
    """

rule star_pe:
  input:
    reference="resources/rnaseq/reference",
    fastq_r1="resources/rnaseq/fastq/{sampleID}_R1.fastq.gz",
    fastq_r2="resources/rnaseq/fastq/{sampleID}_R2.fastq.gz",
    idx="output/rnaseq/genome_loaded.done",
  output:
    bamfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
    sjfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_SJ.out.tab",
    logfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Log.out",
    logprogressfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Log.progress.out",
    logfinalfile="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Log.final.out",
  conda:
    "../envs/star.yaml"
  log:
    "output/rnaseq/bam/{sampleID}/{sampleID}_STAR.log"
  params:
    sampleID="{sampleID}"
  message:
    "--- STAR: paired-ended alignment {params.sampleID} ---"
  threads: 
    10
  shell:
    """
    STAR --runThreadN {threads} \
         --genomeDir {input.reference} \
         --genomeLoad LoadAndKeep \
         --readFilesIn {input.fastq_r1} {input.fastq_r2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 50000000000 \
         --outFileNamePrefix output/rnaseq/bam/{wildcards.sampleID}/{wildcards.sampleID}_pe_ \
         |& tee -a {log}
    """

rule star_genome_unload:
  input:
    reference="resources/rnaseq/reference",
    bam=expand("output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam", sampleID=sampleIDs),
    idx="output/rnaseq/genome_loaded.done",
  output:
    "output/rnaseq/bam/STAR_genome_unload_Log.out"
  conda:
    "../envs/star.yaml"
  log:
    "output/rnaseq/bam/STAR.log"
  message:
    "--- STAR: unload genome ---"
  threads: 
    8
  shell:
    """
    STAR --genomeDir {input.reference} \
         --genomeLoad Remove \
         --outFileNamePrefix output/rnaseq/bam/STAR_genome_unload_ \
         |& tee -a {log}
    rm {input.idx}
    """

rule samtools_filter:
  input:
    bam="output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam",
  output:
    bam_filtered="output/rnaseq/bam/{sampleID}/{sampleID}_filtered.bam",
  conda:
    "../envs/samtools.yaml"
  log:
    "output/rnaseq/bam/{sampleID}/{sampleID}_samtools.log"
  params:
    sampleID="{sampleID}"
  message:
    "--- SAMtools: post-alignment processing {params.sampleID} ---"
  threads: 
    8
  shell:
    """        
    #Remove:
    # - unmapped reads and multiple mappings (4 + 256 = 260)
    # - reads with mapping score < 10
    # - mitochondrial sequences
    
    samtools view -@ {threads} -S -h -F 260 -q 10 {input.bam} | awk '($1 ~ /^@/) || ($3 != "MT") {{ print $0 }}' | samtools view -@ {threads} -b -o {output.bam_filtered} - |& tee -a {log}
    """

rule samtools_index:
  input:
    bam="output/rnaseq/bam/{sampleID}/{sampleID}_filtered.bam",
  output:
    bai="output/rnaseq/bam/{sampleID}/{sampleID}_filtered.bam.bai",
  conda:
    "../envs/samtools.yaml"
  log:
    "output/rnaseq/bam/{sampleID}/{sampleID}_samtools.log"
  params:
    sampleID="{sampleID}"
  message:
    "--- SAMtools: indexing {params.sampleID} ---"
  threads: 
    8
  shell:
    """        
    samtools index -b -@ {threads} {input.bam} |& tee -a {log}
    """
    
rule featurecounts:
  input:
    bams=expand("output/rnaseq/bam/{sampleID}/{sampleID}_filtered.bam", sampleID=sampleIDs),
    annotationfile=config['rnaseq_gtf_file'],
  output:
    "output/rnaseq/counts/counts.txt"
  conda:
    "../envs/featurecounts.yaml"
  log:
    "output/rnaseq/counts/featurecounts.log"
  message:
    "--- Subread::FeatureCounts: Feature counting ---"
  threads: 
    8
  shell:
    """    
    featureCounts -T {threads} \
                  -a {input.annotationfile} \
                  -t exon \
                  -g gene_id \
                  -p \
                  -o {output} \
                  {input.bams}
    """

# DESeq2 analyses

rule deseq2_import:
  input:
    counts_txt="output/rnaseq/counts/counts.txt",
  output:
    dds_rds="output/rnaseq/r/dds.Rds",
	rlog_rds="output/rnaseq/r/rlog.Rds",
  conda:
    "../envs/r.yaml"
  log:
    "output/rnaseq/counts/deseq2_import.log"
  message:
    "--- Import reads into R DESeq2 ---"
  threads: 
    1
  shell:
    """    
    Rscript --vanilla workflow/scripts/rnaseq_preparation/deseq2_import.R {input.counts_txt} {output.dds_rds} {output.rlog_rds} &> {log}
    """

rule deseq2_eda:
  input:
    dds_rds="output/rnaseq/r/dds.Rds",
	rlog_rds="output/rnaseq/r/rlog.Rds",
  output:
   
  conda:
    "../envs/r.yaml"
  log:
    "output/rnaseq/counts/deseq2_import.log"
  message:
    "--- Import reads into R DESeq2 and perform QC ---"
  threads: 
    1
  shell:
    """    
    Rscript --vanilla workflow/scripts/rnaseq_preparation/deseq2_import.R {input.counts_txt} {output.dds_rds} {output.rlog_rds} &> {log}
    """