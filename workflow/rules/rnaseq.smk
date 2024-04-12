# RNA sequencing

rnaseq_files_df = pd.read_excel(config['rnaseq_files'], header=0)

rnaseq_sampleIDs = rnaseq_files_df['SampleID'].tolist()

rnaseq_fastq_df = pd.concat([pd.Series(pd.concat([rnaseq_files_df['Filepath'].astype(str) + '_R1', rnaseq_files_df['Filepath'].astype(str) + '_R2']), name = 'old_path'), 
                             pd.Series(pd.concat([rnaseq_files_df['SampleID'].astype(str) + '_R1', rnaseq_files_df['SampleID'].astype(str) + '_R2']), name = 'new_path')], 
                            axis=1)

#print(rnaseq_fastq_df)

rnaseq_readfiles = rnaseq_fastq_df['new_path'].tolist()

# Preparation

rule rnaseq_copy_fastqs:
  input:
    readfile=lambda w: (rnaseq_fastq_df[rnaseq_fastq_df.new_path == w.readfile].old_path + '.fastq.gz').tolist()
  output:
    "resources/rnaseq/fastq/{readfile}.fastq.gz",
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying {input.readfile} to scratch ---"
  shell:
    """
    rsync -a "{input.readfile}" "{output}"
    """

rule copy_star_reference:
  input:
    config['rnaseq_reference_dir'],
  output:
    directory("resources/rnaseq/reference"),
  conda:
    "../envs/rsync.yaml"
  message:
    "--- Copying reference index genome to scratch ---"
  shell:
    """
    rsync -ar "{input}/" "{output}"
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
    "output/rnaseq/fastqc/{readfile}_fastqc.log",
  message:
    "--- FastQC {wildcards.readfile} ---"
  threads: 
    8
  benchmark:
    "output/rnaseq/fastqc/{readfile}/{readfile}.txt",
  shell:
    """
    fastqc --outdir="{output.fastqcdir}" \
           --threads {threads} \
           "{input}" \
           &> "{log}"
    """
    
rule rnaseq_multiqc:
  input:
    fastqc=expand("output/rnaseq/fastqc/{readfile}/{readfile}_fastqc.html", readfile=rnaseq_readfiles),
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

rule rnaseq_star_genome_load:
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
    STAR --genomeDir "{input.reference}" \
         --genomeLoad LoadAndExit
    """

rule rnaseq_star_pe:
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

rule rnaseq_star_genome_unload:
  input:
    reference="resources/rnaseq/reference",
    bam=expand("output/rnaseq/bam/{sampleID}/{sampleID}_pe_Aligned.sortedByCoord.out.bam", sampleID=rnaseq_sampleIDs),
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

rule rnaseq_samtools_filter:
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

rule rnaseq_samtools_index:
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
    
rule rnaseq_featurecounts:
  input:
    bams=expand("output/rnaseq/bam/{sampleID}/{sampleID}_filtered.bam", sampleID=rnaseq_sampleIDs),
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
                  -a "{input.annotationfile}" \
                  -t exon \
                  -g gene_id \
                  -p \
                  -o "{output}" \
                  {input.bams}
    """

# DE analyses

rule rnaseq_r_deseq2_import:
  input:
    sample_metadata_xlsx=config['sample_metadata'],
    rnaseq_files_xlsx=config['rnaseq_files'],
    counts_txt="output/rnaseq/counts/counts.txt",
  output:
    dds_rds="{basedir}/output/rnaseq/r/dds.Rds",
	  rld_rds="{basedir}/output/rnaseq/r/rld.Rds",
  conda:
    "../envs/r-deseq2.yaml"
  log:
    "{basedir}/output/rnaseq/deseq2_import.log"
  message:
    "--- Import reads into R and store it as a DESeq2 object ---"
  threads: 
    1
  shell:
    """    
    Rscript --vanilla workflow/scripts/rnaseq/deseq2_import.R "{input.sample_metadata_xlsx}" "{input.rnaseq_files_xlsx}" "{input.counts_txt}" "{output.dds_rds}" "{output.rld_rds}" &> "{log}"
    """

rule rnaseq_r_subsetting:
  input:
    dds_rds="{basedir}/output/rnaseq/r/dds.Rds",
  output:
    dds_subset_rds="{basedir}/output/rnaseq/r/dds_{treatment}.Rds",
    rld_subset_rds="{basedir}/output/rnaseq/r/rld_{treatment}.Rds",
  conda:
    "../envs/r-deseq2.yaml"
  log:
    "{basedir}/output/rnaseq/subsetting_{treatment}.log",
  benchmark:
    "{basedir}/output/rnaseq/subsetting_{treatment}_benchmark.txt",
  message:
    "--- Subsetting RNAseq data ---",
  params:
    treatment="{treatment}",
  threads:
    1
  resources:
    mem_mb=189000,
  shell:
    """
    Rscript --vanilla workflow/scripts/rnaseq/subsetting.R "{input.dds_rds}" "{params.treatment}" "{output.dds_subset_rds}" "{output.rld_subset_rds}"  &> "{log}"
    """

# rule rnaseq_r_eda:
#   input:
# 	  rlog_rds="output/rnaseq/r/rlog.Rds",
#   output:
#     heatmap_pairwise_correlation_pdf,
#     scatterplot_pc1_pc2_pdf,
#   conda:
#     "../envs/r-deseq2.yaml"
#   log:
#     "output/rnaseq/deseq2_eda.log"
#   message:
#     "--- Import reads into R DESeq2 and perform QC ---",
#   params:
#     treatment="{treatment}",
#   threads: 
#     1
#   shell:
#     """    
#     Rscript --vanilla workflow/scripts/rnaseq/eda.R "{input.counts_txt}" "{output.dds_rds}" "{output.rlog_rds}" &> "{log}"
#     """
 
rule rnaseq_r_de:
  input:
	  dds_subset_rds="{basedir}/output/rnaseq/r/dds_{treatment}.Rds",
  output:
    degs_csv="{basedir}/output/rnaseq/degs/degs_{treatment}.csv",
  conda:
    "../envs/r-deseq2.yaml"
  log:
    "{basedir}/output/rnaseq/de_{treatment}.log",
  benchmark:
    "{basedir}/output/rnaseq/de_{treatment}_benchmark.txt",
  message:
    "--- Performing DE analyses ---",
  params:
    treatment="{treatment}",
  threads:
    1
  shell:
    """
    Rscript --vanilla workflow/scripts/rnaseq/de_analyses.R "{input.dds_subset_rds}" "{params.treatment}" "{output.degs_csv}" &> "{log}"
    """

rule rnaseq_r_fgsea:
  input:
    degs_csv="{basedir}/output/rnaseq/degs/degs_{treatment}.csv",
  output:
    fgsea_go_list_rds="{basedir}/output/rnaseq/fgsea/fgsea_go_{treatment}_fgseaobj_list.Rds",
    fgsea_go_csv="{basedir}/output/rnaseq/fgsea/fgsea_go_{treatment}.csv",
    fgsea_kegg_list_rds="{basedir}/output/rnaseq/fgsea/fgsea_kegg_{treatment}_fgseaobj_list.Rds",
    fgsea_kegg_csv="{basedir}/output/rnaseq/fgsea/fgsea_kegg_{treatment}.csv",
  conda:
    "../envs/r-fgsea.yaml"
  log:
    "{basedir}/output/rnaseq/fgsea_{treatment}.log",
  benchmark:
    "{basedir}/output/rnaseq/fgsea_{treatment}_benchmark.txt",
  message:
    "--- Performing DE FGSEA analyses ---",
  threads:
    1
  shell:
    """
    Rscript --vanilla workflow/scripts/rnaseq/fgsea_analyses.R "{input.degs_csv}" "{output.fgsea_go_list_rds}" "{output.fgsea_go_csv}" "{output.fgsea_kegg_list_rds}" "{output.fgsea_kegg_csv}" &> "{log}"
    """
