# External data

# rule download_javierre2016:
#   params:
#     javierre2016=config['javierre2016_pchic']
#   output:
#     "{basedir}/resources/pchic/javierre2016_sm.zip"
#   conda:
#     "../envs/wget.yaml"
#   log:
#     "resources/pchic/rule_download_pchic.log",
#   benchmark:
#     "resources/pchic/rule_download_pchic_benchmark.txt",
#   message:
#     "--- Downloading supplementary data from Javierre et al. 2016 (DOI: 10.1016/j.cell.2016.09.037) ---"
#   shell:
#     """
#     wget -O '{output}' -o '{log}' '{params.javierre2016}' 
#     """  

# rule unzip_javierre2016:
#   input:
#     "{basedir}/resources/pchic/javierre2016_sm.zip",
#   output:
#     "{basedir}/resources/pchic/ActivePromoterEnhancerLinks.tsv",
#   conda:
#     "../envs/unzip.yaml"
#   log:
#     "resources/pchic/rule_unzip_javierre2016.log",
#   benchmark:
#     "resources/pchic/rule_unzip_javierre2016_benchmark.txt",
#   message:
#     "--- Downloading supplementary data from Javierre et al. 2016 (DOI: 10.1016/j.cell.2016.09.037) ---"
#   shell:
#     """
#     unzip "{input}" -d "resources/pchic" &> "{log}"
#     unzip "resources/pchic/DATA_S1.zip" -d "resources/pchic" &> "{log}"
#     """  

# Preprocessing

rule reinstall_preprocesscore:
  output:
    preprocesscore_installed="{basedir}/resources/preprocesscore.installed",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/resources/reinstall_preprocesscore.log",
  benchmark:
    "{basedir}/resources/reinstall_preprocesscore_benchmark.txt",
  message:
    "--- Reinstall preprocessCore ---",
  threads: 
    1
  resources:
    mem_mb=15000,
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/reinstall_preprocesscore.r "{output.preprocesscore_installed}" &> "{log}"
    """  

rule r_import:
  input:
    sample_metadata_xlsx=config['sample_metadata'],
    epic_files_xlsx=config['epic_files'],
    preprocesscore_installed="{basedir}/resources/preprocesscore.installed",
  output:
    rgset_rds="{basedir}/output/epic/import/rgset.Rds",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/import/r_import.log",
  benchmark:
    "{basedir}/output/epic/import/r_import_benchmark.txt",
  message:
    "--- Import the Illumina HumanMethylation EPIC data ---",
  threads: 
    1
  resources:
    mem_mb=64000,
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/import.r "{input.sample_metadata_xlsx}" "{input.epic_files_xlsx}" "{output.rgset_rds}" &> "{log}"
    """  

rule r_qc:
  input:
    rgset_rds="{basedir}/output/epic/import/rgset.Rds",
  output:
    rgset_qc_rds="{basedir}/output/epic/qc/rgset_qc.Rds",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/qc/r_qc.log",
  benchmark:
    "{basedir}/output/epic/qc/r_qc_benchmark.txt",
  message:
    "--- QC of the Illumina HumanMethylation EPIC data ---",
  threads: 
    1
  resources:
    mem_mb=150000,
  params:
    basedir="{basedir}"
  shell:
    """
    Rscript -e \"rmarkdown::render(input = 'workflow/scripts/epic/qc.Rmd', output_dir = 'output/epic/qc', knit_root_dir = '{params.basedir}', params = list(rgset_rds = '{input.rgset_rds}', rgset_qc_rds = '{output.rgset_qc_rds}'))\" &> "{log}"
    """

rule r_subsetting_normalization:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
  output:
    gmset="{basedir}/output/epic/subset/gmset_{treatment}.Rds",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/subset/r_subsetting_normalization_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/subset/r_subsetting_normalization_{treatment}_benchmark.txt",
  message:
    "--- Subsetting and normalizing the Illumina HumanMethylation EPIC data ---",
  params:
    treatment="{treatment}",
  threads:
    1
  resources:
    mem_mb=189000,
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/subsetting_normalization.r "{input.rgset_qc}" "{output.gmset}" "{params.treatment}" &> "{log}"
    """

rule r_epic_annotation:
  input:
    activepromoters="{basedir}/resources/pchic/ActivePromoterEnhancerLinks.tsv",
  output:
    epic_annotation_csv="{basedir}/resources/epic_annotations.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/resources/r_epic_annotation.log",
  benchmark:
    "{basedir}/resources/r_epic_annotation_benchmark.txt",
  threads:
    1
  resources:
    mem_mb=50000,
  message:
    "--- Append important annotations to the standard minfi-provided annotations ---"
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/epic_annotation.r "{input.activepromoters}" "{output.epic_annotation_csv}" &> "{log}"
    """
 
# EDA

rule r_eda:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
  output:
    obs1="{basedir}/workflow/scripts/epic/eda_files/figure-html/obs1-1.png",
    obs2="{basedir}/workflow/scripts/epic/eda_files/figure-html/obs2-1.png",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/eda/r_eda.log",
  benchmark:
    "{basedir}/output/epic/eda/r_eda_benchmark.txt",
  resources:
    mem_mb=20000,
  message:
    "--- EDA of the Illumina HumanMethylation EPIC data ---"
  shell:
    """
    Rscript --vanilla -e \"rmarkdown::render(input = 'workflow/scripts/epic/eda.Rmd', knit_root_dir = '{basedir}/')\" &> "{log}"
    """  

# Differential methylation analyses

rule r_dmp_analyses:
  input:
    gmset="{basedir}/output/epic/subset/gmset_{treatment}.Rds",
  output:
    dmps="{basedir}/output/epic/dmp/dmp_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/dmp/r_dmp_analyses_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/dmp/r_dmp_analyses_{treatment}_benchmark.txt",
  resources:
    mem_mb=40000,
  params:
    treatment="{treatment}"
  message:
    "--- Perform differential methylation analysis on the Illumina HumanMethylation EPIC data ---"
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/dmp_analyses.r "{input.gmset}" "{output.dmps}" "{params.treatment}" &> "{log}"
    """
    
rule r_dmr_analyses:
  input:
    gmset="{basedir}/output/epic/subset/gmset_{treatment}.Rds",
  output:
    dmrcate_rds="{basedir}/output/epic/dmr/dmrcate_{comparison}_{treatment}.Rds",
  conda:
    "../envs/r-dmrcate.yaml"
  log:
    "{basedir}/output/epic/dmr/r_dmr_analyses_{comparison}_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/dmr/r_dmr_analyses_{comparison}_{treatment}_benchmark.txt",
  resources:
    mem_mb=50000,
  params:
    treatment="{treatment}",
    comparison="{comparison}",
    dmr_csv="{basedir}/output/epic/dmr/dmrcate_{comparison}_{treatment}.csv",
  message:
    "--- Perform differential methylation region analysis on the Illumina HumanMethylation EPIC data ---"
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/dmr_analyses.r "{input.gmset}" "{output.dmrcate_rds}" "{params.treatment}" "{params.comparison}" "{params.dmr_csv}" &> "{log}"
    """

rule r_dmg_analysis:
  input:
    dmps_anno="{basedir}/output/epic/dmp/dmp_{treatment}_annotated.csv",
    gmset="{basedir}/output/epic/subset/gmset_{treatment}.Rds",
  output:
    dmgs="{basedir}/output/epic/dmg/dmg_{comparison}_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/dmg/r_dmg_analysis_{comparison}_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/dmg/r_dmg_analysis_{comparison}_{treatment}_benchmark.txt",
  resources:
    mem_mb=40000,
  message:
    "--- Perform differential methylation gene analysis on the DMPs ---"
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/dmg_analyses.r "{input.dmps_anno}" "{input.gmset}" "{output.dmgs}" "{params.comparison}" &> "{log}"
    """

rule r_icc_analysis:
  input:
    gmset="{basedir}/output/epic/subset/gmset_{treatment}.Rds",
  output:
    iccs="{basedir}/output/epic/icc/icc_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/icc/r_icc_analysis_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/icc/r_icc_analysis_{treatment}_benchmark.txt",
  resources:
    mem_mb=98000,
  message:
    "--- Perform intraclass correlation analysis ---"
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/icc_analyses.r "{input.gmset}" "{output.iccs}" "{params.treatment}" &> "{log}"
    """

# Interpreting differential methylation analyses

rule r_dmp_annotation:
  input:
    dmps="{basedir}/output/epic/dmp/dmp_{treatment}.csv",
    epic_annotation_csv="{basedir}/resources/epic_annotations.csv",
  output:
    dmps_anno="{basedir}/output/epic/dmp/dmp_{treatment}_annotated.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/dmp/r_dmp_annotation_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/dmp/r_dmp_annotation_{treatment}_benchmark.txt",
  resources:
    mem_mb=40000,
  message:
    "--- Append annotations to the DMPs ---"
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/dmp_annotation.r "{input.dmps}" "{input.epic_annotation_csv}" "{output.dmps_anno}" &> {log}
    """

rule r_dmp_ora:
  input:
    dmps_anno="{basedir}/output/epic/dmp/dmp_{treatment}_annotated.csv",
  output:
    ora_kegg="{basedir}/output/epic/dmp/ora/dmp_{comparison}_{treatment}_ora_kegg.csv",
    ora_go="{basedir}/output/epic/dmp/ora/dmp_{comparison}_{treatment}_ora_go.csv",
    ora_promhyper_kegg="{basedir}/output/epic/dmp/ora/dmp_{comparison}_{treatment}_ora_promhyper_kegg.csv",
    ora_promhyper_go="{basedir}/output/epic/dmp/ora/dmp_{comparison}_{treatment}_ora_promhyper_go.csv",
    ora_promhypo_kegg="{basedir}/output/epic/dmp/ora/dmp_{comparison}_{treatment}_ora_promhypo_kegg.csv",
    ora_promhypo_go="{basedir}/output/epic/dmp/ora/dmp_{comparison}_{treatment}_ora_promhypo_go.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/dmp/r_dmp_ora_{comparison}_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/dmp/r_dmp_ora_{comparison}_{treatment}_benchmark.txt",
  resources:
    mem_mb=40000,
  message:
    "--- Perform gene set overrepresentation analyses on the DMPs ---"
  params:
    comparison="{comparison}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/dmp_ora.r "{input.dmps_anno}" "{output.ora_kegg}" "{output.ora_go}" "{output.ora_promhyper_kegg}" "{output.ora_promhyper_go}" "{output.ora_promhypo_kegg}" "{output.ora_promhypo_go}" "{params.comparison}" &> "{log}"
    """

# Preparation HORAIZON machine learning

rule r_preparation_horaizon_training:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
  output:
    X="{basedir}/output/epic/horaizon/training/X_{treatment}.csv",
    y="{basedir}/output/epic/horaizon/training/metadata_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/horaizon/training/r_preparation_horaizon_training_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/horaizon/training/r_preparation_horaizon_training_{treatment}_benchmark.txt",
  threads: 
    1
  resources:
    mem_mb=50000,
  message:
    "--- Prepare training data for Horaizon for {params.treatment} ---"
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/preparation_ml_training.r "{input.rgset_qc}" "{params.treatment}" "{output.X}" "{output.y}" &> "{log}"
    """

rule r_preparation_horaizon_validation:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
  output:
    X="{basedir}/output/epic/horaizon/validation/X_{treatment}.csv",
    y="{basedir}/output/epic/horaizon/validation/metadata_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/horaizon/validation/r_preparation_horaizon_validation_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/horaizon/validation/r_preparation_horaizon_validation_{treatment}_benchmark.txt",
  threads: 
    1
  resources:
    mem_mb=50000,
  message:
    "--- Prepare validation data for Horaizon for {params.treatment} ---"
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/preparation_ml_validation.r "{input.rgset_qc}" "{params.treatment}" "{output.X}" "{output.y}" &> "{log}"
    """

rule r_preparation_horaizon_timepoint2_validation:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
  output:
    X="{basedir}/output/epic/horaizon/timepoint2_validation/X_{treatment}.csv",
    y="{basedir}/output/epic/horaizon/timepoint2_validation/metadata_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/horaizon/timepoint2_validation/r_preparation_horaizon_timepoint2_validation_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/horaizon/timepoint2_validation/r_preparation_horaizon_timepoint2_validation_{treatment}_benchmark.txt",
  threads: 
    1
  resources:
    mem_mb=50000,
  message:
    "--- Prepare multibiological validation data for Horaizon for {params.treatment} ---"
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/preparation_ml_timepoint2_validation.r "{input.rgset_qc}" "{params.treatment}" "{output.X}" "{output.y}" &> "{log}"
    """

rule r_preparation_horaizon_multibiological_validation:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
  output:
    X="{basedir}/output/epic/horaizon/multibiological_validation/X_{treatment}.csv",
    y="{basedir}/output/epic/horaizon/multibiological_validation/metadata_{treatment}.csv",
  conda:
    "../envs/r-minfi.yaml"
  log:
    "{basedir}/output/epic/horaizon/multibiological_validation/r_preparation_horaizon_multibiological_validation_{treatment}.log",
  benchmark:
    "{basedir}/output/epic/horaizon/multibiological_validation/r_preparation_horaizon_multibiological_validation_{treatment}_benchmark.txt",
  threads: 
    1
  resources:
    mem_mb=50000,
  message:
    "--- Prepare multibiological validation data for Horaizon for {params.treatment} ---"
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/epic/preparation_ml_multibiological_validation.r "{input.rgset_qc}" "{params.treatment}" "{output.X}" "{output.y}" &> "{log}"
    """

# Interpreting HORAIZON machine learning results

# rule r_horaizon_markers_cdpbmethlts:
#   input:
#     horaizon_markers_csv=config["horaizon_markers"],
#     cdpbmethlts_xlsx=config["cdpbmethlts"],
#   output:
#     horaizon_markers_cdpbmethlts_ICC_csv="{basedir}/output/epic/horaizon/cdpbmethlts/horaizon_markers_{treatment}_cdpbmethlts_ICC.csv",
#     horaizon_markers_cdpbmethlts_ICC_pdf="{basedir}/output/epic/horaizon/cdpbmethlts/horaizon_markers_{treatment}_cdpbmethlts_ICC.pdf",
#   conda:
#     "../envs/r-minfi.yaml"
#   log:
#     "{basedir}/output/epic/horaizon/cdpbmethlts/r_horaizon_markers_cdpbmethlts_{treatment}.log",
#   benchmark:
#     "{basedir}/output/epic/horaizon/cdpbmethlts/r_horaizon_markers_cdpbmethlts_{treatment}_benchmark.txt",
#   resources:
#     mem_mb=40000,
#   message:
#     "--- Annotate output from Horaizon for {params.treatment} with long-term stability data ---"
#   params:
#     treatment="{treatment}",
#   shell:
#     """
#     Rscript --vanilla workflow/scripts/epic/ml_lts_interrogation.r "{input.horaizon_markers_csv}" "{input.cdpbmethlts_xlsx}" "{params.treatment}" "{output.horaizon_markers_cdpbmethlts_ICC_csv}" "{output.horaizon_markers_cdpbmethlts_ICC_pdf}" &> "{log}"
#     """
