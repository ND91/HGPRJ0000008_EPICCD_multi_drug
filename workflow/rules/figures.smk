rule fig_rebuttal_R3P6:
  input:
    rgset_qc="{basedir}/output/epic/qc/rgset_qc.Rds",
    validation_predictions_csv="{basedir}/output/epic/metadata/priorantitnf_summary_glmmodel_predictions_validation.csv",
  output:
    fig_effectsizeplot_pdf="{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR.pdf",
    fig_effectsizeplot_png="{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR.png",
  conda:
    "../envs/r-plot.yaml"
  log:
    "{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR.log",
  benchmark:
    "{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR_benchmark.txt",
  resources:
    mem_mb=20000,
  message:
    "--- Effectsize plot T1RvNR vs T2RvNR ---",
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/effectsizeplot_T1RvNR_T2RvNR.r "{input.dmps_anno}" "{input.horaizon_predictors_xlsx}" "{params.treatment}" "{output.fig_effectsizeplot_pdf}" "{output.fig_effectsizeplot_png}" &> "{log}"
    """

rule effectsizeplot_T1RvNR_T2RvNR:
  input:
    dmps_anno="{basedir}/output/epic/dmp/dmp_{treatment}_annotated.csv",
    horaizon_predictors_xlsx=config['horaizon_markers'],
  output:
    fig_effectsizeplot_pdf="{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR.pdf",
    fig_effectsizeplot_png="{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR.png",
  conda:
    "../envs/r-plot.yaml"
  log:
    "{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR.log",
  benchmark:
    "{basedir}/output/figures/effectsizeplot_{treatment}_T1RvNR_T2RvNR_benchmark.txt",
  resources:
    mem_mb=20000,
  message:
    "--- Effectsize plot T1RvNR vs T2RvNR ---",
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/effectsizeplot_T1RvNR_T2RvNR.r "{input.dmps_anno}" "{input.horaizon_predictors_xlsx}" "{params.treatment}" "{output.fig_effectsizeplot_pdf}" "{output.fig_effectsizeplot_png}" &> "{log}"
    """

rule r_volcanoplot_T2vT1:
  input:
    dmps_anno="{basedir}/output/epic/dmp/dmp_{treatment}_annotated.csv",
  output:
    fig_volcanoplot_RT2vT1_pdf="{basedir}/output/figures/volcanoplot_{treatment}_RT2vT1.pdf",
    fig_volcanoplot_NRT2vT1_pdf="{basedir}/output/figures/volcanoplot_{treatment}_NRT2vT1.pdf",
    fig_effectsizeplot_RT2vT1_NRT2vT1_pdf="{basedir}/output/figures/effectsizeplot_{treatment}_RT2vT1_NRT2vT1.pdf",
  conda:
    "../envs/r-plot.yaml"
  log:
    "{basedir}/output/figures/r_volcanoplot_{treatment}_T2vT1.log",
  benchmark:
    "{basedir}/output/figures/r_volcanoplot_{treatment}_T2vT1_benchmark.txt",
  resources:
    mem_mb=20000,
  message:
    "--- Volcanoplot T2vT1 ---",
  params:
    treatment="{treatment}",
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_volcanoplot_T2vT1.r "{input.dmps_anno}" "{input.horaizon_predictors_xlsx}" "{params.treatment}" "{output.fig_volcanoplot_RT2vT1_pdf}" "{output.fig_volcanoplot_NRT2vT1_pdf}" "{output.fig_effectsizeplot_RT2vT1_NRT2vT1_pdf}" &> "{log}"
    """
