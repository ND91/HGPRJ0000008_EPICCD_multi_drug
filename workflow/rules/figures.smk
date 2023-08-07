# BENIMID 2022

rule r_volcanoplot_RT2vT1:
  input:
    dmps_anno_Ada="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Adalimumab_annotated.csv",
    dmps_anno_Ifx="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Infliximab_annotated.csv",
    dmps_anno_Vdz="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Vedolizumab_annotated.csv",
    dmps_anno_Ust="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Ustekinumab_annotated.csv",
  output:
    fig_volcanoplot_png="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/volcanoplot_RT2vT1.png",
    fig_pairplot_png="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/pairplot_RT2vT1.png",
  conda:
    "../envs/r.yaml"
  log:
    "/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/r_volcanoplot_RT2vT1.log",
  benchmark:
    "/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/r_volcanoplot_RT2vT1_benchmark.txt",
  resources:
    mem_mb=20000,
  message:
    "--- Volcanoplot RT2vT1 ---"
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_volcanoplot_RT2vT1.r {input.dmps_anno_Ada} {input.dmps_anno_Ifx} {input.dmps_anno_Vdz} {input.dmps_anno_Ust} {output.fig_volcanoplot_png} {output.fig_pairplot_png} &> {log}
    """  

rule r_volcanoplot_T1RvNR:
  input:
    dmps_anno_Ada="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Adalimumab_annotated.csv",
    dmps_anno_Ifx="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Infliximab_annotated.csv",
    dmps_anno_Vdz="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Vedolizumab_annotated.csv",
    dmps_anno_Ust="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/epic/dmp/dmp_Ustekinumab_annotated.csv",
    predictor_cpgs_path="config/horaizon/predictor_cpgs.xlsx"
  output:
    fig_volcanoplot_png="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/volcanoplot_RT2vT1.png",
    fig_pairplot_png="/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/pairplot_RT2vT1.png",
  conda:
    "../envs/r.yaml"
  log:
    "/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/r_volcanoplot_RT2vT1.log",
  benchmark:
    "/home/andrewliyim/lkgres/projects/PRJ0000008_EPICCD/multi_drug/output/figures/r_volcanoplot_RT2vT1_benchmark.txt",
  resources:
    mem_mb=20000,
  message:
    "--- Volcanoplot RT2vT1 ---"
  shell:
    """
    Rscript --vanilla workflow/scripts/figures/fig_volcanoplot_RT2vT1.r {input.dmps_anno_Ada} {input.dmps_anno_Ifx} {input.dmps_anno_Vdz} {input.dmps_anno_Ust} {input.predictor_cpgs_path} {output.fig_volcanoplot_png} {output.fig_pairplot_png} &> {log}
    """  
