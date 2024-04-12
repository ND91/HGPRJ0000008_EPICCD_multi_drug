This snakemake pipeline includes the preparation and analyses of the data presented in "Peripheral blood DNA methylation signatures predict response to vedolizumab and ustekinumab in adult patients with Crohn's disease: The EPIC-CD study".

Important to realize is that the code here consists of two consecutive steps. Note that step 2 is dependent on data generated in step 1:
1. Data preparation. This part consists of the import, preparing data for machine learning, and other basic analyses. This part is snakemaked.
2. Figure generation. This part consists of summarizing the analyses into the figures. This part is manual. 

To run this pipeline, you will need to have `mamba`/`conda` and `snakemake` installed. easiest way is to setup a `mamba`/`conda` environment as defined in the `environment.yml` file, which will install `snakemake`, `Python`, `Graphviz`, as well as the Python packages `pandas` and `openpyxl`. Per the Snakemake documentation, it is advisable to use mamba for this. Run the following command to create an environment named `snakemake` with aforementioned software installed.

`mamba env create -f environment.yaml`

## Data needed

In terms of data, you will need:
- The `.idat` files as stored at the European Genome-phenome Archive (EGA) under accession IDs: . This data needs to be stored in the `resources` directory. Note that you will need to update `files_metadata.csv`.
- The promoter-capture Hi-C data as published as "Data S1" by [Javierre et al. 2016](https://www.sciencedirect.com/science/article/pii/S0092867416313228?via%3Dihub) for annotation purposes. The contents of this `.zip` needs to be stored in the `resources/pchic`.
- The ICC as published in "Supplementary Table 1" by [Joustra et al. 2022](https://www.sciencedirect.com/science/article/pii/S2352345X22002624?via%3Dihub). This dataset is used as resource to interrogate the stability of the predictor CpGs against. This `.csv` should be stored as `resources/joustra2022_icc.csv`.

## Running the preparation

### 1. Data preparation

The only option we need to have is `use-conda`. 

```
snakemake --use-conda
```

Feel free to add other options you find necessary. My command looked like the following based on the compute resources I had available:

```
snakemake --cores 12 --resources mem_mb=189000 --use-conda --rerun-triggers mtime
```

### 2. Figure generation

This section can only be run after "1. Data preparation" and will need to be run manually. All R code for generating the figures as presented in the manuscript (and more) are stored in `workflow/scripts/figures`. Open the `.Rmd` files and fix the `basepath` variable in the code to reflect your own locale. By default, the script will create a directory in `output/figures`. Important to note is that the figures will be generated separately, the final figures were assembled in Inkscape and Adobe Illustrator.
