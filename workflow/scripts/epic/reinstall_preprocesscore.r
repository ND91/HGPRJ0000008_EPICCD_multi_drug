#!/usr/bin/env R
# Due to some incompatibilities between the default preprocessCore and OpenBLAS, I cannot run the normalization method (preprocessFunnorm or preprocessQuantile for that matter) that I want to. 
# It will return the following error: "return code from pthread_create() is 22"
# A hack to circumvent this is to reinstall preprocessCore. This will need to be done every time the environment is installed
# For more information, see:
# - https://support.bioconductor.org/p/122925/
# - https://support.bioconductor.org/p/117119/#123061

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop(paste0("Script needs 1 argument. Current input is:", args))
}

done_path <- args[1]

BiocManager::install("preprocessCore", configure.args="--disable-threading", force = T)

BiocManager::install("preprocessCore", configure.args="--disable-threading", force = T)

file.create(done_path)

sessionInfo()