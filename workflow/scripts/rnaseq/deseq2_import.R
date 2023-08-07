#!/usr/bin/env Rscript
# The goal of this script is to import and preprocess the RNA-sequencing data counts.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

require(dplyr)
require(readxl)
require(DEseq2)
require(org.Hs.eg.db)
require(AnnotationDbi)

sample_metadata_path <- args[1] #"config/samples/sample_metadata.xlsx"
rnaseq_files_path <- args[2] #"config/samples/rnaseq_files.xlsx"
rnaseq_counts_path <- args[3] #"output/rnaseq/counts/counts.txt"
dds_path <- args[4] #"output/rnaseq/deseq2/dds.Rds"
rlog_path <- args[5] #"output/rnaseq/deseq2/rld.Rds"

# Import the sample metadata
sample_metadata <- readxl::read_excel(sample_metadata_path)
rnaseq_files <- readxl::read_excel(rnaseq_files_path)

rnaseq_sample_metadata <- rnaseq_files %>%
  dplyr::left_join(sample_metadata, by = "SampleID") %>%
  data.frame()

rownames(rnaseq_sample_metadata) <- rnaseq_sample_metadata$SampleID

# Clean counts
rnaseq_counts_raw <- read.csv(rnaseq_counts_path, sep = "\t", skip = 1)
rnaseq_counts <- rnaseq_counts_raw[,-c(1:6)]
rownames(rnaseq_counts) <- rnaseq_counts_raw$Geneid
colnames(rnaseq_counts) <- gsub("^.+(F[0-9]+)(T[0-9]).+$", "\\1_\\2", colnames(rnaseq_counts))
rnaseq_counts <- rnaseq_counts[which(rowSums(rnaseq_counts) != 0),]

# Prepare feature metadata
rnaseq_feature_metadata <- data.frame(rnaseq_counts_raw[,1:6]) %>%
  dplyr::rename(ENSGv = Geneid) %>%
  dplyr::filter(ENSGv %in% rownames(rnaseq_counts)) %>%
  dplyr::mutate(ENSG = gsub("\\.[0-9]+$", "", ENSGv),
                HGNC = mapIds(x = org.Hs.eg.db, keys = ENSG, column = "SYMBOL", keytype = "ENSEMBL"),
                Entrez = mapIds(x = org.Hs.eg.db, keys = ENSG, column = "ENTREZID", keytype = "ENSEMBL"))

# Prepare DESeq2DataSet object

dds <- DESeq2::DESeqDataSetFromMatrix(countData = rnaseq_counts,
                                      colData = rnaseq_sample_metadata[colnames(rnaseq_counts),],
                                      design= ~1)

mcols(dds) <- cbind(mcols(dds), rnaseq_feature_metadata)

saveRDS(object = dds, file = dds_path, compress = "gzip")

rld <- DESeq2::rlog(dds)

saveRDS(object = rld, file = rlog_path, compress = "gzip")

sessionInfo()