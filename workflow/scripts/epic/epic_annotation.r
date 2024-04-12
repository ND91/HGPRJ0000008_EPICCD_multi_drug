#!/usr/bin/env R
# Append updated genetic variant annotation and promoter capture Hi-C enhancer information.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

### Ugly section start
# The following section is super ugly, but I cannot seem to get conda to install the packages (without failing - they appear to conflict with one another). 
# I tried packrat/renv, but that had its own issues. 

# Install the Illumina HumanMethylation EPIC 10b5 (hg38) annotation
if(!"IlluminaHumanMethylationEPICanno.ilm10b5.hg38" %in% rownames(installed.packages())){
  devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest")
  devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
  require(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
} else{
  require(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
}

if(!"TxDb.Hsapiens.UCSC.hg19.knownGene" %in% rownames(installed.packages())){
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
} else{
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
}

if(!"ChIPpeakAnno" %in% rownames(installed.packages())){
  BiocManager::install("ChIPpeakAnno")
  require(ChIPpeakAnno)
} else{
  require(ChIPpeakAnno)
}

if(!"org.Hs.eg.db" %in% rownames(installed.packages())){
  BiocManager::install("org.Hs.eg.db")
  require(org.Hs.eg.db)
} else{
  require(org.Hs.eg.db)
}

### Ugly section end

library(minfi)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPpeakAnno)
library(GenomicRanges)
library(org.Hs.eg.db)

source("workflow/scripts/epic/functions.r")

activepromoters_path <- args[1]
epic_annotation_csv_path <- args[2]

epic_anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
epic_anno$pos_hg38 <- as.numeric(epic_anno$Start_hg38)+1
epic_anno <- data.frame(epic_anno) %>%
  dplyr::filter(!is.na(CHR_hg38))

epic_anno_hg38_gr <- makeGRangesFromDataFrame(epic_anno, keep.extra.columns = T, seqnames.field = "CHR_hg38", start.field = "pos_hg38", end.field = "pos_hg38")

# Append better enhancer annotations

epic_anno <- as.data.frame(epic_anno_hg38_gr) %>%
  dplyr::rename(chr_hg38 = seqnames,
                pos_hg38 = start) %>%
  dplyr::select(-end)
epic_anno_hg19_gr <- makeGRangesFromDataFrame(epic_anno, keep.extra.columns = T, seqnames.field = "chr", start.field = "pos", end.field = "pos")

pchic <- read.csv(activepromoters_path, sep = "\t")
pchic_gr_bait <- makeGRangesFromDataFrame(pchic, keep.extra.columns = T, start.field = "baitSt", end.field = "baitEnd", seqnames.field = "baitChr")
pchic_gr_anno_bait <- annotatePeakInBatch(myPeakList = pchic_gr_bait, AnnotationData=genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
pchic_gr_anno_bait <- addGeneIDs(annotatedPeak = pchic_gr_anno_bait,
                                 orgAnn = "org.Hs.eg.db",
                                 feature_id_type = "entrez_id",
                                 IDs2Add = "symbol")
pchic_gr_anno_oe <- data.frame(pchic_gr_anno_bait) %>%
  dplyr::rename(baitChr = seqnames,
                baitSt = start,
                baitEnd = end) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, seqnames.field = "oeChr", start.field = "oeSt", end.field = "oeEnd")
epic_anno_hg19_gr$Enhancer_gene <- NA
pchic_epic_anno_hg19_gr_overlap <- findOverlaps(pchic_gr_anno_oe, epic_anno_hg19_gr)
epic_anno_hg19_gr$Enhancer_gene[subjectHits(pchic_epic_anno_hg19_gr_overlap)] <- pchic_gr_anno_oe[queryHits(pchic_epic_anno_hg19_gr_overlap),]$symbol

epic_anno <- data.frame(epic_anno_hg19_gr) %>%
  dplyr::rename(chr_hg19 = seqnames,
                pos_hg19 = start) %>%
  dplyr::select(-c(end, width))

# Save data

write.csv(epic_anno, epic_annotation_csv_path)

sessionInfo()