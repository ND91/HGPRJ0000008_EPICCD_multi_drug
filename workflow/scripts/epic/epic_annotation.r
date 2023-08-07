#!/usr/bin/env R
# Append updated genetic variant annotation and promoter capture Hi-C enhancer information.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

# Install the Illumina HumanMethylation EPIC 10b5 (hg38) annotation
if(!"IlluminaHumanMethylationEPICanno.ilm10b5.hg38" %in% rownames(installed.packages())){
  devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest")
  devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
  require(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
} else{
  require(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
}

# # Install dbSNP (version 155)
# if(!"SNPlocs.Hsapiens.dbSNP151.GRCh38" %in% rownames(installed.packages())){
#   BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
#   require(SNPlocs.Hsapiens.dbSNP151.GRCh38)
# } else{
#   require(SNPlocs.Hsapiens.dbSNP151.GRCh38)
# }

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(org.Hs.eg.db))

source("workflow/scripts/epic/functions.r")

activepromoters_path <- args[1]
epic_annotation_csv_path <- args[2]

epic_anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
epic_anno$pos_hg38 <- as.numeric(epic_anno$Start_hg38)+1
epic_anno <- data.frame(epic_anno) %>%
  dplyr::filter(!is.na(CHR_hg38))

epic_anno_hg38_gr <- makeGRangesFromDataFrame(epic_anno, keep.extra.columns = T, seqnames.field = "CHR_hg38", start.field = "pos_hg38", end.field = "pos_hg38")

# Updated genetic variants from dbSNP

# epic_anno_hg38_gr <- epic_anno_hg38_gr+50
# seqlevelsStyle(epic_anno_hg38_gr) <- seqlevelsStyle(SNPlocs.Hsapiens.dbSNP151.GRCh38)
# 
# epic_anno_hg38_gr_snps <- snpsByOverlaps(x = SNPlocs.Hsapiens.dbSNP151.GRCh38, ranges = epic_anno_hg38_gr)
#  
# epic_anno_hg38_snp_overlap <- findOverlaps(epic_anno_hg38_gr, epic_anno_hg38_gr_snps)

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