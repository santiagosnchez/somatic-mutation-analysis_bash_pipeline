#!/bin/Rscript

.libPaths('/hpf/largeprojects/tabori/shared/software/R_libs/4.0.2/')

library(sigminer)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(pals)
library(ggtext)
library(biomaRt)


# get args
args = base::commandArgs(trailingOnly = TRUE)

# define db_type for sigminer
if (args[1] == "wes"){
  db_type="human-exome"
} else {
  db_type="human-genome"
}

##################
# Get Signatures #
##################

# vector of filtered vcf files
somaticvcfpaths <- args[-1]
# filter out non VCF files
somaticvcfpaths = somaticvcfpaths[ grep("vcf$|vcf.gz$", somaticvcfpaths) ]

# get sample names
sample_names = sapply(somaticvcfpaths, function(x) {
  xx = strsplit(x, "/")[[1]]
  xx = xx[ length(xx) ]
  xx = sub("\\..*","",xx)
  return(xx)
})

# make/read vcf as maf
maf.som <- read_vcf(somaticvcfpaths, samples=sample_names)

# make tally of maf objects
# compares each VCF to the hg38 reference and extracts
# (1) a 6-way (SBS-6) matrix of single nucleotide changes
# (2) a 96-way (SBS-96) matrix of trinucleotide changes (default/main matrix)
# ...
# (6) a 6144-way (SBS-6144) matrix of
# useSyn is a flag to include synonymous changes
# add_trans_bias includes categories under transciptional bias
mt_tally <- sig_tally(
  maf.som,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  useSyn = TRUE,
  mode = "SBS",
  add_trans_bias = TRUE,
  cores = 8
)

# get mutation signatures using Bayesian approach coupled with NMF optimization
# Do first just for the SBS-96 matrix
## add more info here
mt_sig_bayes_sbs_96 <- sig_unify_extract(
  mt_tally$all_matrices$SBS_96,
  range = 1:100,
  nrun = 10,
  cores = 8,
  approach = "bayes_nmf"
)

# match signatures using cosine similarity COSMIC v3
matched_mt_sig_sbs_96 <- get_sig_similarity(
  mt_sig_bayes_sbs_96,
  sig_db="latest_SBS_GRCh38",
  db_type=db_type
)

# match signatures using cosine similarity legacy COSMIC v2
matched_mt_sig_legacy_30 <- get_sig_similarity(
  mt_sig_bayes_sbs_96,
  sig_db="legacy",
  db_type=db_type
)

# consolidate bNMF analyses
etio2 = matched_mt_sig_legacy_30$aetiology_db[[1]]
names(etio2) = colnames(matched_mt_sig_legacy_30$similarity)[order(as.numeric(sub("COSMIC_","",colnames(matched_mt_sig_legacy_30$similarity))))]

# cosmic 2
# bnmf_cosmic2  = as.data.frame(matched_mt_sig_legacy_30$similarity) %>%
#   mutate(raw_signature=rownames(matched_mt_sig_legacy_30$similarity)) %>%
#   pivot_longer(cols=-raw_signature, names_to="cosmic_signature", values_to="cosine_similarity") %>%
#   mutate(sample=sub("\\..*","",args[2])) %>%
#   mutate(cosmic_db="v2") %>%
#   mutate(tumor=TUMOR) %>%
#   mutate(normal=NORMAL) %>%
#   mutate(method="bNMF") %>%
#   mutate(etiology=etio2[cosmic_signature]) %>%
#   dplyr::select(tumor,normal,raw_signature,cosmic_signature,cosmic_db,cosine_similarity,etiology,method)

# cosmic 3
# etiologies db
# etio3 = matched_mt_sig_sbs_96$aetiology_db[[1]]
# names(etio3) = colnames(matched_mt_sig_sbs_96$similarity)[order(as.numeric(gsub("[A-Za-z]","",colnames(matched_mt_sig_sbs_96$similarity))))]
# etio3["SBS9"] = sub("Poli","Poly",etio3["SBS9"])

etio3tab = read.table("/hpf/largeprojects/tabori/shared/resources/cosmic_v3.2/cosmic_db_v3.2_signatures_and_etiologies.txt", sep="\t")
etio3 = etio3tab[,2]
names(etio3) = etio3tab[,1]
etio3[ grep("Poli", etio3) ] = sub("Poli","Poly",etio3[ grep("Poli", etio3) ])

# bnmf_cosmic3  = as.data.frame(matched_mt_sig_sbs_96$similarity) %>%
#   mutate(raw_signature=rownames(matched_mt_sig_sbs_96$similarity)) %>%
#   pivot_longer(cols=-raw_signature, names_to="cosmic_signature", values_to="cosine_similarity") %>%
#   mutate(sample=sub("\\..*","",args[2])) %>%
#   mutate(cosmic_db="v3") %>%
#   mutate(tumor=TUMOR) %>%
#   mutate(normal=NORMAL) %>%
#   mutate(method="bNMF") %>%
#   mutate(etiology=etio3[cosmic_signature]) %>%
#   dplyr::select(tumor,normal,raw_signature,cosmic_signature,cosmic_db,cosine_similarity,etiology,method)
#
# # merge
# bnmf_cosmic = bind_rows(bnmf_cosmic2, bnmf_cosmic3)

# save to csv
#write.csv(bnmf_cosmic, file=paste0("analyses/", sample_name, ".bnmf_cosmic.csv"), row.names=F)

# set threshold for matched signatures
sig_threshold = 0.01

# fit linear decomposition on latest COSMIC db v3.2

linear_decomp_mt_sig_sbs_96 <- sig_fit(catalogue_matrix=mt_tally$all_matrices$SBS_96 %>% t(),
                                            sig=mt_sig_bayes_sbs_96$Signature,
                                            sig_index = "ALL",
                                            db_type=db_type,
                                            method="NNLS",
                                            type="relative",
                                            sig_db = "latest_SBS_GRCh38")
# make it proportional
# linear_decomp_mt_sig_sbs_96_norm <- t(t(linear_decomp_mt_sig_sbs_96_raw) / colSums(linear_decomp_mt_sig_sbs_96_raw))
# # ignore matches below sig_threshold
# for (i in 1:dim(linear_decomp_mt_sig_sbs_96_norm)[2]){
#   # find and set to zero
#   set_to_zero = which(linear_decomp_mt_sig_sbs_96_norm[,i] < sig_threshold)
#   linear_decomp_mt_sig_sbs_96_raw[set_to_zero,i] = 0
#   # rescale
#   linear_decomp_mt_sig_sbs_96_norm[,i] = linear_decomp_mt_sig_sbs_96_raw[,i] / sum(linear_decomp_mt_sig_sbs_96_raw[,i])
# }
# # remove non-matches
# keep = apply(linear_decomp_mt_sig_sbs_96_norm, 1, sum) > 0
# linear_decomp_mt_sig_sbs_96_norm = subset(as.data.frame(linear_decomp_mt_sig_sbs_96_norm), keep)
# linear_decomp_mt_sig_sbs_96_raw = subset(as.data.frame(linear_decomp_mt_sig_sbs_96_raw), keep)

# fit lienar decomposition on "legacy" signaures (COSMIC v2)

linear_decomp_mt_sig_legacy_30 <- sig_fit(catalogue_matrix=mt_tally$all_matrices$SBS_96 %>% t(),
                                            sig=mt_sig_bayes_sbs_96$Signature,
                                            sig_index = "ALL",
                                            db_type=db_type,
                                            method="NNLS",
                                            type="relative",
                                            sig_db = "legacy")
# # make it proportional
# linear_decomp_mt_sig_legacy_30_norm <- t(t(linear_decomp_mt_sig_legacy_30_raw) / colSums(linear_decomp_mt_sig_legacy_30_raw))
# # ignore matches below sig_threshold
# for (i in 1:dim(linear_decomp_mt_sig_legacy_30_norm)[2]){
#   # find and set to zero
#   set_to_zero = which(linear_decomp_mt_sig_legacy_30_norm[,i] < sig_threshold)
#   linear_decomp_mt_sig_legacy_30_raw[set_to_zero,i] = 0
#   # rescale
#   linear_decomp_mt_sig_legacy_30_norm[,i] = linear_decomp_mt_sig_legacy_30_raw[,i] / sum(linear_decomp_mt_sig_legacy_30_raw[,i])
# }
# # remove non-matches
# keep = apply(linear_decomp_mt_sig_legacy_30_norm, 1, sum) > 0
# linear_decomp_mt_sig_legacy_30_norm = subset(as.data.frame(linear_decomp_mt_sig_legacy_30_norm), keep)
# linear_decomp_mt_sig_legacy_30_raw = subset(as.data.frame(linear_decomp_mt_sig_legacy_30_raw), keep)

# rename columns and add data to df
# linear_decomp_cosmic = bind_rows(as.data.frame(linear_decomp_mt_sig_legacy_30) %>%
#   mutate(cosmic_db="v2") %>%
#   mutate(cosmic_signature = rownames(linear_decomp_mt_sig_legacy_30)) %>%
#   mutate(etiology=etio2[cosmic_signature]),
linear_decomp_cosmic_v3 = as.data.frame(linear_decomp_mt_sig_sbs_96) %>%
mutate(cosmic_db="v3.2") %>%
mutate(cosmic_signature = rownames(linear_decomp_mt_sig_sbs_96)) %>%
mutate(etiology=etio3[as.character(cosmic_signature)]) %>%
pivot_longer(cols=-c(cosmic_db,cosmic_signature,etiology), names_to="sample", values_to="contribution_proportion") %>%
extract(sample, c("tumor","normal"), "(.*)__(.*)") %>%
dplyr::select(tumor, normal, cosmic_db, cosmic_signature, etiology, contribution_proportion)

linear_decomp_cosmic_v2 = as.data.frame(linear_decomp_mt_sig_legacy_30) %>%
mutate(cosmic_db="v2") %>%
mutate(cosmic_signature = rownames(linear_decomp_mt_sig_legacy_30)) %>%
mutate(etiology=etio2[as.character(cosmic_signature)]) %>%
pivot_longer(cols=-c(cosmic_db,cosmic_signature,etiology), names_to="sample", values_to="contribution_proportion") %>%
extract(sample, c("tumor","normal"), "(.*)__(.*)") %>%
dplyr::select(tumor, normal, cosmic_db, cosmic_signature, etiology, contribution_proportion)

linear_decomp_cosmic = bind_rows(linear_decomp_cosmic_v3, linear_decomp_cosmic_v2)


write.csv(linear_decomp_cosmic, "analyses/COSMIC_signatures_all_samples_as_cohort.csv", quote=F, row.names=F)
