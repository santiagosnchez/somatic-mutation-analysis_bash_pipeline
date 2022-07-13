#!/bin/Rscript

.libPaths('/hpf/largeprojects/tabori/shared/software/R_libs/4.1.2/')

library(sigminer)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)
#library(viridis)
#library(pals)
#library(ggtext)
#library(biomaRt)


# get args
args = base::commandArgs(trailingOnly = TRUE)
# read args
data_type = args[1]
sample_name = args[2]
organism = args[3]
no_ob = args[4]

# define db_type for sigminer
if (organism == "human"){
  if (data_type == "wes"){
    db_type="human-exome"
  } else {
    db_type="human-genome"
  }
  sig_db = "latest_SBS_GRCh38"
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
  cosmic_exdata = "COSMIC_v3.2_SBS_GRCh38.rds"
} else if (organism == "mouse"){
  db_type=""
  sig_db = "latest_SBS_mm10"
  ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
  cosmic_exdata = "COSMIC_v3.2_SBS_mm10.rds"
}

# tumor - normal
TUMOR = sub("__.*","",sample_name)
NORMAL = sub(".*__","",sample_name)

# get working dir
current_dir=getwd()

##################
# Get Signatures #
##################

# vector of filtered vcf files
if (is.na(no_ob)){
  somaticvcfpath <- paste0("mutect2/", sample_name, ".mutect2.selected.",data_type,".vcf.gz")
  outdir="analyses"
} else {
  somaticvcfpath <- paste0("mutect2/", sample_name, ".mutect2.selected_no-obpriors.",data_type,".vcf.gz")
  outdir="analyses/no-obpriors"
}

# make/read vcf as maf
maf.som <- read_vcf(somaticvcfpath)

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
  ref_genome = ref_genome,
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
  range = 2:10,
  nrun = 20,
  cores = 8,
  approach = "bayes_nmf"
)

# match signatures using cosine similarity COSMIC v3
if (db_type == ""){
  matched_mt_sig_sbs_96 <- get_sig_similarity(
    mt_sig_bayes_sbs_96,
    sig_db=sig_db
  )
} else {
  matched_mt_sig_sbs_96 <- get_sig_similarity(
    mt_sig_bayes_sbs_96,
    sig_db=sig_db,
    db_type=db_type
  )
}

# match signatures using cosine similarity legacy COSMIC v2
if (organism == "human"){
  matched_mt_sig_legacy_30 <- get_sig_similarity(
    mt_sig_bayes_sbs_96,
    sig_db="legacy",
    db_type=db_type
  )
  # consolidate bNMF analyses
  etio2 = matched_mt_sig_legacy_30$aetiology_db[[1]]
  names(etio2) = colnames(matched_mt_sig_legacy_30$similarity)[order(as.numeric(sub("COSMIC_","",colnames(matched_mt_sig_legacy_30$similarity))))]
}

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

# etio3tab = read.table("/hpf/largeprojects/tabori/shared/resources/cosmic_v3.2/cosmic_db_v3.2_signatures_and_etiologies.txt", sep="\t")
# etio3 = etio3tab[,2]
# names(etio3) = etio3tab[,1]
# etio3[ grep("Poli", etio3) ] = sub("Poli","Poly",etio3[ grep("Poli", etio3) ])
cosmic_data = readRDS(system.file("extdata", cosmic_exdata, package = "sigminer"))

# set threshold for matched signatures
#sig_threshold = 0.01

# fit linear decomposition on latest COSMIC db v3.2
if (db_type == ""){
  linear_decomp_mt_sig_sbs_96 <- sig_fit(catalogue_matrix=mt_tally$all_matrices$SBS_96 %>% t(),
                                            sig=mt_sig_bayes_sbs_96$Signature,
                                            sig_index = "ALL",
                                            method="NNLS",
                                            type="relative",
                                            auto_reduce=TRUE,
                                            sig_db=sig_db)
  linear_decomp_mt_sig_legacy_30 <- sig_fit(catalogue_matrix=mt_tally$all_matrices$SBS_96 %>% t(),
                                            sig=mt_sig_bayes_sbs_96$Signature,
                                            sig_index = "ALL",
                                            method="NNLS",
                                            type="relative",
                                            sig_db = "legacy")
} else {
  linear_decomp_mt_sig_sbs_96 <- sig_fit(catalogue_matrix=mt_tally$all_matrices$SBS_96 %>% t(),
                                              sig=mt_sig_bayes_sbs_96$Signature,
                                              sig_index = "ALL",
                                              db_type=db_type,
                                              method="NNLS",
                                              type="relative",
                                              sig_db=sig_db)
}

if (organism == "human"){
  linear_decomp_mt_sig_legacy_30 <- sig_fit(catalogue_matrix=mt_tally$all_matrices$SBS_96 %>% t(),
                                            sig=mt_sig_bayes_sbs_96$Signature,
                                            sig_index = "ALL",
                                            db_type=db_type,
                                            method="NNLS",
                                            type="relative",
                                            sig_db = "legacy")
}

# prep for output 3.2
linear_decomp_cosmic3 = as.data.frame(linear_decomp_mt_sig_sbs_96) %>%
mutate(cosmic_db="v3.2") %>%
mutate(cosmic_signature = rownames(linear_decomp_mt_sig_sbs_96)) %>%
mutate(etiology=cosmic_data$aetiology[rownames(linear_decomp_mt_sig_sbs_96),])
colnames(linear_decomp_cosmic3)[1] = "contribution_proportion"
linear_decomp_cosmic3 = linear_decomp_cosmic3 %>%
  mutate(tumor=TUMOR) %>%
  mutate(normal=NORMAL) %>%
  mutate(method="NNLS")

linear_decomp_cosmic2 = as.data.frame(linear_decomp_mt_sig_legacy_30) %>%
mutate(cosmic_db="v2.0") %>%
mutate(cosmic_signature = rownames(linear_decomp_mt_sig_legacy_30)) %>%
mutate(etiology=cosmic_data$aetiology[rownames(linear_decomp_mt_sig_legacy_30),])
colnames(linear_decomp_cosmic2)[1] = "contribution_proportion"
linear_decomp_cosmic2 = linear_decomp_cosmic2 %>%
  mutate(tumor=TUMOR) %>%
  mutate(normal=NORMAL) %>%
  mutate(method="NNLS")

#write.csv(linear_decomp_cosmic, file="analyses/lrDecomp_cosmic.csv")

# read TMB and coverage data
tmb_data = read.csv("analyses/coverage_tmb_and_mmr_sigs.csv")
tmb = tmb_data %>% filter(Tumor==TUMOR & Normal==NORMAL)

# write to file
write.csv(linear_decomp_cosmic3, file=paste0(outdir, "/", sample_name, ".COSMIC_v3.2.signatures.csv"), quote=F, row.names=F)
write.csv(linear_decomp_cosmic2, file=paste0(outdir, "/", sample_name, ".COSMIC_v2.0.signatures.csv"), quote=F, row.names=F)


# print old output
# tmbs
# "Tumor,Normal,Observed_coverage,Expected_coverage,SNV,Indels,TMB_SNV,TMB_indels,"${sig_heads}
if (file.exists("analyses/old_output.tmbs.tsv")){
  cat(sample_name, tmb$SNV, tmb$TMB_SNV, "\n", sep="\t", append=T, file="analyses/old_output.tmbs.tsv")
} else {
  cat("sample", "total_SNV", "TMB", "\n", sep="\t", file="analyses/old_output.tmbs.tsv")
  cat(sample_name, tmb$SNV, tmb$TMB_SNV, "\n", sep="\t", append=T, file="analyses/old_output.tmbs.tsv")
}
# signatures v3.2
if (file.exists("analyses/old_output.v3.sigs.tsv")){
  tmp_sigs = rep(0, length(linear_decomp_cosmic3[,"contribution_proportion"]))
  names(tmp_sigs) = rownames(linear_decomp_cosmic3[,"etiology"])
  tmp_sigs[ rownames(linear_decomp_cosmic3[,"etiology"]) ] = linear_decomp_cosmic3[,"contribution_proportion"]
  cat(sample_name, tmp_sigs,  "\n", sep="\t", append=T, file="analyses/old_output.v3.sigs.tsv")
} else {
  cat("Signature", rownames(linear_decomp_cosmic3[,"etiology"]), "\n", sep="\t", file="analyses/old_output.v3.sigs.tsv")
  cat("Etiology", linear_decomp_cosmic3[,"etiology"], "\n", sep="\t", append=T, file="analyses/old_output.v3.sigs.tsv")
  tmp_sigs = rep(0, length(linear_decomp_cosmic3[,"contribution_proportion"]))
  names(tmp_sigs) = rownames(linear_decomp_cosmic3[,"etiology"])
  tmp_sigs[ rownames(linear_decomp_cosmic3[,"etiology"]) ] = linear_decomp_cosmic3[,"contribution_proportion"]
  cat(sample_name, tmp_sigs,  "\n", sep="\t", append=T, file="analyses/old_output.v3.sigs.tsv")
}
# signatures v2
if (organism == "human"){
  if (file.exists("analyses/old_output.v2.sigs.tsv")){
    tmp_sigs = rep(0, length(etio2))
    names(tmp_sigs) = names(etio2)
    tmp_sigs[ rownames(linear_decomp_mt_sig_legacy_30) ] = linear_decomp_mt_sig_legacy_30[,1]
    cat(sample_name, tmp_sigs,  "\n", sep="\t", append=T, file="analyses/old_output.v2.sigs.tsv")
  } else {
    cat("Signature", gsub("COSMIC_","Signature\\.", names(etio2)), "\n", sep="\t", file="analyses/old_output.v2.sigs.tsv")
    cat("Etiology", etio2, "\n", sep="\t", append=T, file="analyses/old_output.v2.sigs.tsv")
    tmp_sigs = rep(0, 30)
    names(tmp_sigs) = names(etio2)
    tmp_sigs[ rownames(linear_decomp_mt_sig_legacy_30) ] = linear_decomp_mt_sig_legacy_30[,1]
    cat(sample_name, tmp_sigs,  "\n", sep="\t", append=T, file="analyses/old_output.v2.sigs.tsv")
  }
}

print(paste("Done for", sample_name))

# save R objects to disk
# save.image(file="analyses/mutational_signatures_as_R_object.Rdata")
