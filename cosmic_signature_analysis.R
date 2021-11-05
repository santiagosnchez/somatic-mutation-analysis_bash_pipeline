#!/bin/Rscript

.libPaths('/hpf/largeprojects/tabori/software/R_libs/4.1.0/')

library(sigminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(pals)

# function to match signature to etiology
extract_matched_sigs <- function(mt_sig, sim){
  exposure = as.data.frame(mt_sig$Exposure.norm)
  # signatures as column
  exposure$signatures = rownames(exposure)
  # make DF longer
  exposure = exposure %>% pivot_longer(cols=-signatures, values_to="exposure.norm", names_to="samples")
  # add matched COSMIC signature
  # check if object is Signature
  if (class(mt_sig) == "Signature"){
    matched_names = unlist(sapply(sim$best_match[ exposure$signatures ], function(x) x["aetiology"]))
    exposure$ethiology = matched_names
  } else {
    all_names = colnames(sim$rss)
    ethiology = sim$aetiology_db[[1]]
    names(ethiology) = all_names
    exposure$ethiology = ethiology[ exposure$signatures ]
  }
  # get raw exposure values
  exposure$exposure.raw = (pivot_longer(as.data.frame(mt_sig$Exposure), cols=everything(), values_to="exposure.raw", names_to="samples"))$exposure.raw
  # return DF
  return(exposure)
}

produce_bar_plots <- function(df, file_name, legend_names=c(NULL,NULL)){
  # load colors
  cols = c(stepped(20), stepped2(20), stepped(20))
  cols = cols[-2]
  names(cols) = NULL
  # calculate height from samples
  n_samples = length(unique(df$samples))
  # set minimal size of the plot
  if (n_samples < 4){
    n_samples = 4
  }
  # start plotting device
  dev.new(file="tmp_Rplot.pdf", height=n_samples, width=14)
  # get number of legend elements
  # guides(fill=guide_legend(ncol=2))

  # adjust sig labels to match ethiologies
  max_size_of_ethiology_label = max(sapply(df$ethiology, nchar))
  max_size_of_sig_label = max(sapply(df$signatures, nchar))
  # df$dummy_sig_label = paste0(df$signatures, paste(rep(" ", max_size_of_ethiology_label - max_size_of_sig_label), collapse=""))
  # # get sample label size
  max_size_of_samp_label = max(sapply(df$samples, nchar))

  # plot each panel
  sig1 = ggplot(df, aes(y=samples, x=exposure.norm, fill=signatures)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=cols, name=legend_names[1]) +
    labs(x="fraction", y="", title="normalized") +
    guides(fill=guide_legend(ncol=3)) +
    theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      legend.position="right")
  sig2 = ggplot(df, aes(y=samples, x=exposure.raw, fill=signatures)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=cols, name=legend_names[1]) +
    labs(x="SNVs", y="", title="stacked counts") +
    theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      axis.text.y=element_blank(),
      legend.position="none")
  eti3 = ggplot(df, aes(y=samples, x=exposure.norm, fill=ethiology)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=cols, name=legend_names[2]) +
    labs(x="fraction", y="", title="normalized") +
    theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      legend.position="right")
  eti4 = ggplot(df, aes(y=samples, x=exposure.raw, fill=ethiology)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=cols, name=legend_names[2]) +
    labs(x="SNVs", y="", title="stacked counts") +
    theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      axis.text.y=element_blank(),
      legend.position="none")

  # get legends as ggplot objects
  legend_sig1 = get_legend(sig1)
  legend_eti3 = get_legend(eti3)

  # plot grid
  plot_grid(sig1 + theme(legend.position="none", axis.title.x=element_blank()),
            sig2 + theme(axis.title.x=element_blank()),
            legend_sig1,
            eti3 + theme(legend.position="none", plot.title=element_blank()),
            eti4 + theme(plot.title=element_blank()),
            legend_eti3,
            ncol=3,
            rel_heights=c(6,5.55),
            rel_widths=c(max_size_of_samp_label,max_size_of_samp_label*.6,max_size_of_samp_label*1.25)
            )
  # save PNG
  ggsave(paste0("analyses/",file_name,".png"), bg="white", type="cairo", height=n_samples, width=14)
  ggsave(paste0("analyses/",file_name,".pdf"), device="pdf", height=n_samples, width=14)
  # close device
  dev.off()
  # remove tmp plot file
  #file.remove("tmp_Rplot.pdf")
  #
  print("Done")
}

produce_lolipop_plots <- function(df){
  # calculate height from samples
  n_samples = length(unique(df$tumor))
  # set minimal size of the plot
  if (n_samples < 4){
    n_samples = 4
  }
  # start plotting device
  dev.new(file="tmp_Rplot.pdf", height=n_samples, width=14)
  # SNVs
  tmb_snvs = ggplot(tmb_data, aes(x=tumor, y=tmb_snvs, color=tmb_snvs)) +
    geom_point(size=3) +
    scale_y_continuous(expand=c(0,0), limits=c(0,max(pretty(tmb_data$tmb_snvs))+5)) +
    geom_linerange(aes(ymin=0, ymax=tmb_snvs)) +
    scale_color_gradient(high="red", low="blue", name=NULL) +
    labs(y="mutations/mb", title="tumor mutation burden (SNV)") +
    background_grid(major="y") +
    theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      axis.title.x=element_blank()
    )
  # Indels
  tmb_indels = ggplot(tmb_data, aes(x=tumor, y=tmb_indels, color=tmb_indels)) +
    geom_point(size=3) +
    scale_y_continuous(expand=c(0,0), limits=c(0,max(pretty(tmb_data$tmb_indels))+5)) +
    geom_linerange(aes(ymin=0, ymax=tmb_indels)) +
    scale_color_gradient(high="red", low="blue", name=NULL) +
    labs(y="mutations/mb", title="tumor mutation burden (indels)") +
    background_grid(major="y") +
    theme(
      axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
      axis.title.x=element_blank()
    )
  # plot grid
  plot_grid(tmb_snvs, tmb_indels, ncol=1)
  ggsave("analyses/snv_and_indel_tmb.png", bg="white", type="cairo", height=8, width=n_samples)
  ggsave("analyses/snv_and_indel_tmb.pdf", height=8, width=n_samples)
  # close device
  dev.off()
}

# get args and set db_type
args = base::commandArgs(trailingOnly = TRUE)
if (args[1] == "wes"){
  db_type="human-exome"
} else {
  db_type="human-genome"
}

# get working dir
current_dir=getwd()

# vector of filtered vcf files
vcfs <- list.files(paste0(current_dir, "/vcf"), "*.vcf.gz$", full.names = TRUE)

# make/read vcf as maf
maf <- read_vcf(vcfs)

# make tally of maf objects
# compares each VCF to the hg38 reference and extracts
# (1) a 6-way (SBS-6) matrix of single nucleotide changes
# (2) a 96-way (SBS-96) matrix of trinucleotide changes (default/main matrix)
# ...
# (6) a 6144-way (SBS-6144) matrix of
# useSyn is a flag to include synonymous changes
# add_trans_bias includes categories under transciptional bias
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  useSyn = TRUE,
  mode = "ALL",
  add_trans_bias = TRUE,
  cores = 8
)

# get mutation signatures using Bayesian approach coupled with NMF optimization
# Do first just for the SBS-96 matrix
## add more info here
mt_sig_bayes_sbs_96 <- sig_unify_extract(
  mt_tally$SBS_96,
  range = 10,
  nrun = 10,
  cores = 8,
  approach = "bayes_nmf"
)

# match signatures using cosine similarity COSMIC v3
matched_mt_sig_sbs_96 <- get_sig_similarity(
  mt_sig_bayes_sbs_96,
  sig_db="SBS_hg38",
  db_type=db_type
)

# match signatures using cosine similarity legacy COSMIC v2
matched_mt_sig_legacy_30 <- get_sig_similarity(
  mt_sig_bayes_sbs_96,
  sig_db="legacy",
  db_type=db_type
)

# set threshold for matched signatures
sig_threshold = 0.05

# fit linear decomposition on SBS_hg38
linear_decomp_mt_sig_sbs_96 = list()
linear_decomp_mt_sig_sbs_96_raw <- sig_fit(mt_tally$SBS_96 %>% t(), sig_index = "ALL", sig_db = "SBS_hg38")
# make it proportional
linear_decomp_mt_sig_sbs_96_norm <- t(t(linear_decomp_mt_sig_sbs_96_raw) / colSums(linear_decomp_mt_sig_sbs_96_raw))
# ignore matches below sig_threshold
for (i in 1:dim(linear_decomp_mt_sig_sbs_96_norm)[2]){
  # find and set to zero
  set_to_zero = which(linear_decomp_mt_sig_sbs_96_norm[,i] < sig_threshold)
  linear_decomp_mt_sig_sbs_96_raw[set_to_zero,i] = 0
  # rescale
  linear_decomp_mt_sig_sbs_96_norm[,i] = linear_decomp_mt_sig_sbs_96_raw[,i] / sum(linear_decomp_mt_sig_sbs_96_raw[,i])
}
# remove non-matches
keep = apply(linear_decomp_mt_sig_sbs_96_norm, 1, sum) > 0
linear_decomp_mt_sig_sbs_96_norm = subset(as.data.frame(linear_decomp_mt_sig_sbs_96_norm), keep)
linear_decomp_mt_sig_sbs_96_raw = subset(as.data.frame(linear_decomp_mt_sig_sbs_96_raw), keep)

# store in list
linear_decomp_mt_sig_sbs_96$Exposure = linear_decomp_mt_sig_sbs_96_raw
linear_decomp_mt_sig_sbs_96$Exposure.norm = linear_decomp_mt_sig_sbs_96_norm

# fit lienar decomposition on "legacy" signaures (COSMIC v2)
linear_decomp_mt_sig_legacy_30 = list()
linear_decomp_mt_sig_legacy_30_raw <- sig_fit(mt_tally$SBS_96 %>% t(), sig_index = "ALL", sig_db = "legacy")
# make it proportional
linear_decomp_mt_sig_legacy_30_norm <- t(t(linear_decomp_mt_sig_legacy_30_raw) / colSums(linear_decomp_mt_sig_legacy_30_raw))
# ignore matches below sig_threshold
for (i in 1:dim(linear_decomp_mt_sig_legacy_30_norm)[2]){
  # find and set to zero
  set_to_zero = which(linear_decomp_mt_sig_legacy_30_norm[,i] < sig_threshold)
  linear_decomp_mt_sig_legacy_30_raw[set_to_zero,i] = 0
  # rescale
  linear_decomp_mt_sig_legacy_30_norm[,i] = linear_decomp_mt_sig_legacy_30_raw[,i] / sum(linear_decomp_mt_sig_legacy_30_raw[,i])
}
# remove non-matches
keep = apply(linear_decomp_mt_sig_legacy_30_norm, 1, sum) > 0
linear_decomp_mt_sig_legacy_30_norm = subset(as.data.frame(linear_decomp_mt_sig_legacy_30_norm), keep)
linear_decomp_mt_sig_legacy_30_raw = subset(as.data.frame(linear_decomp_mt_sig_legacy_30_raw), keep)
# store in list
linear_decomp_mt_sig_legacy_30$Exposure = linear_decomp_mt_sig_legacy_30_raw
linear_decomp_mt_sig_legacy_30$Exposure.norm = linear_decomp_mt_sig_legacy_30_norm

# match extracted COSMIC signatures to database of known etiologies
# Bayesian NMF
matched_exposures_bNMF_cosmic3 = extract_matched_sigs(mt_sig_bayes_sbs_96, matched_mt_sig_sbs_96)
# adjust names
matched_exposures_bNMF_cosmic3$samples = sub("\\..*","", matched_exposures_bNMF_cosmic3$samples)
# linear regression decomposition on COSMIC v3
matched_exposures_lrDecomp_cosmic3 = extract_matched_sigs(linear_decomp_mt_sig_sbs_96, matched_mt_sig_sbs_96)
# adjust names
matched_exposures_lrDecomp_cosmic3$samples = sub("\\..*","", matched_exposures_lrDecomp_cosmic3$samples)
# linear regression decomposition on "legacy" COSMIC v2
matched_exposures_lrDecomp_cosmic2 = extract_matched_sigs(linear_decomp_mt_sig_legacy_30, matched_mt_sig_legacy_30)
# adjust names
matched_exposures_lrDecomp_cosmic2$samples = sub("\\..*","", matched_exposures_lrDecomp_cosmic2$samples)

# write files
write.csv(matched_exposures_bNMF_cosmic3, file="analyses/matched_exposures_bNMF_cosmic3.csv")
write.csv(matched_exposures_lrDecomp_cosmic3, file="analyses/matched_exposures_lrDecomp_cosmic3.csv")
write.csv(matched_exposures_lrDecomp_cosmic2, file="analyses/matched_exposures_lrDecomp_cosmic2.csv")

# set cowplot theme
theme_set(theme_cowplot())

# make plots
# start with barplots
# second arg is the file name saved to "analyses"
produce_bar_plots(matched_exposures_bNMF_cosmic3, "barplot_matched_signatures_bNMF_cosmicV3", c("Bayesian NMF signatures","matched etiologies (COSMIC v3)"))
produce_bar_plots(matched_exposures_lrDecomp_cosmic3, "barplot_matched_signatures_lrDecomp_cosmicV3", c("COSMIC v3","matched etiologies"))
produce_bar_plots(matched_exposures_lrDecomp_cosmic2, "barplot_matched_signatures_lrDecomp_cosmicV2", c("COSMIC v2","matched etiologies"))

# plot tumor mutation burden
# read data
tmb_data = read.csv("analyses/coverage_and_tmb.csv")
# as lolipop plots
produce_lolipop_plots(tmb_data)


# save R objects to disk
save.image(file="analyses/mutational_signatures_as_R_object.Rdata")
