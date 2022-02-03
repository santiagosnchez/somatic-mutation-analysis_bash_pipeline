#!/bin/Rscript

.libPaths('/hpf/largeprojects/tabori/shared/software/R_libs/4.0.2/')

library(sigminer)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(pals)
library(ggtext)
library(biomaRt)


# get args and set db_type
args = base::commandArgs(trailingOnly = TRUE)
if (args[1] == "wes"){
  db_type="human-exome"
} else {
  db_type="human-genome"
}

# sample name
sample_name = args[2]

# tumor - normal
TUMOR = sub("__.*","",sample_name)
NORMAL = sub(".*__","",sample_name)

# get working dir
current_dir=getwd()

##################
# Get Signatures #
##################

# vector of filtered vcf files
somaticvcfpath <- paste0(current_dir, "/vcf/", args[2], ".mutect2.all.Somatic.annotated-funcotator.",args[1],".vcf.gz")

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

# consolidate bNMF analyses
etio2 = matched_mt_sig_legacy_30$aetiology_db[[1]]
names(etio2) = colnames(matched_mt_sig_legacy_30$similarity)[order(as.numeric(sub("COSMIC_","",colnames(matched_mt_sig_legacy_30$similarity))))]

# cosmic 2
bnmf_cosmic2  = as.data.frame(matched_mt_sig_legacy_30$similarity) %>%
  mutate(raw_signature=rownames(matched_mt_sig_legacy_30$similarity)) %>%
  pivot_longer(cols=-raw_signature, names_to="cosmic_signature", values_to="cosine_similarity") %>%
  mutate(sample=sub("\\..*","",args[2])) %>%
  mutate(cosmic_db="v2") %>%
  mutate(tumor=TUMOR) %>%
  mutate(normal=NORMAL) %>%
  mutate(method="bNMF") %>%
  mutate(etiology=etio2[cosmic_signature]) %>%
  dplyr::select(tumor,normal,raw_signature,cosmic_signature,cosmic_db,cosine_similarity,etiology,method) %>%
  filter(cosine_similarity > 0.5)

# cosmic 3
# etiologies db
etio3 = matched_mt_sig_sbs_96$aetiology_db[[1]]
names(etio3) = colnames(matched_mt_sig_sbs_96$similarity)[order(as.numeric(gsub("[A-Za-z]","",colnames(matched_mt_sig_sbs_96$similarity))))]

bnmf_cosmic3  = as.data.frame(matched_mt_sig_sbs_96$similarity) %>%
  mutate(raw_signature=rownames(matched_mt_sig_sbs_96$similarity)) %>%
  pivot_longer(cols=-raw_signature, names_to="cosmic_signature", values_to="cosine_similarity") %>%
  mutate(sample=sub("\\..*","",args[2])) %>%
  mutate(cosmic_db="v3") %>%
  mutate(tumor=TUMOR) %>%
  mutate(normal=NORMAL) %>%
  mutate(method="bNMF") %>%
  mutate(etiology=etio3[cosmic_signature]) %>%
  dplyr::select(tumor,normal,raw_signature,cosmic_signature,cosmic_db,cosine_similarity,etiology,method) %>%
  filter(cosine_similarity > 0.5)

# merge
bnmf_cosmic = bind_rows(bnmf_cosmic2, bnmf_cosmic3)

# save to csv
#write.csv(bnmf_cosmic, file=paste0("analyses/", sample_name, ".bnmf_cosmic.csv"), row.names=F)

# set threshold for matched signatures
sig_threshold = 0.05

# fit linear decomposition on SBS_hg38

linear_decomp_mt_sig_sbs_96_raw <- sig_fit(mt_tally$all_matrices$SBS_96 %>% t(), sig_index = "ALL", sig_db = "SBS_hg38")
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

# fit lienar decomposition on "legacy" signaures (COSMIC v2)

linear_decomp_mt_sig_legacy_30_raw <- sig_fit(mt_tally$all_matrices$SBS_96 %>% t(), sig_index = "ALL", sig_db = "legacy")
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

# rename columns and add data to df
linear_decomp_cosmic = bind_rows(linear_decomp_mt_sig_legacy_30_norm %>%
  mutate(cosmic_db="v2") %>%
  mutate(cosmic_signature = rownames(linear_decomp_mt_sig_legacy_30_norm)) %>%
  mutate(etiology=etio2[cosmic_signature]),
linear_decomp_mt_sig_sbs_96_norm %>%
mutate(cosmic_db="v3") %>%
mutate(cosmic_signature = rownames(linear_decomp_mt_sig_sbs_96_norm)) %>%
mutate(etiology=etio3[cosmic_signature]))
colnames(linear_decomp_cosmic)[1] = "contribution_proportion"
linear_decomp_cosmic = linear_decomp_cosmic %>%
  mutate(tumor=TUMOR) %>%
  mutate(normal=NORMAL) %>%
  mutate(method="LRD")

#write.csv(linear_decomp_cosmic, file="analyses/lrDecomp_cosmic.csv")

# read TMB and coverage data
tmb_data = read.csv("analyses/coverage_and_tmb.csv")
tmb = tmb_data %>% filter(tumor==TUMOR & normal==NORMAL)

# write merged usable output
# merge all
all_signatures = bind_rows(bnmf_cosmic, linear_decomp_cosmic)
all_data = inner_join(tmb_data, all_signatures %>% dplyr::select(-normal), by="tumor")
write.csv(all_data, file=paste0("analyses/", sample_name, ".signatures.csv"), quote=F, row.names=F)

# print old output
# tmbs
if (file.exists("analyses/old_output.tmbs.tsv")){
  cat(sample_name, tmb$snvs, tmb$tmb_snvs, "\n", sep="\t", append=T, file="analyses/old_output.tmbs.tsv")
} else {
  cat("sample", "total_SNV", "TMB", "\n", sep="\t", file="analyses/old_output.tmbs.tsv")
  cat(sample_name, tmb$snvs, tmb$tmb_snvs, "\n", sep="\t", append=T, file="analyses/old_output.tmbs.tsv")
}
# signatures
if (file.exists("analyses/old_output.sigs.tsv")){
  tmp_sigs = rep(0, 30)
  names(tmp_sigs) = names(etio2)
  tmp_sigs[ rownames(linear_decomp_mt_sig_legacy_30_norm) ] = linear_decomp_mt_sig_legacy_30_norm[,1]
  cat(sample_name, tmp_sigs,  "\n", sep="\t", append=T, file="analyses/old_output.sigs.tsv")
} else {
  cat("Signature", gsub("COSMIC_","Signature\\.", names(etio2)), "\n", sep="\t", file="analyses/old_output.sigs.tsv")
  cat("Etiology", etio2, "\n", sep="\t", append=T, file="analyses/old_output.sigs.tsv")
  tmp_sigs = rep(0, 30)
  names(tmp_sigs) = names(etio2)
  tmp_sigs[ rownames(linear_decomp_mt_sig_legacy_30_norm) ] = linear_decomp_mt_sig_legacy_30_norm[,1]
  cat(sample_name, tmp_sigs,  "\n", sep="\t", append=T, file="analyses/old_output.sigs.tsv")
}

##############################
# Plot Tumor Mutation Burden #
##############################

cols_tmb = rev(stepped(4))

# plot tumor mutation burden

pl_tmb = ggplot(tmb, aes(y=tumor)) +
  geom_point(aes(x=tmb_snvs)) +
  geom_rect(aes(xmin=1,xmax=10,ymin=0, ymax=2), fill=cols_tmb[1]) +
  geom_rect(aes(xmin=10,xmax=100,ymin=0, ymax=2), fill=cols_tmb[2]) +
  geom_rect(aes(xmin=100,xmax=1000,ymin=0, ymax=2), fill=cols_tmb[3]) +
  geom_vline(xintercept=c(1,10,100), linetype=2) +
  annotate(geom="line", x=rep(tmb$tmb_snvs,2), y=c(0,0.8), arrow=arrow(length=unit(0.30,"cm"), ends="first", type = "closed")) +
  annotate(geom="text", y=1.8, x=c(1.1,11,110), label=c("Low","Hypermutant","Ultrahypermutant"), hjust=0 ) +
  geom_label(aes(x=tmb_snvs, label=paste(tmb_snvs,"(mt/mb)"))) +
  scale_x_log10(expand=c(0,0),
    limits=c(0.5, 5000),
    breaks=c(1,10,100,1000),
    labels=c("1","10","100","1,000")) +
  annotation_logticks(sides="b") +
  labs(y=TUMOR, x="mutations per million nucleotides\n[Single Nucleotide Variants (SNV), log10 scale]",
       title="Tumor Mutation Burden") +
  background_grid(major="x") +
  theme(
    panel.background=element_blank(),
    plot.background=element_rect(fill="grey95", color="black"),
    axis.text.y=element_blank(),
    axis.text.x=element_text(size=10),
    axis.ticks.y=element_blank(),
    axis.line.x=element_line(),
    plot.title=element_text(face="bold", size=10),
    plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
  )

##################################
# Plot COSMIC signatures         #
##################################

# set colors
cols_cosmic2 = rep("", 30)
cols_cosmic2[ grep("DNA mismatch", etio2) ] = stepped3(12)[9:12]
cols_cosmic2[ grep("polymerase", etio2) ] = stepped3(8)[5:6]
cols_cosmic2[ grep("eami", etio2) ] = stepped3(4)[1:3]
cols_cosmic2[ grep("Unknown", etio2) ] = stepped3(20)[18]
cols_cosmic2[ grep("unknown", etio2) ] = stepped3(20)[18]
cols_cosmic2[ grep("tobacco", etio2) ] = stepped2(12)[9:10]
cols_cosmic2[ grep("exposures* to a", etio2) ] = stepped2(16)[13:15]
cols_cosmic2[ grep("UV exposure", etio2) ] = stepped3(16)[13]
cols_cosmic2[ grep("DNA-DSB repair by HR", etio2) ] = stepped2(20)[19]
names(cols_cosmic2) = names(etio2)

prep_caption_ld=paste(paste(etio2[(linear_decomp_cosmic %>% filter(cosmic_db == "v2"))$cosmic_signature], " (",(linear_decomp_cosmic %>% filter(cosmic_db == "v2"))$cosmic_signature, ")", sep=""), collapse="\n")

# break COSMIC_1 etiology
if (length(grep("COSMIC_1", prep_caption_ld)) != 0){
  prep_caption_ld = sub(" \\(","\n\\(", prep_caption_ld)
}

pl_ld_v2 = ggplot(linear_decomp_cosmic %>% filter(cosmic_db == "v2"),
    aes(x="", y=round(contribution_proportion * 100,1), fill=cosmic_signature)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(contribution_proportion * 100,0),"%")), position = position_stack(vjust=0.5)) +
  coord_polar("y", start=0) +
  labs(title="Contribution from COSMIC Signature\n(legacy db.v2)",
       caption=prep_caption_ld) +
  scale_fill_manual(values=cols_cosmic2[(linear_decomp_cosmic %>% filter(cosmic_db == "v2"))$cosmic_signature]) +
  theme_void() +
  theme(
    legend.title=element_blank(),
    plot.title=element_text(face="bold", size=10),
    plot.caption=element_text(hjust=0),
    plot.background=element_rect(fill="grey95", color="black"),
    plot.margin=unit(c(0.5,0.5,1,0.5), "cm")
  )

# plot cosine similarity (bNMF)

prep_caption_bnmf=paste(paste(etio2[(bnmf_cosmic2 %>% filter(cosmic_db == "v2"))$cosmic_signature], " (",(bnmf_cosmic2 %>% filter(cosmic_db == "v2"))$cosmic_signature, ")", sep=""), collapse="\n")

# break COSMIC_1 etiology
if (length(grep("COSMIC_1", prep_caption_bnmf)) != 0){
  prep_caption_bnmf = sub(" \\(","\n\\(", prep_caption_bnmf)
}

pl_bnmf_v2 = ggplot(bnmf_cosmic2 %>% mutate(cosmic_signature=factor(cosmic_signature, levels=cosmic_signature[order(cosine_similarity)])),
    aes(y=cosmic_signature, x=round(cosine_similarity * 100,0), fill=cosmic_signature)) +
  geom_bar(stat="identity", show.legend=F) +
  geom_text(aes(label=paste0(round(cosine_similarity * 100,0),"%")), hjust=1.1) +
  scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
  labs(x="cosine similarity (%)", title="Similarity to COSMIC Signature\n(legacy db.v2)",
       caption=prep_caption_bnmf) +
  scale_fill_manual(values=cols_cosmic2[(bnmf_cosmic2 %>% filter(cosmic_db == "v2"))$cosmic_signature]) +
  background_grid(major="xy", color.major="grey85") +
  theme(
    panel.background=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(angle=45, vjust=0.5),
    axis.ticks.y=element_blank(),
    axis.line.x=element_line(),
    legend.title=element_blank(),
    plot.title=element_text(face="bold", size=10),
    plot.caption=element_text(hjust=0),
    plot.background=element_rect(fill="grey95", color="black"),
    plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
  )

###################################################
# Process/plot germline and somatic RRD mutations #
###################################################

# set RRD genes to be tested
RRD_genes = c("MLH1","MSH2","MSH6","PMS2","POLD1","POLD2","POLD3","POLD4","POLE","POLE2")
RRD_genes_num = seq_along(RRD_genes)
names(RRD_genes_num) = RRD_genes

# set ensembl database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
RRD_ensembl = getBM(attributes=c('chromosome_name', 'start_position', 'end_position','exon_chrom_start','exon_chrom_end'),
  filters="hgnc_symbol",
  values=RRD_genes, mart=ensembl)
RRD_ensembl$gene = RRD_genes[ factor(RRD_ensembl$start_position, levels=as.character(unique(RRD_ensembl$start_position))) ]
RRD_ensembl = RRD_ensembl %>% mutate(exon_start = exon_chrom_start - start_position) %>% mutate(exon_end = exon_chrom_end - start_position)
RRD_ensembl$y = seq_along(RRD_genes)[ factor(RRD_ensembl$gene, levels=RRD_genes) ]
RRD_transcript_sizes = RRD_ensembl %>% group_by(gene, y) %>% summarize(start=1, end=sum(exon_end-exon_start))

# plot germline and somatic mutations RRD

germline_file = paste("analyses/",sample_name,".germline_MMR_mutations.genes.csv",sep="")
if (file.size(germline_file) != 0){
  germline_mmr = read.csv(germline_file, head=F, stringsAsFactors=F)
  colnames(germline_mmr) = c("tumor","normal","chromosome","position","reference","alternate","gt","frequency","gene","class","type","transcript","strand","position_tr","nt_change","aa_change")
  # add transcript base position
  germline_mmr$y = RRD_genes_num[germline_mmr$gene]
  germline_mmr$genotype = "het"
  germline_mmr$genotype[ grep("1[\\/\\|]1", germline_mmr$gt) ] = "hom"
  germline_mmr$position_tr = as.numeric(gsub("_.*","", germline_mmr$position_tr))
}

# mutation colors
cols_gt = rev(stepped(5)[c(1,5)])

# plot germline variants

if (length(ls(pattern="germline_mmr")) != 0){
  pl_germl = ggplot(RRD_transcript_sizes, aes(y=y)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=y-0.1, ymax=y+0.1), fill="grey") +
    geom_point(data=germline_mmr %>% filter(!(is.na(position_tr) | class == "SILENT")),
      aes(y=y, x=position_tr, color=factor(genotype, levels=c("het","hom"))), shape="|", size=3, show.legend=T) +
    #geom_text(data=germline_mmr %>% filter(!(is.na(position_tr) | class == "SILENT")), aes(x=position_tr, label=sub("[a-z]\\.","",aa_change), color=genotype), size=2.5, hjust=0, vjust=0, angle=45, nudge_y=0.2, show.legend=F) +
    geom_text_repel(data=germline_mmr %>% filter(!(is.na(position_tr) | class == "SILENT")),
      aes(x=position_tr, label=sub("[a-z]\\.","", aa_change), color=genotype),
      direction="x",
      nudge_x=500,
      nudge_y=0.1,
      size=2.5,
      hjust=0,
      vjust=1,
      angle=90,
      max.overlaps = Inf,
      segment.size = 0.2,
      max.iter = 1e4, max.time = 1,
      show.legend=F) +
    scale_y_continuous(breaks=seq_along(RRD_genes), limits=c(0.5,11.5), labels=RRD_genes) +
    scale_x_continuous(expand=c(0,0), labels=scales:::comma) +
    scale_color_manual(values=cols_gt, limits=c("het","hom")) +
    guides(colour=guide_legend(override.aes=list(shape=15, size=2))) +
    labs(title="Germline RRD variants", x="transcript position (base pair)") +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size=10),
      axis.text.y=element_text(angle=45, vjust=0.5),
      axis.ticks.y=element_blank(),
      axis.line.x=element_line(),
      legend.position=c(0.8,0.95),
      legend.direction="horizontal",
      #legend.position="none",
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      plot.title=element_text(face="bold", size=10),
      plot.background=element_rect(fill="grey95", color="black"),
      plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
    )
} else {
  pl_germl = ggplot(RRD_transcript_sizes, aes(y=y)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=y-0.1, ymax=y+0.1), fill="grey") +
    #geom_point(data=germline_mmr %>% filter(class == "MISSENSE"), aes(y=y, x=position_tr, color=genotype), shape="|", size=3) +
    #geom_text(data=germline_mmr %>% filter(class == "MISSENSE"), aes(x=position_tr, label=sub("[a-z]\\.","",aa_change), color=genotype), size=2.5, hjust=0, vjust=0, angle=45, nudge_y=0.2, show.legend=F) +
    scale_y_continuous(breaks=seq_along(RRD_genes), limits=c(0.5,11.5), labels=RRD_genes) +
    scale_x_continuous(expand=c(0,0), labels=scales:::comma) +
    scale_color_manual(values=cols_gt) +
    labs(title="Germline RRD variants", x="transcript position (base pair)") +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size=10),
      axis.text.y=element_text(angle=45, vjust=0.5),
      axis.ticks.y=element_blank(),
      axis.line.x=element_line(),
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.box.background=element_rect(color="black"),
      plot.title=element_text(face="bold", size=10),
      plot.background=element_rect(fill="grey95", color="black"),
      plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
      legend.position="none"
    )
}

# load somatic variants, potential driver mutations
somatic_file = paste("analyses/",sample_name,".somatic_MMR_mutations.genes.csv",sep="")
if (file.size(somatic_file) != 0){
  somatic_mmr = read.csv(somatic_file, head=F, stringsAsFactors=F)
  colnames(somatic_mmr) = c("tumor","normal","chromosome","position","reference","alternate","gt","frequency","gene","class","type","transcript","strand","position_tr","nt_change","aa_change")
  # add transcript base position
  somatic_mmr$y = RRD_genes_num[somatic_mmr$gene]
  somatic_mmr$genotype = "het"
  somatic_mmr$genotype[ grep("1[\\/\\|]1", somatic_mmr$gt) ] = "hom"
  somatic_mmr$position_tr = as.numeric(gsub("_.*","", somatic_mmr$position_tr))
}

# from MAF
# somatic_maf = maf@data %>% filter(Hugo_Symbol %in% RRD_genes) %>%
#   dplyr::select(Hugo_Symbol, Chromosome, Start_Position, cDNA_Change, Protein_Change, tumor_f)
# colnames(somatic_maf) = c("gene","chromosome","position","base_change","aa_change","tumor_f")
# somatic_maf$position_tr = as.numeric(gsub("\\D+","", somatic_maf$base_change))
# somatic_maf$y = which(RRD_genes %in% somatic_maf$gene)

# plot somatic variants

if (length(ls(pattern="somatic_mmr")) != 0){
  pl_soma = ggplot(RRD_transcript_sizes, aes(y=y)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=y-0.1, ymax=y+0.1), fill="grey") +
    geom_point(data=somatic_mmr %>% filter(!(is.na(position_tr) | class == "SILENT")),
      aes(y=y, x=position_tr, color=frequency), shape="|", size=3, show.legend=T) +
    #geom_text(data=somatic_mmr %>% filter(!(is.na(position_tr) | class == "SILENT")), aes(x=position_tr, label=sub("[a-z]\\.","",aa_change)), size=2.5, hjust=0, vjust=0, angle=45, nudge_y=0.2, show.legend=F) +
    geom_text_repel(data=somatic_mmr %>% filter(!(is.na(position_tr) | class == "SILENT")),
      aes(x=position_tr, label=sub("[a-z]\\.","", aa_change), color=frequency),
      direction="x",
      nudge_x=500,
      nudge_y=0.1,
      size=2.5,
      hjust=0,
      vjust=1,
      angle=90,
      max.overlaps = Inf,
      segment.size = 0.2,
      max.iter = 1e4, max.time = 1,
      show.legend=F) +
    scale_y_continuous(breaks=seq_along(RRD_genes), limits=c(0.5,11.5), labels=RRD_genes) +
    scale_x_continuous(expand=c(0,0), labels=scales:::comma) +
    scale_color_gradient(limits=c(0,1), breaks=c(0,0.5,1), labels=c("0","0.5","1"), name="freq.") +
    labs(title="Somatic RRD variants", x="transcript position (base pair)") +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size=10),
      axis.text.y=element_text(angle=45, vjust=0.5),
      axis.ticks.y=element_blank(),
      axis.line.x=element_line(),
      #legend.position="none",
      legend.position=c(0.8,0.95),
      legend.direction="horizontal",
      #legend.title=element_blank(),
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.box.background=element_blank(),
      plot.title=element_text(face="bold", size=10),
      plot.background=element_rect(fill="grey95", color="black"),
      plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
    )
} else {
  pl_soma = ggplot(RRD_transcript_sizes, aes(y=y)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=y-0.1, ymax=y+0.1), fill="grey") +
    #geom_point(data=somatic_mmr %>% filter(class == "MISSENSE"), aes(y=y, x=position_tr, color=frequency), shape="|", size=3) +
    #geom_text(data=somatic_mmr %>% filter(class == "MISSENSE"), aes(x=position_tr, label=sub("[a-z]\\.","",aa_change)), size=2.5, hjust=0, vjust=0, angle=45, nudge_y=0.2, show.legend=F) +
    scale_y_continuous(breaks=seq_along(RRD_genes), limits=c(0.5,11.5), labels=RRD_genes) +
    scale_x_continuous(expand=c(0,0), labels=scales:::comma) +
    scale_color_gradient(limits=c(0,1)) +
    labs(title="Somatic RRD variants", x="transcript position (base pair)") +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size=10),
      axis.text.y=element_text(angle=45, vjust=0.5),
      axis.ticks.y=element_blank(),
      axis.line.x=element_line(),
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.box.background=element_rect(color="black"),
      plot.title=element_text(face="bold", size=10),
      plot.background=element_rect(fill="grey95", color="black"),
      plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
      legend.position="none"
    )
}
#####################
# Plot sample name  #
#####################

pl_name = ggplot(data.frame(x=1,y=5, z=sample_name), aes(x=x, y=y, label=sample_name)) +
  geom_text(size=7, fontface="bold") +
  theme_void() +
  scale_y_continuous(expand=c(0,0))

#####################
# Render all plots  #
#####################

pdf(file=paste0("analyses/",sample_name,".wes_analysis.pdf"), height=10, width=9)

plot_grid(plot_grid(NULL,pl_name,NULL, nrow=1, rel_widths=c(0.2,5,3)),
  NULL,
    plot_grid(NULL,pl_germl,NULL,pl_soma,NULL, rel_widths=c(1,9,0.5,9,1), nrow=1),
      NULL,
        plot_grid(NULL, pl_tmb, NULL, rel_widths=c(1,8,1), nrow=1),
          NULL,
            plot_grid(NULL, pl_ld_v2, NULL, pl_bnmf_v2, NULL, rel_widths=c(0.1,9,0.05,9,1), nrow=1),
              NULL,
                ncol=1, rel_heights=c(0.5,0,4,0.2,2,0.2,3,0.2))
dev.off()

print(paste("Done for", sample_name))

# save R objects to disk
# save.image(file="analyses/mutational_signatures_as_R_object.Rdata")
