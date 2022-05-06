#!/bin/Rscript

.libPaths('/hpf/largeprojects/tabori/shared/software/R_libs/4.1.2/')

library(MutationalPatterns)
library(NMF)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(pals)

#################
# read/set args #
#################

# get args
args = base::commandArgs(trailingOnly = TRUE)
mode = args[1]
sample_name = args[2]
vcfpath = args[3]
organism = args[4]
reference = args[5]

# ref genome
if (organism == "human"){
  ref_genome = paste0("BSgenome.Hsapiens.UCSC.", reference)
} else if (organism == "mouse"){
  ref_genome = paste0("BSgenome.Mmusculus.UCSC.", reference)
}

# get COSMIC signature database
sig_loc = paste0("/hpf/largeprojects/tabori/shared/resources/cosmic_v3.2/",reference,"/cosmic_db_v3.2_signature_matrix.txt")
signatures3 = read.table(sig_loc, sep="\t", head=T, row.names=1)
# get etiologies
etio_loc = paste0("/hpf/largeprojects/tabori/shared/resources/cosmic_v3.2/",reference,"/cosmic_db_v3.2_signatures_and_etiologies.txt")
etio3tab = read.table(etio_loc, sep="\t")


# load genome ref lib
library(ref_genome, character.only=T)

# tumor - normal
TUMOR = sub("__.*","",sample_name)
NORMAL = sub(".*__","",sample_name)

##################
# Get Signatures #
##################

# read vcf as genomic ranges
message("Reading VCF as GenomicRanges...")
grl = read_vcfs_as_granges(vcfpath, sample_name, ref_genome)

# make tally of maf objects
message("Extracting SBS96 trinucleotide contexts...")
mut_mat = mut_matrix(vcf_list = grl, ref_genome = ref_genome)
pseudo_count_prop = 0.0001
message(paste("Adding a small pseudocount proportion:", pseudo_count_prop))
mut_mat = mut_mat + pseudo_count_prop

# extract de novo signatures
#estimate = nmf(mut_mat, rank = 2:5, method = "brunet",
#                nrun = 10, seed = 123456, .opt = "v-p")

# reorder signatures matrix
signatures = as.matrix(signatures)
signatures = signatures[rownames(mut_mat),]

# minor edits to the table
etio3 = etio3tab[,2]
names(etio3) = etio3tab[,1]
etio3[ grep("Poli", etio3) ] = sub("Poli","Poly",etio3[ grep("Poli", etio3) ])

# signature refitting/fitting to COSMIC v3.2
# fit_res = fit_to_signatures(mut_mat, signatures)

# Decreasing this number will make the refitting less strict,
# while increasing it will make the refitting more strict.
# Trying out different values can sometimes be useful to achieve the best results.
# delta = 0.004
strict_refit = fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
strict_refit$fit_res$contribution_prop = strict_refit$fit_res$contribution / colSums(strict_refit$fit_res$contribution)

# prep table for plotting
cosmic_decomp = as.data.frame(strict_refit$fit_res$contribution_prop) %>%
                               rename(contribution_proportion=all_of(sample_name)) %>%
                               mutate(sample=sample_name, tumor=TUMOR, normal=NORMAL) %>%
                               mutate(contribution_raw=strict_refit$fit_res$contribution[,sample_name]) %>%
                               mutate(cosmic_signature=rownames(.)) %>%
                               mutate(etiology=etio3[cosmic_signature])


##########
## Plots #
##########

# raw/reconstructed signatures

# set classic colors
mutcols = c("#02bced","#010101","#c8332f","#cac8c9", "#a0ce62","#ecc6c5")
# prep data frame
sigs.df = data.frame(raw = mut_mat[,1] - pseudo_count_prop,
                     extracted = strict_refit$fit_res$reconstructed[,1]) %>%
                     mutate(normalized = extracted / sum(extracted)) %>%
                     mutate(mutations = rownames(.))
sigs.df$to = gsub("[ATCG]\\[|\\][ATCG]","",sigs.df$mutations)
sigs.df$from = gsub("\\[[TC]>|\\]","",sigs.df$mutations)

# barplot with proportionls of SBS96 mutational patterns
p1 = ggplot(sigs.df, aes(x=from, y=normalized, fill=to)) +
  geom_bar(stat="identity", width=0.8, show.legend=F) +
  facet_wrap(~to, nrow=1, scale="free_x") +
  scale_fill_manual(values=mutcols) +
  scale_y_continuous(expand=c(0,0)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=6),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.background=element_blank(),
    panel.grid=element_blank(),
    panel.grid.major.y=element_line(size=0.5, color="grey"),
    axis.line.y=element_line()
  )

# set colors for cosmic 3.2
cols_cosmic3 = rep("", length(etio3))
names(cols_cosmic3) = etio3
# unknown data
cols_cosmic3[ grep("^Unknown$", etio3) ] = stepped3(20)[18]
cols_cosmic3[ grep("Possible sequencing artefacts", etio3) ] = stepped3(20)[17]
# MMR
cols_cosmic3[ grep("DNA mismatch", etio3) ] = colorRampPalette(stepped3(12)[9:12][c(1,4)])(7)
cols_cosmic3[ grep("Polymerase|POLD1 proofreading", etio3) ] = colorRampPalette(stepped3(8)[5:8][c(1,4)])(5)
# clocklike
cols_cosmic3[ grep("clock-like signature", etio3) ] = stepped3(4)[1:2]
# UV
cols_cosmic3[ grep("Ultraviolet", etio3) ] = stepped3(16)[13:16]
# tobacco
cols_cosmic3[ grep("Tobacco", etio3) ] = stepped(8)[5:7]
# treatment
cols_cosmic3[ grep("treatment", etio3) ] = colorRampPalette(stepped2(4)[c(1,4)])(7)
# homologous recombination DNA damage repair
cols_cosmic3[ grep("Defective homologous recombination DNA damage repair", etio3) ] = stepped2(12)[10]
# APOBEC
cols_cosmic3[ grep("APOBEC|cytidine deaminase", etio3) ] = stepped2(16)[13:16]
# carcinogens
cols_cosmic3[ grep("exposure", etio3)[5:length(grep("exposure", etio3))] ] = colorRampPalette(stepped2(20)[c(1,4)])(5)
# DNA base excision repair
cols_cosmic3[ grep("DNA base excision", etio3) ] = stepped2(8)[5:6]
# reactive oxigen
cols_cosmic3[ grep("reactive oxygen", etio3) ] = stepped(4)[2]
# indirect UV
cols_cosmic3[ grep("Indirect effect of ultraviolet", etio3) ] = stepped(16)[1]
# rename
cols_cosmic3.2 = cols_cosmic3
names(cols_cosmic3.2) = names(etio3)

# prep captionf or figure
# break long signature names/etiologies
caption = etio3
caption[nchar(etio3) > 60] = sub("bacteria ","bacteria\n   ",sub("and ","and\n   ",sub("of ","of\n   ",etio3[nchar(etio3) > 60])))
caption = sub("^","*", caption)
prep_caption_ld=paste(paste(caption[(cosmic_decomp %>% filter(contribution_proportion > 0))$cosmic_signature], " (",(cosmic_decomp %>%
  filter(contribution_proportion > 0))$cosmic_signature, ")", sep=""), collapse="\n")

p2 = ggplot(cosmic_decomp %>%
                  filter(contribution_proportion > 0) %>%
                  mutate(cosmic_signature=as.character(cosmic_signature)),
    aes(x="", y=round(contribution_proportion * 100,1), fill=cosmic_signature)) +
  geom_bar(stat="identity", show.legend=T, width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label=paste0(round(contribution_proportion * 100,0),"%")), position=position_stack(vjust=0.5), color="white") +
  #scale_y_continuous(expand=c(0,0), breaks=seq(0,100,10), limits=c(0,30)) +
  # scale_y_discrete(limits=(cosmic_decomp %>%
  #                   filter(contribution_proportion > 0) %>%
  #                   mutate(cosmic_signature=as.character(cosmic_signature)) %>%
  #                   arrange(contribution_proportion))$cosmic_signature) +
  # labs(title="COSMIC Signature Analysis\n(SBS-96 v3.2, hg38)",
  #      caption=prep_caption_ld, x="contribution to mutational signature (%)") +
  scale_fill_manual(values=cols_cosmic3.2[(cosmic_decomp %>% filter(contribution_proportion > 0))$cosmic_signature]) +
  background_grid(major="xy", color.major="grey85") +
  theme_void()
  # theme(
  #   panel.background=element_blank(),
  #   axis.title.y=element_blank(),
  #   axis.text.x=element_text(size=10),
  #   axis.text.y=element_text(angle=45, vjust=0.5),
  #   axis.ticks.y=element_blank(),
  #   axis.line.x=element_line(),
  #   legend.title=element_blank(),
  #   plot.title=element_text(face="bold", size=10),
  #   plot.caption=element_text(hjust=0),
  #   plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")
  # )

# plot grid
plot_grid(p1, p2, nrow=1, rel_widths=c(6,3))

# save image
ggsave(file="MD1553_signatures.png", bg="white")
