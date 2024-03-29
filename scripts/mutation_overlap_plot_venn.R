#!/bin/Rscript

# load lib path
.libPaths('/hpf/largeprojects/tabori/shared/software/R_libs/4.1.0/')

library(eulerr)

plot_venn = function(x, title){
  # prepare vector
  overlaps = as.vector(table(x))
  names(overlaps) = names(table(x))
  # fit euler function
  fit = euler(overlaps)
  # plot
  plot(fit,
    fills = list(fill = c("red", "steelblue4"), alpha = 0.5),
    labels = list(col = "black", font = 4),
    quantities=TRUE,
    main=title)
}

# read file
file = base::commandArgs(trailingOnly = TRUE)[1]
dat = read.csv(file)

# get which vars are SNVs and which are indels
snvs = apply(dat[c("REF","ALT")], 1, function(x) nchar(x[1]) == 1 & nchar(x[2]) == 1 )
indels = !snvs

# get sample names
samples = colnames(dat)[5:dim(dat)[2]]

# all var overlaps
overlaps_all = apply(dat[samples], 1, function(x)
  if (length(grep("0/1|0\\|1|1\\|0", x)) == length(samples)){
    return(paste(samples, collapse="&"))
  } else if (length(grep("0/1|0\\|1|1\\|0", x)) == 1 & grep("0/1|0\\|1|1\\|0", x) == 1){
    return(samples[1])
  } else if (length(grep("0/1|0\\|1|1\\|0", x)) == 1 & grep("0/1|0\\|1|1\\|0", x) == 2){
    return(samples[2])
  }
)

overlaps_all = apply(dat[samples], 1, function(x) paste(names(x[grep("0/1|0\\|1|1\\|0", x)]), collapse="&") )

# SNVs only
overlaps_snvs = apply(dat[snvs,samples], 1, function(x)
  if (length(grep("0/1|0\\|1|1\\|0", x)) == length(samples)){
    return(paste(samples, collapse="&"))
  } else if (length(grep("0/1|0\\|1|1\\|0", x)) == 1 & grep("0/1|0\\|1|1\\|0", x) == 1){
    return(samples[1])
  } else if (length(grep("0/1|0\\|1|1\\|0", x)) == 1 & grep("0/1|0\\|1|1\\|0", x) == 2){
    return(samples[2])
  }
)

# indels only
overlaps_indels = apply(dat[indels,samples], 1, function(x)
  if (length(grep("0/1|0\\|1|1\\|0", x)) == length(samples)){
    return(paste(samples, collapse="&"))
  } else if (length(grep("0/1|0\\|1|1\\|0", x)) == 1 & grep("0/1|0\\|1|1\\|0", x) == 1){
    return(samples[1])
  } else if (length(grep("0/1|0\\|1|1\\|0", x)) == 1 & grep("0/1|0\\|1|1\\|0", x) == 2){
    return(samples[2])
  }
)

# plot all 3 venn diagrams
# all
png(file=paste0(paste(samples,collapse="__"), ".all.venn.png"))
plot_venn(overlaps_all, "all")
dev.off()
# snvs
png(file=paste0(paste(samples,collapse="__"), ".snvs.venn.png"))
plot_venn(overlaps_all, "SNVs")
dev.off()
# indels
png(file=paste0(paste(samples,collapse="__"), ".id.venn.png"))
plot_venn(overlaps_all, "ID")
dev.off()
