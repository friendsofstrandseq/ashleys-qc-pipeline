library(edgeR)
library(ggplot2)

# fetch arguments
args = commandArgs(trailingOnly = T)

filter <- ifelse('filter' %in% args, TRUE, FALSE)
chosen_transform <- ifelse('anscombe' %in% args, 'anscombe', 
                           ifelse('laubschner' %in% args, 'laubschner', 'anscombe'))
plot <- ifelse(is.na(args[2]), FALSE, TRUE)

# open count file
counts_raw <- data.table::fread(snakemake@input[["counts"]])

# force cell column to factor
counts_raw$cell <- as.factor(counts_raw$cell)

# fuse bin coordinates
counts_raw$bin <- paste(counts_raw$chrom, counts_raw$start, counts_raw$end, sep='_')

# add tot counts
counts_raw$tot_count <- counts_raw$c + counts_raw$w

# sum of counts per cell
sum_counts <- aggregate(counts_raw$tot_count, by= list(counts_raw$cell), FUN = sum)
colnames(sum_counts) <- c('cell', 'sum')


# centromere read spikes exclusion
exclude_centromeres <- function(counts_raw, exclusion_list) {
  exclusion_list <- GenomicRanges::GRanges(seqnames = exclusion_list$chrom, ranges = IRanges::IRanges(start = exclusion_list$start, end=exclusion_list$end))
  counts_iranges <- GenomicRanges::GRanges(seqnames = counts_raw$chrom, ranges = IRanges::IRanges(start = counts_raw$start, end=counts_raw$end))

  overlaps <- data.table::as.data.table(GenomicRanges::findOverlaps(exclusion_list, counts_iranges))
  counts <- counts_raw
  counts <- counts[-overlaps$subjectHits,]
  return(counts)
}

if (!is.na(args[3])) {
  message('blacklisting...')
  counts <- exclude_centromeres(counts_raw, data.table::fread(args[3]))
} else {
  counts <- counts_raw
}

# convert to bin count matrix
to_matrix <- function(counts) {
  mat_tot <- reshape2::dcast(counts, bin ~ cell, value.var = "tot_count")
  rownames(mat_tot) <- mat_tot$bin
  mat_tot <- mat_tot[,2:ncol(mat_tot)]
  
  return(mat_tot)
}

mat_tot <- to_matrix(counts)

generate_dgelist <- function(mat_tot, filter) {
  
  # create the DGEList object
  # all samples belong to group 1
  y.raw <- DGEList(as.matrix(mat_tot), group = rep(1, ncol(mat_tot)))
  
  # filter out bins with too low counts to be informative
  if (filter) {
    keep.exprs <- edgeR::filterByExpr(y.raw, group=y.raw$samples$group)
    y <- y.raw[keep.exprs, keep.lib.sizes=FALSE]
    message(paste('filtering out', (dim(y.raw)[1]-dim(y)[1]), 'bins out of', dim(y.raw)[1], "due to low count"))
  } else {
    y <- y.raw
    message('no filtering')
  }
  
  # calculate library size and composition normalization
  y <- calcNormFactors(y)
  
  return(y)
}

estimate_dispersion <- function(y) {
  
  # estimate common dispersion
  # default settings for DGEList objects
  message('Estimating common dispersion with edgeR...')
  disp <- estimateDisp(y, design=NULL, prior.df=NULL, trend.method="locfit", tagwise=TRUE,
                       span=NULL, min.row.sum=5, grid.length=21, grid.range=c(-10,10), robust=FALSE, 
                       winsor.tail.p=c(0.05,0.1), tol=1e-06)
  phi <- disp$common.dispersion
  message(paste("phi =",phi))
  write(paste(args[1], phi, sep="\t"), './log.txt', append=TRUE)
  
  return(phi)
}

y <- generate_dgelist(mat_tot, filter)
phi <- estimate_dispersion(y)


# VARIANCE STABILIZING TRANSFORMATION

# TRANSFORMS
# for negative binomial distributions
# Anscombe, 1948
# Laubschner, 1961
anscombe_transform <- function(x, phi) {
  a <- x + (3/8)
  b <- (1/phi) - (3/4)
  c <- sqrt(a/b)
  y <- asinh(c)
  return(y)
}
laubscher_transform <- function(x, phi) {
  a <- sqrt(phi)
  b <- asinh( sqrt( x/phi ) )
  c <- sqrt( phi-1 )
  d <- anscombe_transform(x, phi)
  y <- a*b + c*d
  return(y)
}
transform_list <- list('anscombe'= anscombe_transform, 'laubscher'= laubscher_transform)
transform <- transform_list[[chosen_transform]]

message(paste("Transforming data with", chosen_transform, "VST"))
#ans <- transform(y$counts, phi)
counts$w_corr <- transform(counts$w, phi)
counts$c_corr <- transform(counts$c, phi)
counts$w_corr <- counts$w_corr - min(counts$w_corr)
counts$c_corr <- counts$c_corr - min(counts$c_corr)
counts$tot_count_corr <- counts$w_corr + counts$c_corr

# scaling counts to original range
k <- mean(counts$tot_count)/mean(counts$tot_count_corr)
scaled <- counts$tot_count_corr * k
counts$w <- (counts$w_corr / counts$tot_count_corr) * scaled
counts$c <- (counts$c_corr / counts$tot_count_corr) * scaled
counts$w[is.na(counts$w)] <- 0
counts$c[is.na(counts$c)] <- 0
counts$tot_count <- scaled

df <- data.table::data.table(counts[,c('chrom', 'start', 'end', 'sample', 'cell', 'w', 'c', 'class', 'tot_count')])

data.table::fwrite(df, snakemake@output[["counts_vst"]])

if (plot) {
  library(ggplot2)
  library(ggpubr)
  
  p1 <- ggplot(counts_raw, aes(x=tot_count)) + 
    geom_histogram(bins=256) +
    ggtitle("raw") +
    xlab('read count') +
    ylab('bin count')
  
  p2 <- ggplot(df, aes(x=tot_count)) + 
    geom_histogram(bins=256) +
    ggtitle(paste(chosen_transform, "VST")) +
    xlab('read count') +
    ylab('bin count')
  
  corr_plot <- ggarrange(p1, p2)
  
  ggsave(args[2], corr_plot, width=12, height=6)
}
