# set arguments
input_path <- snakemake@input[["counts_scaled_gc"]]
save_path <- snakemake@output[["counts_scaled_gc_vst"]]
chosen_transform <- "anscombe"
plot <- TRUE


# PRE-PROCESSING DATA
# open count file
counts_raw <- data.table::fread(input_path)

# force cell column to factor
counts_raw$cell <- as.factor(counts_raw$cell)

# fuse bin coordinates
counts_raw$bin <- paste(counts_raw$chrom, counts_raw$start, counts_raw$end, sep = "_")

# add tot counts
counts_raw$tot_count <- counts_raw$c + counts_raw$w

counts <- counts_raw

# convert to bin count matrix
to_matrix <- function(counts) {
  # fuse bin coordinates
  counts$bin <- paste(counts$chrom, counts$start, counts$end, sep = "_")
  mat_tot <- reshape2::dcast(counts, bin ~ cell, value.var = "tot_count")
  rownames(mat_tot) <- mat_tot$bin
  mat_tot <- mat_tot[, 2:ncol(mat_tot)]

  return(mat_tot)
}

# VARIANCE STABILIZING TRANSFORMATION

# TRANSFORMS
# for negative binomial distributions
# Anscombe, 1948
# Laubschner, 1961
anscombe_transform <- function(x, phi) {
  a <- x + (3 / 8)
  b <- (1 / phi) - (3 / 4)
  c <- sqrt(a / b)
  y <- asinh(c)
  return(y)
}
laubscher_transform <- function(x, phi) {
  a <- sqrt(phi)
  b <- asinh(sqrt(x / phi))
  c <- sqrt(phi - 1)
  d <- anscombe_transform(x, phi)
  y <- a * b + c * d
  return(y)
}
transform_data <- function(counts, transform, phi) {
  cols <- colnames(counts)
  counts$w_corr <- transform(counts$w, phi)
  counts$c_corr <- transform(counts$c, phi)
  counts$w_corr <- counts$w_corr - min(counts$w_corr) # shift values to start from zero
  counts$c_corr <- counts$c_corr - min(counts$c_corr)
  counts$tot_count_corr <- counts$w_corr + counts$c_corr
  # rescale around original mean
  k <- mean(counts$tot_count) / mean(counts$tot_count_corr)
  scaled <- counts$tot_count_corr * k
  counts$w <- (counts$w_corr / counts$tot_count_corr) * scaled
  counts$c <- (counts$c_corr / counts$tot_count_corr) * scaled
  counts$w[is.na(counts$w)] <- 0
  counts$c[is.na(counts$c)] <- 0
  counts$tot_count <- scaled
  counts <- counts[, ..cols]
  return(counts)
}
transform_list <- list("anscombe" = anscombe_transform, "laubscher" = laubscher_transform)
transform <- transform_list[[chosen_transform]]

disp_score <- function(counts, transform, phi, design=NULL){
  counts$tot_count <- transform(counts$tot_count, phi)
  mat <- to_matrix(counts)
  # if multiple samples are present design matrix can be used
  if (is.null(design))
    design <- matrix(1, ncol=1, nrow=ncol(mat))
  res <- as.matrix(mat) %*% MASS::Null(design)
  rsd <- sqrt(rowMeans(res*res))
  score <- sd(rsd)/mean(rsd)
  return(score)
}

message(paste("Transforming data with", chosen_transform, "VST"))

# estimate dispersion by residual variance
opt <- optimize(disp_score, counts = counts, transform = transform, interval = c(0.00001,1))
phi <- opt$minimum
message(paste("Estimated dispersion - phi: ", phi))

# correction
corr_counts <- transform_data(counts, transform, phi)
corr_counts <- data.table::data.table(corr_counts[, c("chrom", "start", "end", "sample", "cell", "w", "c", "class", "tot_count")])

message("saving...")
data.table::fwrite(corr_counts, save_path)

if (plot) {
  library(ggplot2)
  library(ggpubr)

  p1 <- ggplot(counts_raw, aes(x = tot_count)) +
    geom_histogram(bins = 256) +
    ggtitle("raw") +
    xlab("read count") +
    ylab("bin count")

  p2 <- ggplot(corr_counts, aes(x = tot_count)) +
    geom_histogram(bins = 256) +
    ggtitle(paste(chosen_transform, "VST")) +
    xlab("read count") +
    ylab("bin count")

  corr_plot <- ggarrange(p1, p2)

  ggsave(snakemake@output[["plot"]], corr_plot, width = 12, height = 6)
}
