# fetch arguments
# args <- commandArgs(trailingOnly = T)

library(data.table)
# args[1] <- '/data/r-workspace/strandseq_utils/H3JNHAFX3_MNIphotoconverted_200000_fixed.norm.txt.gz'
# args[2] <- '/data/r-workspace/strandseq_utils/GC_matrix_200000.txt'
print(snakemake@params[["gc_matrix"]])
# open files
counts <- data.table::fread(snakemake@input[["counts_scaled"]], header = T)
# counts <- data.table::fread(args[1], header = T)
GC_matrix <- data.table::fread(snakemake@params[["gc_matrix"]], header = T)
print(GC_matrix)
# GC_matrix <- data.table::fread(args[2], header = T)
save_path <- snakemake@output[["counts_scaled_gc"]]
# save_path <- args[3]

# check GC plots
# plot <- ifelse(is.na(args[4]), FALSE, TRUE)
plot <- TRUE
if (plot) {
  # import libraries
  library(ggplot2)
  library(ggpubr)
}


# check files
if (!all(c("cell", "chrom", "start", "end", "w", "c") %in% colnames(counts))) {
  message("count file does not contain required columns: 'cell', 'chrom', 'start', 'end', 'w', 'c'")
  message("Usage: Rscript GC_correction.R count-file.txt.gz gc-matrix.txt output.txt.gz")
  stop()
}
if (!all(c("chrom", "start", "end", "GC%") %in% colnames(GC_matrix))) {
  message("GC_matrix file does not contain required columns: 'chrom', 'start', 'end', 'GC%'")
  message("Usage: Rscript GC_correction.R count-file.txt.gz gc-matrix.txt output.txt.gz")
  stop()
}
if (!(all(unique(counts$chrom) %in% unique(GC_matrix$chrom)) &
  all(unique(counts$start) %in% unique(GC_matrix$start)) &
  all(unique(counts$end) %in% unique(GC_matrix$end)))) {
  message("bin features ('crhom', 'start', 'end') do not match between count file and GC matrix")
  message("make sure to choose files with identical bin sizes")
}


# green light message
# message(paste("\ncount file:", args[1]))
# message(paste("GC matrix file:", args[2]))
# message(paste("savepath:", args[3]))
message("preprocessing...\n")


#################
# Preprocessing #
#################
min_reads <- snakemake@params[["gc_min_reads"]]
n_subsample <- snakemake@params[["gc_n_subsample"]]
# min_reads <- ifelse(is.na(args[5]), 5, args[5])
# n_subsample <- ifelse(is.na(args[6]), 1000, args[6])


# force cell column to factor
counts$cell <- as.factor(counts$cell)

# convert strandseq count file to count matrix
counts$tot_count <- counts$c + counts$w

# ######################
# # GC bias correction #
# ######################


counts <- merge(counts, GC_matrix[, c("chrom", "start", "GC%")], by = c("chrom", "start"), all.x = T)

# filter data for subsampling
c <- counts[counts$tot_count >= min_reads]
if (dim(c)[[1]] == 0) {
  stop(paste("there are no bins with more than", min_reads, "reads"))
}
c$`GC%` <- as.numeric(c$`GC%`)
c$log_count_norm <- log(c$tot_count) - log(median(c$tot_count))
not.na <- !is.na(c$`GC%`)
s <- c[not.na]

# subsample from quantiles
s$GC_bin <- cut(s$`GC%`, breaks = c(quantile(s$`GC%`, probs = seq(0, 1, by = 1 / 10))), labels = seq(1, 10, by = 1), include.lowest = TRUE)
# s$GC_bin <- cut(s$`GC%`, breaks = c(quantile(s$`GC%`, probs = seq(0, 1, by = 1 / 10))), labels = seq(0, 1, by = 1 / 10))

subsample <- data.frame()
for (i in seq(10)) {
  sbin <- s[s$GC_bin == i]
  m <- min(dim(sbin)[1], n_subsample)
  sa <- sbin[sample(nrow(sbin), size = m), ]
  subsample <- rbind(subsample, sa)
}

#############################
# lowess fit and correction #
#############################
# lowess fit
z <- lowess(subsample$`GC%`, subsample$log_count_norm)

# ################
# # SAVING PLOTS #
# ################

# adjust tot count to closest predicted GC value
idxs <- sapply(as.numeric(counts$`GC%`), FUN = function(a) {
  which.min(abs(z$x - a))
})
idxs[lapply(idxs, length) == 0] <- NA
counts$pred <- z$y[unlist(idxs)]
counts$tot_count_gc <- log(counts$tot_count / median(counts$tot_count)) - counts$pred
counts$tot_count_gc <- exp(counts$tot_count_gc) * median(counts$tot_count)

if (plot) {
  sidxs <- sapply(as.numeric(subsample$`GC%`), FUN = function(a) {
    which.min(abs(z$x - a))
  })
  sidxs[lapply(sidxs, length) == 0] <- NA


  subsample$pred <- z$y[unlist(sidxs)]
  subsample$tot_count_gc <- log(subsample$tot_count / median(subsample$tot_count)) - subsample$pred
  subsample$tot_count_gc <- exp(subsample$tot_count_gc) * median(subsample$tot_count)

  z$y2 <- exp(z$y) * median(subsample$tot_count)
  ymin <- min(cbind(subsample$tot_count, subsample$tot_count_gc))
  ymax <- max(cbind(subsample$tot_count, subsample$tot_count_gc))

  p1 <- ggplot(subsample, aes(`GC%`, tot_count)) +
    geom_point(size = 1, alpha = .2) +
    ggtitle("raw") +
    ylim(ymin, ymax) +
    xlab("GC_content") +
    ylab("read count") +
    geom_line(data = as.data.frame(z), aes(x, y2), color = "red")

  p2 <- ggplot(subsample, aes(`GC%`, tot_count_gc)) +
    geom_point(size = 1, alpha = .2) +
    ggtitle("gc corrected") +
    ylim(ymin, ymax) +
    xlab("GC content") +
    ylab("read count")

  corr_plot <- ggarrange(p1, p2)


  # save plots
  ggsave(snakemake@output[["plot"]], corr_plot, width = 12, height = 6)
}

# adjust w, c and fill NAs
counts$w <- (counts$w / counts$tot_count * counts$tot_count_gc)
counts$c <- (counts$c / counts$tot_count * counts$tot_count_gc)
counts$w[is.na(counts$w)] <- 0
counts$c[is.na(counts$c)] <- 0
counts$tot_count[is.na(counts$tot_count)] <- 0

output <- counts[, c("chrom", "start", "end", "sample", "cell", "w", "c", "tot_count", "class")]

data.table::fwrite(output, save_path)
