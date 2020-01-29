#!/bin/env Rscript
#
# combines feature associations from multiple sources
#
suppressMessages(library(arrow))
suppressMessages(library(metap))
suppressMessages(library(tidyverse))

# load individual dataset associations
infiles <- unlist(snakemake@input)

wts_list <- lapply(infiles, read_parquet)
names(wts_list) <- infiles

# reorder associations list so that the entry with the most rows appears first
wts_list <- wts_list[order(-rank(unlist(lapply(wts_list, nrow)), ties.method = 'first'))]

# get feature id field (e.g. "symbol" or "gene set");
# should be the same for all datasets
id_field <- colnames(wts_list[[1]])[1]

# sanity check; make sure no duplicated keys exist
no_dups <- sapply(wts_list, function(x) {
  nrow(x) == length(unique(pull(x, id_field)))
})
if (!all(no_dups)) {
  stop("Duplicated feature ids found!")
}

# combine individual associations into a single tibble
wts <- wts_list %>%
  reduce(left_join, by = id_field)

# create a matrix version of the weights without the id column
wts_mat <- wts %>%
  select(-id_field) %>%
  as.matrix()

# count number of missing values for each gene or gene set
num_missing <- apply(wts_mat, 1, function(x) {
  sum(is.na(x))
})
num_present <- ncol(wts_mat) - num_missing

# wrap sumz and sumlog functions to handle missing values and insufficient data cases
sumlog_wrapper <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) {
    NA
  } else {
    as.numeric(sumlog(x)$p)
  }
}

sumz_wrapper <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) {
    NA
  } else {
    as.numeric(sumz(x)$p)
  }
}

# combine p-values across datasets and covariates
res <- data.frame(
  pull(wts, id_field),
  mean_pval   = apply(wts_mat, 1, mean, na.rm = TRUE),
  median_pval = apply(wts_mat, 1, median, na.rm = TRUE),
  min_pval    = apply(wts_mat, 1, min, na.rm = TRUE),
  max_pval    = apply(wts_mat, 1, max, na.rm = TRUE),
  sumlog_pval = suppressWarnings(apply(wts_mat, 1, sumlog_wrapper)),
  sumz_pval   = suppressWarnings(apply(wts_mat, 1, sumz_wrapper)),
  mean_score   = apply(wts_mat, 1, function(x) { mean(-log10(pmax(x, 1E-20)), na.rm=TRUE) }),
  median_score = apply(wts_mat, 1, function(x) { median(-log10(pmax(x, 1E-20)), na.rm=TRUE) }),
  num_present,
  num_missing
)
colnames(res)[1] <- id_field

res <- res %>%
  arrange(sumlog_pval)

# store results
write_feather(wts, snakemake@output[['scores']])
write_feather(res, snakemake@output[['combined_scores']])
