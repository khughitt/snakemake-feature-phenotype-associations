#!/bin/env Rscript
#
# combines feature associations from multiple sources
#
#suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# load individual dataset associations (p-values, correlations, etc.)
coef_infiles <- unlist(snakemake@input[['coefs']])
pval_infiles <- unlist(snakemake@input[['pvals']])
stat_infiles <- unlist(snakemake@input[['stats']])

coefs_list <- lapply(coef_infiles, read_csv)
names(coefs_list) <- coef_infiles

pvals_list <- lapply(pval_infiles, read_csv)
names(pvals_list) <- pval_infiles

stats_list <- lapply(stat_infiles, read_csv)
names(stats_list) <- stat_infiles

# reorder associations list so that the entry with the most rows appears first
num_rows <- unlist(lapply(pvals_list, nrow))
ind <- order(-rank(num_rows, ties.method = 'first'))

coefs_list <- coefs_list[ind]
pvals_list <- pvals_list[ind]
stats_list <- stats_list[ind]

# get feature id field (e.g. "symbol" or "gene set");
# should be the same for all datasets
id_field <- colnames(pvals_list[[1]])[1]

# sanity check; make sure no duplicated keys exist
no_dups <- sapply(pvals_list, function(x) {
  nrow(x) == length(unique(pull(x, id_field)))
})

if (!all(no_dups)) {
  stop("Duplicated feature ids found!")
}

# combine inidividual feature-phenotype associations into a single file
coefs_merged <- coefs_list %>%
  reduce(left_join, by = id_field)

pvals_merged <- pvals_list %>%
  reduce(left_join, by = id_field)

stats_merged <- stats_list %>%
  reduce(left_join, by = id_field)

# normalize column order
ind <- c(id_field, sort(colnames(coefs_merged)[-1]))

coefs_merged <- coefs_merged[, ind]
pvals_merged <- pvals_merged[, ind]
stats_merged <- stats_merged[, ind]

# store results
write_csv(coefs_merged, snakemake@output[['coefs']])
write_csv(pvals_merged, snakemake@output[['pvals']])
write_csv(stats_merged, snakemake@output[['stats']])
