#!/bin/env Rscript
#
# combines feature associations from multiple sources
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# load individual dataset associations (p-values, correlations, etc.)
coef_infiles <- unlist(snakemake@input[['coefs']])
pval_infiles <- unlist(snakemake@input[['pvals']])
stat_infiles <- unlist(snakemake@input[['stats']])

coefs_list <- lapply(coef_infiles, read_feather)
names(coefs_list) <- coef_infiles

pvals_list <- lapply(pval_infiles, read_feather)
names(pvals_list) <- pval_infiles

stats_list <- lapply(stat_infiles, read_feather)
names(stats_list) <- stat_infiles

# drop missing values
PVAL_IND <- 2

for (i in 1:length(pvals_list)) {
  mask <- !is.na(pvals_list[[i]][, PVAL_IND])

  pvals_list[[i]] <- pvals_list[[i]][mask, ]
  coefs_list[[i]] <- coefs_list[[i]][mask, ]
  stats_list[[i]] <- stats_list[[i]][mask, ]
}

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
  reduce(full_join, by = id_field)

pvals_merged <- pvals_list %>%
  reduce(full_join, by = id_field)

stats_merged <- stats_list %>%
  reduce(full_join, by = id_field)

# normalize column order
ind <- c(id_field, sort(colnames(coefs_merged)[-1]))

coefs_merged <- coefs_merged[, ind]
pvals_merged <- pvals_merged[, ind]
stats_merged <- stats_merged[, ind]

# exclude genes which are missing in all or all but one dataset/covariate pair
num_nas <- apply(pvals_merged[, -1], 1, function(x) {
  sum(is.na(x));
})

num_covariates <- ncol(pvals_merged) - 1
max_missing <- num_covariates - 1

mask <- num_nas < max_missing

coefs_merged <- coefs_merged[mask, ]
pvals_merged <- pvals_merged[mask, ]
stats_merged <- stats_merged[mask, ]

# store results
write_feather(coefs_merged, snakemake@output[['coefs']])
write_feather(pvals_merged, snakemake@output[['pvals']])
write_feather(stats_merged, snakemake@output[['stats']])
