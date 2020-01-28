#!/bin/env Rscript
#
# Compute feature-phenotype logit regression model
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

source('src/utils.R')

set.seed(1)

# load feature data
feat_dat <- load_data(snakemake@input$features)

# get feature settings
feat_level <- snakemake@wildcards$feature_level
feat_type <- snakemake@wildcards$feature_type
phenotype <- snakemake@wildcards$phenotype

feat_config <- snakemake@params$config$features[[feat_level]][[feat_type]]

# get phenotype settings
phenotype <- snakemake@wildcards$phenotype
pheno_config <- snakemake@params$config$phenotype[[phenotype]]

# association parameters
assoc_params <- pheno_config$associations[[snakemake@wildcards$covariate]]$params

# load phenotype data
pheno_dat <- read_tsv(snakemake@input$phenotype, col_types = cols())

# first column in the phenotype file corresponds to sample id
sample_id_col <- colnames(pheno_dat)[1]

# ensure that sample order is consistent between feature/pheno data
ind <- match(colnames(feat_dat)[-1], pull(pheno_dat, sample_id_col))
pheno_dat <- pheno_dat[ind, ]

# create a version of the feature data without the id columns, for convenience
feat_mat <- as.matrix(feat_dat[, -1])

if (!all(colnames(feat_mat) == pull(pheno_dat, sample_id_col))) {
  stop("Feature/phenotype sample IDs do not match!")
}

# get covariate data vector
covariate <- factor(pull(pheno_dat, assoc_params$field))

# build feature-covariate logit regression models and record p-values
pvals <- apply(feat_mat, 1, function(x) {
  # ignore infrequent "fitted probabilities numerically 0 or 1" warnings that may arise
  # due to quasi- or perfect separation
  mod <- suppressWarnings(glm(covariate ~ x, family = "binomial"))

  # in some cases (e.g. when x contains mostly 0's), fit will only include an
  # intercept term, making it necessary to check the coef dimensions
  coefs <- coef(summary(mod))

  # "2" and "4" are the expected number of rows and columns in coef() output, and the
  # indices associated with the p-value for the fitted model
  if (all(dim(coefs) == c(2, 4))) {
    coefs[2, 4]
  } else {
    1
  }
})

# store result
feat_id_col <- colnames(feat_dat)[1]
feat_ids <- pull(feat_dat, feat_id_col)

res <- data.frame(feat_ids, pvals)

colnames(res) <- c(feat_id_col, sprintf("%s_pval", snakemake@wildcards$covariate))

# for microarray data which may include multiple gene symbols for a single
# row (e.g. "ABC1 // ABC2 // ETC"), split each such entries into multiple
# rows and then collapse duplicated gene symbols, keeping only the lowest
# p-value for each symbol
if ((feat_id_col == 'symbol') && (any(grepl('//', feat_ids)))) {
  res <- res %>%
    separate_rows(symbol, sep = " ?//+ ?") %>%
    group_by(symbol) %>%
    summarise_all(min) %>%
    ungroup()
}

res <- res %>%
  arrange(get(feat_id_col))

arrow::write_parquet(res, snakemake@output[[1]], compression = 'ZSTD')
