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
pheno_config <- snakemake@params$config$phenotypes[[phenotype]]

# association parameters
assoc_params <- pheno_config$associations[[snakemake@wildcards$covariate]]$params

# load phenotype data
pheno_dat <- load_data(snakemake@input$phenotype)

# determine name of sample id column
if ('sample_id' %in% names(pheno_config)) {
  sample_id_col <- pheno_config$sample_id
} else {
  # if not specified in config, default to first column
  sample_id_col <- colnames(pheno_dat)[1]
}

# ensure that sample order is consistent between feature/pheno data
ind <- match(colnames(feat_dat)[-1], pull(pheno_dat, sample_id_col))
pheno_dat <- pheno_dat[ind, ]

# limit to specific rows, if requested
if ('filter' %in% names(assoc_params)) {
  # determine which samples pass filter
  mask <- pull(pheno_dat, assoc_params$filter$field) %in% assoc_params$filter$values

  if (sum(mask) == 0) {
    stop("Filtering resulted in all rows being removed! Check filter settings.")
  }

  pheno_dat <- pheno_dat[mask, ]

  # apply same mask to the feature data, including the id column
  feat_dat <- feat_dat[, c(TRUE, mask)]
}

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
  # mod <- suppressWarnings(glm(covariate ~ x, family = "binomial"))
  mod <- glm(covariate ~ x, family = "binomial")

  # in some cases (e.g. when x contains mostly 0's), fit will only include an
  # intercept term, making it necessary to check the coef dimensions
  coefs <- coef(summary(mod))

  # check to see if independent variant is included in fit
  if ('x' %in% rownames(coefs)) {
    # get p-value
    coefs['x', 'Pr(>|z|)']
  } else {
    # in cases where gene is not included in fitted model, return "NA" instead
    # of "1" to avoid artificially inflating the tail of the p-value distribution
    NA
  }
})

# store result
feat_id_col <- colnames(feat_dat)[1]
feat_ids <- as.character(pull(feat_dat, feat_id_col))

res <- data.frame(feat_ids, pvals, stringsAsFactors = FALSE)

# update column names; dataset id is added as a prefix to avoid collisions when
# joining results from multiple datasets later on
cname <- sprintf("%s_%s_pval", snakemake@wildcards$dataset, snakemake@wildcards$covariate)
colnames(res) <- c(feat_id_col, cname)

# for microarray data which may include multiple gene symbols for a single
# row (e.g. "ABC1 // ABC2 // ETC"), split each such entries into multiple
# rows
if ((feat_id_col == 'symbol') && (any(grepl('//', feat_ids)))) {
  res <- res %>%
    separate_rows(symbol, sep = " ?//+ ?")
}

feat_ids <- pull(res, feat_id_col)

# for datasets with multi-mapped identifiers, collapse to a single row keeping
# the small p-value
if (length(feat_ids) != length(unique(feat_ids))) {
  res <- res %>%
    group_by_at(feat_id_col) %>%
    summarise_all(min) %>%
    ungroup()
}

# sanity check
if (length(pull(res, feat_id_col)) != length(unique(pull(res, feat_id_col)))) {
  stop("Duplicated identifiers found in result!")
}

res <- res %>%
  arrange(get(feat_id_col))

arrow::write_parquet(res, snakemake@output[[1]], compression = 'ZSTD')
