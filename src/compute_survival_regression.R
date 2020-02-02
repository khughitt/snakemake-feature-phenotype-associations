#!/bin/env Rscript
#
# Compute feature-phenotype survival regression model
#
suppressMessages(library(arrow))
suppressMessages(library(survival))
suppressMessages(library(tidyverse))

source('src/utils.R')

set.seed(1)

# load feature data
feat_dat <- load_data(snakemake@input$feat_infile)

# get feature settings
feat_level <- snakemake@wildcards$feature_level
feat_type <- snakemake@wildcards$feature_type

# get phenotype settings
phenotype <- snakemake@wildcards$phenotype
pheno_config <- snakemake@params$config$phenotypes$associations[[phenotype]]

# load phenotype data
pheno_dat <- load_data(snakemake@input$pheno_infile)

# determine name of sample id column
sample_id_col <- colnames(pheno_dat)[1]

# ensure that sample order is consistent between feature/pheno data
ind <- match(colnames(feat_dat)[-1], pull(pheno_dat, sample_id_col))
pheno_dat <- pheno_dat[ind, ]

# limit to specific rows, if requested
if ('filter' %in% names(pheno_config$params)) {
  # determine which samples pass filter
  mask <- pull(pheno_dat, pheno_config$params$filter$field) %in% pheno_config$params$filter$values

  if (sum(mask) == 0) {
    stop("Filtering resulted in all rows being removed! Check filter settings.")
  }

  pheno_dat <- pheno_dat[mask, ]

  # apply same mask to the feature data, including the id column
  feat_dat <- feat_dat[, c(TRUE, mask)]
}

if (!all(colnames(feat_dat)[-1] == pull(pheno_dat, sample_id_col))) {
  stop("Feature/phenotype sample IDs do not match!")
}

# create a version of the feature data without the id columns, for convenience
feat_mat <- as.matrix(feat_dat[, -1])

# get survival time and event data
surv_time <- as.numeric(pull(pheno_dat, pheno_config$params$time))
surv_event <- as.logical(pull(pheno_dat, pheno_config$params$event))

# build cox regression models and record p-values
surv <- survival::Surv(time = surv_time, event = surv_event)

fits <- apply(feat_mat, 1, function(x) {
  # hide "Loglik converged before variable  1 ; coefficient may be infinite" warnings
  suppressWarnings(summary(survival::coxph(surv ~ x)))
})

# extract cox regression p-values
pvals <- as.numeric(unlist(lapply(fits, function(x) {
  x$waldtest['pvalue']
})))

# store result
feat_id_col <- colnames(feat_dat)[1]
feat_ids <- as.character(pull(feat_dat, feat_id_col))

res <- data.frame(feat_ids, pvals, stringsAsFactors = FALSE)

# update column names; dataset id is added as a prefix to avoid collisions when
# joining results from multiple datasets later on
cname <- sprintf("%s_%s_pval", snakemake@wildcards$dataset, snakemake@wildcards$phenotype)
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
