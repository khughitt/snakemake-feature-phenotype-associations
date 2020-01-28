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

if (!all(colnames(feat_dat)[-1] == pull(pheno_dat, sample_id_col))) {
  stop("Feature/phenotype sample IDs do not match!")
}

# create a version of the feature data without the id columns, for convenience
feat_mat <- as.matrix(feat_dat[, -1])

# get survival time and event data
surv_time <- as.numeric(pull(pheno_dat, assoc_params$time))
surv_event <- as.logical(pull(pheno_dat, assoc_params$event))

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
