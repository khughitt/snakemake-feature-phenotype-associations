#!/bin/env Rscript
#
# Compute feature-phenotype differential expression
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))

source('src/utils.R')

set.seed(1)

# load feature data
feat_dat <- load_data(snakemake@input$feat_infile)

# tmp / sanity check...
if (sum(is.na(feat_dat)) > 0) {
  stop("missing values!")
}

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
sample_ids <- pull(pheno_dat, sample_id_col)

# debugging..
save.image(sprintf('/tmp/%s_%s_log.rda', snakemake@wildcards$dataset, snakemake@wildcards$phenotype))

if (!all(colnames(feat_dat)[-1] %in% sample_ids)) {
  stop("Sample ID mismatch! Check to make sure first column in metadata matches colnames in expression data.")
}

ind <- match(colnames(feat_dat)[-1], sample_ids)
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

# create a version of the feature data without the id columns, for convenience
feat_mat <- as.matrix(feat_dat[, -1])

# DESeq usually expected unprocessed counts; rounding here to allow pre-processed
# counts to be used..
feat_mat <- round(feat_mat)

if (!all(colnames(feat_mat) == pull(pheno_dat, sample_id_col))) {
  stop("Feature/phenotype sample IDs do not match!")
}

#####################

# fit DESeq2 model
dds <- DESeqDataSetFromMatrix(countData = feat_mat,
                              colData = pheno_dat,
                              design= formula(pheno_config$design$full))

dds <- DESeq(dds, test="LRT", reduced=formula(pheno_config$design$reduced))

deseqResult <- results(dds)

coefs <- deseqResult$log2FoldChange
test_stats <- deseqResult$stat
pvals <- deseqResult$pvalue

# store result
feat_id_col <- colnames(feat_dat)[1]
feat_ids <- as.character(pull(feat_dat, feat_id_col))

res <- data.frame(feat_ids, coefs, test_stats, pvals, stringsAsFactors = FALSE)

# update column names; dataset id is added as a prefix to avoid collisions when
# joining results from multiple datasets later on
col_prefix <- sprintf("%s_%s_", snakemake@wildcards$dataset, snakemake@wildcards$phenotype)
colnames(res) <- c(feat_id_col, paste0(col_prefix, c('coef', 'stat', 'pval')))

# for microarray data which may include multiple gene symbols for a single
# row (e.g. "ABC1 // ABC2 // ETC"), split each such entries into multiple
# rows
if ((feat_id_col == 'symbol') && (any(grepl('//', feat_ids)))) {
  res <- res %>%
    separate_rows(symbol, sep = " ?//+ ?")
}

feat_ids <- pull(res, feat_id_col)

pval_field <- colnames(res)[4]

# for datasets with multi-mapped identifiers, collapse to a single row keeping
# the smallest p-value
if (length(feat_ids) != length(unique(feat_ids))) {
  res <- res %>%
    group_by_at(feat_id_col) %>%
    slice(which.min(get(pval_field))) %>%
    ungroup()
}

# sanity check
if (length(pull(res, feat_id_col)) != length(unique(pull(res, feat_id_col)))) {
  stop("Duplicated identifiers found in result!")
}

res <- res %>%
  arrange(get(feat_id_col))

coefs <- res %>%
  select(-ends_with('stat'), -ends_with('pval'))

pvals <- res %>%
  select(-ends_with('stat'), -ends_with('coef'))

stats <- res %>%
  select(-ends_with('pval'), -ends_with('coef'))

# strip "_pval" and "_stat" column names suffixes; no longer needed
colnames(coefs) <- sub("_coef$", "", colnames(coefs))
colnames(pvals) <- sub("_pval$", "", colnames(pvals))
colnames(stats) <- sub("_stat$", "", colnames(stats))

# store coefficients, p-values, and test statistics
write_feather(coefs, snakemake@output[["coefs"]])
write_feather(pvals, snakemake@output[["pvals"]])
write_feather(stats, snakemake@output[["stats"]])
