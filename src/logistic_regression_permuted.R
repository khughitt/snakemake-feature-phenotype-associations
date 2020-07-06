#!/bin/env Rscript
#
# Compute feature-phenotype logit regression model using permuted data
#
suppressMessages(library(arrow))
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

# create a version of the feature data without the id columns, for convenience
feat_mat <- as.matrix(feat_dat[, -1])

if (!all(colnames(feat_mat) == pull(pheno_dat, sample_id_col))) {
  stop("Feature/phenotype sample IDs do not match!")
}

# get phenotype data vector
phenotype <- factor(pull(pheno_dat, pheno_config$params$field))

feat_id_col <- colnames(feat_dat)[1]

# lists to store permutation results in
coefs_list <- list()
pvals_list <- list()
stats_list <- list()

# iterate over permutations and compute test statistics
for (i in 1:snakemake@config$num_permutations) {
  message(sprintf("Permutation %d...", i))

  # permute data
  feat_mat <- feat_mat[, sample(ncol(feat_mat))]

  # build feature-phenotype logit regression models and record p-values
  mod_results <- apply(feat_mat, 1, function(x) {
    # ignore infrequent "fitted probabilities numerically 0 or 1" warnings that may arise
    # due to quasi- or perfect separation
    # mod <- suppressWarnings(glm(phenotype ~ x, family = "binomial"))
    mod <- glm(phenotype ~ x, family = "binomial")

    # in some cases (e.g. when x contains mostly 0's), fit will only include an
    # intercept term, making it necessary to check the coef dimensions
    coefs <- coef(summary(mod))

    # check to see if independent variant is included in fit
    if ('x' %in% rownames(coefs)) {
      # get test statistic and p-value
      coefs['x', c('Estimate', 'z value', 'Pr(>|z|)')]
    } else {
      # in cases where gene is not included in fitted model, return "NA" instead
      # of "1" to avoid artificially inflating the tail of the p-value distribution
      rep(NA, 3)
    }
  })

  # result is in the form of a 2 x <# gene> matrix with the two rows corresponding
  # to the test statistic and p-value.
  coefs <- mod_results[1, ]
  test_stats <- mod_results[2, ]
  pvals <- mod_results[3, ]

  # store result
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
  # the small p-value
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
  colnames(coefs) <- sub("_coef", "", colnames(coefs))
  colnames(pvals) <- sub("_pval", "", colnames(pvals))
  colnames(stats) <- sub("_stat", "", colnames(stats))

  # add permutation number as a suffix to results so that column have different names
  # when combined later
  colnames(coefs)[-1] <- sprintf("%s_%d", colnames(coefs)[-1], i)
  colnames(pvals)[-1] <- sprintf("%s_%d", colnames(pvals)[-1], i)
  colnames(stats)[-1] <- sprintf("%s_%d", colnames(stats)[-1], i)

  # store results
  coefs_list[[i]] <- coefs
  pvals_list[[i]] <- pvals
  stats_list[[i]] <- stats
}

# combine permutation results into single dataframes and store result
coefs <- Reduce(function(...) inner_join(..., by = feat_id_col), coefs_list)
pvals <- Reduce(function(...) inner_join(..., by = feat_id_col), pvals_list)
stats <- Reduce(function(...) inner_join(..., by = feat_id_col), stats_list)

# store p-values and test statistics in two separate files
write_feather(coefs, snakemake@output[["coefs"]])
write_feather(pvals, snakemake@output[["pvals"]])
write_feather(stats, snakemake@output[["stats"]])
