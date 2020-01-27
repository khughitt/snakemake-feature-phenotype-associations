#!/bin/env Rscript
#
# Compute feature-phenotype logic regression model
#
library(arrow)
library(tidyverse)

source('src/utils.R')

set.seed(1)

save.image('~/tmp2.rda')

# load feature data
feat_dat <- load_data(snakemake@input$features)

# get feature config
feat_level <- snakemake@wildcards$feature_level
feat_type <- snakemake@wildcards$feature_type
phenotype <- snakemake@wildcards$phenotype
covariate <- snakemake@wildcards$covariate

feat_config <- snakemake@params$config$features[[feat_level]][[feat_type]]

# get phenotype config
phenotype <- snakemake@wildcards$phenotype
sample_id <- snakemake@params$config$phenotype[[phenotype]]$sample_id

# load phenotype data
pheno_dat <- read_tsv(snakemake@input$phenotype, col_types = cols())

# ensure that sample order is consistent between feature/pheno data
ind <- match(colnames(feat_dat)[-1], pull(pheno_dat, sample_id))
pheno_dat <- pheno_dat[ind, ]

if (!all(colnames(feat_dat)[-1] == pull(pheno_dat, sample_id))) {
  stop("Feature/phenotype sample IDs do not match!")
}

# create a version of the feature data without the id columns, for convenience
feat_mat <- as.matrix(feat_dat[, -1])

# get covariate data vector
covariate <- factor(pull(pheno_dat, snakemake@wildcards$covariate))

# perform multinomial logistic regression
# null <- glm(covariate ~ 1, family = "binomial")

# build feature-covariate logit regression models and record p-values
pvals <- apply(feat_mat, 1, function(x) {
  mod <- glm(covariate ~ x, family = "binomial")

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
res <- data.frame(pull(feat_dat, 1), pvals)

colnames(res) <- c(colnames(feat_dat)[1], 
                   sprintf("%s_pval", snakemake@wildcards$covariate))

arrow::write_parquet(res, snakemake@output[[1]], compression = 'ZSTD')
