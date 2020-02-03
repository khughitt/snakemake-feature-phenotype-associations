#!/bin/env Rscript
#
# combines feature associations from multiple sources
#
suppressMessages(library(arrow))
suppressMessages(library(feather))
suppressMessages(library(tidyverse))

# load individual dataset associations (p-values, correlations, etc.)
infiles <- unlist(snakemake@input)

scores_list <- lapply(infiles, read_parquet)
names(scores_list) <- infiles

# reorder associations list so that the entry with the most rows appears first
num_rows <- unlist(lapply(scores_list, nrow))
scores_list <- scores_list[order(-rank(num_rows, ties.method = 'first'))]

# get feature id field (e.g. "symbol" or "gene set");
# should be the same for all datasets
id_field <- colnames(scores_list[[1]])[1]

# sanity check; make sure no duplicated keys exist
no_dups <- sapply(scores_list, function(x) {
  nrow(x) == length(unique(pull(x, id_field)))
})

if (!all(no_dups)) {
  stop("Duplicated feature ids found!")
}

# combine individual associations into a single tibble
scores <- scores_list %>%
  reduce(left_join, by = id_field)

# store results
write_feather(scores, snakemake@output[[1]])
