#
# Utility functions
#

#'
#' load_data - loads a dataframe
#'
load_df <- function(infile) {
  if (endsWith(infile, '.csv') || endsWith(infile, '.csv.gz')) {
    read_csv(infile, col_types = cols())
  } else if (endsWith(infile, '.tsv') || endsWith(infile, '.tsv.gz')) {
    read_tsv(infile, col_types = cols())
  } else if (endsWith(infile, '.parquet')) {
    arrow::read_parquet(infile)
  } else if (endsWith(infile, '.feather')) {
    arrow::read_feather(infile)
  }
}

#'
#' load_matrix - loads a numeric matrix 
#'
load_matrix <- function(infile) {
  df <- load_df(infile)

  df <- df %>%
    column_to_rownames(colnames(df)[1])

  # fix types
  rnames <- rownames(df)
  df <- data.frame(lapply(df, as.numeric))
  rownames(df) <- rnames

  as.matrix(df)
}
