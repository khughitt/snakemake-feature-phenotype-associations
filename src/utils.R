#
# Utility functions
#

#'
#' load_data - loads a tabular data file stored in one of several supported formats
#'
load_data <- function(infile) {
  if (endsWith(infile, '.tsv') || endsWith(infile, '.tsv.gz')) {
    feat_dat <- read_tsv(infile, col_types = cols())
  } else if (endsWith(infile, '.parquet')) {
    feat_dat <- arrow::read_parquet(infile)
  } else if (endsWith(infile, '.feather')) {
    feat_dat <- feather::read_feather(infile)
  }
}
