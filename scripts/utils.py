"""
Python utility functions
"""
import pandas as pd

def load_data(infile):
    """Loads a tabular data file from disk"""
    if infile.endswith('.csv') or infile.endswith('.csv.gz'):
        dat = pd.read_csv(infile)
    elif infile.endswith('.tsv') or infile.endswith('.tsv.gz'):
        dat = pd.read_csv(infile, sep='\t')
    elif infile.endswith('.parquet'):
        dat = pd.read_parquet(infile)
    elif infile.endswith('.feather'):
        dat = pd.read_feather(infile)

    return dat

def maxabs(x):
    """Returns the maximum absolute value of a set of values"""
    return max(abs(x))
