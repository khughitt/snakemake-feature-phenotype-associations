"""
"
" Compute feature-phenotype pearson correlation matrix
"
" Pearson correlation for each feature (e.g. gene expression profile) and each 
" phenotype (e.g. drug response profile). The maximum correlation observed for each
" feature across all phenotypes is then returned.
"
"""
import numpy as np
import pandas as pd
from utils import load_data, maxabs

"""
Efficient correlation implementation in python
Source: https://github.com/ikizhvatov/efficient-columnwise-correlation/blob/master/columnwise_corrcoef_perf.py
"""
def corrcoeff_einsum_optimized(O, P):
    (n, t) = O.shape      # n traces of t samples
    (n_bis, m) = P.shape  # n predictions for each of m candidates

    DO = O - (np.einsum("nt->t", O, optimize='optimal') / np.double(n)) # compute O - mean(O)
    DP = P - (np.einsum("nm->m", P, optimize='optimal') / np.double(n)) # compute P - mean(P)

    cov = np.einsum("nm,nt->mt", DP, DO, optimize='optimal')

    varP = np.einsum("nm,nm->m", DP, DP, optimize='optimal')
    varO = np.einsum("nt,nt->t", DO, DO, optimize='optimal')
    tmp = np.einsum("m,t->mt", varP, varO, optimize='optimal')

    return cov / np.sqrt(tmp)

# load processed feature / phenotype data
feat_dat = load_data(snakemake.input['feat_infile'])
pheno_dat = load_data(snakemake.input['pheno_infile'])

# determine feature sample indices
feat_type = snakemake.wildcards['feature_type']
feat_level = snakemake.wildcards['feature_level']

# get samples present in both the feature and phenotype data;
# in both cases, sample ids are expected to start in second column.
sample_ids = feat_dat.columns[1:]
sample_ids = sorted(list(set(pheno_dat.columns[1:]).intersection(sample_ids)))

feat_key  = feat_dat.columns[0]
pheno_key = pheno_dat.columns[0]

feat_dat  = feat_dat.loc[:, [feat_key] + sample_ids]
pheno_dat = pheno_dat.loc[:, [pheno_key] + sample_ids]

# generate numpy versions of the datasets without the id columns and transpose
feat_mat = feat_dat.loc[:, sample_ids].to_numpy().T
pheno_mat = pheno_dat.loc[:, sample_ids].to_numpy().T

# compute pearson correlation
cor_mat = corrcoeff_einsum_optimized(feat_mat, pheno_mat)

# convert back to a dataframe and add row/column identifiers
cor_mat = pd.DataFrame(cor_mat.T, columns=pheno_dat.iloc[:, 0])
cor_mat.insert(0, feat_key, feat_dat[feat_key])

# store correlation matrix
cor_mat.to_parquet(snakemake.output['cor_mat'], compression='zstd')

# convert id field to an index 
cor_mat = cor_mat.set_index(cor_mat.columns[0])

# determine function to use to collapse correlation scores
phenotype = snakemake.wildcards['phenotype']

pheno_config = snakemake.params.config['phenotypes']
assoc_params = pheno_config['associations'][phenotype]['params']

collapse_func = assoc_params['collapse_func']

# collapse correlations along features using specified aggregation function
res = cor_mat.apply(eval(collapse_func), axis=1)

# create result data frame and update column name
res = pd.DataFrame(res)

cname = "_".join([snakemake.wildcards['dataset'], snakemake.wildcards['phenotype'], "cor"])
res.columns = [cname]

# convert index to a column
res = res.reset_index()

res.to_parquet(snakemake.output['result'], compression='zstd')
