"""
" Compute feature-phenotype pearson correlation matrix
"""
import numpy as np
import pandas as pd

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
if snakemake.input['features'].endswith('.tsv.gz'):
    feat_dat = pd.read_csv(snakemake.input['features'], sep='\t')
elif snakemake.input['features'].endswith('.parquet'):
    feat_dat = pd.read_parquet(snakemake.input['features'])

pheno_dat = pd.read_csv(snakemake.input['phenotype'], sep='\t')

# determine feature sample indices
feat_type = snakemake.wildcards['feature_type']
feat_level = snakemake.wildcards['feature_level']

# get samples present in both the feature and phenotype data
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
cor_mat = pd.DataFrame(cor_mat.T, columns=pheno_dat.drug)
cor_mat.insert(0, feat_key, feat_dat[feat_key])

# cor_mat = pd.concat([feat_dat.loc[:, feat_keys], cor_mat], axis=1)

# store result
cor_mat.to_parquet(snakemake.output[0], compression='zstd')
