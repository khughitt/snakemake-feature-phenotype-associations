#
# Snakemake Feature Weights Pipeline
# V. Keith Hughitt
#
import os
import glob
import yaml
from os.path import join

out_dir = join(config['output_dir'], config['version'])

# load data source configs
datasets = {}

dataset_cfg_paths = glob.glob('datasets/mm25/{}/*.yml'.format(config['version']))

for infile in dataset_cfg_paths:
    with open(infile, 'r') as fp:
        cfg = yaml.safe_load(fp)
        datasets[cfg['name']] = cfg

# dict of feature input filepaths, indexed by association function
feats = {}

expands = {}

params = {}

# determine which rules calculations should be run for each dataset
for id_ in datasets:
    cfg = datasets[id_]
    params[id_] = {}

    for feat_level in cfg['features']:
        for feat_type in cfg['features'][feat_level]:
            for phenotype in cfg['phenotypes']['associations']:
                # filepath to feature data
                feat_path = cfg['features'][feat_level][feat_type]

                method = cfg['phenotypes']['associations'][phenotype]['method']

                # if this is the first time an association of this type is being
                # added, initialize empty dicts to store relevant info
                if method not in feats:
                    feats[method] = []
                    expands[method] = {}

                # append feat/pheno pair to input lists
                feats[method].append(feat_path)

                # update expands dict with relevant info 
                if 'dataset' not in expands[method]:
                    expands[method]['dataset']       = []
                    expands[method]['feature_level'] = []
                    expands[method]['feature_type']  = []
                    expands[method]['phenotype']     = []

                expands[method]['dataset'].append(id_)
                expands[method]['feature_level'].append(feat_level)
                expands[method]['feature_type'].append(feat_type)
                expands[method]['phenotype'].append(phenotype)

                # store params
                phenotype_cfg = cfg['phenotypes']['associations'][phenotype]

                if 'params' in phenotype_cfg.keys():
                    params[id_][phenotype] = phenotype_cfg['params']
                else:
                    params[id_][phenotype] = {}

# function to determine input files associated with a set of wildcards
def get_inputs(wildcards):
    return {
        'feat_infile': datasets[wildcards.dataset]['features'][wildcards.feature_level][wildcards.feature_type],
        'pheno_infile': datasets[wildcards.dataset]['phenotypes']['path']
    }

# expected outputs
all_outputs = []

for method in feats:
    outfile = '{}.parquet'.format(method)

    outputs = expand(join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}', outfile),
                     zip,
                     dataset=expands[method]['dataset'],
                     feature_level=expands[method]['feature_level'],
                     feature_type=expands[method]['feature_type'],
                     phenotype=expands[method]['phenotype'])

    all_outputs = all_outputs + outputs

# split into gene- and gene set-specific outputs
gene_outputs = [x for x in all_outputs if "/genes/" in x]
pathway_outputs = [x for x in all_outputs if "/gene_sets/" in x]

#
# rules
#
rule all:
    input: 
        join(out_dir, 'summary', 'gene_pvals_indiv.feather'),
        join(out_dir, 'summary', 'gene_pvals.feather'),
        join(out_dir, 'summary', 'gene_scores.feather'),
        join(out_dir, 'summary', 'pathway_pvals_indiv.feather'),
        join(out_dir, 'summary', 'pathway_pvals.feather'),
        join(out_dir, 'summary', 'pathway_scores.feather'),
        join(out_dir, 'summary', 'phenotype_metadata.feather')

rule create_phenotype_metadata_table:
    output:
        join(out_dir, 'summary', 'phenotype_metadata.feather')
    params:
        cfg_paths=dataset_cfg_paths
    script:
        "src/create_phenotype_metadata_table.py"

rule combine_pathway_level_associations:
    input: pathway_outputs
    output:
        pvals_indiv=join(out_dir, 'summary', 'pathway_pvals_indiv.feather'),
        pvals=join(out_dir, 'summary', 'pathway_pvals.feather'),
        scores=join(out_dir, 'summary', 'pathway_scores.feather')
    script: 'src/combine_associations.R'

rule combine_gene_level_associations:
    input: gene_outputs
    output:
        pvals_indiv=join(out_dir, 'summary', 'gene_pvals_indiv.feather'),
        pvals=join(out_dir, 'summary', 'gene_pvals.feather'),
        scores=join(out_dir, 'summary', 'gene_scores.feather')
    script: 'src/combine_associations.R'

if 'logit' in feats:
    rule logistic_regression:
        input: unpack(get_inputs)
        params:
            config=lambda wildcards, output: datasets[wildcards.dataset]
        output:
            join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/logit.parquet')
        script: 'src/compute_logistic_regression.R'

if 'survival' in feats:
    rule survival_regression:
        input: unpack(get_inputs)
        params:
            config=lambda wildcards, output: datasets[wildcards.dataset]
        output:
            join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/survival.parquet')
        script: 'src/compute_survival_regression.R'

if 'pearson_cor' in feats:
    rule pearson_correlation:
        input: unpack(get_inputs)
        params:
            config=lambda wildcards, output: datasets[wildcards.dataset]
        output:
            cor_mat=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/pearson_cor_mat.parquet'),
            result=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/pearson_cor.parquet')
        script: 'src/compute_pearson_correlation.py'

if 'spearman_cor' in feats:
    rule spearman_correlation:
        input: unpack(get_inputs)
        params:
            config=lambda wildcards, output: datasets[wildcards.dataset]
        output:
            cor_mat=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/spearman_cor_mat.parquet'),
            result=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/spearman_cor.parquet')
        script: 'src/compute_spearman_correlation.R'
