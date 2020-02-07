#
# Snakemake Feature Weights Pipeline
# V. Keith Hughitt
#
import os
import glob
import yaml
from os.path import join

# directories
out_dir = join(config['output_dir'], config['name'], config['version'])
data_cfg_dir = join('datasets', config['name'], config['version'])

# load data source configs
datasets = {}

dataset_cfg_paths = glob.glob(os.path.join(data_cfg_dir, '*.yml'))

for infile in dataset_cfg_paths:
    with open(infile, 'r') as fp:
        cfg = yaml.safe_load(fp)
        datasets[cfg['name']] = cfg

# dict of feature input filepaths, indexed by association function
feats = {}

# dict of snakemake rule expand components
expands = {}

# dict of parameters to be passed in to different rules
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
output_prefixes = []

for method in feats:
    outputs = expand(join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}', method),
                     zip,
                     dataset=expands[method]['dataset'],
                     feature_level=expands[method]['feature_level'],
                     feature_type=expands[method]['feature_type'],
                     phenotype=expands[method]['phenotype'])

    output_prefixes = output_prefixes + outputs

# split into gene- and gene set-specific outputs
gene_pval_files = ["{}_pvals.feather".format(x) for x in output_prefixes if "/genes/" in x]
gene_stat_files = ["{}_stats.feather".format(x) for x in output_prefixes if "/genes/" in x]

pathway_pval_files = ["{}_pvals.feather".format(x) for x in output_prefixes if "/gene_sets/" in x]
pathway_stat_files = ["{}_stats.feather".format(x) for x in output_prefixes if "/gene_sets/" in x]

#
# rules
#
rule all:
    input: 
        join(out_dir, 'merged', '{}_gene_association_pvals.feather'.format(config['name'])),
        join(out_dir, 'merged', '{}_gene_association_stats.feather'.format(config['name'])),
        join(out_dir, 'merged', '{}_pathway_association_pvals.feather'.format(config['name'])),
        join(out_dir, 'merged', '{}_pathway_association_stats.feather'.format(config['name'])),
        join(out_dir, 'metadata', 'association_metadata.feather')

rule association_metadata:
    output:
        join(out_dir, 'metadata', 'association_metadata.feather')
    params:
        cfg_paths=dataset_cfg_paths
    script:
        "src/association_metadata.py"

rule combine_pathway_level_associations:
    input:
        pvals=pathway_pval_files,
        stats=pathway_stat_files
    output:
        pvals=join(out_dir, 'merged', '{}_pathway_association_pvals.feather'.format(config['name'])),
        stats=join(out_dir, 'merged', '{}_pathway_association_stats.feather'.format(config['name']))
    script: 'src/combine_associations.R'

rule combine_gene_level_associations:
    input:
        pvals=gene_pval_files,
        stats=gene_stat_files
    output:
        pvals=join(out_dir, 'merged', '{}_gene_association_pvals.feather'.format(config['name'])),
        stats=join(out_dir, 'merged', '{}_gene_association_stats.feather'.format(config['name']))
    script: 'src/combine_associations.R'

if 'logit' in feats:
    rule logistic_regression:
        input: unpack(get_inputs)
        params:
            config=lambda wildcards, output: datasets[wildcards.dataset]
        output:
            pvals=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/logit_pvals.feather'),
            stats=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/logit_stats.feather')
        script: 'src/logistic_regression.R'

if 'survival' in feats:
    rule survival_regression:
        input: unpack(get_inputs)
        params:
            config=lambda wildcards, output: datasets[wildcards.dataset]
        output:
            pvals=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/survival_pvals.feather'),
            stats=join(out_dir, 'datasets/{dataset}/{feature_level}/{feature_type}/{phenotype}/survival_stats.feather')
        script: 'src/survival_regression.R'

