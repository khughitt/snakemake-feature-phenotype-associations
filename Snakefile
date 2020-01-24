#
# feature-phenotype association pipeline
#
import os
import glob
import yaml
from os.path import join

out_dir = join(config['output_dir'], config['version'])

# load data source configs
datasets = {}

for infile in glob.glob('datasets/*.yml'):
    with open(infile, 'r') as fp:
        cfg = yaml.safe_load(fp)
        datasets[cfg['name']] = cfg

# dics of feature and phenotype input filepaths, indexed by association function
xinputs = {}
yinputs = {}

expands = {}

params = {}

# determine which rules calculations should be run for each dataset
for id_ in datasets:
    cfg = datasets[id_]
    params[id_] = {}

    for feat_level in cfg['features']:
        for feat_type in cfg['features'][feat_level]:
            for pheno in cfg['phenotypes']:
                for assoc in cfg['associations']:
                    # feature/phenotype inputs
                    xinput = cfg['features'][feat_level][feat_type]
                    yinput = cfg['phenotypes'][pheno]

                    method = cfg['associations'][assoc]['method']

                    # if this is the first time an association of this type is being
                    # added, initialize empty dicts to store relevant info
                    if method not in xinputs:
                        xinputs[method] = []
                        yinputs[method] = []
                        expands[method] = {}

                    # append feat/pheno pair to input lists
                    xinputs[method].append(xinput)
                    yinputs[method].append(yinput)

                    # update expands dict with relevant info 
                    if 'dataset' not in expands[method]:
                        expands[method]['dataset']       = []
                        expands[method]['feature_level'] = []
                        expands[method]['feature_type']  = []
                        expands[method]['phenotype']     = []
                        expands[method]['assoc']         = []

                    expands[method]['dataset'].append(id_)
                    expands[method]['feature_level'].append(feat_level)
                    expands[method]['feature_type'].append(feat_type)
                    expands[method]['phenotype'].append(pheno)
                    expands[method]['assoc'].append(assoc)

                    # store params
                    params[id_][assoc] = cfg['associations'][assoc]['params']

# function to determine input files associated with a set of wildcards
def get_inputs(wildcards):
    return {
        'features':  datasets[wildcards.dataset]['features'][wildcards.feature_level][wildcards.feature_type],
        'phenotype': datasets[wildcards.dataset]['phenotypes'][wildcards.phenotype]
    }

# expected outputs
all_outputs = []

for method in xinputs:
    outfile = '{}.parquet'.format(method)

    outputs = expand(join(out_dir, '{dataset}/{feature_level}/{feature_type}/{phenotype}/{assoc}', outfile),
                     zip,
                     dataset=expands[method]['dataset'],
                     feature_level=expands[method]['feature_level'],
                     feature_type=expands[method]['feature_type'],
                     phenotype=expands[method]['phenotype'],
                     assoc=expands[method]['assoc'])

    all_outputs = all_outputs + outputs

breakpoint()
#
# rules
#
rule all:
    input: all_outputs

if 'logistic_regression' in xinputs:
    rule logistic_regression:
        input: unpack(get_inputs)
        params:
            lambda wildcards, output: datasets[wildcards.dataset]['associations']['logistic_regression']
        output:
            join(out_dir, '{dataset}/{feature_level}/{feature_type}/{phenotype}/logistic_regression.parquet')
        script: 'src/compute_logistic_regression.R'

if 'survival_regression' in xinputs:
    rule survival_regression:
        input: unpack(get_inputs)
        params:
            lambda wildcards, output: datasets[wildcards.dataset]['associations']['survival_regression']
        output:
            join(out_dir, '{dataset}/{feature_level}/{feature_type}/{phenotype}/survival_regression.parquet')
        script: 'src/compute_survival_regression.R'

if 'pearson_correlation' in xinputs:
    rule pearson_correlation:
        input: unpack(get_inputs)
        params:
            lambda wildcards, output: datasets[wildcards.dataset]['associations']['pearson_correlation']
        output:
            join(out_dir, '{dataset}/{feature_level}/{feature_type}/{phenotype}/pearson_correlation.parquet')
        script: 'src/compute_pearson_correlation.py'

if 'spearman_correlation' in xinputs:
    rule spearman_correlation:
        input: unpack(get_inputs)
        params:
            lambda wildcards, output: datasets[wildcards.dataset]['associations']['spearman_correlation']
        output:
            join(out_dir, '{dataset}/{feature_level}/{feature_type}/{phenotype}/spearman_correlation.parquet')
        script: 'test.R'
