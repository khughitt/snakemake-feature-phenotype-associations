#!/bin/env python
#
# Creates a covariate metadata table using the information in the dataset yaml files
#
import yaml
import pandas as pd

entries = []

for infile in snakemake.params['cfg_paths']:
    with open(infile) as fp:
        cfg = yaml.safe_load(fp)

        for phenotype in cfg['phenotypes']:
            for assoc in cfg['phenotypes'][phenotype]['associations']:
                assoc_cfg = cfg['phenotypes'][phenotype]['associations'][assoc]

                entry = [cfg['name'], assoc, assoc_cfg['category'], assoc_cfg['method']]

                entries.append(entry)
                
df = pd.DataFrame(entries, columns=['dataset', 'covariate', 'category', 'method'])

df.to_feather(snakemake.output[0])

