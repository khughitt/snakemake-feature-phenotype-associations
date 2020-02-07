#!/bin/env python
#
# Creates a phenotype metadata table using the information in the dataset yaml files
#
import yaml
import pandas as pd

entries = []

for infile in snakemake.params['cfg_paths']:
    with open(infile) as fp:
        cfg = yaml.safe_load(fp)

        for phenotype in cfg['phenotypes']['associations']:
            pheno_cfg = cfg['phenotypes']['associations'][phenotype]

            entry = [cfg['name'], phenotype, pheno_cfg['category'], pheno_cfg['method']]

            entries.append(entry)
                
df = pd.DataFrame(entries, columns=['dataset', 'phenotype', 'category', 'method'])

df.to_feather(snakemake.output[0])

