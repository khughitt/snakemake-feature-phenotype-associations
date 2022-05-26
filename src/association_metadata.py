#!/bin/env python
#
# Creates a feature-phenotype association metadata table using the information in the
# dataset yaml files
#
import yaml
import pandas as pd
from utils import load_data

entries = []

for infile in snakemake.params['cfg_paths']:
    with open(infile) as fp:
        cfg = yaml.safe_load(fp)

        # iterate over dataset entries
        for feat_level in cfg['features']: 
            for feat_type in cfg['features'][feat_level]:
                # load feature data 
                feat_path = cfg['features'][feat_level][feat_type]
                feat_dat = load_data(feat_path)

                for phenotype in cfg['phenotypes']['associations']:
                    pheno_cfg = cfg['phenotypes']['associations'][phenotype]


                    entry = {
                        "dataset": cfg['name'],
                        "feature_level": feat_level,
                        "feature_type": feat_type,
                        "phenotype": phenotype,
                        "num_samples": feat_dat.shape[1] - 1,
                        "num_features": feat_dat.shape[0],
                        "feature_path": feat_path,
                        "phenotype_path": cfg['phenotypes']['path'],
                        "category": pheno_cfg['category'],
                        "method": pheno_cfg['method']
                    }

                    entries.append(entry)
                
df = pd.DataFrame(entries)

df.to_csv(snakemake.output[0])

