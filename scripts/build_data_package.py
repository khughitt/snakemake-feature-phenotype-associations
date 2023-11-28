#!/bin/env python
"""
Generate datapackages for finalized datasets
"""
import datetime
import os
import shutil
from frictionless import describe

snek = snakemake

target = snek.params["target"]
version = snek.config["version"]

# metadata to include
mdata = {
    "id": f"mm30-fassoc-{target}-{version}",
    "title": f"MM30: feature association P-values ({target})",
    "rows": "genes",
    "columns": "experiment-covariate pairs",
    "biodat": {
        "datatype": "analysis_result",
        "diseases": ["D009101"],
        "species": [9606]
    },
    "provenance": [{
        "action": "package",
        "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%m:%s"),
        "urls": ["https://github.com/khughitt/snakemake-feature-phenotype-associations"],
        "description": f"MM30 {version}: {target}-level feature association results"
    }],
    "contributors": [{
        "title": "V. Keith Hughitt",
        "email": "keith.hughitt@nih.gov",
        "role": "author"
    }]
}

# switch to data dir
out_dir = os.path.dirname(snek.output[0])
os.chdir(out_dir)

# copy relevant data over
shutil.copy(snek.input[0], out_dir)

pkg = describe("*.feather", stats=True)

# add dataset metadata
for key in mdata.keys():
    pkg.custom[key] = mdata[key]

# write datapackage out
pkg.to_yaml("datapackage.yml")
