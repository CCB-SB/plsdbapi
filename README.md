
# API for the plasmid database PLSDB
PLSDB is a resource of plasmid records that are collected from the NCBI database and subsequently processed to remove incomplete, inconsistent or chromosomal entries, and to add annotations. This package allows to connect to the PLSDB webserver and download plasmid data. Hence, a stable internet connection is required.
To find out more about PLSDB please visit our [web server](https://www.ccb.uni-saarland.de/plsdb). 

All PLSDB tools are provided and hosted by the [Chair for Clinical Bioinformatics at Saarland University](https://www.ccb.uni-saarland.de/).


### Data accessible with plsdbapi
- general data about one or more plasmids using the NCBI accession id:  `summary`
- search nucleotide sequences in PLSDB with *mash screen*, *mash dist*, *blastn* and *tblastn*:  `query_plasmid_sequence`
- filter PLSDB according to e.g. name, location or length and download the found plasmids: `filter_nuccore`, `filter_biosample`, `filter_taxonomy`
- fasta sequences can be downloaded with `summary` and `filter_X` as well as with `download_fasta`

To see some examples on how to use **plsdbapi** in python and R read the [quickstart guides](https://github.com/CCB-SB/plsdbapi/tree/master/examples). 


&nbsp;

# Installation

## Python package
```
git clone git@github.com:CCB-SB/plsdbapi.git
cd plsdbapi
pip install .
```

&nbsp;

## Conda [![Conda][conda-badge]][conda-link]
```
conda install -c ccb-sb plsdbapi
```

[conda-badge]: https://anaconda.org/conda-forge/skidl/badges/installer/conda.svg
[conda-link]: https://anaconda.org/ccb-sb/plsdbapi
