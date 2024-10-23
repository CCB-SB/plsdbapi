# R API Example
The [reticulate](https://github.com/rstudio/reticulate) library allows to use the python package plsdbapi in R code.
```R
library(reticulate)
plsdbapi = import("plasdbapi")
```

&nbsp;

## **Query plasmid ids**
Use `query_plasmid_id` to search **PLSDB** for the NCBI sequence accession IDs given as input. A pandas dataframe is returned. When `fasta` is set to *TRUE* a file containing the corresponding fasta sequences will be downloaded.

When a unique accession was found the plasmid is annotated with the `label` **found**, otherwise with **notfound**. In the latter case a list with all possible matches can be found in the dataframe; sequences will not be downloaded.
```R
ids = c('NZ_CP053191.1', 'NZ_CP05319')
df = plsdbapi$query$summary(ids, fasta=FALSE)
```
    
&nbsp;

## **Query sequences (with mash_screen)**
Use `query_plasmid_sequence` to run `mash_screen`, `mash_dist`, `blastn` or `tbastn` with the input fasta against **PLSDB**. A pandas dataframe with the results is returned. 

Here we run `mash_screen` with the [example fasta from PLSDB](https://www.ccb.uni-saarland.de/plsdb/plasmids/search_form/seq/?example_mix) with a maximum p-value of 0 and a minimum identity of 0.999.
Either a fasta file or a sequence string are accepted for querying.
```R
ifile = 'example_seq_mix.fasta' # downloaded from https://www.ccb.uni-saarland.de/plsdb/plasmids/search_form/seq/
df = plsdbapi$query$query_plasmid_sequence('mash_screen', ifile=ifile, mash_max_v=0, mash_min_i=0.999)
```
| identity | shared_hashes | median_multiplicity | pvalue  |  ACC_NUCCORE | UID_NUCCORE | ... |
|----------|---------------|---------------------|---------|--------------|-------------|-----|
| 1.00000  |    1000  |   1  |   0  |  NZ_MT230195.1 |  1884912098 | ...  |
| 1.00000  |     734  |   1  |   0  |  NZ_MT230189.1 |  1884921125 | ...  |
| 1.00000  |     532  |   1  |   0  |  NZ_AJ223173.1 |  1864210711 | ...  |
| 1.00000  |     422  |   1  |   0  |  NZ_MT230196.1 |  1884911723 | ...  |
| 1.00000  |     331  |   1  |   0  |  NZ_MT230192.1 |  1918495211 | ...  |
| 1.00000  |     313  |   1  |   0  |  NZ_MT230193.1 |  1918154458 | ...  |
| 0.99916  |     281  |   1  |   0  |  NZ_MT230306.1 |  1918154850 | ...  |

&nbsp;

For `mash_screen` the maximum p-value, the minimum identity and the winner-takes-all strategy can be set. For `mash dist` the maximum p-value, the maximum distance and individually strategy are available.
For `blastn` the minimum identity and the minimum coverage can be adjusted; for `tblastn` only the minimum coverage is available.  

&nbsp;

## **Download fasta sequence for plasmids**
Run `download_fasta` with a list of plasmid NCBI accessions to download the fasta sequences. *TRUE* is returned if no errors occur; *FALSE* otherwise.
```R
isdownloaded = plsdbapi$query$download_fasta(c('NC_011102.1', 'NC_01110'), opath='plsdbapi.fasta')
```

&nbsp;

## **Download plasmids based on a filter**
Run `query_plasmid_filter` to filter the records in **PLSDB**. Keeping the default parameters stores the whole plasmid table. The default strategy for strings is *contains* .Alternatively *is*, *begins* or *ends* can be used.
In this example all plasmids containing *NZ_CP0110* and a location beginning with *China* are filtered and returned in a pandas dataframe. Additionally, the fastas of these plasmids can be downloaded (`fasta=TRUE`).
```R
plsdbapi$query$filter_nuccore(NUCCORE_Source="RefSeq", NUCCORE_Topology="circular", NUCCORE_has_identical='yes', AMR_genes="espP,toxB,ehxA,katP")
plsdbapi$query$filter_biosample(ECOSYSTEM_tags="fecal", ECOSYSTEM_taxid=9606, DISEASE_ontid_name='Aspiration pneumonia')
plsdbapi$query$filter_taxonomy(TAXONOMY_strain_id=340184)
```