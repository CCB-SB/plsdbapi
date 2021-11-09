# Python API Example
```python
from plsdbapi import query
```

&nbsp;

## **Query plasmid ids**
Use `query_plasmid_id` to search **PLSDB** for the NCBI sequence accession IDs given as input. A pandas dataframe is returned. When `fasta` is set to *True* a file containing the corresponding fasta sequences will be downloaded.

When a unique accession was found the plasmid is annotated with the `label` **found**, otherwise with **notfound**. In the latter case a list with all possible matches can be found in the dataframe; sequences will not be downloaded.
```python
ids = ['NZ_CP053191.1', 'NZ_CP05319']
df = query.query_plasmid_id(ids)
```
| label | searched  | UID_NUCCORE   | ACC_NUCCORE   | ...   |   pmlst   | count | matches   |
|-------|-----------|---------------|---------------|-------|-----------|-------|-----------|
| found | NZ_CP053191.1 | 1841300391 | NZ_CP053191.1 | ...  |IncHI2 DLST(1): smr0018(1);smr0199(1)| NaN | NaN |
| notfound | NZ_CP05319 | NaN | NaN | ... | NaN | 4 | NZ_CP053194.1, NZ_CP053193.1, NZ_CP053192.1, NZ_CP053192.1 |

    
&nbsp;

## **Query sequence with mash_screen**
Use `query_plasmid_sequence` to run `mash_screen`, `mash_dist`, `blastn` or `tbastn` with the input fasta against **PLSDB**. A pandas dataframe with the results is returned.

Here we run `mash_screen` with the [example fasta from PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/search_form/seq/?example_mix) with a maximum p-value of 0 and a minimum identity of 0.999.
Either a fasta file or a sequence string are accepted for querying.
```python
ifile = 'example_seq_mix.fasta'
df = query.query_plasmid_sequence('mash_screen', ifile=ifile, mash_max_v=0, mash_min_i=0.999)
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

For `mash_screen` the maximum p-value, the minimum identity and the winner-takes-all strategy can be set. For `mash dist` the maximum p-value, the maximum distance and the individually strategy are available.
For `blastn` the minimum identity and the minimum coverage can be adjusted; for `tblastn` only the minimum coverage is available.  

&nbsp;

## **Download fasta sequence for plasmids**
Run `download_fasta` with a list of plasmid NCBI accessions to download the fasta sequences. *True* is returned if no errors occur; *False* otherwise.
```python
isdownloaded = query.download_fasta(['NC_011102.1', 'NC_01110'])
```

&nbsp;

## **Download plasmids based on a filter**
Run `query_plasmid_filter` to filter the records in **PLSDB**. Keeping the default parameters stores the whole plasmid table. The default strategy for strings is *contains*. Alternatively *is*, *begins* or *ends* can be used.
In this example all plasmids containing *NZ_CP0110* and a location beginning with *China* are filtered and returned in a pandas dataframe. Additionally, the fasta sequences of these plasmids can be downloaded (`fasta=True`).
```python
df = query.query_plasmid_filter(fasta=True, opath='plsdbapi.fasta', plasmid='NZ_CP0110', location='China', location_strategy='begins')
```
| UID_NUCCORE  |  ACC_NUCCORE     |     Description_NUCCORE |  ...  |  pmlst   |
|--------------|------------------|-------------------------|-------|----------|
| 1016070225 | NZ_CP011067.1 | Escherichia coli str. Sanji plasmid pSJ_2, ... | ... |                                              None |
| 1016070224 | NZ_CP011066.1 | Escherichia coli str. Sanji plasmid pSJ_3, ... | ... |                                              None |
| 1016070223 | NZ_CP011065.1 | Escherichia coli str. Sanji plasmid pSJ_82, ... | ... | IncF RST(F2:A-:B-): FIA(-); FIB(-); FIC(-); FII(2)... |
| 1016070222 | NZ_CP011064.1 | Escherichia coli str. Sanji plasmid pSJ_94, ... | ... | IncF RST(F36:A-:B-): FIA(4,20); FIB(-); FIC(-); FI... |
| 1016070221 | NZ_CP011063.1 | Escherichia coli str. Sanji plasmid pSJ_98, ... | ... |                                              None |
| 1016070220 | NZ_CP011062.1 | Escherichia coli str. Sanji plasmid pSJ_255, ... | ... | IncHI2 DLST(3): smr0018(3); smr0199(2) IncN MLST... |

