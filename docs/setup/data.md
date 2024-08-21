# Input Data
`tribal` requires the input data to be first clustered into clonotypes, i.e., groups of cells that descend from the same naive B cell receptor.  We recommend [Dandelion](https://sc-dandelion.readthedocs.io/en/latest/) to assist with this preprocessing step. 

`tribal` requires two input files for preprocessing. The  [single-cell RNA sequencing data](#sequencing-data) and the [germline root sequences](#germline-clonotype-roots). 

## Sequencing data

Sequencing data should be provided in a csv file with the following columns:   

| Column Name         | Type     | Description | Required |
| :------------------ | :------: | :---- | :------:|
| cell             |  str or int | unqiue id or barcode of the sequnced B cell  | True |
| clonotype           |   str       | unique clonotype id to which that cell belongs | True |
| heavy_chain_isotype   |  str   | the isotype of the constant region of the heavy chain  | True |
| heavy_chain_seq |  str   |the variable region sequence of the heavy chain | True |
|heavy_chain_allele | str | the v allele of the heavy chain | True |
| heavy_chain_isotype   |  str   | the isotype of the constant region of the heavy chain  | True |
| light_chain_seq |  str   |the variable region sequence of the light chain | False |
|light_chain_allele | str | the v allele of the light chain | False, if light_chain_seq not provided |


See [data.csv](https://github.com/elkebir-group/TRIBAL/blob/main/example/data.csv) for an example.

## Germline clonotype roots

Additionally, the germline root sequences by clonotype should be provided in a csv file containing the heavy chain sequence and optionally, the light chain sequence. 

| Column Name         | Type     | Description | Required |
| :------------------ | :------: | :---- | :------:|
| clonotype           |   str       | unique clonotype id of the germline root (naive BCR) | True |
| heavy_chain_root  |  str   | the heavy chain variable region germline root sequence  | True |
| light_chain_root  |  str   | the light chain variable region germline root sequence  | False |

See [roots.csv](https://github.com/elkebir-group/TRIBAL/blob/main/example/roots.csv) for an example.

!!! note
    All light chain columns may be omitted if the `use_light_chain` argument in `preprocess` is `False`.  In other words, `tribal` may be used with only the heavy chain BCR sequences.



