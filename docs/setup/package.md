# Package overview


The `tribal` package can be imported as a package into a python script or jupyter notebookd or it can be used as a command line tool. 

## The `tribal` package 

The `tribal` package provides the following:  
- [example](#example) data to help you properly format your [input](#input) data and familarize yourself with the package.    
- [preprocess](#preprocess) functionality to preare the data for input to tribal, i.e., filter clonotypes and cells, obtain a multiple sequence alignment (MSA) per clonotype and infer a parsimony forest for each clonotype.   
- [Clonotype class](#clonotype) to store the formatted input data, including the MSA, parsimony forest, and isotype labels of the sequenced B cells. 
-  [fit](#fit) functionality to infer both a B cell lineage forest and isotype transition probabilites for a given dataset.    
- [LineageTree class](#lineagetree) to interact with and visualize the inferred B cell lineage trees.     

The API provides additional details on each of these items. 


### Input


`tribal` requires two input files:  
    - sequencing data saved in csv file with the following columns:  
        * cell : the id or barcode of the sequnced B cell   
        * clonotype: the clonotype id to which that cell belongs   
        * heavy_chain_isotype: the isotype of the constant region of the heavy chain  
        * heavy_chain_seq: the variable region sequence of the heavy chain
        * heavy_chain_allele: the v allele of the heavy chain
        * light_chain_seq:  the variable region sequence of the light chain
        * light_chain_allele: the v allele of the light chain  
    - gerlmine roots saved in csv with the folling columns:  
        * clonotype: the clonotype id of the germline root  
        * heavy_chain_root: the heavy chain variable region germline root sequence
        * light_chain_root: the light chain variable region germline root sequence

!!! note
    The light chain columns may be omitted if the `use_light_chain` argument in `preprocess` is `False`. 




```py
import tribal

```





### Example

The `tribal` package includes some example data in order help you get started.


The sequencing data is stored in a [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) named `df`. 

```python 
from tribal import df

print(df.columns)
print(df.head())

```

The corresponding germline roots for both the heavy and light are stored in a [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) named roots. 
```python 
from tribal import roots

print(roots.columns)
print(roots.head())
```

## Preprocess 

The preprocessing function will:  
    -  filter out clonotypes that are below the minimum number of cells  
    -  filter out cells which have v alleles that differ from the majority of the clonotype  
    -  perform a multiple sequence alignment (MSA) for each valid clonotype using [mafft](https://mafft.cbrc.jp/alignment/software/)  
    -  infer a parsimony forest for each clonotype given the MSA using [dnapars](https://phylipweb.github.io/phylip/)  


See [preprocess](../api/preprocess.md) for more details. 


### Usage
```python


   isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']

   clonotypes, df_filt = preprocess(df, roots,isotypes=isotypes, min_size=4, use_light_chain=True, cores=3, verbose=True )


```


## LineageTree 