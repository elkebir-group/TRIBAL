# Package overview


The `tribal` package can be imported as a package into a python script or jupyter notebookd or it can be used as a [command line tool](cli.md). 

The following functions and classes are accessible via the `tribal` package:  

| Name               | Description | Type |  
| ------------------ |------|---- |  
| [preprocess](../api/preprocess.md) | preprocess the input data to the correct format for tribal by finding a multiple sequence aligment and parsimony forest for each clonotype | function |
| [BaseTree](../api/base_tree.md) | a class to model the lineage tree topology | class |
| [Clonotype](../api/clonotype.md)| a dataclass to structure the input data for `tribal`| class |
| [Tribal](../api/tribal.md) | the main class to run the `tribal` algorithm and fit the input data |  class |
| [LineageTree](../api/lineagetree.md) | a class to model an inferred B cell lineage tree| class |
| [LineageTreeList](../api/lineagetree.md) | an extensions of class list to contain a list of B cell lineage trees | class |

The API provides additional details on each of these items. 

## Example data
In addition to the above functions and class, the following example data can be imported to help users better understand the data formatting and package use. 

| Name | Description | Type |
|------|-------------|------|
| df   | input sequencing data | `pandas.DataFrame` |
| roots | input germline roots for sequencing data | `pandas.DataFrame` |
| probabilities | example isotype transition probability matrix |`numpy.ndarray` | 
| clonotypes | dictionary of [Clonotype](../api/clonotype.md) objects | `dict` |
| lineage_tree | an example inferred [B cell lineage tree](../api/lineagetree.md)| [LineageTree](../api/lineagetree.md) | 
| lineage_tree_list | an example inferred [B cell lineage tree list](../api/lineagetreelist.md)| [LineageTreeList](../api/lineagetree.md) | 

See [Data](data.md) for more details on the input format for the data. 

Load and view the example input data:

```python
from tribal import df, roots

print(df.head())
print(roots.head())

```

Load and view the example output data from `preprocess`:

```python
from tribal import clonotypes
for key, clonotype in clonptypes:
    print(key)
    print(clonotype)
```

Load and view the example output data from the `tribal` algorithm:

```python
from tribal import probabilities, lineage_tree, lineage_tree_list
print(probabilities)
print(lineage_tree)
print(lineage_tree_list)
```



## Using the package

Here is a brief walkthrough of how to utilize the functionality of the `tribal` package.
First, load the package:

```python
import tribal

```

or, alternatively load specific functions, classes or example data.

```python
from tribal import preprocess, df, roots

```


### Preprocessing

The [preprocess](../api/preprocess.md) function will:  
    -  filter out clonotypes that are below the minimum number of cells .
    -  filter out cells which have v alleles that differ from the majority of the clonotype  
    -  perform a multiple sequence alignment (MSA) for each valid clonotype using [mafft](https://mafft.cbrc.jp/alignment/software/)  
    -  infer a parsimony forest for each clonotype given the MSA using [dnapars](https://phylipweb.github.io/phylip/)  


See [preprocess](../api/preprocess.md) for more details. 



```python
    from tribal import preprocess, df, roots
    isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']
    clonotypes, df_filt = preprocess(df, roots,isotypes=isotypes, 
                                    min_size=4, use_light_chain=True, 
                                    cores=3, verbose=True)
```

The output dictionary `clonotypes` is the formatted input to `tribal`.  To view
the formatted example data without running the preprocessing step, run the following.

```python
from tribal import clonotypes
for key, clonotype in clonotypes:
    print(clonotype)
```


### Running TRIBAL

Tribal takes the dictionary of [Clonotype](../api/clonotype.md) objects as input and can be run in two modes.   
1. `refinement` (recommended) : the full algorithm where the CSR likelihood is optimized by solving the most parsimonious tree refinement problem.  
2. `score` : the input parsimony lineage trees are not refined and isotypes of the internal nodes are inferred using weighted parsimony via the Sankoff algorithm, with the isotype transition probabilities as weights.   

```python
from tribal import Tribal, clonotypes

#the clonotype data contains the following isotypes encoded from 0 to 7
isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']
tr = Tribal(n_isotypes=len(isotypes), restarts=2, niter=15, verbose=True)
        
#run in refinement mode (recommended)
shm_score, csr_likelihood, best_trees, probabilities = tr.fit(clonotypes=clonotypes, 
                                                                mode="refinement", cores=3)

#run in score mode to infer isotypes using weighted parsimony (Sankoff algorithm) w/o tree refinement
shm_score, csr_likelihood, best_trees, probabilities = tr.fit(clonotypes=clonotypes, 
                                                                mode="score", cores=3)
```
`shm_score` and `csr_likelihood` are floats representing the corresponding SHM or CSR objective values. 

`probabilities` is a numpy array of shape `(n_isotypes, n_isotypes)` containing the inferred isotype transition probabilites. 

`best_trees` is a dictionary with clonotype id as key and the value containing a [LineageTreeList](../api/lineagetreelist.md) with all inferred optimal B cell lineage trees for a given clonotype. 

Additionally, `Tribal` can be `fit` with a user-provided isotype transition probability matrix:

```python
from tribal import Tribal, clonotypes, probabilities

#the clonotype data contains the following isotypes encoded from 0 to 7
isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']
tr = Tribal(n_isotypes=len(isotypes), restarts=2, niter=15, verbose=True)

#specifying the transmat argument will skip the step of inferring isotype transition probabilites
shm_score, csr_likelihood, best_trees, probabilities = tr.fit(clonotypes=clonotypes,
                                                                mode="refinement", transmat=probabilites, 
                                                                cores=3)

```


### Exploring and visualizing the inferred B cell lineage trees

`tribal fit` returns a list of all optimal B cell lineage trees for each clonotype. 
Specifically, in the above examples `best_trees` is a dictionary, with clonotype as key, of [LineageTreeLists](../api/lineagetreelist.md). 

A B cell lineage tree for tribal is a rooted tree with nodes labeled by BCR sequences (concentated heavy and optional light chains) and by isotypes. The [LineageTree](../api/lineagetree.md) class also holds the current  SHM parsimony score (`shm_obj`) and CSR likelihood (`csr_obj`). 