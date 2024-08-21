# Command line tool (cli)

`tribal` can also be run as a command line tool.

```bash
❯ tribal -h
usage: tribal [-h] {preprocess,fit} ...

Tribal CLI Tool

positional arguments:
  {preprocess,fit}  Sub-commands
    preprocess      Preprocess data
    fit             B cell lineage tree inference

optional arguments:
  -h, --help        show this help message and exit
```

## Overview

The cli has two sub-commands:  
  1. [preprocess](#preprocess) - filter the data and find a multiple sequence alignment and parsimony forsest for each clonotype.  
  2. [fit](#fit) - infer a set of optimal B cell lineage trees per clonotype an a shared isotype transition probability matrix.  

!!! tip
    It is recommended to use the preprocessing tool to prepare the input data to the proper format for `tribal`.  

## Preprocess

The preprocessing command will:
    1. filter out clonotypes that are below the minimum number of cells 
    2. filter out cells which have v alleles that differ from the majority of the clonotype
    3. perform a multiple sequence alignment (MSA) for each valid clonotype using [mafft](https://mafft.cbrc.jp/alignment/software/)
    4. infer a parsimony forest for each clonotype given the MSA using [dnapars](https://phylipweb.github.io/phylip/)



### Usage

```bash
❯ tribal preprocess -h
usage: tribal preprocess [-h] -d DATA -r ROOTS -e ENCODING [--min-size MIN_SIZE] [--dataframe DATAFRAME] [-o OUT]
                         [-j CORES] [--heavy] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -d DATA, --data DATA  filename of csv file with the sequencing data
  -r ROOTS, --roots ROOTS
                        filename of csv file with the root sequences
  -e ENCODING, --encoding ENCODING
                        filename isotype encodings
  --min-size MIN_SIZE   minimum clonotype size (default 4)
  --dataframe DATAFRAME
                        path to where the filtered dataframe with additional sequences and isotype encodings should be
                        saved.
  -o OUT, --out OUT     path to where pickled clonotype dictionary input should be saved
  -j CORES, --cores CORES
                        number of cores to use (default 1)
  --heavy               only use the heavy chain and ignore the light chain
  -v, --verbose         print additional messages
```

### Input
The `--data`, `--roots` and `--encoding` are required arguments.  See [data description](data.md) for more details. The `--encoding` argument should be that path to a text file that lists the correct isotype ordering as well as the isotype labels that are present in the input data. 



```
IGHM
IGHG3
IGHG1
IGHA1
IGHG2
IGHG4
IGHE
IGHA2
```

!!! note "isotype labeling"
    Be sure that the labels used in the encoding file exactly match the labeling syntax in the input data. There is no standard convention for isotype labels,  e.g., IgM, M, IghM and IGHM, and therefore the convention must be provided by the user. 

### Output
The main output from the `preprocess` sub-command is the a pickled dictionary storing `Clonotype` objects for each retained clonotype in the data. This output file will be used as the input to the `fit` sub-command. 


### Example

 Here is an example of how to run the preprocess sub-command.

```bash 
 $ tribal preprocess -d data.csv -r roots.csv -e isotypes.txt -j 4 --min-size 4 --dataframe filtered.csv -o clonotypes.pkl -v

```


## fit

The `fit` sub-command will infer a set of B cell lineage trees for each clonotype and a shared isotype transition probability matrix. 

### Usage 

```bash 
tribal fit -h
usage: tribal fit [-h] -c CLONOTYPES --n_isotypes N_ISOTYPES [--stay-prob STAY_PROB] [-t TRANSMAT] [--niter NITER]
                  [--thresh THRESH] [-j CORES] [--max-cand MAX_CAND] [-s SEED] [--restarts RESTARTS]
                  [--mode {score,refinement}] [--score SCORE] [--transmat-infer TRANSMAT_INFER] [--verbose]
                  [--pickle PICKLE] [--write-results WRITE_RESULTS]

optional arguments:
  -h, --help            show this help message and exit
  -c CLONOTYPES, --clonotypes CLONOTYPES
                        path to pickled clonotypes dictionary of parsimony forests, alignments, and isotypes
  --n_isotypes N_ISOTYPES
                        the number of isotypes states to use
  --stay-prob STAY_PROB
                        the lower and upper bound of not class switching, example: 0.55,0.95
  -t TRANSMAT, --transmat TRANSMAT
                        optional filename of isotype transition probabilities
  --niter NITER         max number of iterations during fitting
  --thresh THRESH       theshold for convergence in during fitting
  -j CORES, --cores CORES
                        number of cores to use
  --max-cand MAX_CAND   max candidate tree size per clonotype
  -s SEED, --seed SEED  random number seed
  --restarts RESTARTS   number of restarts
  --mode {score,refinement}
                        mode for fitting B cell lineage trees, one of 'refinment' or 'score'
  --score SCORE         filename where the objective values file should be saved
  --transmat-infer TRANSMAT_INFER
                        filename where the inferred transition matrix should be saved
  --verbose             print additional messages.
  --pickle PICKLE       path where the output dictionary of LineageTree lists should be pickled
  --write-results WRITE_RESULTS
                        path where all optimal solution results are saved

```

### Input

The pickled dictionary of clonotypes (`clonotypes.pkl`) that is output from `tribal preprocess` will be the input to `tribal fit`.  See [Clonotype](../api/clonotype.md) for details. 

### Example

 Assuming `clonotypes.pkl` is in the working directory, here is an example of how to run
 `tribal fit`. 

```python
tribal fit -c clonotypes.pkl -j 3 --transmat-infer transmat.txt --pickle lineage_trees.pkl
--write-results results --score objective.csv

```

### Output

There are four optional outputs from `tribal fit`:  
  1.  `--transmat-infer`: the inferred isotype transition probability matrix.  
  2.  `--pickle` : a dictionary with clonotype id as key and value is a [LineageTreeList](../api/lineagetreelist.md) containing the optimal lineage trees for each clonotype.  
  3.  `--score` : a  csv file containing the SHM parsimony scores and CSR likelihood.  
  4.  `--write-results` : a directory where the lineage trees files will be saved including:    
  + a fasta file containing the inferred BCR sequences,    
  + a csv file containing the inferred isotypes,    
  + a text file containing the  edge list of the lineage tree,    
  + a png file containing a visualization of the lineage tree.    
  