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
    1. [preprocess](#preprocess) 
    2. [fit](#fit)

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

#### Input
The `--data`, `--roots` and `--encoding` are required arguments.  See [input](./package.md#input) for more details. The `--encoding` argument should be that path to a text file that lists the correct isotype ordering as well as the isotype labels that are present in the input data. 



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

```

### Input

