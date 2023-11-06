## TRIBAL: Tree inference of B cell clonal lineages 




![Overview of TRIBAL](Figure1.png)
TRIBAL is a method to infer B cell clonal lineages and isotype transition probabilities for a set of n B cells clustered into k clonotypes. 



## Contents
- [Installation](#installation)
  - [Using github](#using-github)
- [Phases](#phases)
- [IO Formats](#io-formats)
- [Usage](#usage)
  - [Isotype Transition Probability Inference](#isotype-transition-probability-inference)
  - [B cell lineage tree inference](#b-cell-lineage-tree-inference)
  - [Isotype transition probability inference example](#isotype-transition-probability-inference-example)
  - [B cell lineage tree inference example](#b-cell-lineage-tree-inference-example)
- [Snakemake](#snakemake)

<a name="install"></a>

## Installation

  
<a name="compilation"></a> 
### Using github
   1. Clone the repository
      ```bash
            $ git clone https://github.com/elkebir-group/TRIBAL.git
       ```
       <a name="pre-requisites"></a> 
   2. Install dependencies 
      + python3 >=3.9
      + [numpy](https://numpy.org/doc/)
      + [pandas](https://pandas.pydata.org)
      + [ete3](http://etetoolkit.org) >=3.1.2
      + [networkx](https://networkx.org)
      + [pygraphviz](https://pygraphviz.github.io)  
      + [seaborn](https://seaborn.pydata.org)
      + [matplotlib](https://matplotlib.org)
      + [gurobipy](https://pypi.org/project/gurobipy/) ( **Requires the installation of a free academic license.**)
      
      Optional: 
         + [snakemake](https://snakemake.readthedocs.io/en/stable/)
         + [phylip](https://anaconda.org/bioconda/phylip) >3.697
      
      A conda environment named `tribal` with required dependencies can also be created from the provided `tribal_env.yml` file as follows:  
        ``` $ conda env -f tribal_env.yml```
<a name="phases"></a>
## Phases
TRIBAL is run in two phases. 
  1. infers the isotype transition probabilities for a set of k clonotypes. 
  2. uses these probabilities to find the most parsimonious refinement of each inout tree 

Two phases are required because the size of the input sets may be large for some clonotypes. TRIBAL will 
downsample the input set in each iteration to size `--max_cand` to speed up inference. It always retains the best tree found so far in its sample
to ensure convergence of the coordinate ascent appraoch. 
  In the second phase, TRIBAL solves the most parsiminonious refinement for every
input tree given the isotype transition probabilities inferred in phase 1.  Note that if downsampling is not used, then the second step is not required.



<a name="io"></a>
## IO Formats

 
1. *Isotype transition probability inference:* 
    + **Input**:  
        - A text file containing a list of clonotype subdirectory names (see example below) which are to be included in the inference
            ```
             clonotype_1
             clontoype_2
             clonotype_3
            ```
        - A text file containing the ordered list of isotypes found in the data. These will be encoded from $0$ to $r-1$, where $r$ is the number of isotypes listed in the file (see example below). As naming conventions of isotype states varies, the listed isotypes must match the isotype names in the input files. 
             ```
             IgM/D
             IgG3
             IgA
            ```

        - A specified id of the root or germline representing the naive BCR 
        - Each clonotype subdirectory should contain the following two files:  
            1. fasta file for the MSE of the concatenated heavy and light chain variable region sequences  
            2. fasta or csv file with the isotype expression of the heavy chain for each cell (ids should correspond to sequence ids)  
    + **Output**:  
        - a text file containing the inferred isotype transition probabilties   
        - a text file containing the inferred isotype proportions   
        - Optional outputs include heatmaps of the inferred isotype transition probabilites 
2. *Tree inference:* 
    + **Input**:  
        - fasta file for the concatenated and aligned heavy and light chain variable region sequences  
        - a fasta or csv file with the isotype expression of the heavy chain for each cell (ids should correspond to sequence ids)      
        - A text file containing the ordered encoding of isotypes  
        - istoype transition probabilities inferred during the previous phase     
        - tree inference mode (score, refine, search)  
    + **Output**:  
        - a text file containing the tree encoded with each row in the file encoding an edge in the tree as  child, parent    
        - a fasta for csv file containing the inferred sequences of the heavy and light chain   
        - a fasta or csv file containing the inferred isotype states     
        - Optional outputs include png or pdf visualizations of the inferred tree  




 <a name="usage"></a>
## Usage
### Isotype Transition Probability Inference
<a name="probabilities"></a>
```

usage: tribal.py [-h] [-f FOREST] [-p PATH] [-c CLONOTYPES] [-e ENCODING] [--n_isotypes N_ISOTYPES] [--fasta FASTA] [-i ISOTYPES] [-j JUMP_PROB] [-t TRANSMAT]
                 [-r ROOT] [--tree_path TREE_PATH] [--candidates CANDIDATES] [--niter NITER] [--thresh THRESH] [--nworkers NWORKERS] [--max_cand MAX_CAND] [-s SEED]
                 [--restarts RESTARTS] [--mode {score,refine,refine_ilp,search}] [--score SCORE] [--transmat_infer TRANSMAT_INFER] [--state_probs STATE_PROBS]
                 [--heatmap HEATMAP] [--propmap PROPMAP]

optional arguments:
    -h, --help            show this help message and exit
    -f FOREST, --forest FOREST
                          path to pickled clonotypes dictionary of lineage forests
    -p PATH, --path PATH  path to the directory containing input files
    -c CLONOTYPES, --clonotypes CLONOTYPES
                          filename with list of clonotype subdirectories that should be included in the inference. If not provided, scans provided path for all
                          subdirectory names
    -e ENCODING, --encoding ENCODING
                          text file isotype states listed in germline order
    --n_isotypes N_ISOTYPES
                          the number of isotypes states to use if isotype encoding file is not provided and input isotypes are encoded numerically
    --fasta FASTA         filename of input MSA in fasta file
    -i ISOTYPES, --isotypes ISOTYPES
                          filename of isotype fasta file within each clonotype directory
    -j JUMP_PROB, --jump_prob JUMP_PROB
                          for inititalization of transition matrix if not provided
    -t TRANSMAT, --transmat TRANSMAT
                          optional filename of input transition matrix for initialization
    -r ROOT, --root ROOT  the common id of the root in all clonotypes
    --tree_path TREE_PATH
                          path to directory where candidate trees are saved
    --candidates CANDIDATES
                          filename containing newick strings for candidate trees
    --niter NITER         max number of iterations in the fitting phase
    --thresh THRESH       theshold for convergence in fitting phase
    --nworkers NWORKERS   number of workers to use in the event in multiple restarts
    --max_cand MAX_CAND   max candidate tree size per clonotype
    -s SEED, --seed SEED
    --restarts RESTARTS   number of restarts
    --mode {score,refine,refine_ilp,search}
    --score SCORE         filename where the score file should be saved
    --transmat_infer TRANSMAT_INFER
                          filename where the inferred transition matrix should be saved
    --state_probs STATE_PROBS
                          filename where the inferred state probabilities should be saved
    --heatmap HEATMAP     filename where the {png,pdf} of transition matrix should be saved
    --propmap PROPMAP     filename where the {pdf,png} of isotype proportions should be saved

```


### B cell lineage tree inference
<a name="tree"></a>
```

usage: tribal_tree.py [-h] [-f FULL_FOREST] [-c CLONOTYPE] [-a ALIGNMENT] [-i ISOTYPES] [-t TRANSMAT] -r ROOT [--timeout TIMEOUT] [-l LINEAGE] [--forest]
                      [--candidates CANDIDATES] [--mode {score,refine,refine_ilp,search}] -e ENCODING [--alpha ALPHA] [-j JUMP_PROB] [--ntrees NTREES] [-o OUTPUT]
                      [--tree TREE] [--fasta FASTA] [--png PNG] [--all_pngs] [--sequences SEQUENCES] [--score SCORE] [--reversible] [--iso_infer ISO_INFER]
                      [--all_optimal_sol ALL_OPTIMAL_SOL] [--nworkers NWORKERS] [--seed SEED] [--best_tree_diff BEST_TREE_DIFF] [--pickle_best PICKLE_BEST]
                      [--pickle_all PICKLE_ALL]

optional arguments:
  -h, --help            show this help message and exit
  -f FULL_FOREST, --full-forest FULL_FOREST
                        path to pickled clonotypes dictionary of lineeage forests
  -c CLONOTYPE, --clonotype CLONOTYPE
                        name of clonotype lineage to refine from the full forest
  -a ALIGNMENT, --alignment ALIGNMENT
                        filename of input fasta file containing the alignment
  -i ISOTYPES, --isotypes ISOTYPES
                        filename of input file containing the isotype labels
  -t TRANSMAT, --transmat TRANSMAT
                        filename of input transition matrix
  -r ROOT, --root ROOT  the id of the root sequence in the alignment
  --timeout TIMEOUT     max number of hours to let tribal search per tree
  -l LINEAGE, --lineage LINEAGE
                        pickle file of lineage tree/forest returned from tribal.py
  --forest
  --candidates CANDIDATES
                        filename containing newick strings for candidate tree(s)
  --mode {score,refine,refine_ilp,search}
  -e ENCODING, --encoding ENCODING
  --alpha ALPHA
  -j JUMP_PROB, --jump-prob JUMP_PROB
  --ntrees NTREES       number of top scoring trees to return
  -o OUTPUT, --output OUTPUT
                        outputfile of all best trees
  --tree TREE           outputfile of best tree
  --fasta FASTA         filename where reconstructed ancestral sequences should be saved as fasta file
  --png PNG             filename where to save a png of the optimal tree
  --all_pngs
  --sequences SEQUENCES
                        filename where reconstructed ancestral sequences should be saved as csv file
  --score SCORE         filename of the objective function value objective function value
  --reversible          a flag to indicate the standard 0/1 cost function is used (the number of isotype changes is minimized and irreversiblility is ignored)
  --iso_infer ISO_INFER
                        filename of the inferred isotypes for the internal nodes
  --all_optimal_sol ALL_OPTIMAL_SOL
                        path where all optimal solution results are saved
  --nworkers NWORKERS   number of workers to use in the event of multiple input candidate trees
  --seed SEED           random seed for picking a single best tree among all tied trees
  --best_tree_diff BEST_TREE_DIFF
                        best tree RF distances
  --pickle_best PICKLE_BEST
                        filename to pickle the best results
  --pickle_all PICKLE_ALL
                        filename to pickle the best results

```

<!-- ### Isotype transition probability inference example


Here we show an example of how to run `TRIBAL` to infer isotype transition probabilities. We will use experimental dataset GCB_NP_2.  To run TRIBAL for different datasets, replace `GCB_NP_2` with either `day_14` or `GCB_NP_1`.
The input files are located in the `experimental_data/GCB_NP_2`:


    $   python src/tribal.py -c experimental_data/GCB_NP_2/clonotypes.txt \
        -p experimental_data/GCB_NP_2/input \
        -r naive \
        -e experimental_data/mouse_isotype_encoding.txt \
        --tree_path experimental_data/GCB_NP_2/dnapars \
        --alpha 0.75 --niter 10  --thresh 0.1  \
         -j 0.25 --mu 0.075 --sigma 0.05 \
        --nworkers 5 --restarts 5 \
        --transmat_infer experimental_data/GCB_NP_2/transmat.txt  \
        --state_probs experimental_data/GCB_NP_2/proportions.txt \
        --heatmap experimental_data/GCB_NP_2/transmat.png  \
        --propmap experimental_data/GCB_NP_2/proportions.png 

### B cell lineage tree inference example

Here we show an example of how to run `TRIBAL` to infer a B cell lineage tree for `B_12_1_5_24_1_5` for experimental data `GCB_NP_2`. The input files are located in the `experimental_data/GCB_NP_2/B_12_1_5_24_1_5`. First, we narrow down candidate trees to the top 10 by ranking them after tree refinment.

    python src/tribal_tree.py \
        -r naive \
        -a experimental_data/GCB_NP_2/input/B_12_1_5_24_1_5/concat.aln.fasta \
       -t experimental_data/GCB_NP_2/transmat.txt \
       -e experimental_data/mouse_isotype_encoding.txt \
       -i experimental_data/GCB_NP_2/input/B_12_1_5_24_1_5/isotype.fasta \
       --candidates experimental_data/GCB_NP_2/dnapars/B_12_1_5_24_1_5/outtree \
       --alpha 0.75 \
       --mode refine \
        --ntrees 10 \
       --nworkers 5 \
       -o experimental_data/GCB_NP_2/example/forest.pickle -->




