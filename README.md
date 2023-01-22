## TRIBAL: Tree inference of B cell clonal lineages 




![Overview of TRIBAL](overview.png)
TRIBAL is a method to infer isotype transition probabilities, isotype proportions and B cell clonal lineages for a set of n B cells clustering into k clonotypes. 



## Contents

  1. [Installation](#install)
     * [Using github](#compilation)
     * [Dependencies](#pre-requisites)
  2. [Phases](#phases) 
  3. [I/O formats](#io)
  4. [Usage](#usage)
      + [Inferring isotype transition probabilities and proportions](#probabilities)
      + [Inferring B cel clonal lineages](#trees)

<a name="install"></a>

## Installation

  
<a name="compilation"></a> 
### Using github
   1. Clone the repository
      ```bash
            $ git clone https://github.com/elkebir-group/tribal.git
       ```
       <a name="pre-requisites"></a> 
   2. Install dependencies 
      + python3 >=3.9
      + [numpy](https://numpy.org/doc/)
      + [ete3](http://etetoolkit.org) >=3.1.2
      + [networkx](https://networkx.org)
      + [pydot](https://pygraphviz.github.io)  
      
      Optional: 
         + [snakemake](https://snakemake.readthedocs.io/en/stable/)
         + [phylip](https://anaconda.org/bioconda/phylip) >3.697


<a name="phases"></a>
## Phases
TRIBAL is run in two phases. 
  1. It infers the isotype transition probabilities for a set of k clonotypes. 
  2. It uses the inferred isotype transition probabilities to aid in B cell lineage tree inference. 


<a name="io"></a>
## IO Formats

 
 See `example/input` for examples of all input files.  


See `example/output` for examples of all output files


1. *Isotype transition probability inference:* 
    + **Input**:  
        - A text file containing a list of clonotype subdirectory names (see example below) which are to be included in the inference
            ```
             clonotype_1
             clontoype_2
             clonotype_3
            ```
        - A text file containing the ordered encoding of isotypes (see example below). The names of isotypes should match the names in the input files.    
             ```
             IgM/D
             IgG3
             IgA
            ```
        - A specified id of the root sequence  
        - Each clonotype subdirectory should contain the following two files:  
            1. fasta file for the concatenated and aligned heavy and light chain variable region sequences  
            2. fasta or csv file with the isotype expression of the heavy chain for each cell (ids should correpond to sequence ids)  
    + **Output**:  
        - a text file containing the inferred isotype transition probabilties   
        - a text file containing the inferred isotype proportions   
        - Optional outputs include heatmaps and state diagrams of the isotype transition probabilities  
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

        usage: tribal.py [-h] -p PATH [-c CLONOTYPES] [-e ENCODING] [--n_isotypes N_ISOTYPES] [--fasta FASTA] [-i ISOTYPES] [-j JUMP_PROB] [-t TRANSMAT] [-r ROOT] --tree_path TREE_PATH
                        [--candidates CANDIDATES] [--niter NITER] [--thresh THRESH] [--mu MU] [--sigma SIGMA] [--nworkers NWORKERS] [--max_cand MAX_CAND] [-s SEED] [--alpha ALPHA]
                        [--restarts RESTARTS] [--mode {score,refine,search}] [-o OUTPATH] [--score SCORE] [--transmat_infer TRANSMAT_INFER] [--state_probs STATE_PROBS]
                        [--diagram DIAGRAM] [--diagram_pdf DIAGRAM_PDF] [--heatmap HEATMAP] [--save_all_restarts SAVE_ALL_RESTARTS]

        optional arguments:
        -h, --help            show this help message and exit
        -p PATH, --path PATH  path to the directory containing input files
        -c CLONOTYPES, --clonotypes CLONOTYPES
                                filename with list of clonotype subdirectories that should be included in the inference. If not provided, scans provided path for all subdirectory names
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
        --mu MU               mean of gaussian white noise to add for distortion
        --sigma SIGMA         std of gaussian white noise to add for distortion
        --nworkers NWORKERS   number of workers to use in the event in multiple restarts
        --max_cand MAX_CAND   max candidate tree size per clonotype
        -s SEED, --seed SEED
        --alpha ALPHA
        --restarts RESTARTS   number of restarts
        --mode {score,refine,search}
        -o OUTPATH, --outpath OUTPATH
                                path to directory where output files should be saved
        --score SCORE         filename where the score file should be saved
        --transmat_infer TRANSMAT_INFER
                                filename where the inferred transition matrix should be saved
        --state_probs STATE_PROBS
                                filename where the inferred state probabilities should be saved
        --diagram DIAGRAM     filename where the png of transition matrix should be saved
        --diagram_pdf DIAGRAM_PDF
                                filename where the pdf of transition matrix should be saved
        --heatmap HEATMAP     filename where the heatmap pdf of transition matrix should be saved
        --save_all_restarts SAVE_ALL_RESTARTS
                                path where all restarts should be saved

<!-- <a name="cna-mode-example"></a>
### CNA mode example


Here we show an example of how to run `Phertilizer` in CNA Mode.
The input files are located in the `example/input` directory.


    $ phertilizer -f example/input/variant_counts.tsv  \
      --bin_count_data example/input/binned_read_counts.csv \
      --bin_count_normal example/input/normal_cells.tsv --snv_bin_mapping example/input/snv_bin_mapping.csv \
      --min_cells 100 --min_snvs 100 -d 14 --tree example/cna_mode_output/tree.png \
      -n example/cna_mode_output/cell_clusters.csv \
      -m example/cna_mode_output/SNV_clusters.csv -e example/cna_mode_output/CNA_genotypes.csv 

This command generates output files `tree.png`, `cell_clusters.csv`, `SNV_clsuters.csv` and `CNA_genotypes.csv` in directory `example\cna_mode_output`.



<a name="snv-mode-example"></a>
### SNV mode example

Here we show an example of how to run `Phertilizer` in SNV Mode.
The input files are located in the `example/input` directory.


    $ phertilizer snv_mode/run_phertilizer.py -f example/input/variant_counts.tsv \
    --bin_count_data example/input/binned_read_counts.csv  --min_cells 100 --min_snvs 100 -d 14 \
    --tree example snv_mode_output/tree.png -n example/snv_mode_output /cell_clusters.csv \
    -m example/snv_mode_output/SNV_clusters.csv 

This command generates output files `tree.png`, `cell_clusters.csv`, and `SNV_clsuters.csv` in directory `example\snv_mode_output`.
 --> -->
