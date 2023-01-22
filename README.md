## TRIBAL: Tree inference of B cell clonal lineages 




![Overview of TRIBAL](overview.png)
TRIBAL is a method to infer isotype transition probabilities, isotype proportions and B cell clonal lineages for a set of n B cells clustering into k clonotypes. 



## Contents

  1. [Installation](#install)
     * [Using github](#compilation)
     * [Dependencies](#pre-requisites)
  2. [I/O formats](#io)
  3. [Modes](#modes)  
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
   2. Install dependencies
      + python3 (>=3.9)
      + [numpy](https://numpy.org/doc/)
      + [ete3 >=3.1.2] (http://etetoolkit.org)
      + [networkx](https://networkx.org)
      + [pydot](https://pygraphviz.github.io)
      Optional:
       + [phylip >3.697](https://anaconda.org/bioconda/phylip)








<a name="io"></a>
## IO Formats

 
 See `example/input` for examples of all input files.  

<!-- The ouput file options include:   -->
  <!-- 1. A png of the clonal tree with maximum posterior probability
  2. A dataframe mapping cells to nodes
  3. A dataframe mappping SNVs to nodes
  4. A dataframe containing the CNA genotypes
  5. A pickle file of the ouput clonal tree
  6. A pickle file containing a ClonalTreeList of all enumerated clonal trees
  7. A dataframe mapping the internal Phertilizer cell indices to cell labels
  8. A  dataframe mapping the internal Phertilizer SNV indices to cell labels -->

See `example/output` for examples of all output files

<a name="modes"></a>
## Phases
TRIBAL is run in two phases:
 1. *Isotype transition probability inference* 
    + Input: 
      a. A text file containing a list of clonotype subdirectories in the working directory which are to be included in the inference
      b. A text file containing the ordred encoding of isotypes 
      d. A specified id of the root sequence
      c. Each clonotype subdirectory should contain the following two files: 
          1. fasta file for the concatenated and aligned heavy and light chain variable region sequences, 
          2. a fasta or csv file with the isotype expression of the heavy chain for each cell (ids should correpond to sequence ids)
    + Output: 
        a. a text file containing the inferred isotype transition probabilties
        b. a text file containing the inferred isotype proportions
        c. Optional outputs include heatmaps and state diagrams of the isotype transition probabilities
 2. *Tree inference 
    + Input: 
        a. fasta file for the concatenated and aligned heavy and light chain variable region sequences, 
        b. a fasta or csv file with the isotype expression of the heavy chain for each cell (ids should correspond to sequence ids)
        c. A text file containing the ordered encoding of isotypes 
        d. istoype transition probabilities inferred during the previous phase
        e. tree inference mode (score, refine, search)
    +Output:
        a. a text file containing the tree encoded with each row in the file encoding an edge in the tree as  child, parent
        b. a fasta for csv file containing the inferred sequences of the heavy and light chain
        c. a fasta or csv file containing the inferred isotype states
        d. Optional outputs include png or pdf visualizations of the inferred tree




<!-- <a name="usage"></a>
## Usage

      usage: phertilizer [-h] -f FILE --bin_count_data BIN_COUNT_DATA [--bin_count_normal BIN_COUNT_NORMAL] [--snv_bin_mapping SNV_BIN_MAPPING] [-a ALPHA]
                        [--min_cells MIN_CELLS] [--min_snvs MIN_SNVS] [--min_frac MIN_FRAC] [-j ITERATIONS] [-s STARTS] [-d SEED] [--npass NPASS] [--radius RADIUS]
                        [-c COPIES] [--neutral_mean NEUTRAL_MEAN] [--neutral_eps NEUTRAL_EPS] [-m PRED_MUT] [-n PRED_CELL] [-e PRED_EVENT] [--tree TREE]
                        [--tree_pickle TREE_PICKLE] [--tree_path TREE_PATH] [--tree_list TREE_LIST] [--cell_lookup CELL_LOOKUP] [--mut_lookup MUT_LOOKUP]

      optional arguments:
      -h, --help            show this help message and exit
      -f FILE, --file FILE  input file for variant and total read counts with unlabled columns: [chr snv cell base var total]
      --bin_count_data BIN_COUNT_DATA
                              input binned read counts with headers containing bin ids
      --bin_count_normal BIN_COUNT_NORMAL
                              input binned read counts for normal cells with identical bins as the bin count data
      --snv_bin_mapping SNV_BIN_MAPPING
                              a comma delimited file with unlabeled columns: [snv chr bin]
      -a ALPHA, --alpha ALPHA
                              per base read error rate
      --min_cells MIN_CELLS
                              smallest number of cells required to form a clone
      --min_snvs MIN_SNVS   smallest number of SNVs required to form a cluster
      --min_frac MIN_FRAC   smallest proportion of total cells(snvs) needed to form a cluster, if min_cells or min_snvs are given, min_frac is ignored
      -j ITERATIONS, --iterations ITERATIONS
                              maximum number of iterations
      -s STARTS, --starts STARTS
                              number of restarts
      -d SEED, --seed SEED  seed
      --npass NPASS
      --radius RADIUS
      -c COPIES, --copies COPIES
                              max number of copies
      --neutral_mean NEUTRAL_MEAN
                              center of neutral RDR distribution
      --neutral_eps NEUTRAL_EPS
                              cutoff of neutral RDR distribution
      -m PRED_MUT, --pred-mut PRED_MUT
                              output file for mutation clusters
      -n PRED_CELL, --pred_cell PRED_CELL
                              output file cell clusters
      -e PRED_EVENT, --pred_event PRED_EVENT
                              output file cna genotypes
      --tree TREE           output file for png (dot) of Phertilizer tree
      --tree_pickle TREE_PICKLE
                              output pickle of Phertilizer tree
      --tree_path TREE_PATH
                              path to directory where pngs of all trees are saved
      --tree_list TREE_LIST
                              pickle file to save a ClonalTreeList of all generated trees
      --cell_lookup CELL_LOOKUP
                              output file that maps internal cell index to the input cell label
      --mut_lookup MUT_LOOKUP
                              output file that maps internal mutation index to the input mutation label

<a name="cna-mode-example"></a>
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
 -->
