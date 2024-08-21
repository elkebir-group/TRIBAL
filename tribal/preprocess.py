"""A module to preprocess data for TRIBAL."""

import os
import re
import tempfile
from io import StringIO
import subprocess
import multiprocessing as mp
import pandas as pd
from Bio import AlignIO
import networkx as nx
from ete3 import Tree
from .clonotype import Clonotype
from .utils import read_fasta


CONFIG = """
J
1
10
O
1
4
5
.
Y
"""


def _process_clonotype(j, df, roots, use_light_chain,isotype_encoding, verbose):
    """Process each clonotype individually.

    Finds a multiple sequence alignment and parsimony forest for a single clonotype.

    """
    try:
        if verbose:
            print(f"Preprocessing clonotype {j}")


        if j not in roots["clonotype"].values:
            raise ValueError(f"Clonotype {j} is not included in germline root dataframe")

        heavy_root = roots.loc[roots["clonotype"] == j, "heavy_chain_root"].values[0]


        # Further processing
        clono = df[df["clonotype"] == j]
        isotypes = dict(zip(clono["cellid"], clono["isotype"]))
        isotypes["naive"] = 0
        heavy_seqs = dict(zip(clono["seq"], clono["heavy_chain_seq"]))


        # Align heavy chain
        heavy_chain_align = _align(heavy_seqs, heavy_root)

        if use_light_chain:
            light_seqs = dict(zip(clono["seq"], clono["light_chain_seq"]))
            light_root = roots.loc[roots["clonotype"] == j, "light_chain_root"].values[0]
            light_chain_align = _align(light_seqs, light_root)
        else:
            light_chain_align = {}

        # Combine alignments
        alignment = {}
        for key, heavy in heavy_chain_align.items():
            if use_light_chain:
                alignment[key] = heavy.upper() + light_chain_align.get(key, "").upper()
            else:
                alignment[key] = heavy.upper()


        # print(f"Generated {len(tree_list)} trees for clonotype {j}")
        mapping = dict(zip(clono["seq"].values,clono['cellid'].values ))

        # Run DNAPARS
        phylip_str = _convert_alignment_to_phylip(_convert_to_string(alignment))
        outtrees = _run_dnapars(phylip_str)

        # Create trees
        tree_list = _create_trees(outtrees, mapping)


        alignment = _update_labels(alignment, mapping)

        linforest = Clonotype(id=j,
                                alignment=alignment,
                                isotypes=isotypes,
                            #   mapping=mapping,
                                isotype_encoding=isotype_encoding)
        if verbose:
            print(f"Clonotype {j} has {len(tree_list)} max parsimony trees.")
        linforest.generate_from_list(tree_list, root="naive")
        return j,linforest


    except RuntimeError as e:
        print(f"An error occurred: {e}")
        return j, None



def preprocess( df: pd.DataFrame,
            roots: pd.DataFrame,
            isotypes:list,
            min_size=4,
            use_light_chain=True,
            cores:int =1,
            verbose=False
            ):
    """
    Preprocess the input data to prepare data for TRIBAL.

    The clonotypes will first be filtered to ensure each clonotype has at least `min_size` cells.
    Each retained clonotype is aligned to the root sequence and then maximum parsimony forest is
    enumerated for the B cells with that clonotype.

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe with columns including 'clonotype', 'heavy_chain_seq', 'light_chain_seq',
        'heavy_chain_v_allele', 'light_chain_v_allele', and 'heavy_chain_isotype'.
        The 'Light Chain' columns are optional.
    roots : pd.DataFrame
        A dataframe containing the root sequences.
    isotypes : list
        A list of the ordered isotype labels, i.e., ['IghM', ...,'IghA']. These labels should
        match the isotype labels in the input dataframe. See notes below.  
    min_size : int, optional
        The minimum number of B cells needed to form a clonotype (default is 4).  
    use_light_chain : bool, optional
        Should the light chain be included in the BCR sequences (default is True).  
    cores : int, optional
        The number of CPU cores to use (default is 1).  
    verbose : bool, optional
        Should verbose output be printed (default is False).  

    Returns
    -------
    dict 
        A dictionary of clonotype objects formatted for input to tribal with clonotype id as key
    pd.DataFrame
        A dataframe that is filtered to contain only unfiltered B cells
    
    Notes
    -----
    Ensure that the isotype labels in `isotypes` match the labels in the input dataframe.
    """
    isotype = {iso: i for i, iso in enumerate(isotypes)}
    if verbose:
        print("\nPreprocessing input data for tribal...")
        print("\nIsotype ordering:")
        for iso in isotypes:
            print(iso)
        print("\nParameter settings:")
        print(f"minimum clonotype size: {min_size}")
        print(f"include light chain: {use_light_chain}")
        print(f"cores: {cores}")
        print(f"verbose: {verbose}")


    #first filter out clonotypes smaller than min size
    if verbose:

        print("\nPrior to filtering...")
        print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['clonotype'].nunique()}.")
    
    df = df.groupby("clonotype").filter(lambda x: len(x) >= min_size)

    df = df[df["heavy_chain_isotype"].isin(isotype.keys())]

    if verbose:
        print(f"\nFiltering clonotypes with fewer than {min_size} cells...")
        print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['clonotype'].nunique()}.")

    df = _filter_alleles(df, "heavy_chain_v_allele")
    if use_light_chain:
        df = _filter_alleles(df, "light_chain_v_allele")
    df = df[ df['heavy_chain_isotype'].notna()]

    df = df.groupby("clonotype").filter(lambda x: len(x) >= min_size)

    if verbose:
        print(f"\nFiltering cells based on v_alleles {min_size}...")

        print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['clonotype'].nunique()}.")

    if verbose:
        print(f"\nAfter all filtering, the number of cells is {df.shape[0]} and the number of clonotypes is {df['clonotype'].nunique()}.\n")
    
    #prep dnapars sequence ids
    df['seq'] = df.groupby('clonotype').cumcount() + 1
    df['seq'] = 'seq' + df['seq'].astype(str)

    df.columns.values[0] = "cellid"
    roots.columns.values[0] = "clonotype"
    df["isotype"] = df['heavy_chain_isotype'].map(isotype)

    clonodict = {}    
    instances  = [ (j, df.copy(), roots.copy(), use_light_chain, isotypes, verbose)
                   for j in df["clonotype"].unique()]
    if cores > 1:
        with mp.Pool(cores) as pool:
            results = pool.starmap(_process_clonotype, instances)
    else:
        results = []
        for inst in instances:
            results.append(_process_clonotype(*inst))
    
    for j, linforest in results:
        clonodict[j] = linforest

    if verbose:
        print("\nPreprocessing complete!")
    return clonodict, df
    

    

def _filter_alleles( df, col):
    grouped = df.groupby(['clonotype', col]).size().reset_index(name='count')

    # # Identify the allele with the highest count for each Clonotype
    idx = grouped.groupby('clonotype')['count'].idxmax()
    max_alleles = grouped.loc[idx]

    # # Merge with the original dataframe to filter the desired rows
    filtered_df = pd.merge(df, max_alleles[['clonotype', col]], on=['clonotype', col])

    return filtered_df


  


def _run_dnapars(phylip_content):
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Write PHYLIP content to a temporary "infile"
        infile_path = os.path.join(temp_dir, "infile")
        with open(infile_path, "w") as infile:
            infile.write(phylip_content)
        
   
        # if dnapars_path is None:
        def get_dnapars_path():
            # Determine the path to the dnapars executable within the package directory
            return os.path.join(os.path.dirname(__file__), 'dnapars', 'dnapars')
 
        # dnapars_cline = [get_dnapars_path()]
        dnapars_cline = ["dnapars"]
        # else:
        #     dnapars_cline = [dnapars_path]
        # process = subprocess.Popen(dnapars_cline, cwd=temp_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        process = subprocess.Popen(dnapars_cline, cwd=temp_dir, 
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        stdout = process.communicate(input=CONFIG)
        
        # Read the contents of outtree file
        outtree_path = os.path.join(temp_dir, "outtree")
        with open(outtree_path, "r") as outtree:
            outtree_content = outtree.read()
    
    return outtree_content


def _convert_to_string(seqs, root=None):
    if root is not None:
        mystr = f">naive\n{root}\n"
    else:
        mystr = ""
    for key, val in seqs.items():
        mystr += f">{key}\n{val}\n"
    return mystr


def _align(seqs, root):
    #!cat seqs.fasta mafft --quiet - > aligned.fasta
    seq_str = _convert_to_string(seqs, root)
    child = subprocess.Popen(["mafft", "--quiet", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child.stdin.write(seq_str.encode())
    child_out = child.communicate()[0].decode("utf8")

    child.stdin.close()
    aligned_seqs = read_fasta(child_out)

    return aligned_seqs


    

def _convert_alignment_to_phylip(alignment_str, input_format="fasta"):
    # Create a file-like object from the string
    input_handle = StringIO(alignment_str)
    
    # Read the alignment
    alignment = AlignIO.read(input_handle, input_format)
    
    # Write the alignment in PHYLIP format to a string buffer
    output_handle = StringIO()
    AlignIO.write(alignment, output_handle, "phylip")
    
    # Get the content of the PHYLIP formatted alignment
    phylip_str = output_handle.getvalue()
    
    return phylip_str



def _align_clonotypes(df, roots):
    clonotype_alignments= {}
    for j in df["clonotype"].unique():
        clono = df[df["clonotype"]== j]
        heavy_seqs = dict(zip(clono["seq"],clono["heavy_chain_seq"]))
    
        heavy_root = roots[roots["clonotype"]==j]["heavy_chain_root"].values[0]

        heavy_chain_align = _align(heavy_seqs, heavy_root)
        light_seqs = dict(zip(clono["seq"],clono["light_chain_seq"]))
    
        light_root = roots[roots["clonotype"]==j]["light_chain_root"].values[0]
        light_chain_align = _align(light_seqs, light_root)
        concat_seqs = {}
        for key in heavy_chain_align:
            concat_seqs[key] = heavy_chain_align[key].upper() + light_chain_align[key].upper()
        clonotype_alignments[j] = concat_seqs
    return clonotype_alignments



def _convert_to_nx( ete_tree, root, mapping):
    nx_tree = nx.DiGraph()
    internal_node = 1
    internal_node_count = 0
    for node in ete_tree.traverse("preorder"):

        if node.name == "":
            node.name = internal_node
            internal_node_count += 1
            internal_node += 1
        if node.is_root():
            root_name =node.name


        for c in node.children:
            if c.name == "":
                c.name = str(internal_node)
                internal_node += 1
                internal_node_count += 1
    

            nx_tree.add_edge(node.name, c.name)

    if len(list(nx_tree.neighbors(root))) == 0:

        nx_tree.remove_edge(root_name, root)
        nx_tree.add_edge(root, root_name)
    
    nx_tree = nx.relabel_nodes(nx_tree, mapping=mapping)
    return nx_tree
    

def _create_trees(outtrees, mapping):
    cand_trees = []
    exp = '\[.*\]'
    file = StringIO(outtrees)
  
    nw_strings = []
    nw_string = ""
    for nw in file:
            line = nw.strip()
            nw_string += line
            if ";" in line:
                
                nw_strings.append(nw_string)
                nw_string = ""

    for nw in nw_strings:

        nw = re.sub(exp, '', nw)
        

        ete_tree = Tree(nw, format=0)

        nx_tree= _convert_to_nx(ete_tree, "naive", mapping)
        
        cand_trees.append(nx_tree)
    return cand_trees





def _update_labels(mydict, mapping):
        # labels = {key : "".join(value) for key, value in self.sequences.items()}
        remapped_labels = {}
        for key, seq in mydict.items():
            if key in mapping:
          
                    remapped_labels[mapping[key]] = seq
           
        
            else:
        
                    remapped_labels[key] = seq
        
        return remapped_labels







