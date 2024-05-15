import pandas as pd
import argparse
import utils as ut  
from io import StringIO
import subprocess
from Bio import AlignIO
import tempfile
import os, re
from lineage_tree import LineageForest
import networkx as nx
from ete3 import Tree


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

def run_dnapars(phylip_content):
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Write PHYLIP content to a temporary "infile"
        infile_path = os.path.join(temp_dir, "infile")
        with open(infile_path, "w") as infile:
            infile.write(phylip_content)
        
        # Write config content to a temporary "config" file
        # config_path = os.path.join(temp_dir, "dnapars.cfg")
        # with open(config_path, "w") as config:
        #     config.write(config_content)
        # outfile = os.path.join(temp_dir, "outfile")
        # Run dnapars from the temporary directory
        dnapars_cline = ["dnapars"]
        # process = subprocess.Popen(dnapars_cline, cwd=temp_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        process = subprocess.Popen(dnapars_cline, cwd=temp_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        stdout = process.communicate(input=CONFIG)
        
        # if process.returncode != 0:
        #     raise RuntimeError(f"dnapars failed with error code {process.returncode}: {stderr}")
        
        # Read the contents of outtree file
        outtree_path = os.path.join(temp_dir, "outtree")
        with open(outtree_path, "r") as outtree:
            outtree_content = outtree.read()
    
    return outtree_content

def convert_to_string(seqs, root=None):
    if root is not None:
        mystr = f">naive\n{root}\n"
    else:
        mystr = ""
    for key, val in seqs.items():
        mystr += f">{key}\n{val}\n"
    return mystr

def align(seqs, root):

    #!cat seqs.fasta mafft --quiet - > aligned.fasta
    seq_str = convert_to_string(seqs,root)
    child = subprocess.Popen(["mafft", "--quiet", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child.stdin.write(seq_str.encode())
    child_out = child.communicate()[0].decode("utf8")

    child.stdin.close()
    aligned_seqs = ut.read_fasta(child_out)
 
    return aligned_seqs
  

def convert_alignment_to_phylip(alignment_str, input_format="fasta"):
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

    # return stdout

def convert_alignment_to_phylip(alignment_str, input_format="fasta"):
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


def align_clonotypes(df, roots):
    clonotype_alignments= {}
    for j in df["Clonotype"].unique():
        clono = df[df["Clonotype"]== j]
        heavy_seqs = dict(zip(clono["seq"],clono["Heavy Chain Variable Seq"]))
    
        heavy_root = roots[roots["Clonotype"]==j]["Heavy Chain Root"].values[0]
  
        heavy_chain_align = align(heavy_seqs, heavy_root)
        light_seqs = dict(zip(clono["seq"],clono["Light Chain Variable Seq"]))
    
        light_root = roots[roots["Clonotype"]==j]["Light Chain Root"].values[0]
        light_chain_align = align(light_seqs, light_root)
        concat_seqs = {}
        for key in heavy_chain_align:
            concat_seqs[key] = heavy_chain_align[key].upper() + light_chain_align[key].upper()
        clonotype_alignments[j] = concat_seqs
    return clonotype_alignments



def convert_to_nx(ete_tree, root):
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
      

    return nx_tree
def create_trees(outtrees):
    cand_trees = []
    exp = '\[.*\]'
    file = StringIO(outtrees)
    # with open(cand_fname, 'r') as file:
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

        nx_tree= convert_to_nx(ete_tree, "naive")
        
        cand_trees.append(nx_tree)
        return cand_trees
    
def create_isotype_encoding(fname):

    iso_encoding = {}
    counter = 0

    with open(fname, "r") as file:
        for line in file:
            isotype = line.strip()
            if counter ==0:
                start_iso = isotype 
            iso_encoding[isotype] = counter
            counter += 1
    return iso_encoding

def preprocess(df, roots, isotype_encoding, min_size=4, verbose=True):
    #first filter out clonotypes smaller than min size
    if verbose:
        print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
    df = df.groupby("Clonotype").filter(lambda x: len(x) >= min_size)
    df = df[ df['Heavy Chain Isotype'].notna()]
    # print(df[df['Heavy Chain Isotype']=="IGHD"])
    if verbose:
        print(f"After filtering, the number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
    
    #prep dnapars sequence ids
    df['seq'] = df.groupby('Clonotype').cumcount() + 1
    df['seq'] = 'seq' + df['seq'].astype(str)

    df.columns.values[0] = "cellid"
    roots.columns.values[0] = "Clonotype"
    df["isotype"] = df['Heavy Chain Isotype'].map(isotype_encoding)
    # print(df["Clonotype"].unique())
    
    
    #create the dictionary of alignments 
    # alignments = align_clonotypes(df[df["Clonotype"].isin(["Clonotype_44","Clonotype_9239","Clonotype_7747"])], roots)
    clonodict = {}
    #TODO: add a parallel version 
    for j in df["Clonotype"].unique():
        clono = df[df["Clonotype"]== j]
        isotypes = dict(zip(clono["seq"], clono["isotype"]))
        isotypes["naive"] = 0
        heavy_seqs = dict(zip(clono["seq"],clono["Heavy Chain Variable Seq"]))
    
        heavy_root = roots[roots["Clonotype"]==j]["Heavy Chain Root"].values[0]
  
        heavy_chain_align = align(heavy_seqs, heavy_root)
        light_seqs = dict(zip(clono["seq"],clono["Light Chain Variable Seq"]))
    
        light_root = roots[roots["Clonotype"]==j]["Light Chain Root"].values[0]
        light_chain_align = align(light_seqs, light_root)
        alignment = {}
        for key in heavy_chain_align:
            alignment[key] = heavy_chain_align[key].upper() + light_chain_align[key].upper()
    
        
 
        mapping = dict(zip(clono["seq"].values,clono['cellid'].values ))
        if verbose:
            print(f"Running dnapars for clonotype {j} with {len(alignment)} sequences...")
        align_str = convert_to_string(alignment)
        phylip_str = convert_alignment_to_phylip(alignment_str=align_str, input_format="fasta")
  
        outtrees = run_dnapars(phylip_str)
        tree_list = create_trees(outtrees)


            
        linforest = LineageForest(alignment=alignment, isotypes=isotypes, mapping=mapping)
        print(f"Clonotype {j} has {len(tree_list)} max parsimony trees.")
        linforest.generate_from_list(tree_list, root="naive")
        clonodict[j] = linforest
    
    
    return clonodict


        #convert the alignment back to a fasta file
    

def main(args):
    df  =pd.read_csv(args.file)
    roots =pd.read_csv(args.roots)
    iso_encoding = create_isotype_encoding(args.encoding)
    clonodict = preprocess(df, roots, isotype_encoding = iso_encoding, min_size=args.min_size)
    if args.pickle is not None:
        pd.to_pickle(clonodict, args.pickle)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file", type=str,
        help="filename of csv file with the clone data")
    parser.add_argument("-r", "--roots", type=str,
        help="filename of csv file with the root sequences")
    parser.add_argument("-e", "--encoding", type=str,
        help="filename isotype encodings")
    
    parser.add_argument( "--min-size", type=int, default=4,
        help="minimum clonotype size")
    parser.add_argument("-P", "--pickle", type=str,
        help="path to where pickled clonotype dictionary input should be saved")


    args= parser.parse_args()
    args = parser.parse_args([
        "-f", "Human/human_data.csv",
        "-r", "Human/human_data_root_seq.csv",
        "-e", "human_encoding.txt",
        "-P", "Human/clonotypes.pkl"

    ])

    main(args)


  





