
import networkx as nx
import numpy as np
import sys, os, re
import argparse 
import pickle
from copy import deepcopy
import utils as ut
from ete3 import Tree
from alignment import Alignment
from score_class import ScoreList
from lineage_tree import LineageForest



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
        # path = nx.shortest_path(nx_tree, source=root_name, target=root)

        G = nx_tree.to_undirected()
        H = nx.dfs_tree(G,source=root)
        # print("Edges of the re-rooted tree:")
        # print(list(H.edges()))

        if H.out_degree[root_name]==0:
            H.remove(root_name)

   
        # nx_tree.remove_edge(root_name, root)
        # nx_tree.add_edge(root, root_name)
      

    return H
        
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
    return iso_encoding, start_iso, counter


def create_input( path,  tree_path, clonotype, root, seq_fasta_fname, 
                 trees_fname, iso_fasta_fname, iso_encoding=None, start_iso=None):

    tree_fname =f"{tree_path}/dnapars/{clonotype}/outtree"
    align_fname = f"{path}/{clonotype}/{seq_fasta_fname}"
    iso_fname =f"{path}/{clonotype}/{iso_fasta_fname}"
    tree_list = create_trees(tree_fname)

    #simplified alignment 
    alignment = Alignment(align_fname,root=args.root).simplify()
    # alignment = ut.read_fasta(align_fname)
    # alignment = {key: list(value.strip()) for key,value in alignment.items()}


    #only treat isotype file as a fasta file if .fasta is in the name, otherwise, we assume it is a csv file dictionary
    if ".fasta" in iso_fname:
        isotypes = ut.read_fasta(iso_fname)
    else:
        isotypes = ut.read_dict(iso_fname)

    if iso_encoding is not None and start_iso is not None:
        isotypes_filt = {}
        for i in alignment:
                iso = isotypes[i]
                if iso not in iso_encoding:
                    iso = isotypes[i].lower()
                    if 'm' in iso or 'd' in iso:
                        iso = start_iso
                isotypes_filt[i] = iso_encoding[iso]
        isotypes = isotypes_filt
    
    linforest = LineageForest(alignment=alignment, isotypes=isotypes)
    linforest.generate_from_list(tree_list, root)

    return linforest

  

def save_results(outpath, lin_tree_dict, pngs=False, isotype_mapping=None):
   
    for clono, res in lin_tree_dict.items():
        clono_path = f"{outpath}/{clono}"
        os.makedirs(clono_path, exist_ok=True)
        tree = res["tree"]
        seq =  ut.update_labels(res["labels"])
        iso = res["isotypes"]

      
  
        tree.save_tree(f"{clono_path}/tree.txt")
        if pngs:
            tree.save_png(f"{clono_path}/tree.png", iso, isotype_mapping)
     
        if isotype_mapping is not None:
            iso_labs = {key: isotype_mapping[val] for key,val in iso.items()}
        else:
            iso_labs =iso 
        ut.write_fasta(f"{clono_path}/seq.fasta", seq)
        ut.write_fasta(f"{clono_path}/isotypes.fasta", iso_labs)
        ut.save_dict(f"{clono_path}/seq.csv", seq)
        ut.save_dict(f"{clono_path}/isotypes.csv", iso_labs)



def create_trees(cand_fname):
    cand_trees = []
      
    exp = '\[.*\]'
       
    with open(cand_fname, 'r') as file:
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

            nx_tree= convert_to_nx(ete_tree, args.root)
          
            cand_trees.append(nx_tree)
        return cand_trees


def pickle_save(obj, fname):
        with open(fname, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--path", type=str, required=True, help="path to the directory containing input files")
    parser.add_argument("-c", "--clonotypes", required=False, type=str,
        help="filename with list of clonotype subdirectories that should be included in the inference. If not provided, scans provided path for all subdirectory names")
    parser.add_argument("-e", "--encoding", type=str, help="text file isotype states listed in germline order")
    parser.add_argument("--n_isotypes", type=int, default=7, help="the number of isotypes states to use if isotype encoding file is not provided and input isotypes are encoded numerically")
    parser.add_argument( "--fasta", type=str, default= "heavy.aln.fasta", help="filename of input MSA in fasta file")
    parser.add_argument("-i", "--isotypes",  type=str, default= "isotype.fasta",
        help="filename of isotype fasta file within each clonotype directory")
    parser.add_argument("-r", "--root", required=False, default="naive",
        help="the common id of the root in all clonotypes")
    parser.add_argument( "--tree_path", type=str, required=True, help="path to directory where candidate trees are saved")
    parser.add_argument("--candidates", type=str, default="outtree", help="filename containing newick strings for candidate trees")
    parser.add_argument( "-o", "--outfile", type=str, help="filename where clonodict pickle object should be saved")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
  
    # dataset = "day_14"
    # fpath = f"/scratch/projects/tribal/experimental_data"
    # args = parser.parse_args([
    #     "--clonotypes", f"{fpath}/{dataset}/clonotypes_igphyml.txt", 
    #     "--encoding", f"{fpath}/mouse_isotype_encoding.txt",
    #     "-p", f"{fpath}/{dataset}/recomb_input",
    #     "--root", "naive",
    #     "--tree_path", f"{fpath}/{dataset}",
    #     "--fasta", "heavy.aln.fasta",
    #     "--candidates", "dnapars/outtree",
    #     "-o", f"{fpath}/{dataset}/tribal_forest.pickle"
    
    #     ])

    if args.encoding is not None:
        iso_encoding, start_iso, n_isotypes = create_isotype_encoding(args.encoding)
        rev_encoding = {val: key for key, val in iso_encoding.items()}
    else:
        n_isotypes = args.n_isotypes
        iso_encoding = None
        start_iso= None 
        rev_encoding = None
    


    
    if args.clonotypes is not None:
        clonotypes = []
        with open(args.clonotypes, 'r+') as file:
            for line in file:
                clonotypes.append(line.strip())

    else:
         clonotypes = [it.name for it in os.scandir(args.path) if it.is_dir()]

    clonodict = {}
    for c in clonotypes:
        print(f"reading input for clonotype {c}")
        clonodict[c] = create_input(args.path, args.tree_path, c, args.root, args.fasta, 
                        args.candidates, args.isotypes, iso_encoding, start_iso)
    
    if args.outfile is not None:
        ut.pickle_save(clonodict, args.outfile)