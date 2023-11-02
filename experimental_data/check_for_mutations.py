from Bio.Seq import Seq
import sys
sys.path.append("../src")
import score_class as sc 
import utils as ut
from dataclasses import dataclass
import networkx as nx 
import numpy as np
import pandas as pd 
import argparse 
import lineage_tree as lt 


def is_ancestral(tree, u,v):
    try:
        path = nx.shortest_path(tree, source=u, target=v)
    except:
        path = []    
    return len(path) > 0


def dfs_mut(tree, mut_list):
    mut_dict = {}
    for v in nx.dfs_preorder_nodes(tree.T, tree.root):
            if v in mut_list:
                to_delete = []
                mut_dict[v] = []
                for u in mut_list:
                    if nx.has_path(tree.T, v,u):
                        
                        mut_dict[v].append(u)
                        to_delete.append(u)
                mut_list = [m for m in mut_list if m not in to_delete]
    return mut_dict



                  


                


        

# def find_ancestral_descendants(tree, nodes):
#     # Step 1: Find the most ancestral node(s)
#     most_ancestral_nodes = []
#     for node in nodes:
#         is_ancestral = all(nx.has_path(tree, node, other_node) for other_node in nodes if other_node != node)
#         if is_ancestral:
#             most_ancestral_nodes.append(node)

#     # Step 2: For each most ancestral node, find its descendants
#     ancestral_descendants = {}
#     for ancestral_node in most_ancestral_nodes:
#         descendants = [node for node in nodes if nx.has_path(tree, ancestral_node, node)]
#         ancestral_descendants[ancestral_node] = descendants
#     return most_ancestral_nodes, ancestral_descendants


def find_mutations(reference, query):

    mutations = []


    # Check for substitutions at each position
    for i in range(len(reference)):
        if len(query) == 0:
            continue
        if reference[i] ==  'X' or query[i] =='X':
            continue
        if reference[i] != query[i]:
          
            m =f"{reference[i]}{i+1}{query[i]}"

            mutations.append(m)


    
    return mutations

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--trees", type=str, help="path to list of pickled best scores" )
    parser.add_argument("-a", "--alignment", type=str,
        help="filename of input fasta file containing the alignment")
    parser.add_argument("--method", choices=["tribal", "igphyml"])
    parser.add_argument("-o", "--output", type=str,
                        help= "CSV of mutation presence, level and number of subtrees present in" )
    parser.add_argument("-p", "--pairwise", type=str,
                        help= "CSV of pairwise mutation relationship summary" )

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    # clono = "B_82_9_8_148_1_41"
    # pth = "/scratch/projects/tribal/experimental_data/day_14/igphyml/results"
    # args = parser.parse_args([
    #     "--method", "igphyml",
    #     "-t", f"{pth}/{clono}.pickle",
    #     # "-s", f"{pth}/{clono}.sequence.csv"
    #     "-o", "test/mut_summary.csv",
    #     "-p", "test/pairwise_summary.csv"
    # ])

    # pth = "/scratch/projects/tribal/experimental_data/day_14/tribal_recomb/ml/refine_ilp"
    # args = parser.parse_args([
    #     "--method", "tribal",
    #     "-t", f"{pth}/{clono}/forest.pickle",
    #     "-a", f"/scratch/projects/tribal/experimental_data/day_14/recomb_input/{clono}/heavy.aln.fasta",
    #     "-o", "test/mut_summary.csv",
    #     "-p", "test/pairwise_summary.csv"
    # ])


    # forest_file = f"/scratch/projects/tribal/experimental_data/day_14/tribal_recomb/marginal/refine_ilp/{clono}/forest.pickle"
    # fasta = f"/scratch/projects/tribal/experimental_data/day_14/recomb_input/{clono}/heavy.aln.fasta"
    if args.method=="tribal":
        seqs = ut.read_fasta(args.alignment)
        scores = sc.load(args.trees)
        i = 0
        lin_tree = scores[i].tree
        lin_tree.save_png("test/tribal.png", scores[i].isotypes)
        source_node = 'naive'


        pars_score , pars_seqs = lin_tree.sequence_parismony(seqs, compute=True)
        anc_seqs = {key: "".join(val) for key,val in pars_seqs.items()}

    if args.method == "igphyml":
        source_node = "naive"
        forest = lt.load(args.trees)
        alignment = forest.alignment
        anc_seqs = {key: "".join(val) for key,val in alignment.items()}
        lin_tree = forest.forest[0]
        pruned = lin_tree.prune_unifurcations()
        for p in pruned:
            del anc_seqs[p]
        lin_tree.save_png("test/igphyml.png", forest.isotypes, hide_underscore=False)


     
    
    
    dna   = {key: Seq(val) for key, val in anc_seqs.items()}
    protein = {}
    for key,val in dna.items():
        try:
            protein[key] = val.translate()
        except:
            print(f"waring: sequence {key} not translated!")
            protein[key] = ""
    # {key: val.translate() for key,val in dna.items()}

    germline= dna[source_node].translate()
    node_levels = lin_tree.get_node_levels()

# Define the reference and query DNA sequences
# reference_sequence_dna = Seq("ATC")
# query_sequence_dna = Seq("ATG")



    prot_sub = {}
    for key, query in protein.items():
        prot_sub[key] = find_mutations(germline, query)


    muts = ['S66N', 'K59R', 'W33L']
    min_level_dict = {}
    node_mut_mapping = {s:[] for s in muts}
    mut_dict = {}
  
    for s in muts:

        min_level = np.Inf
        for key, mut_list in prot_sub.items():
            if key ==lin_tree.root:
                continue
            if s in mut_list:
                node_mut_mapping[s].append(key)
                if node_levels[key] < min_level:
                    min_level = node_levels[key] -1
        min_level_dict[s] = min_level

    ancestral_dict = {m : dfs_mut(lin_tree, node_mut_mapping[m]) for m in muts}

    # ancestral_dict ={}
    # for s, nodes in node_mut_mapping.items():

    #     most_anc, desc_sets = find_ancestral_descendants(lin_tree.T, nodes)
    #     ancestral_dict[s] = most_anc
    leafs = [k for k in prot_sub if 'seq' in k]

    summary = []
    if 'naive' in prot_sub:
        del prot_sub['naive']
    if 'Germline' in prot_sub:
        del prot_sub['Germline']
    for m, lev in min_level_dict.items():
        leafs_present = [l for l in leafs if m in prot_sub[l]]
        nodes_present = [n for n in prot_sub if m in prot_sub[n]]

        if lev == np.Inf:
            x = np.NAN
        else:
            x = lev 

        nnodes = len(prot_sub) 
        summary.append([m, lev < np.Inf, x, len(ancestral_dict[m]), len(leafs_present), len(leafs), len(nodes_present), nnodes])
    cols =["mut", "is_present", "path_length_from_mrca", "num_subtrees","cells_present", "nleafs", "node_present", "nnodes"]
    if len(summary) == 0:
        data = {c: [] for c in cols }
        summary_df = pd.DataFrame(data)
   
    else:
        summary_df  = pd.DataFrame(summary, columns=cols)

    if args.output is not None:
        print(f"saving {args.output}......")
        summary_df.to_csv(args.output, index=False)

    results = []

    for s in muts:

    
        for m in muts:
            coocur=  [l for l in leafs if s in prot_sub[l] and m in prot_sub[l]]
            num_occur = len(coocur)
            if s == m:
                continue

            total_anc = 0
            total_incomp = 0
            s_anc = ancestral_dict[s]
            m_anc = ancestral_dict[m]
            pairs = 0
            both_present = min_level_dict[s] < np.Inf and min_level_dict[m] < np.Inf
            total_same = 0
            for u in s_anc:
                for v in m_anc:
                    pairs += 1
                    if u ==v:
                        total_same += 1
                    else:
                        x = is_ancestral(lin_tree.T, u,v)
                        total_anc += x 
                        if not x:
                            y = is_ancestral(lin_tree.T, v,u)
                            if not x and not y:
                                total_incomp += 1
            results.append([s, m, total_anc, total_incomp,pairs, both_present, total_same, num_occur])
    cols2 = ["mut1", "mut2", "total_ancestral", "total_incomp", "pairs", "both_present", "total_same", "num_seq_cooccur"]
    if len(results) ==0:
        data = {c: [] for c in cols2}
        relationship = pd.DataFrame(data)
    else:
        relationship = pd.DataFrame(results, columns=cols2)

    if args.pairwise is not None:
        print(f"saving {args.pairwise}......")
        relationship.to_csv(args.pairwise, index=False)

    



    print("done")



