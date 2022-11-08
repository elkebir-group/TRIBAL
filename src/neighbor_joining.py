

import numpy as np
import networkx as nx
import argparse 
import utils as ut
from tree_utils import BaseTree

def compute_q_matrix(dist_mat):
    q_matrix = np.zeros_like(dist_mat)
    n = dist_mat.shape[0]
    for i in range(n):
        for j in range(n):
            if i ==j:
                q_matrix[i,j] = np.Inf
            else:
                q_matrix[i,j] = (n-2)*dist_mat[i,j] - np.sum(dist_mat[i,:].sum()) - np.sum(dist_mat[:,j])
    return q_matrix


def update_distance_matrix(dmat, a,b, ids, next_node):

    new_row = np.zeros(len(ids)-2)
    k = 0
    for i, node in enumerate(ids):
        if node not in [ids[a], ids[b]]:

 
            new_row[k] = 0.5*(dmat[a,i] + dmat[b,i] - dmat[a,b])
            k += 1
    
    dmat2 = np.delete(dmat, [a,b], axis=0)
    dmat2 = np.delete(dmat2, [a,b], axis=1)
    dmat2 = np.vstack([dmat2, new_row.reshape(1,-1)])
 
    new_col = np.zeros(len(ids)-1)
    new_col[0:len(ids)-2] = new_row
    new_col[len(ids)-2] = 0
    dmat2 = np.hstack([dmat2, new_col.reshape(-1,1)])


    return dmat2


def neighbor_joining(dist_mat, ids):

    n = dist_mat.shape[0]
    next_node = str(n)
    tree = nx.Graph()
    center = 2*n -3 
    for i in ids:
        tree.add_edge(i, center)


    for i in range(n-3):
        q_matrix = compute_q_matrix(dist_mat)
        a, b = (np.unravel_index(q_matrix.argmin(), q_matrix.shape))
    
        tree.remove_edge(ids[a],center )
        tree.remove_edge(ids[b], center)
        tree.add_edge(next_node,ids[a])
        tree.add_edge(next_node, ids[b])
        tree.add_edge(center, next_node)
        dist_mat = update_distance_matrix(dist_mat, a,b, ids, next_node)

        to_remove = [ids[a], ids[b]]
        for val in to_remove:
            ids.remove(val)

        ids.append(next_node)
        # print(ids)
        nn = int(next_node) + 1
        next_node = str(nn)


    return tree

def root_tree(tree, root):
        tree = nx.dfs_tree(tree,root)
        # if len(list(tree.neighbors(root)))==1:
        #         tree.remove_node(root)
       
     
        # root_children = list(tree.successors(root))
        # for c in root_children:
        #     grandchildren = tree.successors(c)
        #     for g in grandchildren:
        #         tree.add_edge(root, g)

        # tree.remove_node(c)

        return tree

def hamming_distance(s1, s2):

    return (np.array(s1) != np.array(s2)).sum()

def isotype_distance(i1, i2, metric="euclidean"):
    if metric == "euclidean":
        return np.sqrt( (i1-i2)**2)
    else:
        return np.abs(i1-i2)

def create_distance_matrix(alignment, ids):

    
    dist_mat = np.array([[hamming_distance(alignment[k1], alignment[k2]) for k2 in ids] for k1 in ids])
    # print(dist_mat)
    return dist_mat.astype(float)

def weighted_distance_matrix(alignment, isotypes, ids,alpha):

    
    dist_mat = np.array([[alpha*hamming_distance(alignment[k1], alignment[k2]) + 
                    (1-alpha)*(isotype_distance(isotypes[k1], isotypes[k2])) for k2 in ids] for k1 in ids])
    # print(dist_mat)
    return dist_mat.astype(float)


if __name__ =="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", required=True, type=str,
        help="filename of input alignment fasta file")
    parser.add_argument("-r", "--root", required=True, type=str,
        help="the id of the root sequence in the alignment")
    parser.add_argument("-o", "--output", required=True, type=str,
        help="the output file of the tree")
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("-n", "--newick", type=str, help="filename where newick string should be saved")
    args= parser.parse_args()
    # path = "/scratch/projects/tribal/benchmark_pipeline/sim_data/shm_sim/2.0/0.365/1"
    # args =parser.parse_args([
    #     "-a", f"{path}/GCsim_dedup.fasta",
    #     "-r", "naive",
    #     "-o", "/scratch/projects/tribal/src/test/nj.tree",
    #     "--sequences", "/scratch/projects/tribal/src/test/nj.fasta",
    #     "-n", "/scratch/projects/tribal/src/test/nj.newick"
    # ])

    alignment = ut.read_fasta(args.alignment)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    ids = list(alignment.keys())
    rng = np.random.default_rng(args.seed)
    rng.shuffle(ids)
 
    dmat = create_distance_matrix(alignment, ids)
    tree = neighbor_joining(dmat, ids.copy())
    bt= BaseTree(tree, args.root, is_rooted=False)
    rooted_tree = bt.get_rooted_tree()
    score, labels = bt.label_nodes(alignment)
    if args.fasta is not None:
        ut.write_fasta(args.fasta, labels)
    if args.sequences is not None:
        ut.save_dict(labels, args.sequences)

    if args.output is not None:
        parents = bt.get_parents()
        ut.save_dict(parents, args.output)
    if args.newick is not None:
        newick = bt.tree_to_newick(rooted=False)
        with open(args.newick, "w+") as file:
            file.write(newick)
    # bt.save_tree_to_text(args.output)







# s1 = ["a", "b", "c"]
# s2 = ["c", "b", "c"]
# s3 = ["a", "c", "c"]
# s4 = ["a", "a", "a"]
# alignment = {1: s1, 2: s2, 3 : s3, 4: s4}
# ids = list(alignment.keys())
# print(hamming_distance(s1,s2))
# dmat = create_distance_matrix(alignment)
# tree = neighbor_joining(dmat, ids)
# print(list(tree.edges()))
# print(dmat)
# dist = np.array([[0,5,9,9,8],  [5,0,10,10,9],  [9,10,0,8,7], [9,10,8,0,3],  [8,9,7,3,0]], dtype=float)

# ids = [i for i in range(dist.shape[0])]

# tree = neighbor_joining(dist, ids)
# print(list(tree.edges))