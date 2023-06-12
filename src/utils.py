
import sys, re
import numpy as np
import networkx as nx
import pickle 

def read_fasta(fname):
    '''read a fasta file into a dictionary 
    https://coding4medicine.com/backup/Python/reading-fasta-files.html '''
    seq = {}
    f=open(fname,'r')
    lines=f.readlines()

    hre=re.compile('>(\S+)')
    lre=re.compile('^(\S+)$')


    for line in lines:
            outh = hre.search(line)
            if outh:
                    id=outh.group(1)
            else:
                    outl=lre.search(line)
                    if(id in seq.keys()):
                            seq[id] += outl.group(1)
                    else:
                            seq[id]  =outl.group(1)
    return seq

def write_fasta(fname, mydict):
        with open(fname, 'w') as file:
                for key in mydict:
                        file.write(f">{key}\n")
                        file.write(f"{mydict[key]}\n")


# def root_tree(tree, root):
#         tree = nx.dfs_tree(tree,root)
#         if len(list(tree.neighbors(root)))==1:
#                 tree.remove_node(root)
       
     
#         # root_children = list(tree.successors(root))
#         # for c in root_children:
#         #     grandchildren = tree.successors(c)
#         #     for g in grandchildren:
#         #         tree.add_edge(root, g)

#         # tree.remove_node(c)

#         return tree

def save_dict(  mydict, fname):
        with open(fname, "w+") as file:
            for key, value in mydict.items():
                file.write(f"{key},{value}\n")

def read_trans_mat(fname):
        return np.loadtxt(fname)

def get_parent(tree, node):
        return list(tree.predecessors(node))[0]
def save_tree_to_text(fname, tree, root):

    with open(fname, "w+") as file:
    

        for e in tree.edges:
                file.write(f"{e[0]} {e[1]}\n")


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




def convert_transmat_to_weights(transmat):

        log_transmat = -1*np.log(transmat)
        states = np.arange(transmat.shape[0])
        weights = {}
        for s in states:
            for t in states:
                weights[s,t] = log_transmat[s,t]
        return weights, states

def read_dict(fname):
    mydict= {}
    with open(fname, 'r') as file:
        for line in file:
            line_list =line.strip().split(",")
            mydict[line_list[0]] = line_list[1]
    return mydict

def update_labels(labels):
   labels = {key : "".join(value) for key, value in labels.items()}
   return labels 


def pickle_save(obj, fname):
        with open(fname, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def pickle_load(fname):
    with open(fname, 'rb') as handle:
        obj= pickle.load(handle)
    
    return obj
   

def check_triangle_inequality(distance_matrix):
    n = distance_matrix.shape[0]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if distance_matrix[i, k] > distance_matrix[i, j] + distance_matrix[j, k]:
                    return False
    return True
# tm = fit_dirichlet(2,6)
# print(tm.sum(axis=1))
