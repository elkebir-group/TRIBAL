
import sys, re
import numpy as np
import networkx as nx

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
        with open(fname, 'w+') as file:
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

def save_dict(mydict, fname):
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

# tm = fit_dirichlet(2,6)
# print(tm.sum(axis=1))
