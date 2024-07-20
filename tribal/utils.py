
import sys, re
import numpy as np
import networkx as nx
import pickle 
from io import StringIO

def read_fasta(fasta_input ):
    '''read a fasta file into a dictionary 
    https://coding4medicine.com/backup/Python/reading-fasta-files.html '''
    seq = {}
    if isinstance(fasta_input, str) and "\n" in fasta_input:
       f = StringIO(fasta_input)
    else:
        f = open(fasta_input, 'r+')
    # f=open(fname,'r+')
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


def recode_isotypes(fname, iso_dict):
    iso_encoding = {}
    counter = 0

    with open(fname, "r") as file:
        for line in file:
            isotype = line.strip()
        
            if counter ==0:
                start_iso = isotype 
         
            
            iso_encoding[isotype] = counter
            
            counter += 1
    
    isotype_encoding = {val: key for key,val in iso_encoding.items()}

    #node name as keys and isotype vales 
    iso_dict_encoded = {}
    for key, val in iso_dict.items():
         if val in iso_encoding:
            iso_dict_encoded[key] = iso_encoding[val]
         else:
            iso_dict_encoded[key] = 0

    return iso_dict_encoded, isotype_encoding
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
def read_edge_list(fname):
    edge_list = []
    with open(fname, 'r') as file:
        for line in file:
            line_list =line.strip().split(",")
            edge_list.append((line_list[0], line_list[1]))
    
    return edge_list

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

def get_alignment(fname):
    alignment = ut.read_fasta(fname)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    return alignment

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

def save_graph(fname, G, node_attribute=None):

        color_encoding =  {
                0 : "#f0f0f0",
                1 : "#FFEDA0",
                2 : "#FD8D3C",
                3 : "#E31A1C",
                4 : "#800026",
                5 : "mediumseagreen",
                6 : "#74C476",
                7 : "#6A51A3",
                8 : "darkgoldenrod",
                9 : "thistle1"
            }

        color_encoding =  {
                '0' : "#f0f0f0",
                'ABC' : "#FFEDA0",
                'Prolif. ABC/PB' : "#FD8D3C",
                'GC' : "#E31A1C",
                'FO' : "#800026",
                'MZ' : "mediumseagreen",
                'PB' : "#74C476",
                'Pre-ABC?' : "#6A51A3",
                # 8 : "darkgoldenrod",
                # 9 : "thistle1"
            }
        

        # Generate layout
        # pos = nx.spring_layout(self.G)

        # Draw the graph
        if '.pdf' in fname:
             ext = 'pdf'
        else:
            ext ='png'
        pgv_graph = nx.nx_agraph.to_agraph(G)

        if node_attribute is not None:
            node_colors = {n: color_encoding[val] for n,val in node_attribute.items()}
        else:
            node_colors  ={n: color_encoding[0] for n in G.nodes}

        for node in pgv_graph.nodes():
                if node in node_attribute:
                    pgv_graph.get_node(node).attr["fillcolor"] = node_colors[node]
                else:
                    pgv_graph.get_node(node).attr["fillcolor"] = color_encoding['0']
                pgv_graph.get_node(node).attr["style"] = "filled"
    # Draw pygraphviz graph and save as PDF
        pgv_graph.draw(fname, prog="dot", format=ext)
# tm = fit_dirichlet(2,6)
# print(tm.sum(axis=1))
