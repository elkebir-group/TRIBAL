from dataclasses import dataclass, field
import networkx as nx
from small_parsimony import SmallParsimony
import numpy as np 
from draw_tree import DrawTree
from utils import save_dict
import os 
import pickle 
from pandas import DataFrame
from utils import hamming_distance, read_fasta
from collections import Counter
from ete3 import Tree 



@dataclass
class LineageTree:
    """B cell lineage tree class"""
    T: nx.DiGraph
    root: str
    id: int = 0
    name: str = None
    

    def postorder_traversal(self) -> list:
        return list(nx.dfs_postorder_nodes(self.T, source=self.root))
    

    def preorder_traversal(self) -> list:
        return list(nx.dfs_preorder_nodes(self.T, source=self.root))

    def parent(self,n):
        preds = list(self.T.predecessors(n))
        if len(preds) ==0:
            return None
    
        
        return preds[0]

    def relabel(self, dict):
       

        self.T = nx.relabel_nodes(self.T, dict)
        if self.root in dict:
            self.root = dict[self.root]
    
    def children(self,n):
        return list(self.T.neighbors(n))
    
    def nodes(self):
        return list(self.T.nodes())

    def is_leaf(self, node):
        return self.T.out_degree(node) ==0
    
    def get_leafs(self):
        return [n for n in self.T if self.is_leaf(n)]
    
    def is_root(self, node):
        if type(node) == int:
            node =str(node)
        return node == self.root
    
    def get_parents(self):
        parents = {}
        for n in self.T:
            parent = self.parent(n)
            if parent is not None:
                parents[n] = parent
        return parents
    
    def set_id(self,id):
        self.id =id 
    
    def set_name(self, name):
        self.name= name
    
    def set_tree(self, tree):
        self.T= tree
    
    def root_to_tip(self, new_root):
        kids = self.children(self.root)
  
        if len(kids) ==1:
            for k in kids:
                self.T.add_edge(new_root, k)
                self.T.remove_edge(self.root, k)
            self.T.add_edge(new_root, self.root)
        else:
            for k in kids:
                self.T.add_edge(new_root, k)
                self.T.remove_edge(self.root, k)
            self.T.add_edge(new_root, self.root)
            # raise ValueError("Expecting a root node with a single child!")

    def root_outgroup(self, outgroup):
        self.T.remove_node(outgroup)
        self.T.add_edge(outgroup, self.root)
        self.root = outgroup
    
    def prune_unifurcations(self):
        unifur = [n for n in self.T if self.T.out_degree[n]==1 and n != self.root]
        if len(unifur) > 0:
            print(f"Removing unifurcations {unifur}")
        for n in self.preorder_traversal():
            if n in unifur:
                gparent = list(self.T.predecessors(n))[0]
                child = list(self.T.neighbors(n))[0]
                self.T.remove_edge(n, child)
                self.T.remove_node(n)
            
                self.T.add_edge(gparent, child)
        return unifur

    def to_newick(self, fname=None, format=9):
        T_prev = self.T.copy()
        _ =self.prune_unifurcations()
        self.root_to_tip(new_root=-1)
    
        def get_node(nodename):
            if nodename in nodes_by_name:
                return nodes_by_name[nodename]
            else:
                nodes_by_name[nodename] = Tree(name=nodename)
                return nodes_by_name[nodename]

        nodes_by_name = {}
        for parent, child in self.T.edges:
            parent = get_node(parent)
            parent.add_child(get_node(child))
        self.T = T_prev
        tree= get_node(-1)
        if fname is not None:
            tree.write(format=format, outfile=fname)
        return tree.write(format=format)
    
    @staticmethod
    def find_leaf_descendants(node, graph):
        leaf_descendants = set()

        # Helper function to traverse the graph
        def dfs(current_node):
            nonlocal leaf_descendants
            # If the current node is a leaf, add it to the set
            if graph.out_degree(current_node) == 0:
                leaf_descendants.add(current_node)
            else:
                # Traverse all child nodes recursively
                for child_node in graph.successors(current_node):
                    dfs(child_node)

        # Start the depth-first search from the specified node
        dfs(node)
        return leaf_descendants

    def rf_distance(self, lintree):
        if type(lintree) == LineageTree:
            t1 = self.get_clade_set(self.T)
            t2 = self.get_clade_set(lintree.T)
            return (0.5*len(t1.symmetric_difference(t2)))


    def get_state_changes(self, node, labels, nisotypes):
        counts =np.zeros(shape=(nisotypes, nisotypes))
        path = nx.shortest_path(self.T, source=self.root, target=node)
        for i,n in enumerate(path):
            if i ==len(path)-1:
                break
            iso_par = labels[n]    
            iso_child = labels[path[i+1]] 
            counts[iso_par, iso_child] += (iso_child != iso_par)
        return counts 



    def pickle_tree(self, fname):
        with open(fname, 'wb') as file:
                pickle.dump(self, file)

    def get_clade_set(self, tree):
        clade_set = []
        for node in tree:
            clade_set.append(self.find_leaf_descendants(node, tree))
    
        return(set(map(frozenset, clade_set)))
    
    def get_node_levels(self):
        source_node = self.root


        node_levels = {source_node: 0}

        # Perform a BFS traversal and calculate the level of each node
        queue = [(source_node, 0)]  # (node, level)
        while queue:
            node, level = queue.pop(0)
            for neighbor in self.T.neighbors(node):
                if neighbor not in node_levels:
                    node_levels[neighbor] = level + 1
                    queue.append((neighbor, level + 1))
        return node_levels 

    def sequence_parismony(self, alignment, 
                           alphabet=None, 
                           cost_function=None,
                           compute=False):
        if not compute:
            
            return 0, alignment
        
        else:
            sp = SmallParsimony(self.T, 
                                self.root,
                                alphabet= alphabet,
                                cost = cost_function)
            seq_score, labels = sp.sankoff(alignment)
            return seq_score, labels
    
    
    def parsimony(self, alignment, iso_leaves, transMat, alphabet=None, cost=None, convert=False):
      
        states =np.arange(transMat.shape[0]).tolist()

        iso_score, iso_labels = self.isotype_parsimony_polytomy(iso_leaves, transMat,states, convert=convert)
   
        seq_score, anc_labels = self.sequence_parismony(alignment, alphabet, cost)
    
        return seq_score, iso_score, anc_labels, iso_labels

    def number_of_changes(self, labels):

        nodes = set(self.T.nodes) 
        keys = set(labels.keys())
        #check that every node has a label
        if nodes <= keys:
            return sum([hamming_distance(labels[u], labels[v]) for u,v in self.T.edges])
        else:
            print("Warning: labels were not provided for all nodes!")
            return np.NAN
    
    def isotype_parsimony(self, iso_leaves, weights, states):
   
        
        
        sp = SmallParsimony(self.T, self.root, alphabet=states, cost=weights)
        iso_score, iso_labels = sp.sankoff(iso_leaves)
        return iso_score, iso_labels 

    def isotype_parsimony_polytomy(self, iso_leaves, weights, states):

    

        sp = SmallParsimony(self.T, self.root)
        iso_score, labels, tree = sp.polytomy_resolver(iso_leaves,weights, states)

       
        return iso_score, labels 
    

    def save_png(self,fname, isotypes=None, iso_encoding=None, 
                 show_legend=False, show_labels=True, hide_underscore=True,
                 color_encoding = None):
    
        parents = self.get_parents()
        dt = DrawTree(parents, isotypes, show_legend=show_legend, root=self.root,
                       isotype_encoding=iso_encoding, show_labels=show_labels,
                       hide_underscore=hide_underscore,
                       color_encoding=color_encoding)
        dt.save(fname)

    def save_tree(self,fname):
    
        parents = self.get_parents()
        save_dict( parents, fname)
    
    def save_edges(self, fname):
        with open(fname, "w+") as file:
            for u,v in self.T.edges:
                file.write(f"{u},{v}\n")
    
    def get_edge_df(self):
        
        u_list =[u for u,v in self.T.edges]
        v_list = [v for u,v in self.T.edges]
        return DataFrame({'parent': u_list, 'child': v_list})
    
    @staticmethod
    def compute_purity(labels):
        label_counts = Counter(labels)
        majority_label = label_counts.most_common(1)[0][0]
        total_samples = len(labels)
        purity = label_counts[majority_label] / total_samples
        return purity

    @staticmethod
    def compute_entropy(labels):
        label_counts = Counter(labels)
        total_samples = len(labels)
        probabilities = np.array(list(label_counts.values())) / total_samples
        entropy = -np.sum(probabilities * np.log2(probabilities))
        return entropy


    def get_clade_nodes(self, node):
        clade_nodes = [node]  

        successors = list(self.T.successors(node))
        while successors:
            clade_nodes.extend(successors)
            successors = [child for successor in successors for child in self.T.successors(successor)]

        return clade_nodes
    
    def avg_node_score(self, func,labels):
        node_score =[]
        clade_score = {}
        for n in self.T:

            nodes = self.find_leaf_descendants(n, self.T)
            clade_labels = [labels[n] for n in nodes]
            score = func(clade_labels)
            clade_score[n] = score

            if n == self.root or self.T.out_degree[n]==1:
                continue
 
            node_score.append(score)
        
        return np.mean(node_score),clade_score
    
    def avg_purity(self, labels):
        return self.avg_node_score(self.compute_purity, labels)


    def avg_entropy(self, labels):
        return self.avg_node_score(self.compute_entropy, labels)

    
    def get_node_degrees(self):
        return {n: self.T.out_degree[n] for n in self.T}
    
    def collapse(self, labels, ignore=[]):
        leaves = [n for n in self.T if self.T.out_degree[n]==0]
        ignore = ignore + leaves
        if set(self.T.nodes) <= set(labels.keys()):

            nodes = self.preorder_traversal()
            for n in nodes:
                if n==self.root or n in ignore:
                    continue
                children = list(self.T.neighbors(n))
                
                dist = sum([hamming_distance(labels[n], labels[c]) for c  in children])
          
                if dist ==0:
                    print(f"collapsing node {n}")
                    gp= list(self.T.predecessors(n))[0]
                    self.T.remove_node(n)
                    for c in children:
              
                        self.T.add_edge(gp, c)

          
    


@dataclass
class LineageForest:
    alignment: dict = None
    isotypes: dict = None
    forest: list = field(default_factory=list)
 
 
    def generate_from_list(self, tree_list, root=None):

        for i,t in enumerate(tree_list):
            if type(t) == LineageTree:
                t.set_id(i)
                self.add(t)
                
            else:
                self.add(LineageTree(t,root,i))

    def add(self, tree):
        self.forest.append(tree)


    def __getitem__(self, key):
        return self.forest[key]
    
    
    def size(self):
        return len(self.forest)
    
    def get_trees(self):
        return self.forest
    
  
    def save_forest(self, fname):
        with open(fname, 'wb') as file:
            pickle.dump(self, file)

    def save_trees(self, path):
        for tree in self.forest:
            if tree.name != "":
                folder = tree.name
            else:
                folder = tree.id
            clono_path = f"{path}/{folder}"
            os.makedirs(clono_path, exist_ok=True)
            fname = f"{clono_path}/fit_tree.pickle"
            with open(fname, 'wb') as file:
                pickle.dump(tree, file)


#####load scripts 
from ete3 import Tree 
import re
def lintrees_from_newick(fname,root):
    '''
    takes in a filename containing line separated newick strings and returns a list of LineageTrees
    '''
    trees = read_trees(fname, root)
    lin_trees = []
    for i,t in enumerate(trees):
        lin_trees.append(LineageTree(t,root,id=i))
    
    return lin_trees 

def is_ancestral(self, u,v):
    path = nx.shortest_path(self.T, source=u, target=v)
    return len(path) > 0

def load(fname):
    with open(fname, "rb") as file:
        lt = pickle.load(file)
    
    return lt 

def read_trees(fname,root):
    '''
    Takes in a fname and a name of the root
    returns a list of networkx trees 
    '''
    print(f"\nreading trees from {fname}....")

    exp = '\[.*\]'
    trees = []
       
    with open(fname, 'r+') as file:
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

            nx_tree= convert_to_nx(ete_tree, root)
          
            trees.append(nx_tree)
        print(f"\n{len(trees)} read from {fname}!")
        return trees
    
def get_alignment(fname):
    alignment = read_fasta(fname)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    return alignment

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


        G = nx_tree.to_undirected()
        H = nx.dfs_tree(G,source=root)
   

        if H.out_degree[root_name]==0:
            H.remove(root_name)

    return H

