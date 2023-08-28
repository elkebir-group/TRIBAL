from dataclasses import dataclass, field
import networkx as nx
from small_parsimony import SmallParsimony
import numpy as np 
from draw_tree import DrawTree
from utils import save_dict
import os 
import pickle 
import pygraphviz as pgv
from pandas import DataFrame
from utils import hamming_distance
from collections import Counter



@dataclass
class LineageTree:
    """Class for keeping track of lineage trees"""
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
    
    def children(self,n):
        return list(self.T.neighbors(n))
    
    def nodes(self):
        return list(self.T.nodes())

    def is_leaf(self, node):
        return self.T.out_degree(node) ==0
    
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
    
    def root_to_tip(self):
        kids = self.children(self.root)
  
        if len(kids) ==1:
            new_root = "root_0"
            for k in kids:
                self.T.add_edge(new_root, k)
                self.T.remove_edge(self.root, k)
            self.T.add_edge(new_root, self.root)

    def root_outgroup(self, outgroup):
        self.T.remove_node(outgroup)
        self.T.add_edge(outgroup, self.root)
        self.root = outgroup
    
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

    def sequence_parismony(self, alignment, alphabet=None, cost_function=None):

        

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
    # @staticmethod
    # def convert_transmat_to_weights(transmat):

    #     log_transmat = -1*np.log(transmat)
    #     states = np.arange(transmat.shape[0])
    #     weights = {}
    #     for s in states:
    #         for t in states:
    #             weights[s,t] = log_transmat[s,t]
    #     return weights, states
    
    def isotype_parsimony(self, iso_leaves, weights, states):
   
        
        
        sp = SmallParsimony(self.T, self.root, alphabet=states, cost=weights)
        iso_score, iso_labels = sp.sankoff(iso_leaves)
        return iso_score, iso_labels 

    def isotype_parsimony_polytomy(self, iso_leaves, weights, states):

    

        sp = SmallParsimony(self.T, self.root)
        iso_score, labels, tree = sp.polytomy_resolver(iso_leaves,weights, states)

       
        return iso_score, labels 
    

    def save_png(self,fname, isotypes, iso_encoding=None, show_legend=False):
    
        parents = self.get_parents()
        dt = DrawTree(parents, isotypes, show_legend=show_legend, isotype_encoding=iso_encoding)
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
        clade_nodes = [node]  # Start with the provided node

        successors = list(self.T.successors(node))
        while successors:
            clade_nodes.extend(successors)
            successors = [child for successor in successors for child in self.T.successors(successor)]

        return clade_nodes
    
    def avg_node_score(self, func,labels):
        node_score =[]

        for n in self.T:
            #skip the germline unifurication
            if n == self.root and self.T.out_degree[n]==1:
                continue
            # nodes = self.get_clade_nodes(n)
            nodes = self.find_leaf_descendants(n, self.T)
            clade_labels = [labels[n] for n in nodes]
            purity = func(clade_labels)
 
            node_score.append(purity)
        
        return np.mean(node_score)
    
    def avg_purity(self, labels):
        return self.avg_node_score(self.compute_purity, labels)

    def avg_entropy(self, labels):
        return self.avg_node_score(self.compute_entropy, labels)

    
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
 
 

    # def __post_init__(self):
     
    #     for i,tree in enumerate(self.forest):
    #         tree.set_id(i)
    
    def generate_from_list(self, tree_list, root=None):

        for i,t in enumerate(tree_list):
            if type(t) == LineageTree:
                t.set_id(i)
                self.add(t)
                
            else:
                self.add(LineageTree(t,root,i))

    def add(self, tree):
        self.forest.append(tree)
        # tree.set_id(len(self.forest))

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


