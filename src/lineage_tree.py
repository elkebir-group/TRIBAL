from dataclasses import dataclass, field
import networkx as nx
from small_parsimony import SmallParsimony
import numpy as np 
from draw_tree import DrawTree
from utils import save_dict
import os 
import pickle 

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
        if n == self.root:
            return None
        return list(self.T.predecessors(n))[0]
    
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
    

    def save_png(self,fname, isotypes, iso_encoding=None):
    
        parents = self.get_parents()
        dt = DrawTree(parents, isotypes, show_legend=False, isotype_encoding=iso_encoding)
        dt.save(fname)

    def save_tree(self,fname):
    
        parents = self.get_parents()
        save_dict( parents, fname)
    


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