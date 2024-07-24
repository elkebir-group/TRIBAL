from dataclasses import dataclass, field
import networkx as nx

import numpy as np 
from .draw_tree import DrawTree
import os 
import pickle 
from pandas import DataFrame
from .utils import hamming_distance, save_dict
from ete3 import Tree 



@dataclass
class BaseTree:
    """
    B cell lineage tree class
    """
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
        if type(lintree) == Tree:
            t1 = self.get_clade_set(self.T)
            t2 = self.get_clade_set(lintree.T)
            return (0.5*len(t1.symmetric_difference(t2)))





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



    def number_of_changes(self, labels):

        nodes = set(self.T.nodes) 
        keys = set(labels.keys())
        #check that every node has a label
        if nodes <= keys:
            return sum([hamming_distance(labels[u], labels[v]) for u,v in self.T.edges])
        else:
            print("Warning: labels were not provided for all nodes!")
            return np.NAN
    


    # def isotype_parsimony_polytomy(self, iso_leaves, weights, states):

    

    #     sp = SmallParsimony(self.T, self.root)
    #     iso_score, labels, tree = sp.polytomy_resolver(iso_leaves,weights, states)

       
    #     return iso_score, labels 
    

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
    


   
          
    


@dataclass
class Clonotype:
    """
    Clonotype data 

    Attributes
    ----------

    id: str
        the label of the clonotype
    
    alignment: dict
        the multiple sequence alignment for the BCR sequences

    isotypes: dict
        the encoded isotypes of the sequenced B cells
    
    forest: list
        a list of networkx.Digraphs containing the maximum parsimony forest for a 
        the multiple sequence alignment
    
    mapping: dict
        a mapping of the internal ids to the original cell ids

    """
    id: str   #id label of the clonotype
    alignment: dict = None #dictiona
    isotypes: dict = None
    forest: list = field(default_factory=list)
    mapping: dict = field(default_factory=dict)
 


    def generate_from_list(self, tree_list, root=None):

        for i,t in enumerate(tree_list):
            if type(t) == BaseTree:
                t.set_id(i)
                self.add(t)
                
            else:
                self.add(BaseTree(t,root,i))

    def add(self, tree):
        self.forest.append(tree)


    def __getitem__(self, key):
        return self.forest[key]
    
    
    def size(self):
        return len(self.forest)
    
    def get_forest(self):
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




