
import networkx as nx
import numpy as np 
from ete3 import Tree
from small_parsimony import SmallParsimony
import utils as ut
class BaseTree:
    def __init__(self, tree, root, is_rooted=False) -> None:
        
        
        self.root = root
        if is_rooted:
            self.rooted_T=tree
            self.T = tree
        else:
            self.T = tree
            self.root_tree()

    

    def root_tree(self):
        self.rooted_T = nx.dfs_tree(self.T,self.root)

        #remove unifurcation of root node by removing its child and adding edges from root to grandkids
        root_kids = self.get_children(self.root)
        if len(root_kids) == 1:
            child = root_kids[0]
            grandkids = self.get_children(child)
            self.rooted_T.remove_node(child)
            for g in grandkids:
                self.rooted_T.add_edge(self.root, g)

    
    def get_rooted_tree(self):
        return self.rooted_T

    def get_parent(self, node):
        parents = list(self.rooted_T.predecessors(node))
        if len(parents) > 0:
            return parents[0]
        return None
    
    def get_children(self,node):
        return list(self.rooted_T.neighbors(node))
    
    def save_tree_to_text(self, fname):

        with open(fname, "w+") as file:
    
            for e in self.rooted_T.edges:
                file.write(f"{e[0]} {e[1]}\n")
    
    def tree_to_newick(self, rooted=False):
        if rooted:
            return  "(" + self.root + "," + self.tree_to_newick_helper(self.rooted_T, self.root) + ");"
        else:
            return self.tree_to_newick_helper(self.rooted_T, self.root) + ";"

   

    #https://stackoverflow.com/questions/46444454/save-networkx-tree-in-newick-format
 
    def tree_to_newick_helper(self,g, root):
    # if root is None:
    #     roots = list(filter(lambda p: p[1] == 0, g.in_degree()))
    #     assert 1 == len(roots)
    #     root = roots[0][0]
        subgs = []
        for child in g[root]:
            if len(g[child]) > 0:
                subgs.append(self.tree_to_newick_helper(g, root=child))
            else:
                subgs.append(str(child))
        return "(" + ','.join(subgs) + ")"  



    def label_nodes(self, alignment, fname=None, cost=None):
        sp = SmallParsimony(self.rooted_T, self.root)
        score, labels = sp.sankoff(alignment)

        for key in labels:
            labels[key] = "".join(labels[key])

        if fname is not None:
            ut.save_dict(labels, fname)

        return score, labels
    
    def networkx_to_ete3(self):
        pass 

        # if len(list(tree.neighbors(root)))==1:
        #         tree.remove_node(root)
       
     
        # root_children = list(tree.successors(root))
        # for c in root_children:
        #     grandchildren = tree.successors(c)
        #     for g in grandchildren:
        #         tree.add_edge(root, g)

        # tree.remove_node(c)

        # return tree

    def get_parents(self):
        parents = {}
        for n in self.rooted_T.nodes:
            parent = self.get_parent(n)
            if parent is not None:
                parents[n] = parent
        return parents
