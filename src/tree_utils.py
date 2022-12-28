
import networkx as nx
import numpy as np 
from ete3 import Tree
from small_parsimony import SmallParsimony
import utils as ut
class BaseTree:
    def __init__(self, tree=None, root="naive", 
                is_rooted=False, ids=None, n=7,
                rng = None) -> None:
        
        self.root = root
        if rng is None:
            self.rng = np.random.default_rng()
        else:
            self.rng  = rng
        if tree is None:
            if ids is None:
                ids = [str(i) for i in range(n)]
                ids.append(root)
            self.T = self.random_tree(ids)
            self.root_tree()
        else:
            if is_rooted:
                self.rooted_T=tree
                self.unroot_tree()
            else:
                self.T = tree
                self.root_tree()

    

    def root_tree(self):
        self.rooted_T = nx.dfs_tree(self.T,self.root)

        #remove unifurcation of root node by removing its child and adding edges from root to grandkids
        # root_kids = self.get_children(self.root)
        # if len(root_kids) == 1:
        #     child = root_kids[0]
        #     grandkids = self.get_children(child)
        #     self.rooted_T.remove_node(child)
        #     for g in grandkids:
        #         self.rooted_T.add_edge(self.root, g)

    
    def get_rooted_tree(self):
        return self.rooted_T

    def get_parent(self, node):
        parents = list(self.rooted_T.predecessors(node))
        if len(parents) > 0:
            return parents[0]
        return None
    
    def get_children(self,node):
        return list(self.rooted_T.neighbors(node))
    
    def is_leaf(self, node):
        return len(self.get_children(node)) ==0
    
    def save_tree_to_text(self, fname):

        with open(fname, "w+") as file:
    
            for e in self.rooted_T.edges:
                file.write(f"{e[0]} {e[1]}\n")
    
    def tree_to_newick(self, rooted=False):
        if rooted:
            return  "(" + self.root + "," + self.tree_to_newick_helper(self.rooted_T, self.root) + ");"
        else:
            return self.tree_to_newick_helper(self.rooted_T, self.root) + ";"

    def random_tree(self, ids):
        tree = nx.Graph()
     
   
        n = len(ids)
        center = str(2*n -3 )
        for i in ids:
            tree.add_edge(i, center)
      
        next_node = n
       
        for i in range(n-3):
            pair = self.rng.choice(ids,2, replace=False)
    

            for p in pair:
                tree.add_edge(p, str(next_node))
                tree.remove_edge(p, center)
                tree.add_edge(center, str(next_node))
                ids.remove(p)
            ids.append(str(next_node))
            next_node += 1
         
  
        return tree
    
    def has_polytomies(self):
        for n in self.rooted_T.nodes:
            children = self.get_children(n)
            if len(children) > 2:
                return True 
        return False 
    
    def get_next_node(self):
        next_node = len(self.T.nodes) +1
        while True:
            if str(next_node) not in self.T.nodes:
                return str(next_node)
            else:
                next_node += 1

    
    def resolve_polytomies_randomly(self):
        polytomies = [n for n in self.rooted_T.nodes if len(self.get_children(n)) > 2 ]
        for p in polytomies:
         
            while True:
                children = self.get_children(p)
                if len(children) <=2:
                    break
                siblings = np.random.choice(children, 2, replace=False )
                next_node = self.get_next_node()
                self.rooted_T.add_edge(p, next_node)
                for s in siblings:
                    self.rooted_T.remove_edge(p, s)
                  
                    self.rooted_T.add_edge(next_node, s)
               
     
        self.unroot_tree()
            



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

    def set_rooted_tree(self, rooted_tree):
        self.rooted_T = rooted_tree
        self.unroot_tree()
    
    def set_unrooted_tree(self, tree):
        self.T =tree 
        self.root_tree()
    
    def get_unrooted_tree(self):
        return self.T.copy()
    
    def unroot_tree(self):
        '''change tree from digraph to undirected graph
        '''
        self.T = nx.to_undirected(self.rooted_T)
    
    def num_nodes(self):
        return len(list(self.rooted_T.nodes))
        # n = len(self.T.nodes) + 1
        # while True:
        #     if n in self.T:
        #         n+= 1
        #     else:
        #         break
        # mapping = {self.root : n}
        # self.T= nx.relabel_nodes(self.T, mapping)
        # self.T.add_edge(self.root, n)

class TribalTree(BaseTree):
    def __init__(self, tree=None, root="naive", is_rooted=False, ids=None, n=7, rng=None) -> None:
    
        super().__init__( tree, root, is_rooted, ids, n, rng)
        


        self.seq_score = np.Inf 
        self.iso_score = np.Inf
  
        # self.nodes = list(self.rooted_T.nodes)
    
        # self.leaves = [n for n in self.rooted_T.nodes if self.rooted_T.out_degree(n)==0]
        # self.ntaxa = len(self.leaves)

    

    def __str__(self):
            mystr = str(self.ntaxa) + " taxa\n"
            mystr+=str(list(self.tree.edges()))
            return mystr

    
    def sequence_parismony(self, alignment, alphabet=None, cost_function=None):

     

        sp = SmallParsimony(self.rooted_T, 
                            self.root,
                            alphabet= alphabet,
                            cost = cost_function)
        self.seq_score, labels = sp.sankoff(alignment)
        return labels
    
    
    def parsimony(self, alignment, iso_leaves, transMat, alphabet=None, cost=None, ighd=None):
      
        states =np.arange(transMat.shape[0]).tolist()
        # if self.has_polytomies():
        iso_labels = self.isotype_parsimony_polytomy(iso_leaves, transMat,states, ighd=ighd)
        # else:
        #     #perform regular parsimony for binary trees
        #     pass 
        anc_labels = self.sequence_parismony(alignment, alphabet, cost)
        # iso_labels = self.isotype_ml(iso_leaves, transMat)
        # print(self.isotypes)
        return self.seq_score, self.iso_score, anc_labels, iso_labels

    
    def isotype_parsimony_polytomy(self, iso_leaves, transmat, states, ighd=None):
        log_transmat = -1*np.log(transmat)
        weights = {}
        for s in states:
            for t in states:
                weights[s,t] = log_transmat[s,t]
       
        best_iso = np.Inf
        best_labels = None
        while True:
            tree = self.rooted_T.copy()
            sp = SmallParsimony(tree, self.root)
            iso_score, labels, tree = sp.polytomy_resolver(iso_leaves,weights, states, ighd=ighd)
            # if iso_score >= best_iso:
            #     break 
            # else:
            best_labels = labels 
            self.iso_score = iso_score 
            best_iso = iso_score 
            self.set_rooted_tree(tree)
            break
        return best_labels 

        

    def isotype_parsimony(self, transMat):
      
            sp = SmallParsimony(self.rooted_T, self.root)
            self.iso_score, self.isotypes = sp.fitch(self.iso_leaves, transMat)
    
    def isotype_ml(self, iso_leaves, transMat):
   
        sp = SmallParsimony(self.rooted_T, self.root)
        self.iso_score, iso_labels = sp.fastml(iso_leaves, transMat)
        return iso_labels
        

    def tree_score(self, alpha):
        return alpha*self.get_score() + (1-alpha)*self.get_iso_score()
        
    def get_score(self):
        return self.seq_score

    def get_iso_score(self):
        return self.iso_score
    
    def ntaxa(self):
        return self.ntaxa

 


        
    # def get_output(self, alignment, isotypes, transMat, alphabet=None, cost=None):
    #     labels = self.parsimony(alignment, isotypes, transMat, alphabet, cost)

    #     return labels, isotypes, self.get_score(), self.get_iso_score() 