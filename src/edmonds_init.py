
import networkx as nx
from networkx.algorithms.tree.branchings import Edmonds
from networkx.algorithms.isomorphism import rooted_tree_isomorphism
import numpy as np
from itertools import combinations

class EdmondsInit:
    def __init__(self, dmat, isotypes, seed, root=0) -> None:
        
        
        self.rng = np.random.default_rng(seed)
        self.seed = seed 
        self.tree = nx.DiGraph()
 
        for n1 in isotypes:
            for n2 in isotypes:
                if n2 == root:
                    continue
                if isotypes[n1] <= isotypes[n2] and n1 != n2:
                    self.tree.add_edge(n1,n2, weight=dmat[n1,n2])
        
        self.n_edges=  len(list(self.tree.edges))
        

    
    def random_edge_removal(self, prop):
        tree = self.tree.copy()
        num_edges_to_remove = int(prop*self.n_edges)
        for i in range(num_edges_to_remove):
            cand_nodes = [n for n in self.tree.nodes if self.tree.out_degree[n] > 1]
            u = self.rng.choice(cand_nodes, 1)[0]
            cand_v = list(self.tree.successors(u))
            v = self.rng.choice(cand_v, 1)[0]
        
            tree.remove_edge(u,v)
        return tree
        
            


    def generate(self,ntrees=1, prop=0.1):
        trees = []

        for i in range(ntrees):
            if i ==0:
                in_tree = self.tree   

            else:
                self.seed += 1
                in_tree = self.random_edge_removal(prop)
     
            out_tree  =Edmonds(in_tree).find_optimum(kind="min", style="arborescence", seed=self.seed)

            trees.append(out_tree)
        
        # combos = combinations(trees, 2)
        # to_delete =[]
        # for t1, t2 in combos:
        #     if rooted_tree_isomorphism(t1,0,t2,0):
        #         print(list(t1.edges))
        #         print(list(t2.edges))
        #         to_delete.append(t1)

        # trees = [t for t in trees if t not in to_delete]
        return trees


# d_mat = np.array([[0,2,3,2], [2,0,5,1], [3,5,0,3], [2,1,3,0]])
# isos = {i: 0 for i in range(d_mat.shape[0])}

# ed = EdmondsInit(d_mat, isos, 1)

# trees = ed.generate(ntrees=5,prop=0.3)
# for t in trees:
#     print(list(t.edges))
