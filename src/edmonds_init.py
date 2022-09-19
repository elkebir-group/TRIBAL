
import networkx as nx
from networkx.algorithms.tree.branchings import Edmonds
import numpy as np

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
            u = rng.choice(cand_nodes, 1)[0]
            cand_v = list(self.tree.successors(u))
            v = rng.choice(cand_v, 1)[0]
        
            tree.remove_edge(u,v)
        
            


    def generate(self,ntrees=1, prop=0.1):
        trees = []

        for i in range(ntrees):
            if i ==0:
                in_tree = self.tree   

            else:
                self.seed += 1
                in_tree = self.random_edge_removal(prop)
     
            out_tree  =Edmonds(in_tree).find_optimum(kind="min", style="arborescence", seed=self.seed)



        return trees


d_mat = np.array([[0,2,3], [2,0,5], [3,5,0]])
isos = {i: 0 for i in range(3)}

ed = EdmondsInit(d_mat, isos, 1)
ed.generate()
