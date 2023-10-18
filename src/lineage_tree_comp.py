import numpy as np
import networkx as nx 
    
class LinTreeComp:
    def __init__(self, gt_score, inf_score) -> None:
        self.gt_iso = gt_score.isotypes
        self.inf_iso = inf_score.isotypes 

        self.gt_tree  = gt_score.tree
        self.inf_tree = inf_score.tree

        gt_leafs = self.gt_tree.get_leafs()
        mapping = {}
        for l in self.gt_tree.T:
            if 'seq' in l:
                seqs = l.split("_")
                if len(seqs) > 0:
                    mapping[l] = seqs[0]
        
        self.gt_tree.relabel(mapping)
        for key, val in mapping.items():
            old = self.gt_iso[key]
            del self.gt_iso[key]
            self.gt_iso[val] = old 




    @staticmethod 
    def compute_iso_changes(lin_tree, iso, tgt):
        path = nx.shortest_path(lin_tree.T, source=lin_tree.root, target=tgt)
        total =0
        for i in range(1,len(path)):
            u = path[i-1]
            v = path[i]
            total += iso[u] != iso[v]
        return total
            


    def num_isotype_changes(self, l):

        leafs = self.inf_tree.get_leafs()
        if l not in leafs:
            raise ValueError(f"Leaf {l} is not in inferred tree.")


        gt_total = self.compute_iso_changes(self.gt_tree, self.gt_iso, l)
        inf_total =self.compute_iso_changes(self.inf_tree, self.inf_iso, l)
        return gt_total, inf_total, np.abs(gt_total- inf_total)
   

            