from dataclasses import dataclass, field
from lineage_tree import LineageTree

@dataclass
class Score:
    objective: float
    seq_obj: float 
    iso_obj: float 
    labels: dict 
    isotypes: dict 
    tree: LineageTree= None



    def improvement(self, cand_score ):

        return self.objective < cand_score.objective
    
    # def __lt__(self, cand_score):
    #     return isinstance(cand_score, Score) and self.
    def __str__(self):
        return f"Objective: {self.objective} Sequence: {self.seq_obj} Isotype: {self.iso_obj}"
    
    def save_score(self, fname, alpha):
         with open(fname, 'w+') as file:
            file.write(f"{alpha},{self.objective},{self.seq_obj},{self.iso_obj}\n")
    
    def check_score(self, weights):
       
        iso = self.isotypes
        score = 0

        nodes = self.tree.preorder_traversal()
        for n in nodes:
            t= iso[n]
            for c in self.tree.children(n):
                s = iso[c]
                score += weights[t,s]
        assert round(score,3) == round(self.iso_obj,3)
        print(f"Score: {score} Iso Obj: {self.iso_obj}")
        return score


