from dataclasses import dataclass, field
from lineage_tree import LineageTree
from functools import total_ordering
import pickle

@dataclass
@total_ordering
class Score:
    objective: float
    seq_obj: float 
    iso_obj: float 
    labels: dict 
    isotypes: dict 
    tree: LineageTree= None


    def _validate_item(self, item):
        if isinstance(item, (Score)):
            return item
        raise TypeError(
            f"Score expected, item got {type(item).__name__}"
        )
    def __eq__(self, __value: object) -> bool:
        item = self._validate_item(__value)
        self.objective ==item.objective
    
    def __lt__(self, __value: object) -> bool:
        item = self._validate_item(__value)
        self.objective <item.objective

    def improvement(self, cand_score ):

        return self.objective < cand_score.objective
    

    def __str__(self):
        return f"Objective: {self.objective} Sequence: {self.seq_obj} Isotype: {self.iso_obj}"
    
    def strline(self, *args):
        sep= ","
        dat = [self.tree.id, self.objective, self.seq_obj, self.iso_obj]
        for arg in args:
            dat.append(arg)
   
       
        dat = [str(d) for d in dat]
        all_data = sep.join(dat)
        return f"{all_data}\n"

    
    def pickle_score(self, fname):
        with open(fname, 'wb') as file:
                pickle.dump(self, file)
                
    def save_score(self, fname, alpha):
         with open(fname, 'w+') as file:
            file.write(f"{alpha},{self.objective},{self.seq_obj},{self.iso_obj}\n")
    
    def seq_len(self):
        for key, val in self.labels.items():
            return len(val)
    
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
        # print(f"Score: {score} Iso Obj: {self.iso_obj}")
        return score


class ScoreList(list):
    def append(self, item):
        super().append(self._validate_item(item))
    
    #only Score objects are allowd in
    def _validate_item(self, item):
        if isinstance(item, (Score)):
            return item
        raise TypeError(
            f"Score expected, item got {type(item).__name__}"
        )

    def write(self, fname, sep=","):
        with open(fname,'w+') as file:
            file.write(f"id{sep}objective{sep}seq_score{sep}iso_score\n")

            for score in self:
            
                file.write(f"{score.tree.id}{sep}{score.objective}{sep}{score.seq_score}{sep}{score.iso_score}{sep}\n")
    
