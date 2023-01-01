from dataclasses import dataclass, field

@dataclass
class Score:
    objective: float
    seq_obj: float 
    iso_obj: float 
    labels: dict 
    isotypes: dict 



    def improvement(self, cand_score ):
        return self.objective < cand_score.objective
    

    def __str__(self):
        return f"Objective: {self.objective} Sequence: {self.seq_obj} Isotype: {self.iso_obj}"
    
    def save_score(self, fname, alpha):
         with open(fname, 'w+') as file:
            file.write(f"{alpha},{self.objective},{self.seq_obj},{self.iso_obj}\n")
