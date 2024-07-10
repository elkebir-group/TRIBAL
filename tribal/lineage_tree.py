from dataclasses import dataclass, field
from base_tree import BaseTree
from functools import total_ordering
import pickle
import numpy as np
from collections import Counter
from small_parsimony import SmallParsimony
from steiner_tree import ConstructGraph, SteinerTree

@dataclass
@total_ordering
class LineageTree:
    
    clonotype: str
    tree: BaseTree= None
    csr_obj: float = 0
    isotypes: dict = field(default_factory=dict)
    shm_obj: float = 0
    sequences: dict = field(default_factory=dict)

  


    def __post_init__(self):
        self.objective = (self.shm_obj,self.csr_obj)
        self.root = self.tree.root

    def _validate_item(self, item):
        if isinstance(item, (LineageTree)):
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
        return f"Objectives\tSHM: {self.shm_obj} CSR: {self.csr_obj}"
    
    def strline(self, *args):
        sep= ","
        dat = [self.tree.id, self.shm_obj, self.csr_obj]
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
            file.write(f"{self.shm_obj},{self.csr_obj}\n")
    
    def seq_len(self):
        for key, val in self.sequences.items():
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
        assert round(score,3) == round(self.csr_obj,3)
        # print(f"Score: {score} Iso Obj: {self.iso_obj}")
        return score
    
    def compute_score(self, weights):
       
        iso = self.isotypes
        score = 0

        nodes = self.tree.preorder_traversal()
        for n in nodes:
            t= iso[n]
            for c in self.tree.children(n):
                s = iso[c]
                score += weights[t,s]
        
    
        return score
    


    def ancestral_sequence_reconstruction(self, alignment, 
                           alphabet=None, 
                           cost_function=None):

        sp = SmallParsimony(self.tree, 
                            self.root,
                            alphabet= alphabet,
                            cost = cost_function)
        self.shm_obj, self.sequences = sp.sankoff(alignment)
    


    def isotype_parsimony(self, iso_leaves, cost, states):

        #TODO: add a check that the iso_leaves keys is identical to the tree leafset
   
        
    
        sp = SmallParsimony(self.tree, alphabet=states, cost=cost)
        self.csr_obj, self.isotypes = sp.sankoff(iso_leaves)
    

    def refinement(self, isotype_labels, cost, threads=1):
        #TODO: Refactor ConstructGraph & SteinerTree
        cg = ConstructGraph(cost, isotype_labels, root_identifier=self.root)
        fg = cg.build(self.tree)
        st = SteinerTree(fg.G, self.tree.T, fg.find_terminals(), fg.seq_weights, 
                             fg.iso_weights,fg.node_mapping, fg.tree_to_graph,
                               fg.node_out_degree, root=self.root, threads=threads )
        self.csr_obj, tree = st.run()

        tree, self.isotypes = cg.decodeTree(tree)
        self.tree = BaseTree(tree, self.root, self.tree.id, self.tree.name)

class LineageTreeList(list):
    def append(self, item):
        super().append(self._validate_item(item))
    
    #only Score objects are allowd in
    def _validate_item(self, item):
        if isinstance(item, (LineageTree)):
            return item
        raise TypeError(
            f"Score expected, item got {type(item).__name__}"
        )

    def write(self, fname, sep=","):
        with open(fname,'w+') as file:
            file.write(f"id{sep}objective{sep}seq_score{sep}iso_score\n")

            for score in self:
            
                file.write(f"{score.tree.id}{sep}{score.objective}{sep}{score.seq_score}{sep}{score.iso_score}{sep}\n")
    
    def find_best_scores(self):
        min_score = min(self, key=lambda x: x.iso_obj).iso_obj

        # Filter objects with the minimum score
        min_score_object = [obj for obj in self if round(obj.iso_obj, 5) == round(min_score,5)][0]
        return min_score, [min_score_object]
    
    def find_all_best_scores(self):
        min_score = min(self, key=lambda x: x.iso_obj).iso_obj

        # Filter objects with the minimum score
        min_score_object = [obj for obj in self if round(obj.iso_obj, 5) == round(min_score,5)]
        return min_score, min_score_object
    
    def get_all_trees(self):

        return [x.tree for x in self]
    
    def get_ids(self):
        return [x.tree.id for x in self]
    
    def pickle_scores(self, fname):
        with open(fname, 'wb') as file:
            pickle.dump(self, file)
    
 
    
    def sample_best_scores(self, rng=None, seed=1016):
        if rng is None:
            rng = np.random.default_rng(seed)
        _, best_scores = self.find_all_best_scores()
        sampled_index = np.random.choice(len(best_scores))
        return best_scores[sampled_index]
    

    @staticmethod
    def compute_entropy(labels):
        label_counts = Counter(labels)
        total_samples = len(labels)
        probabilities = np.array(list(label_counts.values())) / total_samples
        entropy = -np.sum(probabilities * np.log2(probabilities))
        return entropy


    def get_clade_nodes(self, node):
        clade_nodes = [node]  

        successors = list(self.T.successors(node))
        while successors:
            clade_nodes.extend(successors)
            successors = [child for successor in successors for child in self.T.successors(successor)]

        return clade_nodes
    
    def avg_node_score(self, func,labels):
        node_score =[]
        clade_score = {}
        for n in self.T:

            nodes = self.find_leaf_descendants(n, self.T)
            clade_labels = [labels[n] for n in nodes]
            score = func(clade_labels)
            clade_score[n] = score

            if n == self.root or self.T.out_degree[n]==1:
                continue
 
            node_score.append(score)
        
        return np.mean(node_score),clade_score

    def avg_entropy(self, labels):
        return self.avg_node_score(self.compute_entropy, labels)

    

    def entropy_permutation_test(self, labels, reps=1000, seed=10):
        ent_vals = []
        rng = np.random.default_rng(seed)
        labels = {key: val for key,val in labels.items()if self.is_leaf(key)}
        # print(labels)
        for i in range(reps):
            vals = [val for key, val in labels.items() if key != self.root]
     
            perm = rng.permutation(vals)
            new_labels = {}

            assert perm.shape[0] == len(labels)
            for key, lab in zip(list(labels.keys()),perm):
                new_labels[key] = lab 
       
            avg_score, clade_scores = self.avg_entropy(new_labels)
            cvals = np.array([v for k, v in clade_scores.items() if not self.is_leaf(k) and k != self.root])
     
            ent_vals.append(cvals.mean())
        return ent_vals
            




    def get_node_degrees(self):
        return {n: self.T.out_degree[n] for n in self.T}
    
    def collapse(self, labels, ignore=[]):
        leaves = [n for n in self.T if self.T.out_degree[n]==0]
        ignore = ignore + leaves
        if set(self.T.nodes) <= set(labels.keys()):

            nodes = self.preorder_traversal()
            for n in nodes:
                if n==self.root or n in ignore:
                    continue
                children = list(self.T.neighbors(n))
                
                dist = sum([hamming_distance(labels[n], labels[c]) for c  in children])
          
                if dist ==0:
                    print(f"collapsing node {n}")
                    gp= list(self.T.predecessors(n))[0]
                    self.T.remove_node(n)
                    for c in children:
              
                        self.T.add_edge(gp, c)
    
    def get_state_changes(self, node, labels, nisotypes):
        counts =np.zeros(shape=(nisotypes, nisotypes))
        path = nx.shortest_path(self.T, source=self.root, target=node)
        for i,n in enumerate(path):
            if i ==len(path)-1:
                break
            iso_par = labels[n]    
            iso_child = labels[path[i+1]] 
            counts[iso_par, iso_child] += (iso_child != iso_par)
        return counts 

    
    
    def parsimony(self, alignment, iso_leaves, transMat, alphabet=None, cost=None, convert=False):
      
        states =np.arange(transMat.shape[0]).tolist()

        iso_score, iso_labels = self.isotype_parsimony_polytomy(iso_leaves, transMat,states, convert=convert)
   
        seq_score, anc_labels = self.sequence_parismony(alignment, alphabet, cost)
    
        return seq_score, iso_score, anc_labels, iso_labels


        
# def load( fname):
#     with open(fname,'rb') as file:
#             obj = pickle.load(file)
#     return obj
    
    
