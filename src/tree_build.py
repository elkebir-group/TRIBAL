
from os import stat
from re import I
import networkx as nx
from neighbor_joining import create_distance_matrix
from small_parsimony import SmallParsimony
from tree_moves import TreeMoves
import numpy as np
import neighbor_joining as nj

class Tribal:
    def __init__(self, alignment, isotype_labels, seed=1026, root=0):
        
        self.rng = np.random.default_rng(seed)
        self.alignment = alignment
        self.n = len(self.alignment)
        self.root = root
        self.compute_character_sets()
        self.isotype_labels = isotype_labels
    # def __init__(self, sequences, root, switch_graph=None, cost_bases=None, costs_isotypes=None, seed=1234) -> None:
    #     self.sequences = {key : [lower(i) for i in sequences] for key in sequences}
    #     self.switch_graph = switch_graph 

    #     self.T = self.initialize(seed)
    #     self.bases = set("a","c", "g", "t")
    #     self.isotypes = set(switch_graph.nodes())
        self.ISO= "isotype"
        self.SEQ = "sequence"
    
        self.alphabet = ("A", "C", "G", "T", "-")
        self.cost_function = {}

        for i in self.alphabet:
            for j in self.alphabet:
                if (i,j) not in self.cost_function:
                    if i == "-" and j == "-":
                        self.cost_function[i,j] =0
                    if i == "-" or j == "-":
                        self.cost_function[i,j] = 1
                    elif  i ==j:
                        self.cost_function[i,j] = 0
                    
                    else:
                        self.cost_function[i,j] = 1


        self.root = 0

    @staticmethod
    def print_tree(res):
        tree = res[1]
        print(f" OPT Seq: {res[0]} OPT Iso: {res[2]}")
        alignment = nx.get_node_attributes(tree, "sequence")
        isotype_labels =  nx.get_node_attributes(tree, "isotype")
        print(list(tree.edges()))
        for n in tree.nodes():
            seq  ="".join([i for i in alignment[n]])
            print(f"Node {n}: Seq: {seq}, isotype: {isotype_labels[n]}")

    @staticmethod
    def jaccard_distance(a,b):
        """
        computes the jaccard distance between to numpy arrays a and b

        """
        return (np.intersect1d(a,b).shape[0])/np.union1d(a,b).shape[0]

    def compute_character_sets(self):
        self.char_sets = {self.root : np.empty(shape=0, dtype=int)}
        root = np.array(self.alignment[self.root])
        for i in self.alignment:
            if i != self.root:
                self.char_sets[i] = (np.array(self.alignment[i]) != root).nonzero()
        
        self.dmat= np.zeros((self.n, self.n))
        for i in self.alignment:
            for j in self.alignment:
                self.dmat[i,j] = self.jaccard_distance(self.char_sets[i], self.char_sets[j])
        
        

        


    def assign_sequences(self, seqs, parents):
        tribals = []
        #subset distance matrix
        dmat_cluster = np._ix((seqs, parents))
        assign_index = dmat_cluster.argmin(axis=1)
        seq_assign = parents[assign_index]
        for p in parents:
            alignment = {p : self.alignment[p]}
            for i,s in seqs:
                if seq_assign[i] == p:
                
                    alignment[s] = self.alignment[s]
            
            tribals.append(Tribal(alignment, root=p))
        return tribals 
    

    

        









    def hill_climbing(self, tree):
        
            best_score, best_tree = SmallParsimony(tree,  self.root, self.SEQ, self.alphabet,self.cost_function).sankoff()
            best_iso_score, best_tree = SmallParsimony(best_tree, self.root, self.ISO).fitch()

            while True:
                one_spr = TreeMoves(best_tree, self.root).spr()
                improvement = False
           
                for t in one_spr:
              
                    cand_score, t = SmallParsimony(t, self.root, self.SEQ, self.alphabet,self.cost_function).sankoff()
                 
                
                    if cand_score < best_score: #or (cand_score == best_score and cand_iso_score < best_iso_score):
                        # print(list(t.edges))
                        seq_labels = nx.get_node_attributes(t, self.SEQ)
                        cand_iso_score, t = SmallParsimony(t, self.root, self.ISO).fitch()
                        improvement = True
                        best_score = cand_score
                        best_tree = t
                        best_iso_score = cand_iso_score

                
                if not improvement:
                    break
                # else:
                #     iso_score, best_tree = SmallParsimony(best_tree, 0, "state").fitch()

            return best_score, best_tree, best_iso_score

    @staticmethod
    def root_tree(tree, root):
        tree = nx.dfs_tree(tree,root)
     
        root_children = list(tree.successors(root))
        for c in root_children:
            grandchildren = tree.successors(c)
            for g in grandchildren:
                tree.add_edge(root, g)

        tree.remove_node(c)

        return tree

    def initialize_trees(self, states=None, seed=None, num_trees=None):

        dmat = nj.create_distance_matrix(self.alignment)
        tree = nj.neighbor_joining(dmat, list(self.alignment.keys()))
        #test 
        # mytree = nx.Graph()
        # mytree.add_edges_from([(0,5), (5,1), (4,5), (2,4), (3,4)])
        # mytree = self.root_tree(mytree, self.root)

        tree= self.root_tree(tree, self.root)
        nx.set_node_attributes(tree, self.isotype_labels, self.ISO)
        nx.set_node_attributes(tree, self.alignment, self.SEQ)
    
  
        return tree

    def initialize_random(self,n):
        return [self.random_tree()  for i in range(n)]


    def random_tree(self):
        tree = nx.Graph()
        ids = list(self.alignment.keys())
        n = len(ids)
        center = 2*n -3 
        for i in range(n):
            tree.add_edge(i, center)
      
        next_node = n
       
        for i in range(n-3):
            pair = self.rng.choice(ids,2, replace=False)
    

            for p in pair:
                tree.add_edge(p, next_node)
                tree.remove_edge(p, center)
                tree.add_edge(p, next_node)
                tree.add_edge(center, next_node)
                ids.remove(p)
            ids.append(next_node)
            next_node += 1
        

        tree= self.root_tree(tree, self.root)
        nx.set_node_attributes(tree, self.isotype_labels, self.ISO)
        nx.set_node_attributes(tree, self.alignment, self.SEQ)
     
    
  
        return tree

    

    def run(self, seed=None,n=None):
        best_score = np.Inf
        local_opt = []
        # start_tree = self.initialize_trees()

        start_trees= self.initialize_random(5)
    
        
        for s in start_trees:
            local_opt.append(self.hill_climbing(s))
        
        best_res = []
        for l in local_opt:
            if l[0] == best_score:
                best_res.append(l)

            elif l[0] < best_score: 
                best_res = [l]
                best_score = l[0]
         
        
        for l in best_res:
            self.print_tree(l)
      
        

        return best_res
            
# s1 = ["a", "b", "c"]
# s2 = ["c", "b", "c"]
# s3 = ["a", "c", "c"]
# s4 = ["a", "a", "a"]
# alignment = {1: s1, 2: s2, 3 : s3, 4: s4}

seq0 = [ "A", "C", "G", "G"]
seq1 = [ "C", "T", "A", "G"]
seq2 = [ "C", "T", "G", "G"]
seq3 = [ "T", "T", "G", "A"]
seq4 = [ "A", "C", "A", "C"]
seq5 = [ "A", "C", "T", "C"]
isotype_labels = {1: 1, 2: 3, 3:2, 4:0, 5:1}
alignment = {0: seq0, 1: seq1, 2: seq2, 3: seq3, 4: seq4, 5: seq5}
ids = list(alignment.keys())
tr  = Tribal(alignment, isotype_labels)
opt_trees = tr.run()
print(len(opt_trees))
# print(f" OPT Seq: {opt_seq_score} OPT Iso: {opt_iso_score}")
# tr.print_tree()







