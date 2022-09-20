
import networkx as nx
from neighbor_joining import create_distance_matrix
from small_parsimony import SmallParsimony
from tree_moves import TreeMoves
import numpy as np
import neighbor_joining as nj
from edmonds_init import EdmondsInit
import argparse 
import pandas as pd 

class Tribal:
    def __init__(self, alignment, root, seed, isotype_labels, n_init=1, init="random"):
        
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.alignment = alignment
        self.n = len(self.alignment)
        self.root = root
        # self.compute_character_sets()
        self.isotype_labels = isotype_labels

        self.ids = list(alignment.keys())
        self.ids.sort()
        

    
        self.n_init = n_init 
        
        if init=="nj":
            self.init = self.initialize_nj
        elif init == "edmonds":
            self.init = self.initialize_edmonds
        else:
            self.init = self.initialize_random

        self.ISO= "ISO"
        self.SEQ = "SEQ"
    
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

        

    def hill_climbing(self, tree):
        
            best_score, best_tree = SmallParsimony(tree,  self.root,
                                                    attribute  =self.SEQ, 
                                                    sequences = self.alignment, 
                                                    alphabet= self.alphabet,
                                                    cost = self.cost_function).sankoff()
            best_iso_score, best_tree = SmallParsimony(best_tree, self.root, 
                                                        attribute= self.ISO,
                                                        isotypes=self.isotype_labels).fitch()

            while True:
                one_spr = TreeMoves(best_tree, self.root, iso_att_name=self.ISO).spr()
                improvement = False
           
                for t in one_spr:
              
                    cand_score, t = SmallParsimony(t, 
                                                    self.root, 
                                                    attribute=self.SEQ,
                                                    sequences = self.alignment, 
                                                    alphabet= self.alphabet,
                                                    cost=self.cost_function).sankoff()
                 
                
                    if cand_score < best_score: #or (cand_score == best_score and cand_iso_score < best_iso_score):
                        # print(list(t.edges))
                        seq_labels = nx.get_node_attributes(t, self.SEQ)
                        cand_iso_score, t = SmallParsimony(t, self.root,
                                                            attribute =self.ISO,
                                                            isotypes=self.isotype_labels).fitch()
                        improvement = True
                        best_score = cand_score
                        best_tree = t
                        best_iso_score = cand_iso_score

                
                if not improvement:
                    break
 

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



    def initialize_random(self, **kwargs):
        
        return [self.random_tree()  for i in range(self.n_init)]
    

    def initialize_nj(self,  **kwargs):

        dmat= kwargs['dmat']

        tree = nj.neighbor_joining(dmat, self.ids)

        tree= self.root_tree(tree, self.root)
        return [tree]

    def initialize_edmonds(self, **kwargs):
        dmat= kwargs['dmat']
        trees = EdmondsInit(dmat, self.isotype_labels, seed=1026, root=self.root).generate(self.n_init, prop=0.2)

       #post process Edmonds Trees so all sequences are leaves and not internal nodes

        for t in trees:
            print(list(t.edges))
            next_label = max(t.nodes) +1
            nodes = list(t.nodes)
            for n in nodes:
                if t.out_degree[n] >0 and n !=self.root:
                  child_nodes=  t.successors(n)
                  parent = list(t.predecessors(n))[0]
                  t.remove_node(n)
                  t.add_edge(parent, next_label)
                  t.add_edge(next_label, n)
                  for c in child_nodes:
                    t.add_edge(next_label, c)
                  print(list(t.edges))
                  next_label += 1
            print(list(t.edges))
        return trees
                    
    
    def random_tree(self):
        tree = nx.Graph()
        ids = self.ids.copy()
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

    

    def run(self):
        
        dmat = nj.create_distance_matrix(self.alignment, self.ids)
        
        start_trees = self.init(dmat=dmat)

        best_score = np.Inf
        local_opt = []
      
    
        
        for s in start_trees:

            nx.set_node_attributes(s, self.isotype_labels, self.ISO)
            nx.set_node_attributes(s, self.alignment, self.SEQ)
            local_opt.append(self.hill_climbing(s))
        
        best_res = []
        for l in local_opt:
            if l[0] == best_score:
                best_res.append(l)

            elif l[0] < best_score: 
                best_res = [l]
                best_score = l[0]
         
        
        # for l in best_res:
        #     self.print_tree(l)
      
        
        print(best_score)
        return best_res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", required=True, type=str,
        help="filename of input alignment and isotype labels")
    parser.add_argument("-r", "--root", required=True, type=int,
        help="the id of the root sequence in the alignment")
    parser.add_argument("--init", choices=["random", "nj", "edmonds"], default="random")
    parser.add_argument("--n_init", type=int, default=1)
    parser.add_argument("-s", "--seed", type=int, default=1026)

    inpth = "/scratch/projects/tribal/simulator/simulations/B_130_9_9_3_1_43/n7/s1_hot0.25_j0.15_mu0.9_b3"
    # args= parser.parse_args()
    args= parser.parse_args([
        "-a", f"{inpth}/alignment.csv",
        "-r", "0",
        "--init", "edmonds"
    ])

    seqs = {}
    isotypes = {}

    with open(args.alignment, "r") as file:
        index = 0
        for line in file:
            if index ==0:
                index +=1
                continue
            node = [x.strip() for x in line.split(',')]
            n = int(node[0])
            seq = node[1]
            iso = node[2]
            seqs[n] = list(seq)
            isotypes[n] = int(iso)

    
    
    tr = Tribal(
        alignment=seqs,
        root= args.root,
        seed = args.seed,
        isotype_labels= isotypes,
        n_init= args.n_init,
        init = args.init

    )

    tr.run()




  



# s1 = ["a", "b", "c"]
# s2 = ["c", "b", "c"]
# s3 = ["a", "c", "c"]
# s4 = ["a", "a", "a"]
# alignment = {1: s1, 2: s2, 3 : s3, 4: s4}

# seq0 = [ "A", "C", "G", "G"]
# seq1 = [ "C", "T", "A", "G"]
# seq2 = [ "C", "T", "G", "G"]
# seq3 = [ "T", "T", "G", "A"]
# seq4 = [ "A", "C", "A", "C"]
# seq5 = [ "A", "C", "T", "C"]
# isotype_labels = {1: 1, 2: 3, 3:2, 4:0, 5:1}
# alignment = {0: seq0, 1: seq1, 2: seq2, 3: seq3, 4: seq4, 5: seq5}
# ids = list(alignment.keys())
# tr  = TreeBuild(alignment, isotype_labels)
# opt_trees = tr.run()
# print(len(opt_trees))
# print(f" OPT Seq: {opt_seq_score} OPT Iso: {opt_iso_score}")
# tr.print_tree()







