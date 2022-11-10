
import networkx as nx
import numpy as np
import neighbor_joining as nj
from edmonds_tree import EdmondsTree
import argparse 
import pandas as pd 
import pickle
from spr import SPR
from copy import deepcopy
from trans_matrix import TransMat
import utils as ut
from tree_utils import TribalTree


class Tribal:
    def __init__(self, alignment, root, isotype_labels, seed,  n_init=1, 
                init="random", transmat=None, alpha=0.5, jump_prob=0.5, n_isotypes=7,
                cost_function=None, alphabet= ("A", "C", "G", "T", "-")):
        
        self.alignment = alignment
        self.n = len(self.alignment)
        self.root = root
        self.isotype_labels = isotype_labels
        self.seed = seed
        self.rng = np.random.default_rng(seed)

        self.n_init = n_init 
        
        if init=="nj":
            self.init = self.initialize_nj
        elif init == "edmonds":
            self.init = self.initialize_edmonds
        else:
            self.init = self.initialize_random


        if transmat is  None:
            self.n_isotypes =n_isotypes
            self.Q_isotype = TransMat(self.rng, n_isotypes).fit_dirichlet()
        else:
            self.Q_isotype = transmat
            self.n_isotypes = self.Q_isotype.shape[0]
  
        self.ids = list(alignment.keys())
        self.ids.sort()
        self.alpha=alpha
    
    
        self.alphabet = alphabet
        self.cost_function = cost_function

        self.dmat = ut.create_distance_matrix(self.alignment, self.ids)

   

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

        

    def hill_climbing(self, TTree):

            iterations = 0
            best_trees = []
     
            best_score, best_iso_score = TTree.parsimony(self.Q_isotype,self.alphabet, self.cost_function, )
            best_tree = deepcopy(TTree)
            best_trees.append(best_tree)
            # best_score = self.alpha*(best_score) + (1-self.alpha)*(best_iso_score)
            # best_score = best_score
            # best_iso_score = best_iso_score
            # while True:
            #     if best_score <= np.NINF:
            #         spr_trees = SPR(best_tree)
            #         for t in spr_trees:
            #             best_tree.set_tree(t)

            #             cand_score, iso_score  =best_tree.parsimony(self.alphabet, self.cost_function, self.Q_isotype)
            #             best_score = self.alpha*(cand_score) + (1-self.alpha)*(iso_score)
            #             break
                 

            #     else:
            #         break
            count = 0
            
         
            
            while True:
                iterations += 1
                cand_tribal = deepcopy(best_tree)
               
                spr_trees = SPR(best_tree, self.Q_isotype)
                for t in spr_trees:
                    count += 1
                  
                    improvement = False
                    # print("best tree")
                    # print(best_tree)
                    cand_tribal.set_tree(t)
                    # print(cand_tribal)
                    cand_score, iso_score =cand_tribal.parsimony(self.Q_isotype,self.alphabet, self.cost_function)
                 
                    lin_comb_score = self.alpha*(cand_score) + (1-self.alpha)*(-1*iso_score)
                    
                    if cand_score < best_score or (cand_score==best_score and iso_score > best_iso_score): #or (cand_score == best_score and cand_iso_score < best_iso_score):

                        # print(list(t.edges))
                        print(f"Candidate score: {cand_score} Best Score: {best_score} Cand Iso: {iso_score} Best Iso: {best_iso_score}")
                        improvement = True
                        best_score = cand_score
                        best_tree = cand_tribal
                        best_iso_score = iso_score
                        best_trees = [deepcopy(cand_tribal)]
                        break
                    elif cand_score == best_score and (best_iso_score == iso_score):

                        # print(cand_tribal)
                        best_trees.append(deepcopy(cand_tribal))
                    

                
                if not improvement:
                    print("Total Trees Examined: " + str(count))
                    break
                
                if iterations > 2:
                    break
                    # if self.rng.random() < 0.99/iterations:
                    #     best_score = cand_score
                    #     best_tree = cand_tribal
                    #     best_trees = [deepcopy(cand_tribal)]
                    # else:
                    #     break
 

            return best_trees

    # @staticmethod
    # def root_tree(tree, root):
    #     tree = nx.dfs_tree(tree,root)
     
    #     root_children = list(tree.successors(root))
    #     for c in root_children:
    #         grandchildren = tree.successors(c)
    #         for g in grandchildren:
    #             tree.add_edge(root, g)

    #     tree.remove_node(c)

    #     return tree

    def random_tree_old(self):
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
        # nx.set_node_attributes(tree, self.isotype_labels, self.ISO)
        # nx.set_node_attributes(tree, self.alignment, self.SEQ)
     
    
  
        return tree

    def random_tree(self):
        tree = nx.Graph()
  
        ids = []
        for s in range(self.n_isotypes):
            for i in self.ids:
                if self.isotype_labels[i] ==s:
                    ids.append(i)
        
        n = len(ids)
        center = 2*n -3 
        for i in ids:
            tree.add_edge(i, center)
      
        next_node = n
       
        for i in range(n-3):
            pair = [ids[0], ids[1]]
    

            for p in pair:
                tree.add_edge(p, next_node)
                tree.remove_edge(p, center)
                tree.add_edge(center, next_node)
                ids.remove(p)
            ids.append(next_node)
            next_node += 1
        
        # nx.set_node_attributes(tree, self.isotype_labels, self.ISO)
        # nx.set_node_attributes(tree, self.alignment, self.SEQ)
  
        return tree

    def initialize_random(self):

        iso_labs = {k: [v] for k,v in self.isotype_labels.items()}

        tree = self.random_tree()
        ttree =TribalTree(tree, root=self.root, is_rooted=False, 
                        sequences=self.alignment.copy(), isotypes=iso_labs)
        
        return ttree
    

    def initialize_nj(self):


        tree = nj.neighbor_joining(self.dmat, self.ids)
        iso_labs = {k: [v] for k,v in self.isotype_labels.items()}
        ttree =TribalTree(tree, self.root, is_rooted=False , 
                    sequences =self.alignment.copy(), 
                    isotypes=iso_labs)

        return ttree

    # def initialize_edmonds(self, **kwargs):
    #     dmat= kwargs['dmat']
    #     iso_labs = {k: [v] for k,v in self.isotype_labels.items()}
    #     tribal_tree_list = []
    #     trees = EdmondsInit(dmat, self.isotype_labels, seed=1026, transMat=self.Q_isotype, root=self.root, ).generate(1, prop=0.2)
    #     for t in trees:
    #         ttree =TribalTree(t, self.alignment.copy(), 
    #                    iso_labs, 
    #                    root=self.root)
    #         tribal_tree_list.append(ttree)
    #     return tribal_tree_list
                    
    

    

    def run(self):
        best_score = np.Inf
        opt = []
        for i in range(self.n_init):

            #get an intiial TribalTree to start each run
            start_tree = self.init()

            opt.append(self.hill_climbing(start_tree))
        
        # best_trees = None
        # best_score = np.Inf
        # for trees in local_opt:
        #     t = trees[0]
        #     score = t.tree_score(self.alpha)
        #     if score < best_score:
        #         best_trees = trees
        best_tree = None
        best_score =np.Inf 
        for trees in opt:
            t = trees[0]
            score = t.tree_score(self.alpha)
            if score < best_score:
                best_tree = t
      
        return best_tree
        
        # best_res = []
        # for l in local_opt:
        #     if l[0] == best_score:
        #         best_res.append(l)

        #     elif l[0] < best_score: 
        #         best_res = [l]
        #         best_score = l[0]
         
        
        # for l in best_res:
        #     self.print_tree(l)
      
        




def pickle_save(obj, fname):
        with open(fname, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", required=True, type=str,
        help="filename of input fasta file containing the alignment")
    parser.add_argument("-i", "--isotypes", required=True, type=str,
        help="filename of input file containing the isotype labels")
    parser.add_argument("-t", "--transmat", required=False, type=str,
        help="filename of input transition matrix")
    parser.add_argument("-r", "--root", required=True,
        help="the id of the root sequence in the alignment")
    parser.add_argument("--init", choices=["random", "nj", "edmonds"], default="random")
    parser.add_argument("--n_init", type=int, default=1)
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument("-o", "--output", type=str, help="outputfile of all best trees")
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("-j", "--jump-prob", type=float, default=0.5)
    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("-n", "--newick", type=str, help="filename where newick string should be saved")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")


    # args= parser.parse_args()
    path = "/scratch/projects/tribal/benchmark_pipeline/sim_data/shm_sim/2.0/0.365/1"
    args =parser.parse_args([
        "-a", f"{path}/GCsim_dedup.fasta",
        "-r", "naive",
        "-t", f"{path}/GCsim_transition_matrix_geometric.txt",
        "-i", f"{path}/GCsim_isotypes_geometric.csv",
        "-o", "/scratch/projects/tribal/src/test/tribal.tree",
        "--init", "nj",
        "-s", "1",
        "--sequences", "/scratch/projects/tribal/src/test/tribal.seq",
        "--fasta","/scratch/projects/tribal/src/test/tribal.fasta",
        "-n", "/scratch/projects/tribal/src/test/ednmonds.newick",
        "--score","/scratch/projects/tribal/src/test/tribal.score",
        "--iso_infer", "/scratch/projects/tribal/src/test/tribal.isotypes"
    ])

    
    # args= parser.parse_args([
    #     "-a", f"{inpth}/alignment.csv",
    #     # "-a", "/scratch/projects/tribal/mouse_data/input/B_25_5_3_3_1_43.csv",
    #     "-r", "0",
    #     "--init", "edmonds",
    #     "-o", f"{outpath}/tree.pickle",
    #     "--alpha", "0",
    #     "-j", "0.95"

    # ])

    alignment = ut.read_fasta(args.alignment)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    # ids = list(alignment.keys())

    isotypes = ut.read_dict(args.isotypes)
    isotypes[args.root] =0
    isotypes_filt = {i: int(isotypes[i]) for i in alignment}

    transMat = np.loadtxt(args.transmat)
    

    tr = Tribal(
        alignment=alignment,
        root= args.root,
        isotype_labels= isotypes_filt,
        seed = args.seed,
        n_init= args.n_init,
        init = args.init,
        transmat = transMat,
        alpha= args.alpha,
        jump_prob = args.jump_prob

    )

    best_tree = tr.run()
    labels, isotypes, seq_score, iso_score = best_tree.get_output()



    if args.fasta is not None:
        ut.write_fasta(args.fasta, labels)
    if args.sequences is not None:
        ut.save_dict(labels, args.sequences)

    if args.output is not None:
        parents = best_tree.get_parents()
        ut.save_dict(parents, args.output)
    if args.newick is not None:
        newick = best_tree.tree_to_newick(rooted=False)
        with open(args.newick, "w+") as file:
            file.write(newick)
    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write(f"{seq_score},{iso_score}\n")
    if args.iso_infer is not None:
        ut.save_dict(isotypes, args.iso_infer)


  



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







