
import networkx as nx
import numpy as np
import neighbor_joining as nj
from edmonds_tree import EdmondsTree
import argparse 
import pandas as pd 
import pickle
from spr import SPR
from unrooted_spr import USPR
from copy import deepcopy
from trans_matrix import TransMat
import utils as ut
from tree_utils import TribalTree


class Tribal:
    def __init__(self, alignment, root, isotype_labels, seed,  n_init=1, 
                init="random", transmat=None, alpha=0.9, jump_prob=0.5, n_isotypes=7,
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

    @staticmethod
    def prob_state_change(e_new, e, T):
        return np.exp(-1*(e_new - e)/T)
    

    def parsimony(self, ttree):
        curr_pars, curr_iso, labels, isotypes = ttree.parsimony(
                                self.alignment,
                                self.isotype_labels,
                                self.Q_isotype,
                                self.alphabet, 
                                self.cost_function )
        curr_obj= self.alpha*(curr_pars) + (1-self.alpha)*(-1*curr_iso)
        return curr_obj, curr_pars, curr_iso, labels, isotypes

    def simulated_annealing(self, curr_state, restart, temp=50, k_max=50):
    
      

        curr_obj, curr_pars, curr_iso, _, _ = self.parsimony(curr_state)
      
        best_obj = curr_obj
        best_state = deepcopy(curr_state)

        #initialize neighbors
        neighbors = iter(USPR(curr_state.get_unrooted_tree(), self.rng))

        for k in range(k_max):
            if (curr_obj - best_obj)/best_obj >= 0.1:
                curr_obj = best_obj 
                curr_state = best_state
            cand_state = deepcopy(curr_state)
              
            #get the next candidate spr tree 
            try:
                t = next(neighbors)
            except:
                # neighbors =iter(USPR(cand_state.get_unrooted_tree(), self.rng,min_radius=3, max_radius=3))
                # try:
                #     t = next(neighbors)
                break

            cand_state.set_unrooted_tree(t)
            curr_temp = temp*np.power(0.99,k)

            #compute the objective value
            cand_obj, cand_pars, cand_iso, _, _ = self.parsimony(cand_state)
            if k % 25 ==0:
                print(f"{restart},{k},{cand_obj},{curr_obj},{best_obj}")     

            #check if solution is better than current state or randomly jump to new state with probability P(cand_obj, curr_obj,curr_temp)
            if cand_obj < curr_obj or  self.rng.random() <= self.prob_state_change(cand_obj, curr_obj, curr_temp ):
       
                curr_state = deepcopy(cand_state)
                neighbors = iter(USPR(curr_state.get_unrooted_tree(), self.rng))
                curr_obj = cand_obj
                if curr_obj < best_obj:
                    best_obj = curr_obj
                    best_state = deepcopy(curr_state)
                    
           
              
        return best_obj, best_state


        #compute the initial solution



    def greedy_hill_climbing(self, TTree):

            iterations = 0
            best_trees = []
     
            best_score, best_iso_score = TTree.parsimony(self.Q_isotype,self.alphabet, self.cost_function )
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
    
               
                # spr_trees = SPR(best_tree, self.Q_isotype)
                spr_trees = USPR(cand_tribal.get_unrooted_tree())
                for t in spr_trees:
                    cand_tribal.set_unrooted_tree(t)
                    
                    count += 1
                  
                    improvement = False
                    # print("best tree")
                    # print(best_tree)
    
                    # print(cand_tribal)
                    cand_score, iso_score =cand_tribal.parsimony(self.Q_isotype,
                                                                self.alphabet, 
                                                                self.cost_function)
                 
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
                
                # if iterations > 2:
                #     break
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

    # def random_tree_old(self):
    #     tree = nx.Graph()
    #     ids = self.ids.copy()
   
    #     n = len(ids)
    #     center = str(2*n -3 )
    #     for i in ids:
    #         tree.add_edge(i, center)
      
    #     next_node = n
       
    #     for i in range(n-3):
    #         pair = self.rng.choice(ids,2, replace=False)
    

    #         for p in pair:
    #             tree.add_edge(p, str(next_node))
    #             tree.remove_edge(p, center)
    #             tree.add_edge(center, str(next_node))
    #             ids.remove(p)
    #         ids.append(str(next_node))
    #         next_node += 1
          
    #     return tree

    # def random_tree(self):
    #     tree = nx.Graph()
  
    #     ids = []
    #     for s in range(self.n_isotypes):
    #         for i in self.ids:
    #             if self.isotype_labels[i] ==s:
    #                 ids.append(i)
        
    #     n = len(ids)
    #     center = 2*n -3 
    #     for i in ids:
    #         tree.add_edge(i, center)
      
    #     next_node = n
       
    #     for i in range(n-3):
    #         pair = [ids[0], ids[1]]
    

    #         for p in pair:
    #             tree.add_edge(p, next_node)
    #             tree.remove_edge(p, center)
    #             tree.add_edge(center, next_node)
    #             ids.remove(p)
    #         ids.append(next_node)
    #         next_node += 1
        
    
    #     return tree

    def initialize_random(self, n_init):

        ttrees = []
        for i in range(n_init):
            ttrees.append(TribalTree(root=self.root, ids = self.ids.copy(), rng=self.rng))
        
        return ttrees
    

    def initialize_nj(self):


        tree = nj.neighbor_joining(self.dmat, self.ids)
        ttree =TribalTree(tree, self.root, is_rooted=False)
               
        #TODO: initialize with random perturbation of NJ tree
        return [ttree]

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
                    
    

    

    def run(self, n_init=1, alpha=0.9, k_max=1000, temp=50):
        best_tree = None
        best_score =np.Inf 
        self.alpha  = alpha

        print("start,k,cand_obj,curr_obj,best_obj")
        start_trees = self.init(n_init)
        for i,start_tree in enumerate(start_trees):

    
            # opt.append(self.greedy_hill_climbing(start_tree))
            score, tree = self.simulated_annealing(start_tree, i, temp=temp, k_max=k_max)

            if score < best_score:
                best_tree = tree
                best_score = score
        obj, pars_obj, iso_obj, labels, isotypes= self.parsimony(best_tree)
        labels = {key : "".join(value) for key, value in labels.items()}
        isotypes = {key: value[0] for key,value in isotypes.items()}
        return best_tree, obj, pars_obj, iso_obj, labels, isotypes
        



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
    parser.add_argument("-k", "--kmax", type=int, default=1000)
    parser.add_argument("--temp", type=float, default=50)
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument("-o", "--output", type=str, help="outputfile of all best trees")
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("-j", "--jump-prob", type=float, default=0.5)
    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("-n", "--newick", type=str, help="filename where newick string should be saved")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")

    args= parser.parse_args()
    # path = "/scratch/projects/tribal/benchmark_pipeline/sim_data/shm_sim/2.0/0.365/1"
    # args =parser.parse_args([
    #     "-a", f"{path}/GCsim_dedup.fasta",
    #     "-r", "naive",
    #     "-t", f"{path}/GCsim_transition_matrix_geometric.txt",
    #     "-i", f"{path}/GCsim_isotypes_geometric.csv",
    #     "-o", "/scratch/projects/tribal/src/test/tribal.tree",
    #     "--init", "random",
    #     "--n_init", "1",
    #     "-k", "26",
    #     "--temp", "55",
    #     "-s", "1",
    #     "--sequences", "/scratch/projects/tribal/src/test/tribal.seq",
    #     "--fasta","/scratch/projects/tribal/src/test/tribal.fasta",
    #     "-n", "/scratch/projects/tribal/src/test/ednmonds.newick",
    #     "--score","/scratch/projects/tribal/src/test/tribal.score",
    #     "--iso_infer", "/scratch/projects/tribal/src/test/tribal.isotypes"
    # ])

    
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
        init = args.init,
        transmat = transMat,
        jump_prob = args.jump_prob

    )

    best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.run(args.n_init, args.alpha, 
                                                                args.kmax, args.temp)



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
            file.write(f"{args.alpha},{obj},{par_obj},{iso_obj}\n")
    if args.iso_infer is not None:
        ut.save_dict(isotypes, args.iso_infer)


  









