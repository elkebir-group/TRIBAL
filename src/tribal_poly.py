
import networkx as nx
import numpy as np
import neighbor_joining as nj
from edmonds_tree import EdmondsTree
import argparse 
import pickle
from spr import SPR
from unrooted_spr import USPR
from copy import deepcopy
from trans_matrix import TransMat
import utils as ut
from tree_utils import TribalTree
from ete3 import Tree, TreeNode
from draw_tree import DrawTree
from multi_spr import MultSPR
from em_weight_matrix import EMProbs
from lineage_tree import LineageTree, LineageForest

class TribalPoly:
    def __init__(self, transmat=None, alpha=0.9, n_isotypes=7,
                cost_function=None, 
                alphabet= ("A", "C", "G", "T","N", "-"), 
                isotype_encoding=None):
        
        # self.alignment = alignment
        # self.isotype_labels =isotype_labels
        # self.isotype_encoding = isotype_encoding

        if transmat is None:
            self.Q_isotype = np.ones((n_isotypes, n_isotypes))
            for i in range(n_isotypes):
                self.Q_isotype[i,i] = 0
            for s in range(n_isotypes):
                for t in range(n_isotypes):
                    if s > t:
                        self.Q_isotype[s,t] = 1e6
        else:
            self.Q_isotype = transmat
            self.n_isotypes = self.Q_isotype.shape[0]
  
        self.alpha=alpha
        self.alphabet = alphabet
        self.cost_function = cost_function
        self.states = [i for i in range(self.Q_isotype.shape[0])]

   


    

    # def parsimony(self, ttree):
    #     curr_pars, curr_iso, labels, isotypes = ttree.parsimony(
    #                             self.alignment,
    #                             self.isotype_labels,
    #                             self.Q_isotype,
    #                             self.alphabet, 
    #                             self.cost_function )
    #     curr_obj= self.alpha*(curr_pars) + (1-self.alpha)*(curr_iso)
    #     return curr_obj, curr_pars, curr_iso, labels, isotypes

    def simulated_annealing(self, curr_state, restart, temp=50, k_max=1000):
    
      

        curr_obj, curr_pars, curr_iso, _, _ = self.parsimony(curr_state)
      
        best_obj = curr_obj
        best_state = deepcopy(curr_state)

        #initialize neighbors
        neighbors = iter(USPR(curr_state.get_unrooted_tree(), self.rng))

        for k in range(k_max):
            if (curr_obj - best_obj)/best_obj >= 0.50:
                curr_obj = best_obj 
                curr_state = best_state
            cand_state = deepcopy(curr_state)
              
            #get the next candidate spr tree 
            try:
                t = next(neighbors)
            except StopIteration:
                # print("no other trees in neighborhood")
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


    def save_tree(self, ttree, isotypes, fname):
        pth = "/scratch/projects/tribal/src/test"
        fname = f"{pth}/{fname}"
        parents = ttree.get_parents()
        dt = DrawTree(parents, isotypes, show_legend=True, isotype_encoding=self.isotype_encoding)
        dt.save(fname)




    def greedy_hill_climbing(self, lin_tree, alignment, isotype_labels, search="multSPR"):

            iterations = 0
            best_trees = []
            best_tree = deepcopy(lin_tree)
            best_trees.append(best_tree)

     
            best_pars_score, best_iso_score, best_labels,best_isotypes = best_tree.parsimony(alignment,
                                                                                    isotype_labels,
                                                                                    self.Q_isotype,
                                                                                    self.alphabet,
                                                                                    self.cost_function)
      

            best_score = self.alpha*best_pars_score - (1-self.alpha)*best_iso_score
            # self.save_tree(best_tree, isotypes, "start_tree.png")
            # self.save_tree(best_tree, isotypes, "best_tree.png")

 
            count = 0
            
         
            
            while True:
                iterations += 1
                cand_tribal = deepcopy(best_tree)
    
                if search == "multSPR":
                    spr_trees = MultSPR(cand_tribal, isotype_labels, self.states)
                # spr_trees = SPR(best_tree, self.Q_isotype)
                # else:
                #     spr_trees = USPR(cand_tribal.get_unrooted_tree(), self.rng)
                for t in spr_trees:
              
      
                    cand_tribal.set_tree(t)
                    
                    count += 1
                  
                    improvement = False
            
                    pars_score, iso_score, labels,isotypes = cand_tribal.parsimony(alignment,
                                                                                    isotype_labels,
                                                                                    self.Q_isotype,
                                                                                    self.alphabet,
                                                                                    self.cost_function)

                 
                    lin_comb_score = self.alpha*pars_score - (1-self.alpha)*iso_score                    
                    if lin_comb_score < best_score:
                   
                      
                        print(f"Candidate score: {lin_comb_score} Best Score: {best_score}\nCand Iso: {iso_score} Best Iso: {best_iso_score}\nCand Pars {pars_score} Best Pars: {best_pars_score}")
                        improvement = True
                        best_score = lin_comb_score
                        best_pars_score =  pars_score
                        best_tree = deepcopy(cand_tribal)
                        best_labels = labels
                        best_isotypes = isotypes 
                        # self.save_tree(best_tree, isotypes,"best_tree.png")
                        best_iso_score = iso_score
                        # best_trees = [deepcopy(cand_tribal)]
                        break
             

                
                if not improvement:
                    print("Total Trees Examined: " + str(count))
                    break
                
 

            return best_score, best_tree, best_labels, best_isotypes

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

    def initialize_candidates(self, n_init=None, random_resolution=False):
        ttrees = []
        for tree in self.candidates:
            ttree = TribalTree(tree, self.root, is_rooted=True)
            if random_resolution:
                tree.resolve_polytomies_randomly()
            ttrees.append(ttree)
        return ttrees

    def initialize_random(self, n_init):

        ttrees = []
        for i in range(n_init):
            ttrees.append(TribalTree(root=self.root, ids = self.ids.copy(), rng=self.rng))
        
        return ttrees
    

    def initialize_nj(self, n_init):

  
        tree = nj.neighbor_joining(self.dmat, self.ids)
        start_trees =  [TribalTree(tree, self.root, is_rooted=False)]

        spr = USPR(tree, self.rng)
        for i, t in enumerate(spr):
            if i >= n_init:
                break
            ttree = TribalTree(t, self.root, is_rooted=False)

            start_trees.append(ttree)    
        return start_trees

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
                    
    
    def update_labels(self, labels):
        # seq = self.full_alignment[self.root]
        
        # for key in labels:
        #     seq = self.full_alignment[self.root].copy()
        #     new_seq= labels[key]
        #     for j, index in enumerate(self.diff):
        #         seq[index] = new_seq[j]
        #     labels[key] = seq

            
    
        labels = {key : "".join(value) for key, value in labels.items()}
        return labels 

    def run(self, n_init=1, alpha=0.9, k_max=1000, temp=50):
        best_tree = None
        best_score =np.Inf 
        self.alpha  = alpha

        # print("start,k,cand_obj,curr_obj,best_obj")
        start_trees = self.init(n_init)
        for i,start_tree in enumerate(start_trees):

    
            # opt.append(self.greedy_hill_climbing(start_tree))
            # score, tree = self.simulated_annealing(start_tree, i, temp=temp, k_max=k_max)
            score, tree =self.greedy_hill_climbing(start_tree, search="multSPR")

            if score < best_score:
                print("best tree:"+ str(i))
                best_tree = tree
                best_score = score
             
        obj, pars_obj, iso_obj, labels, isotypes= self.parsimony(best_tree)
        labels = self.update_labels(labels)
       
        # isotypes = {key: value[0] for key,value in isotypes.items()}
        return best_tree, obj, pars_obj, iso_obj, labels, isotypes
        
    def fit(self, lin_forest, convert=False):
        best_tree = None
        best_score =np.Inf 
        score_dict = {}
        
        for tree in lin_forest.get_trees():
            score_dict[tree.id] = {}
            curr_pars, curr_iso, labels, isotypes =tree.parsimony(
                                lin_forest.alignment,
                                lin_forest.isotypes,
                                self.Q_isotype,
                                self.alphabet, 
                                self.cost_function,
                                convert=False )
            curr_obj= self.alpha*(curr_pars) + (1-self.alpha)*(curr_iso)
            score_dict[tree.id]["obj"] = curr_obj 
            score_dict[tree.id]["par_obj"] = curr_pars
            score_dict[tree.id]["iso_obj"] = curr_iso
            score_dict[tree.id]["labels"] = labels
            score_dict[tree.id]["isotypes"] = isotypes 
       
            
            if curr_obj < best_score:
                    best_tree = tree
                    best_score = curr_obj
            
        return  best_score, best_tree, score_dict
    
    def search(self, lin_tree, alignment, isotypes):
        self.greedy_hill_climbing(lin_tree, alignment, isotypes)

def convert_to_nx(ete_tree, root):
    nx_tree = nx.DiGraph()
    internal_node = 1
    internal_node_count = 0
    for node in ete_tree.traverse("preorder"):

        if node.name == "":
            node.name = internal_node
            internal_node_count += 1
            internal_node += 1
        if node.is_root():
            root_name =node.name

        # print(node.name)
        for c in node.children:
            if c.name == "":
                c.name = str(internal_node)
                internal_node += 1
                internal_node_count += 1
            # else:
            #     print(c.name)
            # if isinstance(c.name, str):
            #     print(c.name)
       

            nx_tree.add_edge(node.name, c.name)
    
    if len(list(nx_tree.neighbors(root))) == 0:
        # print("root is outgroup")
        nx_tree.remove_edge(root_name, root)
        nx_tree.add_edge(root, root_name)
        # children = list(nx_tree.neighbors(root_name))
        # for child in children:
        #     nx_tree.remove_edge(root_name, child)
        #     nx_tree.add_edge(root, child)
        # nx_tree.remove_node(root_name)
        # print(list(nx_tree.edges))
        # print(len(list(nx_tree.nodes)))


    



    return nx_tree
        


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
    parser.add_argument("--init", choices=["random", "nj", "edmonds", "candidates"], default="random")
    parser.add_argument("--n_init", type=int, default=1)
    parser.add_argument("-k", "--kmax", type=int, default=1000)
    parser.add_argument("--temp", type=float, default=50)
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument("-e", "--encoding", type=str)
    parser.add_argument("-o", "--output", type=str, help="outputfile of all best trees")
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("-j", "--jump-prob", type=float, default=0.5)
    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--png", type=str, help="filename where to save a png of the optimal tree")

    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("-n", "--newick", type=str, help="filename where newick string should be saved")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")
    parser.add_argument("--candidates", type=str, help="filename containing newick strings for candidate trees")
    parser.add_argument("--save_candidates", type=str, help="directory where to save data for candidate trees")
    parser.add_argument("--fit", action="store_true", 
                help="if given an input set of candidate trees, fit will only compute the objective value of each tree and not search tree space")
    # args= parser.parse_args()


    # path = "/scratch/projects/tribal/benchmark_pipeline/sim_data/shm_sim2/2.0/0.365/2"
    # args =parser.parse_args([
    #     "-a", f"{path}/GCsim_dedup.fasta",
    #     "-r", "naive",
    #     "-t", f"/scratch/projects/tribal/benchmark_pipeline/sim_data//transmat_horns_etal.txt",
    #     "-i", f"{path}/GCsim_dedup.isotype.fasta",
    #     "-o", "/scratch/projects/tribal/src/test/tribal.tree",
    #     "--init", "candidates",
    #     # "--n_init", "1",
    #     # "-k", "500",
    #     # "--temp", "55",
    #     "-s", "1",
    #     "--alpha", "0.9",
    #     "--sequences", "/scratch/projects/tribal/src/test/tribal.seq",
    #     "--fasta","/scratch/projects/tribal/src/test/tribal.fasta",
    #     "-n", "/scratch/projects/tribal/src/test/ednmonds.newick",
    #     "--score","/scratch/projects/tribal/src/test/tribal.score",
    #     "--iso_infer", "/scratch/projects/tribal/src/test/tribal.isotypes",
    #     "--candidate", f"{path}/dnapars/outtree",
    #     "--save_candidates", "/scratch/projects/tribal/src/test",
    #     "--fit"
    # ])

    path = "/scratch/projects/tribal/real_data"
    dataset = "day_14"
    clonotype = "B_120_2_8_210_1_13"

    args =parser.parse_args([
        "-a", f"{path}/{dataset}/input/{clonotype}/concat.aln.fasta",
        "-r", "naive",
        "-t", f"{path}/mouse_transmat2.txt",
        "-e", f"{path}/mouse_isotype_encoding.txt",
        "-i", f"{path}/{dataset}/input/{clonotype}/isotype.fasta",
        "-o", f"{path}/test/tree.txt",
        "--init", "candidates",
        "-s", "1",
        "--alpha", "0.75",
        "--sequences", f"{path}/test/triba.seq",
        "--fasta",f"{path}/test/triba.fasta",
        "--score",f"{path}/test/tribal.score",
        "--iso_infer", f"{path}/test/tribal.isotypes",
        "--candidate", f"{path}/{dataset}/dnapars/{clonotype}/outtree",
        "--save_candidates",f"{path}/test",
        "--png", f"{path}/test/best_tree.png",
        # "--fit"
    ])

    if args.candidates is not None:
        args.init = "candidates"


    alignment = ut.read_fasta(args.alignment)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    ids = list(alignment.keys())

    isotypes = ut.read_fasta(args.isotypes)
    if args.encoding is not None:
        iso_encoding = {}
        counter = 0
    
        with open(args.encoding, "r") as file:
            for line in file:
                isotype = line.strip()
                if counter ==0:
                    start_iso = isotype 
                iso_encoding[isotype] = counter
                counter += 1
        
        isotypes_filt = {}
        for i in alignment:
            iso = isotypes[i]
            if iso not in iso_encoding:
                iso = isotypes[i].lower()
                if 'm' in iso or 'd' in iso:
                    iso = start_iso
            isotypes_filt[i] = iso_encoding[iso]
            
        # isotypes_filt = {i: iso_encoding[isotypes[i]] for i in alignment}
        isotype_encoding = {val: key for key, val in iso_encoding.items()}
    else:
        isotypes[args.root] =0
        isotypes_filt = {i: isotypes[i] for i in alignment}
        isotype_encoding = {}

    transMat = np.loadtxt(args.transmat)
  
    if args.candidates is None:
        tr = Tribal(
            alignment=alignment,
            root= args.root,
            isotype_labels= isotypes_filt,
            seed = args.seed,
            init = args.init,
            transmat = transMat,
            jump_prob = args.jump_prob,
            isotype_encoding= isotype_encoding

        )

        best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.run(args.n_init, args.alpha, 
                                                                    args.kmax, args.temp)
    else:
    
        cand_trees = []
        import re
        exp = '\[.*\]'
       
        with open(args.candidates, 'r') as file:
            nw_strings = []
            nw_string = ""
            for nw in file:
                    line = nw.strip()
                    nw_string += line
                    if ";" in line:
                        
                        nw_strings.append(nw_string)
                        nw_string = ""

            for nw in nw_strings:
    
                nw = re.sub(exp, '', nw)
             

                ete_tree = Tree(nw, format=0)
                # ete_tree.name = args.root
                # print(ete_tree)
                nx_tree= convert_to_nx(ete_tree, args.root)
                # print(list(nx_tree.edges))
                # ttree = TribalTree(nx_tree, root=args.root, is_rooted=True)
                cand_trees.append(nx_tree)
        
        

        tr = Tribal(
            alignment=alignment,
            root= args.root,
            isotype_labels= isotypes_filt,
            seed = args.seed,
            init = args.init,
            transmat = transMat,
            jump_prob = args.jump_prob,
            candidates= cand_trees,
            isotype_encoding=isotype_encoding

        )

        #test em weight matrix
        # states = np.arange(transMat.shape[0])
        # forest = LineageForest()
        # for c in cand_trees:
        #     lt = LineageTree(c, "naive")
        #     forest.add(lt)
   
        # em = EMProbs(forest,  transMat )
        # exp_log_like, state_probs, transMat = em.fit(isotypes_filt)
   

        if args.fit:
            best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.fit(cand_trees, args.alpha, save_candidates=args.save_candidates)
        else:
            best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.run(len(cand_trees), args.alpha, args.kmax, args.temp)






        # for k in alignment:
        #     print(k)
        #     seq = "".join(alignment[k])
        #     if seq != labels[k]:
        #         print(best_tree.is_leaf(k))

 



    if args.fasta is not None:
        ut.write_fasta(args.fasta, labels)
    if args.sequences is not None:
        ut.save_dict(labels, args.sequences)

    if args.output is not None:
        parents = best_tree.get_parents()
        ut.save_dict(parents, args.output)
        if args.png is not None:
            dt = DrawTree(parents, isotypes, show_legend=True, isotype_encoding=isotype_encoding)
            dt.save(args.png)
    if args.newick is not None:
        newick = best_tree.tree_to_newick(rooted=False)
        with open(args.newick, "w+") as file:
            file.write(newick)
    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write(f"{args.alpha},{obj},{par_obj},{iso_obj}\n")
    if args.iso_infer is not None:
        ut.save_dict(isotypes, args.iso_infer)
    


  









