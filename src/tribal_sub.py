
import networkx as nx
import numpy as np
import neighbor_joining as nj
import argparse 
import re
from spr import SPR
from unrooted_spr import USPR
from copy import deepcopy
import utils as ut

from ete3 import Tree

from multi_spr import MultSPR
from score_class import Score
from lineage_tree import LineageForest


class TribalSub:
    def __init__(self, transmat=None, alpha=0.9, n_isotypes=7,
                cost_function=None, 
                alphabet= ("A", "C", "G", "T","N", "-")):
        


        if transmat is None:
            self.Q_isotype = np.ones((n_isotypes, n_isotypes))
            for i in range(n_isotypes):
                self.Q_isotype[i,i] = 0
            for s in range(n_isotypes):
                for t in range(n_isotypes):
                    if s > t:
                        self.Q_isotype[s,t] = 1e8
        else:
            self.Q_isotype = transmat
            self.n_isotypes = self.Q_isotype.shape[0]
  
    
        self.alpha=alpha
        self.alphabet = alphabet
        self.cost_function = cost_function
        self.states = [i for i in range(self.Q_isotype.shape[0])]

   


    
    def score(self, lin_tree, alignment, isotype_labels):
        
        seq_score, seq_labels = lin_tree.sequence_parismony(alignment, 
                                                    alphabet=self.alphabet, 
                                                    cost_function=self.cost_function)
        iso_score, iso_labels = lin_tree.isotype_parsimony(isotype_labels, transmat= self.Q_isotype)
        obj = self.compute_score(seq_score, iso_score)

        return   Score(obj, seq_score, iso_score, seq_labels, iso_labels), lin_tree


    def refine(self, lin_tree, alignment, isotype_labels):
        iso_score, iso_labels = lin_tree.isotype_parismony_polytomy(isotype_labels, transmat= self.transmat )

        seq_score, seq_labels = lin_tree.sequence_parismony(alignment, 
                                                    alphabet=self.alphabet, 
                                                    cost_function=self.cost_function)
        obj = self.compute_score(seq_score, iso_score)
        return Score(obj, seq_score, iso_score, seq_labels, iso_labels), lin_tree
    
    def search(self, lin_tree, alignment, isotype_labels, mode="multSPR"):

            iterations = 0
            best_tree = lin_tree
     
            best_result = self.refine(best_tree, alignment, isotype_labels)
    
            count = 0
            improvement = False
            while True:
                iterations += 1
                cand_tribal = deepcopy(best_tree)
    
                if mode == "multSPR":
                    spr_trees = MultSPR(cand_tribal, best_result.isotype_labels, self.states)
                # spr_trees = SPR(best_tree, self.Q_isotype)
                # else:
                #     spr_trees = USPR(cand_tribal.get_unrooted_tree(), self.rng)
                for t in spr_trees:
                    improvement = False
                    cand_tribal.set_tree(t)            
                    count += 1
            
                    current_result = self.refine(cand_tribal, alignment, isotype_labels)

                 
                    if current_result.improvement(best_result):
                        print("Candidate Result:")
                        print(current_result)
                        print("\nBest Result:")
                        print(best_result)
                        improvement = True
                        best_result = current_result
                        best_tree = deepcopy(cand_tribal)
                
                        break
             
                if not improvement:
                    print("Total Trees Examined: " + str(count))
                    break
                
 
            
            return best_result, best_tree


    
    def compute_score(self, pars, iso):
        return self.alpha*pars + (1-self.alpha)*iso

    
   
    def forest_mode(self, lin_forest, alignment=None, isotypes=None, mode="score"):
            best_result = None 
            best_tree = None

            all_results = {}
            
            if alignment is None:
                alignment = lin_forest.alignment
            if isotypes is None:
                isotype_labels = lin_forest.isotypes 

            if mode =="refine":
                mode_func = self.refine
            elif mode == "search":
                mode_func = self.search 
            else:
                mode_func = self.score 

            for lin_tree in lin_forest.get_trees():
                result, out_tree = mode_func(lin_tree,alignment, isotype_labels)
                all_results[lin_tree.id] = result

                if best_result is None:
                    best_result = result 
                    print(best_result)
                else:
                    if result.improvement(best_result):
                        best_result = result 
                        best_tree = out_tree
                        print(best_result)

      
            return best_result, best_tree, all_results
      
  
def create_trees(cand_fname):
    cand_trees = []
      
    exp = '\[.*\]'
       
    with open(cand_fname, 'r') as file:
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
        return cand_trees

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
  


    return nx_tree
        


def get_alignment(fname):
    alignment = ut.read_fasta(fname)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    return alignment

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", type=str,
        help="filename of input fasta file containing the alignment")
    parser.add_argument("-i", "--isotypes",  type=str,
        help="filename of input file containing the isotype labels")
    parser.add_argument("-t", "--transmat", required=False, type=str,
        help="filename of input transition matrix")
    parser.add_argument("-r", "--root", required=True,
        help="the id of the root sequence in the alignment")
    parser.add_argument("-l", "--lineage", type=str, help="pickle file of lineage tree/forest")
    parser.add_argument("--candidates", type=str, help="filename containing newick strings for candidate tree(s)")
    parser.add_argument("--mode", choices=["score", "refine", "search"], default="score")
    parser.add_argument("-e", "--encoding", type=str, required=True)
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("-j", "--jump-prob", type=float, default=0.25)
    parser.add_argument("--forest",  action="store_true")


    
    
    
    parser.add_argument("-o", "--output", type=str, help="outputfile of all best trees")

    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--png", type=str, help="filename where to save a png of the optimal tree")
    parser.add_argument( "--all_pngs", action="store_true")

    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("-n", "--newick", type=str, help="filename where newick string should be saved")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")
    parser.add_argument("--save_candidates", type=str, help="directory where to save data for candidate trees")
    parser.add_argument("--save_all_scores", type=str, help="file where to save data for candidate trees")

    parser.add_argument("--fit", action="store_true", 
                help="if given an input set of candidate trees, fit will only compute the objective value of each tree and not search tree space")
    args= parser.parse_args()




    # path = "/scratch/projects/tribal/real_data"
    # dataset = "day_14"
    # folder = "tribal"
    # clonotype = "B_120_2_8_210_1_13"

    # args =parser.parse_args([
    #     "-a", f"{path}/{dataset}/input/{clonotype}/concat.aln.fasta",
    #     "-r", "naive",
    #     "-t", f"{path}/{dataset}/{folder}/transmat.txt",
    #     # "-l", f"{path}/{dataset}/{folder}/{clonotype}/fit_tree.pickle",
    #     "-e", f"{path}/mouse_isotype_encoding.txt",
    #     "-i", f"{path}/{dataset}/input/{clonotype}/isotype.fasta",
    #     "-o", f"{path}/test/tree.txt",
    #     # "--init", "candidates",
    #     # "-s", "1",
    #     "--alpha", "0.75",
    #     "--sequences", f"{path}/test/triba.seq",
    #     "--fasta",f"{path}/test/triba.fasta",
    #     "--score",f"{path}/test/tribal.score",
    #     "--iso_infer", f"{path}/test/tribal.isotypes",
    #     "--candidate", f"{path}/{dataset}/dnapars/{clonotype}/outtree",
    #     "--save_candidates",f"{path}/test",
    #     "--save_all_scores",f"{path}/test/all_scores.csv",
    #     "--png", f"{path}/test/best_tree.png",
    #     "--all_pngs"
    #     # "--start_png", f"{path}/test/start_tree.png",
    #     # "--fit"
    # ])

    # if args.candidates is not None:
    #     args.init = "candidates"


    
    iso_encoding = {}
    counter = 0

    with open(args.encoding, "r") as file:
        for line in file:
            isotype = line.strip()
            if counter ==0:
                start_iso = isotype 
            iso_encoding[isotype] = counter
            counter += 1
    
    isotype_encoding = {val: key for key,val in iso_encoding.items()}


                
    if args.alignment is not None:
           alignment = get_alignment(args.alignment)
    else:
            alignment = None
    
    if args.isotypes is not None:

        isotypes = ut.read_fasta(args.isotypes)
        
            
        isotypes_filt = {}
        for i in isotypes:
            iso = isotypes[i]
            if iso not in iso_encoding:
                iso = isotypes[i].lower()
                if 'm' in iso or 'd' in iso:
                    iso = start_iso
            isotypes_filt[i] = iso_encoding[iso]
                
    
    transmat= None 
    if args.transmat is not None:
       transmat = np.loadtxt(args.transmat)



    if args.candidates is not None:
        cand = create_trees(args.candidates)
        lin_forest = LineageForest(alignment=alignment, isotypes=isotypes_filt)
        lin_forest.generate_from_list(cand, args.root)
      

    else:
        if args.lineage is not None:
            lin= ut.pickle_load(args.lineage)
            if args.forest:
                lin_forest = lin
                if lin_forest.alignment is None:
                    lin_forest.alignment = alignment 
                if lin_forest.isotypes is None:
                    lin_forest.isotypes = isotypes_filt 
               
            else:
                lin_forest = LineageForest([lin], alignment, isotypes_filt)

    tr = TribalSub(transmat, args.alpha)
    print(f"Input:\nncells: {lin_forest.ncells()}\nforest size: {lin_forest.size()}\nmode: {args.mode}\n")

    best_result, best_tree, all_results =  tr.forest_mode(lin_forest, mode =args.mode)
      

    print(f"{args.mode} complete! \nBest Tree: {best_tree.id}")
    print(best_result)
    print("saving results......")


    if args.output is not None:
        best_tree.save_tree(args.output)
        if args.png is not None:
            best_tree.save_png(args.png, best_result.isotypes, isotype_encoding)
        
    
    best_labels = ut.update_labels(best_result.labels)

    if args.fasta is not None:
        
        ut.write_fasta(args.fasta, best_labels)
    if args.sequences is not None:
        ut.save_dict(best_labels, args.sequences)
    

    if args.score is not None:
        best_result.save_score(args.score, args.alpha)
       
    if args.iso_infer is not None:
        ut.save_dict(best_result.isotypes, args.iso_infer)
    
    if args.save_all_scores:
        with open(args.save_all_scores, 'w+') as file:
            for key, val in all_results.items():
                file.write(f"{key},{args.alpha},{val.objective},{val.seq_obj},{val.iso_obj}\n")
    
    if args.save_candidates is not None:
        pth =args.save_candidates
        for tree in lin_forest:
            res= all_results[tree.id]
            parents = tree.get_parents()
            labs = ut.update_labels(res.labels)
            ut.write_fasta(f"{pth}/tree{tree.id}.seq.fasta", labs)
            ut.save_dict(labs, f"{pth}/tree{tree.id}.seq.csv")
            ut.save_dict(res.isotypes, f"{pth}/tree{tree.id}.isotypes.csv")
            tree.save_tree(f"{pth}/tree{tree.id}.txt")
            if args.all_pngs:
                tree.save_png(f"{pth}/tree{tree.id}.png", res.isotypes, isotype_encoding)
        







            # dt = DrawTree(parents, isotypes, show_legend=True, isotype_encoding=isotype_encoding)
            # dt.save(args.png)
    # if args.newick is not None:
    #     newick = best_tree.tree_to_newick(rooted=False)
    #     with open(args.newick, "w+") as file:
    #         file.write(newick)

    



    # if args.candidates is None:
    #     tr = Tribal(
    #         alignment=alignment,
    #         root= args.root,
    #         isotype_labels= isotypes_filt,
    #         seed = args.seed,
    #         init = args.init,
    #         transmat = transMat,
    #         jump_prob = args.jump_prob,
    #         isotype_encoding= isotype_encoding

    #     )

    #     best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.run(args.n_init, args.alpha, 
    #                                                                 args.kmax, args.temp)
    # else:
    
    #     cand_trees = []
    #     import re
    #     exp = '\[.*\]'
       
    #     with open(args.candidates, 'r') as file:
    #         nw_strings = []
    #         nw_string = ""
    #         for nw in file:
    #                 line = nw.strip()
    #                 nw_string += line
    #                 if ";" in line:
                        
    #                     nw_strings.append(nw_string)
    #                     nw_string = ""

    #         for nw in nw_strings:
    
    #             nw = re.sub(exp, '', nw)
             

    #             ete_tree = Tree(nw, format=0)
    #             # ete_tree.name = args.root
    #             # print(ete_tree)
    #             nx_tree= convert_to_nx(ete_tree, args.root)
    #             # print(list(nx_tree.edges))
    #             # ttree = TribalTree(nx_tree, root=args.root, is_rooted=True)
    #             cand_trees.append(nx_tree)
        
        

    #     tr = Tribal(
    #         alignment=alignment,
    #         root= args.root,
    #         isotype_labels= isotypes_filt,
    #         seed = args.seed,
    #         init = args.init,
    #         transmat = transMat,
    #         jump_prob = args.jump_prob,
    #         candidates= cand_trees,
    #         isotype_encoding=isotype_encoding

    #     )

    #     #test em weight matrix
    #     # states = np.arange(transMat.shape[0])
    #     # forest = LineageForest()
    #     # for c in cand_trees:
    #     #     lt = LineageTree(c, "naive")
    #     #     forest.add(lt)
   
    #     # em = EMProbs(forest,  transMat )
    #     # exp_log_like, state_probs, transMat = em.fit(isotypes_filt)
   

    #     if args.fit:
    #         best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.fit(cand_trees, args.alpha, save_candidates=args.save_candidates)
    #     else:
    #         best_tree, obj, par_obj, iso_obj, labels, isotypes = tr.run(len(cand_trees), args.alpha, args.kmax, args.temp)






        # for k in alignment:
        #     print(k)
        #     seq = "".join(alignment[k])
        #     if seq != labels[k]:
        #         print(best_tree.is_leaf(k))

 





  










####### SCRATCH##############

# def simulated_annealing(self, curr_state, restart, temp=50, k_max=1000):
    
      

#         curr_obj, curr_pars, curr_iso, _, _ = self.parsimony(curr_state)
      
#         best_obj = curr_obj
#         best_state = deepcopy(curr_state)

#         #initialize neighbors
#         neighbors = iter(USPR(curr_state.get_unrooted_tree(), self.rng))

#         for k in range(k_max):
#             if (curr_obj - best_obj)/best_obj >= 0.50:
#                 curr_obj = best_obj 
#                 curr_state = best_state
#             cand_state = deepcopy(curr_state)
              
#             #get the next candidate spr tree 
#             try:
#                 t = next(neighbors)
#             except StopIteration:
#                 # print("no other trees in neighborhood")
#                 # neighbors =iter(USPR(cand_state.get_unrooted_tree(), self.rng,min_radius=3, max_radius=3))
#                 # try:
#                 #     t = next(neighbors)
#                 break

#             cand_state.set_unrooted_tree(t)
#             curr_temp = temp*np.power(0.99,k)

#             #compute the objective value
#             cand_obj, cand_pars, cand_iso, _, _ = self.parsimony(cand_state)
#             if k % 25 ==0:
#                 print(f"{restart},{k},{cand_obj},{curr_obj},{best_obj}")     

#             #check if solution is better than current state or randomly jump to new state with probability P(cand_obj, curr_obj,curr_temp)
#             if cand_obj < curr_obj or  self.rng.random() <= self.prob_state_change(cand_obj, curr_obj, curr_temp ):
       
#                 curr_state = deepcopy(cand_state)
#                 neighbors = iter(USPR(curr_state.get_unrooted_tree(), self.rng))
#                 curr_obj = cand_obj
#                 if curr_obj < best_obj:
#                     best_obj = curr_obj
#                     best_state = deepcopy(curr_state)
                    
           
              
#         return best_obj, best_state