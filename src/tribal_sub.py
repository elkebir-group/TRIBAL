
import networkx as nx
import numpy as np
import neighbor_joining as nj
import argparse 
import re
from spr import SPR
from unrooted_spr import USPR
from copy import deepcopy
import utils as ut
import time
from ete3 import Tree
from itertools import repeat
from multi_spr import MultSPR
from score_class import Score
from lineage_tree import LineageForest
from multiprocessing import Pool


class TribalSub:
    def __init__(self, isotype_weights=None, alpha=0.9, n_isotypes=7,
                cost_function=None, 
                alphabet= ("A", "C", "G", "T","N", "-"), timeout=2, nworkers=1 ):

        
        #abort search after timeout hours
        self.abort_after = timeout*60*60
        


        if isotype_weights is None:

            #use unweighted sankoff cost function
            self.states = [i for i in range(n_isotypes)]
            self.iso_weights = {}

            for s in range(n_isotypes):
                for t in range(n_isotypes):
                    if s > t:
                        self.iso_weights[s,t] = np.Inf
                    elif s < t:
                        self.iso_weights[s,t] = 1
                    else:
                        self.iso_weights[s,t] =0



        else:
            self.iso_weights = isotype_weights
            self.states = list(set([s for s,t in self.iso_weights]))

  
    
        self.alpha=alpha
        self.alphabet = alphabet

        self.cost_function = cost_function
        self.nworkers = nworkers
      



   


    
    def score(self, lin_tree, alignment, isotype_labels):
        
        seq_score, seq_labels = lin_tree.sequence_parismony(alignment, 
                                                    alphabet=self.alphabet, 
                                                    cost_function=self.cost_function)
        iso_score, iso_labels = lin_tree.isotype_parsimony(isotype_labels,  weights=self.iso_weights, states=self.states)
        obj = self.compute_score(seq_score, iso_score)

        return   Score(obj, seq_score, iso_score, seq_labels, iso_labels, lin_tree)


    def refine(self, lin_tree, alignment, isotype_labels):
        iso_score, iso_labels = lin_tree.isotype_parsimony_polytomy(isotype_labels, weights=self.iso_weights, states=self.states )

        seq_score, seq_labels = lin_tree.sequence_parismony(alignment, 
                                                    alphabet=self.alphabet, 
                                                    cost_function=self.cost_function)
        obj = self.compute_score(seq_score, iso_score)
        return Score(obj, seq_score, iso_score, seq_labels, iso_labels, lin_tree)
    
    def search(self, lin_tree, alignment, isotype_labels, mode="multSPR"):

            iterations = 0
            timeout_start = time.time()
     
            best_result = self.refine(lin_tree, alignment, isotype_labels)
    
            count = 0
            improvement = False
            while time.time() < timeout_start + self.abort_after:
                time.sleep(0.25)
                iterations += 1
                cand_tribal = deepcopy(best_result.tree)
    
                if mode == "multSPR":
                    spr_trees = MultSPR(cand_tribal, best_result.isotypes, self.states)
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
                
                        break

                    if time.time() > timeout_start + self.abort_after:
                        break
             
                if not improvement:
                    print("Total Trees Examined: " + str(count))
                    break
                
 
            
            return best_result


    
    def compute_score(self, pars, iso):
        return self.alpha*pars + (1-self.alpha)*iso

    @staticmethod
    def update_best_results(new_result, results, ntrees):
        #find index where new result should be inserted
        new_score = new_result.objective
       
        for index, res in enumerate(results):
            if new_score < res.objective:
                 #insert new result in its proper place
                results.insert(index, new_result)
                break 
            if index == len(results) -1 and len(results) <= ntrees:
                results.append(new_result)
        
        if len(results) > ntrees:
            #remove the highest scoring tree in the list
            del results[-1]

           
   
    def forest_mode(self, lin_forest, alignment=None, isotypes=None, mode="score", ntrees=1):
            
            best_results = []
            # best_result = None 
            # best_tree = None

            all_results = {}
            
            if alignment is None:
                alignment = lin_forest.alignment
            if isotypes is None:
                isotype_labels = lin_forest.isotypes 
            else:
                isotype_labels = isotypes

            if mode =="refine":
                mode_func = self.refine
            elif mode == "search":
                mode_func = self.search 
            else:
                mode_func = self.score 

                #    M = pool.starmap(func, zip(a_args, repeat(second_arg)))
            with Pool(self.nworkers) as pool:
                all_results = pool.starmap(mode_func,zip(lin_forest.get_trees(), repeat(alignment), repeat(isotype_labels)))     
            # for lin_tree in lin_forest.get_trees():
         
            #     result = mode_func(lin_tree,alignment, isotype_labels)
            #     all_results[lin_tree.id] = result
            #scan through results and find the top ntrees results 
            for result in all_results:
                if len(best_results) ==0:
                    best_results.append( result)
                  
                    # print(best_result)
                else:
                    if result.improvement(best_results[-1]) or len(best_results) < ntrees:
                        self.update_best_results(result, best_results, ntrees)
                       
                        # print(best_result)

      
            return best_results, all_results
      
  
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
    parser.add_argument("--timeout", type = float, help="max number of hours to let tribal search per tree", default=8)
    parser.add_argument("-l", "--lineage", type=str, help="pickle file of lineage tree/forest")
    parser.add_argument("--candidates", type=str, help="filename containing newick strings for candidate tree(s)")
    parser.add_argument("--mode", choices=["score", "refine", "search"], default="score")
    parser.add_argument("-e", "--encoding", type=str, required=True)
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("-j", "--jump-prob", type=float, default=0.25)
    parser.add_argument("--forest",  action="store_true")
    parser.add_argument("--ntrees", type=int, help="number of top scoring trees to return", default=1)


    parser.add_argument("-o", "--output", type=str, help="outputfile of all best trees")
    parser.add_argument("--tree",  type=str, help="outputfile of best tree")

    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--png", type=str, help="filename where to save a png of the optimal tree")
    parser.add_argument( "--all_pngs", action="store_true")

    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("-n", "--newick", type=str, help="filename where newick string should be saved")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")
    parser.add_argument("--save_candidates", type=str, help="directory where to save data for candidate trees")
    parser.add_argument("--save_all_scores", type=str, help="file where to save data for candidate trees")
    parser.add_argument("--nworkers", type=int, default=1, help="number of workers to use in the event in multiple restarts")



    args= parser.parse_args()




    # path = "/scratch/projects/tribal/real_data"
    # dataset = "day_14"
    # folder = "tribal"
    # clonotype = "B_147_6_76_148_1_41"

    # args =parser.parse_args([
    #     "-a", f"{path}/{dataset}/input/{clonotype}/concat.aln.fasta",
    #     "-r", "naive",
    #     "-t", f"{path}/{dataset}/{folder}/0.25/transmat.txt",
    #     # "-l", "/scratch/projects/tribal/real_data/test/best_forest.pickle",
    #     # "--forest",
    #     "-e", f"{path}/mouse_isotype_encoding.txt",
    #     "-i", f"{path}/{dataset}/input/{clonotype}/isotype.fasta",
    #     "-o", f"{path}/test/best_forest.pickle",
    #     # "--init", "candidates",
    #     # "-s", "1",
    #     "--alpha", "0.75",
    #     "--sequences", f"{path}/test/tribal.seq",
    #     "--fasta",f"{path}/test/tribal.fasta",
    #     "--score",f"{path}/test/tribal.score",
    #     "--iso_infer", f"{path}/test/tribal.isotypes",
    #     "--candidates", f"{path}/{dataset}/dnapars/{clonotype}/outtree",
    #     # "--save_candidates",f"{path}/test",
    #     # "--save_all_scores",f"{path}/test/all_scores.csv",
    #     "--png", f"{path}/test/best_tree.png",
    #     "--mode", "score",
    #     # "--timeout", "0.1"
    #     "--ntrees", "10"
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
        if ".fasta" in args.isotypes:

            isotypes = ut.read_fasta(args.isotypes)
        else:
            isotypes = ut.read_dict(args.isotypes)
        
            
        isotypes_filt = {}
        for i in isotypes:
            iso = isotypes[i]
            if iso not in iso_encoding:
                iso = isotypes[i].lower()
                if 'm' in iso or 'd' in iso:
                    iso = start_iso
            isotypes_filt[i] = iso_encoding[iso]
                
    
    isotype_weights= None 
    if args.transmat is not None:
       isotype_weights, states = ut.convert_transmat_to_weights(np.loadtxt(args.transmat))



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
                lin_forest = LineageForest( alignment, isotypes_filt, [lin])
              


    tr = TribalSub(isotype_weights, args.alpha, timeout=args.timeout, nworkers=args.nworkers)
    ncells = len(lin_forest.alignment)
    print(f"\nInput:\nncells: {ncells}\nforest size: {lin_forest.size()}\nmode: {args.mode}\n")

    best_results, all_results =  tr.forest_mode(lin_forest, mode =args.mode, ntrees=args.ntrees)

    
    lin_forest_out = LineageForest(lin_forest.alignment, lin_forest.isotypes, [res.tree for res in best_results])
    

    best_tree = best_results[0].tree
  
    print(f"{args.mode} complete! \nBest Tree: {best_results[0].tree.id}")
    print(best_results[0])
    print("\nsaving results......")


    if args.output is not None:
        lin_forest_out.save_forest(args.output)
        if args.png is not None:
            best_tree.save_png(args.png, best_results[0].isotypes, isotype_encoding)
        
    
    best_labels = ut.update_labels(best_results[0].labels)

    if args.fasta is not None:
        
        ut.write_fasta(args.fasta, best_labels)
    
    if args.sequences is not None:
        ut.save_dict(best_labels, args.sequences)
    

    if args.tree is not None:
        best_results[0].tree.save_tree(args.tree)
       
    if args.iso_infer is not None:
        ut.save_dict(best_results[0].isotypes, args.iso_infer)
    
    if args.score:
        with open(args.score, 'w+') as file:
            for result in all_results:
                file.write(f"{result.tree.id},{args.alpha},{result.objective},{result.seq_obj},{result.iso_obj}\n")
    
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