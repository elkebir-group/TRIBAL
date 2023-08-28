
import networkx as nx
import numpy as np
import argparse 
import sys, re, os

from copy import deepcopy
import utils as ut
import time
from ete3 import Tree
from itertools import repeat
from multi_spr import MultSPR
from score_class import Score
from lineage_tree import LineageForest, LineageTree
from multiprocessing import Pool
from steiner_tree import ConstructGraph, SteinerTree



class TribalSub:
    def __init__(self, 
                isotype_weights=None, 
                alpha=0.9, 
                n_isotypes=7,
                cost_function=None, 
                alphabet= ("A", "C", "G", "T","N", "-"), 
                timeout=2, 
                nworkers=1, 
                root_id="naive",
                reversible=False ):

        
        #abort search after timeout hours
        self.abort_after = timeout*60*60
        


        if isotype_weights is None:

            if not reversible:
                rev_val = np.Inf
            else:
                rev_val = 0

         
            print("isotype weights is none, using standard cost function")
                #use unweighted sankoff cost function
            self.states = [i for i in range(n_isotypes)]
            self.iso_weights = {}

            for s in range(n_isotypes):
                for t in range(n_isotypes):
                    if s > t:
                        self.iso_weights[s,t] = rev_val
                    elif s < t:
                        self.iso_weights[s,t] = 1
                    else:
                        self.iso_weights[s,t] =0



        else:
            if type(isotype_weights) is np.ndarray:
                self.iso_weights, self.states = ut.convert_transmat_to_weights(isotype_weights)
            else:
                self.iso_weights = isotype_weights
                self.states = list(set([s for s,t in self.iso_weights]))


        self.alpha=alpha
        self.alphabet = alphabet

        self.cost_function = cost_function
        self.nworkers = nworkers
        self.root_id = root_id
      



   


    
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

    def refine_ilp(self, lin_tree, alignment, isotype_labels):
          
            # lin_tree.save_png(f"test/init_tree{lin_tree.id}.png", isotype_labels)
  
            cg = ConstructGraph(self.iso_weights, isotype_labels, root_identifier=self.root_id)
       
                # lin_tree.save_png("curr_tree.png", isotype_labels, isotype_encoding)
            seq_score_prior, seq_labels = lin_tree.sequence_parismony(alignment)
            fg = cg.build(lin_tree, seq_labels)
            # fg.save_graph("test/flow_graph.png")
            # fg.save_graph(f"test/graphs/G{lin_tree.id}.png")
            st = SteinerTree(fg.G, fg.find_terminals(), fg.seq_weights, fg.iso_weights,root=self.root_id, lamb=self.alpha, threads=1 )
            obj, tree = st.run()
            # print(f"tree: {lin_tree.id} obj: {iso_score}")
            out_tree, out_iso = cg.decodeTree(tree)
            

            out_lt = LineageTree(out_tree, "naive", lin_tree.id, lin_tree.name)
            # out_lt.save_png(f"test/out_tree{lin_tree.id}.png", out_iso)

            #TODO: Speed up by not recomputing the sequences 

            seq_score, seq_labels = out_lt.sequence_parismony(alignment)
            assert seq_score_prior ==seq_score

            iso_score =(obj - self.alpha*seq_score)/ (1-self.alpha)

            sc = Score(obj, seq_score, iso_score, seq_labels, out_iso, out_lt)
            sc.check_score(self.iso_weights)

            return sc 

      
    
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
            
                    current_result = self.refine_ilp(cand_tribal, alignment, isotype_labels)

                 
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




    # def forest_infer(self, lin_forest, alignment=None, isotypes=None):
    #         if alignment is None:
    #             alignment = lin_forest.alignment
    #         if isotypes is None:
    #             isotype_labels = lin_forest.isotypes 
    #         else:
    #             isotype_labels = isotypes
            
    #         cg = ConstructGraph(isotype_weights, isotype_labels, root_identifier=self.root_id)
    #         refine_scores = {}
    #         best_score = np.Inf
    #         best_tree = None
    #         for lin_tree in lin_forest.get_trees():
    #             # lin_tree.save_png("curr_tree.png", isotype_labels, isotype_encoding)
    #             score, seq_labels = lin_tree.sequence_parismony(alignment)
    #             fg = cg.build(lin_tree, seq_labels)
    #             st = SteinerTree(fg.G, fg.seq_weights, fg.iso_weights,root=self.root_id, lamb=self.alpha )
    #             tree_score, tree = st.run()
    #             if tree_score < best_score:
    #                 best_tree, best_iso = cg.decodeTree(tree)
    #                 best_score = tree_score
    #             refine_scores[lin_tree.id] = tree_score
            
    #         best_lt = LineageTree(best_tree, "naive")
    #         _, seq_labels = best_lt.sequence_parismony(alignment)
            
      
    #         return refine_scores, best_lt, best_iso, seq_labels
            # combined_fg = cg.combineGraphs()
            # st =SteinerTree(combined_fg.G, combined_fg.seq_weights, combined_fg.iso_weights, root=self.root_id, lamb=self.alpha)
            # combined_score, combined_tree = st.run()
            # combined_tree, all_isotypes = cg.decodeTree(combined_tree)
            # combined_lt =LineageTree(combined_tree, "naive")
       

          
         
            # return refine_scores, combined_score, combined_lt, seq_labels, all_isotypes, best_lt, best_iso
        
            


    def forest_mode_loop(self, lin_forest, alignment=None, isotypes=None, mode="score"):
            
            # best_results = []
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
            elif mode == "refine_ilp":
                mode_func = self.refine_ilp
            elif mode == "search":
                mode_func = self.search 
            else:
                mode_func = self.score 
     
            all_results = [] 
            for lin_tree in lin_forest.get_trees():
         
                result = mode_func(lin_tree,alignment, isotype_labels)
                all_results.append(  result)
            #scan through results and find the top ntrees results 
            # for result in all_results:
            #     if len(best_results) ==0:
            #         best_results.append( result)
                  
            #         # print(best_result)
            #     else:
            #         if result.improvement(best_results[-1]) or len(best_results) < ntrees:
            #             self.update_best_results(result, best_results, ntrees)
                       
                      

      
            return  all_results



    def forest_mode(self, lin_forest, alignment=None, isotypes=None, mode="score", ntrees=1):
            
            best_results = []
        

            all_results = {}
            
            if alignment is None:
                alignment = lin_forest.alignment
            if isotypes is None:
                isotype_labels = lin_forest.isotypes 
            else:
                isotype_labels = isotypes

            if mode =="refine":
                mode_func = self.refine
            elif mode == "refine_ilp":
                mode_func = self.refine_ilp
            elif mode == "search":
                mode_func = self.search 
            else:
                mode_func = self.score 

        
            with Pool(self.nworkers) as pool:
                all_results = pool.starmap(mode_func,zip(lin_forest.get_trees(), repeat(alignment), repeat(isotype_labels)))     

         
         
            #scan through results and find the top ntrees results 
            # for result in all_results:
            #     if len(best_results) ==0:
            #         best_results.append( result)
                  
        
            #     else:
            #         if result.improvement(best_results[-1]) or len(best_results) < ntrees:
            #             self.update_best_results(result, best_results, ntrees)
                       
            

    
            return all_results
      
  
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
      
            nx_tree= convert_to_nx(ete_tree, args.root)
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


        for c in node.children:
            if c.name == "":
                c.name = str(internal_node)
                internal_node += 1
                internal_node_count += 1
    
       

            nx_tree.add_edge(node.name, c.name)
    
    if len(list(nx_tree.neighbors(root))) == 0:
   
        nx_tree.remove_edge(root_name, root)
        nx_tree.add_edge(root, root_name)
  


    return nx_tree
        
def update_best_results(new_result, results, ntrees):
        #find index where new result should be inserted
        new_score = new_result.objective
        added = False 
        for index, res in enumerate(results):
            if new_score < res.objective:
                 #insert new result in its proper place
                results.insert(index, new_result)
                added = True
                break 
        if not added and len(results) < ntrees:
            results.append(new_result)
        
        if len(results) > ntrees:
            #remove the highest scoring tree in the list
            del results[-1]

# def find_leaf_descendants(node, graph):
#     leaf_descendants = set()

#     # Helper function to traverse the graph
#     def dfs(current_node):
#         nonlocal leaf_descendants
#         # If the current node is a leaf, add it to the set
#         if graph.out_degree(current_node) == 0:
#             leaf_descendants.add(current_node)
#         else:
#             # Traverse all child nodes recursively
#             for child_node in graph.successors(current_node):
#                 dfs(child_node)

#     # Start the depth-first search from the specified node
#     dfs(node)
#     return leaf_descendants

# def get_clade_set(tree):
#     clade_set = []
#     for node in tree:
#         clade_set.append(find_leaf_descendants(node, tree))
    
#     return(set(map(frozenset, clade_set)))

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
    parser.add_argument("-l", "--lineage", type=str, help="pickle file of lineage tree/forest returned from tribal.py")
    parser.add_argument("--forest",  action="store_true")
    parser.add_argument("--candidates", type=str, help="filename containing newick strings for candidate tree(s)")
    parser.add_argument("--mode", choices=["score", "refine", "refine_ilp", "search"], default="score")
    parser.add_argument("-e", "--encoding", type=str, required=True)
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("-j", "--jump-prob", type=float, default=0.25)

    parser.add_argument("--ntrees", type=int, help="number of top scoring trees to return", default=1)
    parser.add_argument("-o", "--output", type=str, help="outputfile of all best trees")
    parser.add_argument("--tree",  type=str, help="outputfile of best tree")

    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--png", type=str, help="filename where to save a png of the optimal tree")
    parser.add_argument( "--all_pngs", action="store_true")
    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--reversible",  action="store_true", 
                        help="a flag to indicate the standard 0/1 cost function is used (the number of isotype changes is minimized and irreversiblility is ignored)")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")
    # parser.add_argument( "--ilp_sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    # parser.add_argument("--ilp_tree",  type=str, help="outputfile of best tree")
    # parser.add_argument("--ilp_iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")

    parser.add_argument("--all_optimal_sol",  help="path where all optimal solution results are saved"  )
    parser.add_argument("--nworkers", type=int, default=1, help="number of workers to use in the event of multiple input candidate trees")
    parser.add_argument("--seed", type=int, default=1026, help="random seed for picking a single best tree among all tied trees")
    # parser.add_argument("--all_obj", type=str, help="comparison of ilp heuristic versus min heuristc")
    parser.add_argument("--best_tree_diff", type=str, help="best tree RF distances")
    parser.add_argument("--pickle_best", type=str, help="filename to pickle the best results")
    # parser.add_argument("--combined_tree", type=str, help="png for combined tree")
    # parser.add_argument("--best_ilp", type=str, help="png for best ilp heuristic tree")


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    # path = "/scratch/projects/tribal/benchmark_pipeline"
    # n = 35
    # k = 25
    # r = 1
    # ttype = "seq"
    # clonotype = 15
    # alignment= f"{path}/sim_data/tmat_inf/{ttype}/cells{n}/size{k}/rep{r}/2.0/0.365/{clonotype}/GCsim_dedup.fasta"
    # isotypes = f"{path}/sim_data/tmat_inf/{ttype}/cells{n}/size{k}/rep{r}/2.0/0.365/{clonotype}/GCsim.isotypes"
    # transmat =   f"{path}/sim_data/tmat_inf/{ttype}/cells{n}/size{k}/rep{r}/2.0/0.365/tribal/transmat.txt"
    # candidates = f"{path}/sim_data/tmat_inf/{ttype}/cells{n}/size{k}/rep{r}/2.0/0.365/{clonotype}/dnapars/outtree"
    # lf =  f"{path}/sim_data/tmat_inf/{ttype}/cells{n}/size{k}/rep{r}/2.0/0.365/tribal_refine_ilp/20/forest.pickle"
    # encoding = f"{path}/sim_encoding.txt"
    # mode= "refine_ilp"
    
    # args = parser.parse_args([
    #     "-a", alignment,
    #     "-i", isotypes,
    #     "-t", transmat,
    #     "--nworkers", "1",
    #     "-r", "naive",
    #     "--candidates" ,candidates,
    #     # "-l", lf,
    #     # "--forest",
    #     "-e", encoding,
    #     "--mode", mode,
    #     # "--all_obj", f"test/{mode}.csv",
    #     "--png", f"test/{mode}.png",
    #     "--ntrees", "1",
    #     # "--all_optimal_sol", "test/opt_tree",
    #      "--tree", f"test/{mode}.txt",
    #     "--sequences", f"test/{mode}_seq.csv",
    #     "--iso_infer", f"test/{mode}_iso.csv",
    #     "--best_tree_diff", f"test/best_tree_rf.csv",
    #     # "--pickle_best", f"test/{mode}.pickle",
    #     "--score", f"test/{mode}.scores.csv"
    # ])

   

    # path = "/scratch/projects/tribal/experimental_data/day_14"
    # clonotype = "B_120_2_8_210_1_13"
    # alignment= f"{path}/input/{clonotype}/concat.aln.fasta"
    # candidates = f"{path}/dnapars/{clonotype}/outtree" 
    # isotypes = f"{path}/input/{clonotype}/isotype.fasta"
    # encoding = "/scratch/projects/tribal/experimental_data/mouse_isotype_encoding.txt"
    # transmat = "/scratch/projects/tribal/experimental_data/day_14/tribal/0.1/transmat.txt"
    # # transmat = "/scratch/projects/tribal/benchmark_pipeline/sim_data/tmat_inf/direct/transmats/transmat2.txt"
    # mode= "refine_ilp"
    # args = parser.parse_args([
    #     "-a", alignment,
    #     "-i", isotypes,
    #     "-t", transmat,
    #     "--nworkers", "10",
    #     "-r", "naive",
    #     "--candidates" ,candidates,
    #     "-e", encoding,
    #     "--mode", mode,
    #     # "--all_obj", f"test/{mode}.csv",
    #     "--png", f"test/{mode}3.png",
    #     "--ntrees", "2",
    #     "--all_optimal_sol", "test/opt_tree",
    #      "--tree", f"test/{mode}.txt",
    #     "--sequences", f"test/{mode}_seq.csv",
    #     "--iso_infer", f"test/{mode}_iso.csv",
    #     "--best_tree_diff", f"test/best_tree_rf.csv",
    #     "--pickle_best", f"test/{mode}3.pickle"
    # ])
  

  

    
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
              


    tr = TribalSub(isotype_weights, args.alpha, timeout=args.timeout, nworkers=args.nworkers, root_id=args.root, reversible=args.reversible)
    ncells = len(lin_forest.alignment)
    print(f"\nInput:\nncells: {ncells}\nforest size: {lin_forest.size()}\nmode: {args.mode}\n")

    # flow_scores, combined_score,  combined_lt, seq_labels, all_isotypes, best_ilp_lt, best_iso_ilp = tr.forest_infer(lin_forest)
    # flow_scores, best_ilp_lt, best_iso_ilp, ilp_seq_labels =tr.forest_infer(lin_forest)
    # if args.combined_tree is not None:
    #     combined_lt.save_png(args.combined_tree, all_isotypes, isotype_encoding)

    # if args.best_ilp is not None:
    #     best_ilp_lt.save_png(args.best_ilp, best_iso_ilp, isotype_encoding)
    
    # ilp_best_labels = ut.update_labels(ilp_seq_labels)

    # if args.ilp_sequences is not None:
    #     ut.save_dict(ilp_best_labels, args.ilp_sequences)
    

    # if args.ilp_tree is not None:
    #     best_ilp_lt.save_tree(args.ilp_tree)
       
    # if args.ilp_iso_infer is not None:
    #     ut.save_dict(best_iso_ilp, args.ilp_iso_infer)
    # for lin in lin_forest:

    #     if lin.id ==3:
    #         clade_sets_11 = []
    #         for node in lin.T:
    #             clade_sets_11.append(find_leaf_descendants(node, lin.T))
    #         t11 = set(map(frozenset, clade_sets_11))
    #     if lin.id ==21:
    #         clade_sets_21 = []
    #         for node in lin.T:
    #             clade_sets_21.append(find_leaf_descendants(node, lin.T))
    #         t21 = set(map(frozenset, clade_sets_21))
    # symmetric_diff = t11.symmetric_difference(t21)
    # print(len(symmetric_diff))
    # rf_dist = 0.5*len(symmetric_diff)
    # print(rf_dist)
      
            
            # lin.save_png(f"test/tree_{lin.id}.png", lin_forest.isotypes, isotype_encoding)
    all_results =  tr.forest_mode_loop(lin_forest, mode =args.mode)

    print(len(all_results))
    for a in all_results:
        print(a)
         #scan through results and find the top ntrees results 
    best_results = []
    top_ntrees = []
    best_obj = np.Inf
    for res in all_results:
        if res.objective < best_obj:
            best_obj = res.objective
       
    
    for res in all_results:
        if round(res.objective,5) == round(best_obj,5):
            best_results.append(res)

        if len(top_ntrees) ==0:
            top_ntrees.append(res)
            
            # print(best_result)
        else:
            if res.improvement(top_ntrees[-1]) or len(top_ntrees) < args.ntrees:
                update_best_results(res, top_ntrees, args.ntrees)

    if args.best_tree_diff is not None:       
        with open(args.best_tree_diff, "w+") as outfile:
            outfile.write("tree1,tree2,rf\n")
            for i,res1 in enumerate(best_results):
                

                for j, res2 in enumerate(best_results):
                    if i < j:
                        rf = res1.tree.rf_distance(res2.tree)
                        outfile.write(f"{res1.tree.id},{res2.tree.id},{rf}\n")
    



    
    lin_forest_out = LineageForest(lin_forest.alignment, lin_forest.isotypes, [res.tree for res in top_ntrees])
    
    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write("tree,alpha,objective,sequence,isotype,\n")
            
            # file.write(f"combined,{combined_score},NA\n")    
            for res in all_results:
                file.write(f"{res.tree.id},{args.alpha},{res.objective},{res.seq_obj},{res.iso_obj}\n")
    #randomly pick the best result if there are multiple trees with the same optimal objective
    # top_score = best_results[0].objective
    # tied_trees = []
    # for res in best_results:
    #     if res.objective ==top_score:
    #         tied_trees.append(res)
    # rng = np.random.default_rng(args.seed)
    # best_result = rng.choice(tied_trees, 1)[0]
    best_result = best_results[0]
    
    best_tree= best_result.tree
    print(f"{args.mode} complete! \nBest Tree: {best_result.tree.id}")
    print(best_results[0])
    print("\nsaving results......")
    # if args.all_obj is not None:
    #     with open(args.all_obj, 'w+') as file:
    #         file.write("tree,ilp,heuristic\n")
            
    #         # file.write(f"combined,{combined_score},NA\n")    
    #         for key,val in flow_scores.items():
    #             heur = all_results[key]
    #             file.write(f"{key},{val},{heur.objective}\n")


    if args.output is not None:
        lin_forest_out.save_forest(args.output)
    

    
    if args.png is not None:
        best_tree.save_png(args.png, best_result.isotypes, isotype_encoding)
    

    best_labels = ut.update_labels(best_result.labels)

    if args.fasta is not None:
        
        ut.write_fasta(args.fasta, best_labels)
    
    if args.sequences is not None:
        ut.save_dict(best_labels, args.sequences)
    

    if args.tree is not None:
        best_tree.save_tree(args.tree)
    
    if args.iso_infer is not None:
        ut.save_dict(best_result.isotypes, args.iso_infer)
    

    # for i,res in enumerate(best_results):
    #     if i ==0:
    #         png_fname = args.png 
    #         fasta_fname =args.fasta
    #         seq_fname = args.sequences 
    #         tree_fname = args.tree
    #         iso_fname =args.iso_infer
    #     elif args.all_optimal_sol:
    #         if args.png is not None:
    #             png_fname =  f"test/{mode}.{i}.png"
            
         
    #         if args.fasta is not None:
    #             fasta_fname =args.fasta + "." + str(i)
    #         if args.sequences is not None:
    #             seq_fname = args.sequences + "." + str(i)
    #         if args.tree is not None:
    #             tree_fname = args.tree + "." + str(i)
    #         if args.iso_infer is not None:
    #             iso_fname =args.iso_infer + "." + str(i)
    #     else:
    #         break
   
    
        # if args.score:
        #     with open(args.score, 'w+') as file:
        #         for result in all_results:
        #             file.write(f"{result.tree.id},{args.alpha},{result.objective},{result.seq_obj},{result.iso_obj}\n")
        
    if args.all_optimal_sol is not None:
        pth =args.all_optimal_sol
        if not os.path.exists(pth):
            # Create the pth
            os.makedirs(pth)
            print("Directory created:", pth)
        else:
            print("Directory already exists:", pth)

        for res in best_results:
            labs = ut.update_labels(res.labels)

            ut.write_fasta(f"{pth}/tree_{res.tree.id}.seq.fasta", labs)
            ut.save_dict(labs, f"{pth}/tree_{res.tree.id}.seq.csv")
            ut.save_dict(res.isotypes, f"{pth}/tree_{res.tree.id}.isotypes.csv")
            res.tree.save_tree(f"{pth}/tree_{res.tree.id}.txt")
            res.tree.save_png(f"{pth}/tree_{res.tree.id}.png", res.isotypes, isotype_encoding)

    if args.pickle_best is not None:
        ut.pickle_save(best_results, args.pickle_best)
    