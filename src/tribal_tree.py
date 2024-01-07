
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
from score_class import Score, ScoreList
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

  
            cg = ConstructGraph(self.iso_weights, isotype_labels, root_identifier=self.root_id)
           

            seq_score_prior, seq_labels = lin_tree.sequence_parismony(alignment)
            fg = cg.build(lin_tree, seq_labels)

            st = SteinerTree(fg.G, lin_tree.T, fg.find_terminals(), fg.seq_weights, 
                             fg.iso_weights,fg.node_mapping, fg.tree_to_graph,
                               fg.node_out_degree, pars_score = seq_score_prior, root=self.root_id, lamb=self.alpha, threads=1 )
            obj, tree = st.run()

            out_tree, out_iso = cg.decodeTree(tree)
            

            out_lt = LineageTree(out_tree, "naive", lin_tree.id, lin_tree.name)
    

            seq_score, seq_labels = out_lt.sequence_parismony(alignment)
            assert seq_score_prior ==seq_score



            sc = Score(obj, seq_score, obj, seq_labels, out_iso, out_lt)
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




        
            


    def forest_mode_loop(self, lin_forest, alignment=None, isotypes=None, mode="score"):
            
    

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
                       
                      

      
            return  all_results



    def forest_mode(self, lin_forest, alignment=None, isotypes=None, mode="score", ntrees=1):
            
         

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

        new_score = new_result.objective
        added = False 
        for index, res in enumerate(results):
            if new_score < res.objective:

                results.insert(index, new_result)
                added = True
                break 
        if not added and len(results) < ntrees:
            results.append(new_result)
        
        if len(results) > ntrees:

            del results[-1]



def get_alignment(fname):
    alignment = ut.read_fasta(fname)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}
    return alignment

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--input-forest", type=str, help="path to pickled clonotypes dictionary of lineeage forests" )
    parser.add_argument("-c", "--clonotype", type=str, help="name of clonotype lineage to refine from the full forest" )
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
    parser.add_argument("--tree",  type=str, help="outputfile of best tree")
    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--png", type=str, help="filename where to save a png of the optimal tree")
    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--reversible",  action="store_true", 
                        help="a flag to indicate the standard 0/1 cost function is used (the number of isotype changes is minimized and irreversiblility is ignored)")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")
    parser.add_argument("--all_optimal_sol",  help="path where all optimal solution results are saved"  )
    parser.add_argument("--nworkers", type=int, default=1, help="number of workers to use in the event of multiple input candidate trees")
    parser.add_argument("--pickle_best", type=str, help="filename to pickle the best results")
    parser.add_argument("--pickle_all", type=str, help="filename to pickle the best results")


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])


  
    
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
        if args.clonotype is not None and args.input_forest is not None:
            full_forest = ut.pickle_load(args.input_forest)
            lin_forest = full_forest[args.clonotype]

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
              
 

    ncells = len(lin_forest.alignment)
    print(f"\nInput:\nncells: {ncells}\nforest size: {lin_forest.size()}\nmode: {args.mode}\n")

    tr = TribalSub(isotype_weights, 0.9, timeout=args.timeout, nworkers=args.nworkers, root_id=args.root, reversible=args.reversible)

    all_results =  tr.forest_mode(lin_forest, mode =args.mode)

    
    for a in all_results:
        print(a)

    score_list = ScoreList(all_results)
    best_score, best_results = score_list.find_all_best_scores()


        
    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write("tree,alpha,objective,sequence,isotype,\n")
            
            for res in all_results:
                file.write(f"{res.tree.id},{res.objective},{res.seq_obj},{res.iso_obj}\n")
    best_result = best_results[0]
    
    best_tree= best_result.tree
    print(f"{args.mode} complete! \nBest Tree: {best_result.tree.id}")
    print(best_results[0])
    print("\nsaving results......")



    

    
    if args.png is not None:
        best_tree.save_png(args.png, best_result.isotypes, isotype_encoding, show_labels=True)
    

    best_labels = ut.update_labels(best_result.labels)

    if args.fasta is not None:
        
        ut.write_fasta(args.fasta, best_labels)
    
    if args.sequences is not None:
        ut.save_dict(best_labels, args.sequences)
    

    if args.tree is not None:
        best_tree.save_tree(args.tree)
    
    if args.iso_infer is not None:
        ut.save_dict(best_result.isotypes, args.iso_infer)
    

        
    if args.all_optimal_sol is not None:
        pth =args.all_optimal_sol
        if not os.path.exists(pth):
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
    
    if args.pickle_all is not None:
        ut.pickle_save(all_results, args.pickle_all)
    