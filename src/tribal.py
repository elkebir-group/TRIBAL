
import networkx as nx
import numpy as np
import sys, os, re
import argparse 
import pickle
from copy import deepcopy
import utils as ut
from ete3 import Tree
from em_weight_matrix import EMProbs
from max_likelihood_trans_probs import MaxLike
from lineage_tree import LineageForest
import init_transmat as tm
from alignment import Alignment

from tribal_tree import TribalSub
from draw_state_diagram import DrawStateDiag
from score_class import ScoreList

class Tribal:
    def __init__(self, 
                clonotypes,  # a dictionary of lineage forests (set of candidate trees, alignment, isotypes for observed cells)
                transmat=None, 
                alpha=0.9, 
                alphabet= ("A", "C", "G", "T","N", "-"), 
                isotype_encoding=None, 
                seed= 1026, 
                max_cand=50, 
                niter=10,
                threshold=0.5, 
                restarts=5,
                n_isotypes = 7, 
                not_trans_prob=0.65, 
                mode="refine_ilp" ):
        
        self.mode= mode
        self.clonotypes = clonotypes
        self.isotype_encoding = isotype_encoding
        self.alphabet = alphabet
        self.alpha = alpha
        self.not_trans_prob = not_trans_prob
        self.min_iterations = min(4, niter)

        if transmat is None:
            if isotype_encoding is not None:
                self.n_isotypes = len(isotype_encoding)
            else:
                self.n_isotypes = n_isotypes 
            print(f"generating transitiom matrix with stay probability {self.not_trans_prob}")
            self.transmat = tm.gen_trans_mat(self.not_trans_prob, self.n_isotypes)
      
        else:
            self.transmat = transmat
        
        self.init_transmat = self.transmat.copy()

    
        # self.states = np.arange(self.transmat.shape[0])
        self.states = np.arange(n_isotypes)

        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.candidates = {}
        self.max_cand = max_cand
        self.n_isotypes = self.transmat.shape[0]

        self.obs_states = {key: val.isotypes for key,val in self.clonotypes.items()}
        self.threshold = threshold
        self.niterations = niter
        self.restarts = restarts 



     
    

    def intialize_candidates(self, best_tree_ids=None):
        '''
        randomly initialize a set of candidate trees up to max_cands 
        for each clonotype, including the best trees found so far is dictionary
        is given.'''
        candidates = {}
        for c in self.clonotypes:
            # print(f"clontoype {c}: size: {self.clonotypes[c].size()}")
            if self.clonotypes[c].size() > self.max_cand:
                cand = self.rng.choice(self.clonotypes[c].size(), self.max_cand, replace=False).tolist()
            else:
                cand = [i for i in range(self.clonotypes[c].size())]
            
            candidates[c]= LineageForest( alignment=self.clonotypes[c].alignment, isotypes=self.clonotypes[c].isotypes)
            for index in cand:
                candidates[c].add(deepcopy(self.clonotypes[c][index]))
            if best_tree_ids is not None:
                for i in best_tree_ids[c]:
                    if i not in cand:
                        candidates[c].add(deepcopy(self.clonotypes[c][i]))

        return candidates 
    




    def score_candidates(self, candidates, transmat=None, mode=None, nproc=1):
            if mode is None:
                mode = self.mode 
            total_likelihood = 0
      
            best_scores = []
            best_tree_ids = {}
            all_best_scores = []
            refined_cands = {}
            for i,c in enumerate(self.clonotypes):
                    ts = TribalSub( isotype_weights=transmat,alpha=self.alpha, nworkers=nproc)
                    if mode=="refine_ilp":
                        all_scores = ts.forest_mode(candidates[c], mode=mode)
                    else:
                        all_scores = ts.forest_mode_loop(candidates[c], mode=mode)
                    sl = ScoreList(all_scores)
                    best_obj, best_scores = sl.find_best_scores()
                    total_likelihood += best_obj
                    all_best_scores.append(best_scores)
           

                    best_tree_ids[c] = [score.tree.id for score in best_scores]
         
                    bs = best_scores[0]
                    lin_tree = bs.tree
                    leaf_isotypes = {l: bs.isotypes[l] for l in lin_tree.get_leafs() }
                    msa = {l: bs.labels[l] for l in lin_tree.get_leafs() }
                    leaf_isotypes[lin_tree.root] = bs.isotypes[lin_tree.root]
                    msa[lin_tree.root] = bs.labels[lin_tree.root]
                    refined_cands[c] = LineageForest(msa, leaf_isotypes, [lin_tree])
               

            return total_likelihood,all_best_scores, best_tree_ids, refined_cands 
    
    @staticmethod
    def check_convergence(old, new, threshold):
        return np.abs(old -new) < threshold


    
    def fit(self, nproc=1):
  
     
        ''' 
            jointly infer isotype transition probabilities for a set of clonotypes
            and a B cell lineage tree that maximizes the CSR likelihood
            
        '''
       
        cand_tmat = []
        cand_scores = []
        best_trees = []
        stay_probs = np.linspace(0.55,0.95, self.restarts)
        for i in range(self.restarts):
            print(f"\nStarting Cycle {i}...")

            transmat =  tm.gen_trans_mat(stay_probs[i], self.n_isotypes)
            best_tree_ids = None 
            old_score = np.Inf
            
            for j in range(self.niterations):
 
                candidates = self.intialize_candidates(best_tree_ids)
                current_score, best_scores, best_tree_ids, _ = self.score_candidates(candidates, transmat=transmat, nproc=nproc)
                print(f"iteration: {j} old score: {old_score} current score: {current_score}")

                if self.check_convergence(current_score, old_score, self.threshold):
                        cand_tmat.append((transmat, state_probs))
                        cand_scores.append(current_score)
                        best_trees.append(best_scores)
                        break
                else:
                    old_score = current_score
                    transmat, state_probs = MaxLike(self.n_isotypes).infer(best_scores) 
        
        min_value = min(cand_scores)
        min_index = cand_scores.index(min_value)
        transmat, state_probs = cand_tmat[min_index]
        best_scores = best_trees[min_index]
        best_scores = ScoreList([b[0] for b in best_scores])
        
        return min_value, transmat, state_probs, best_scores

         



  
      

##Helper functions to preproces input files


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
        
def create_isotype_encoding(fname):

    iso_encoding = {}
    counter = 0

    with open(fname, "r") as file:
        for line in file:
            isotype = line.strip()
            if counter ==0:
                start_iso = isotype 
            iso_encoding[isotype] = counter
            counter += 1
    return iso_encoding, start_iso, counter


def create_input( path,  tree_path, clonotype, root, seq_fasta_fname, 
                 trees_fname, iso_fasta_fname, iso_encoding=None, start_iso=None):

    tree_fname =f"{tree_path}/{clonotype}/{trees_fname}"
    align_fname = f"{path}/{clonotype}/{seq_fasta_fname}"
    iso_fname =f"{path}/{clonotype}/{iso_fasta_fname}"
    tree_list = create_trees(tree_fname)

    #simplified alignment 
    alignment = Alignment(align_fname,root=args.root).simplify()
    if ".fasta" in iso_fname:
        isotypes = ut.read_fasta(iso_fname)
    else:
        isotypes = ut.read_dict(iso_fname)

    if iso_encoding is not None and start_iso is not None:
        isotypes_filt = {}
        for i in alignment:
                iso = isotypes[i]
                if iso not in iso_encoding:
                    iso = isotypes[i].lower()
                    if 'm' in iso or 'd' in iso:
                        iso = start_iso
                isotypes_filt[i] = iso_encoding[iso]
        isotypes = isotypes_filt
    
    linforest = LineageForest(alignment=alignment, isotypes=isotypes)
    linforest.generate_from_list(tree_list, root)

    return linforest

  

def save_results(outpath, lin_tree_dict, pngs=False, isotype_mapping=None):
   
    for clono, res in lin_tree_dict.items():
        clono_path = f"{outpath}/{clono}"
        os.makedirs(clono_path, exist_ok=True)
        tree = res["tree"]
        seq =  ut.update_labels(res["labels"])
        iso = res["isotypes"]

      
  
        tree.save_tree(f"{clono_path}/tree.txt")
        if pngs:
            tree.save_png(f"{clono_path}/tree.png", iso, isotype_mapping)
     
        if isotype_mapping is not None:
            iso_labs = {key: isotype_mapping[val] for key,val in iso.items()}
        else:
            iso_labs =iso 
        ut.write_fasta(f"{clono_path}/seq.fasta", seq)
        ut.write_fasta(f"{clono_path}/isotypes.fasta", iso_labs)
        ut.save_dict(f"{clono_path}/seq.csv", seq)
        ut.save_dict(f"{clono_path}/isotypes.csv", iso_labs)



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


def pickle_save(obj, fname):
        with open(fname, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--forest", type=str, help="path to pickled clonotypes dictionary of lineeage forests" )
    parser.add_argument("-p", "--path", type=str, required=False, help="path to the directory containing input files")
    parser.add_argument("-c", "--clonotypes", required=False, type=str,
        help="filename with list of clonotype subdirectories that should be included in the inference. If not provided, scans provided path for all subdirectory names")
    parser.add_argument("-e", "--encoding", type=str, help="text file isotype states listed in germline order")
    parser.add_argument("--n_isotypes", type=int, default=7, help="the number of isotypes states to use if isotype encoding file is not provided and input isotypes are encoded numerically")
    parser.add_argument( "--fasta", type=str, default= "concat.aln.fasta", help="filename of input MSA in fasta file")
    parser.add_argument("-i", "--isotypes",  type=str, default= "isotype.fasta",
        help="filename of isotype fasta file within each clonotype directory")
    parser.add_argument("-j", "--jump_prob", type=float, default=0.25, help="for inititalization of transition matrix if not provided")
    parser.add_argument("-t", "--transmat", required=False, type=str,
        help="optional filename of input transition matrix for initialization")
    parser.add_argument("-r", "--root", required=False, default="naive",
        help="the common id of the root in all clonotypes")
    parser.add_argument( "--tree_path", type=str, required=False, help="path to directory where candidate trees are saved")
    parser.add_argument("--candidates", type=str, default="outtree", help="filename containing newick strings for candidate trees")
    parser.add_argument("--niter", type=int, help="max number of iterations in the fitting phase", default=10)
    parser.add_argument("--thresh", type=float, help="theshold for convergence in fitting phase" ,default=0.1)
    parser.add_argument("--nworkers", type=int, default=2, help="number of workers to use in the event in multiple restarts")
    parser.add_argument("--max_cand", type=int, default = 20,  help="max candidate tree size per clonotype")
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument("--restarts",  type=int, default=1, help="number of restarts")
    parser.add_argument("--mode", choices=["score", "refine", "refine_ilp", "search"], default="score")
    parser.add_argument("--score", type=str, help="filename where the score file should be saved")
    parser.add_argument("--transmat_infer", type=str, help="filename where the inferred transition matrix should be saved")
    parser.add_argument("--state_probs", type=str, help="filename where the inferred state probabilities should be saved")
    parser.add_argument("--heatmap", type=str, help="filename where the {png,pdf} of transition matrix should be saved")
    parser.add_argument("--propmap", type=str, help="filename where the {pdf,png} of isotype proportions should be saved")


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
  

    if args.encoding is not None:
        iso_encoding, start_iso, n_isotypes = create_isotype_encoding(args.encoding)
        rev_encoding = {val: key for key, val in iso_encoding.items()}
    else:
        n_isotypes = args.n_isotypes
        iso_encoding = None
        start_iso= None 
        rev_encoding = None
    
    
    if args.transmat is not None:

        transmat = np.loadtxt(args.transmat)
    else:
        transmat= None


    
    if args.forest is not None:
        clonodict = ut.pickle_load(args.forest)
    
    else:
        if args.clonotypes is not None:
            clonotypes = []
            with open(args.clonotypes, 'r+') as file:
                for line in file:
                    clonotypes.append(line.strip())

        else:
            clonotypes = [it.name for it in os.scandir(args.path) if it.is_dir()]

        clonodict = {}
        for c in clonotypes:
            print(f"reading input for clonotype {c}")
            clonodict[c] = create_input(args.path, args.tree_path, c, args.root, args.fasta, 
                            args.candidates, args.isotypes, iso_encoding, start_iso)
    

    tr= Tribal(clonodict, 
                transmat, 
                seed = args.seed,
                isotype_encoding= iso_encoding,
                max_cand= args.max_cand,
                niter = args.niter,
                threshold=args.thresh,
                not_trans_prob= 1-args.jump_prob,
                restarts=args.restarts,
                mode = args.mode
                )
    

  
    obj_score, transmat, state_probs,  best_trees= tr.fit(args.nworkers)


    print("\nTRIBAL Complete!, saving results...")

 

    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write(str(obj_score))

    if args.transmat_infer is not None:
        np.savetxt(args.transmat_infer, transmat)
    if args.state_probs is not None:
        np.savetxt(args.state_probs, state_probs)
    
    if args.heatmap is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).heatmap(args.heatmap)
   
    if args.propmap is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).state_heatmap(args.propmap)



        
    


  









