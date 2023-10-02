
import networkx as nx
import numpy as np
import sys, os, re
import argparse 
import pickle
from multiprocessing import Pool
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
                threshold=0.1, 
                restarts=5,
                n_isotypes = 7, 
                not_trans_prob=0.65, 
                mu=0.07, 
                sigma=0.05, 
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
        self.mu = mu
        self.sigma= sigma


     
    

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
            # if best_tree_ids is not None:
            #     for i in best_tree_ids[c]:
            #         if i not in cand:
            #             candidates[c].add(deepcopy(self.clonotypes[c][i]))

        return candidates 
    
    # @staticmethod
    # def find_best_scores(scores):
    #     min_obj = np.Inf
    #     best_trees = []
    #     best_ids = []
    #     for s in scores:
    #         if s.objective < min_obj:
    #             min_obj = s.objective
    #             best_trees = [s.tree]
    #             best_ids = [s.tree.id]
            
    #         elif round(min_obj, 5) == round(s.objective, 5):
    #             best_trees.append(s.tree)
    #             best_ids.append(s.tree.id)
    #         else:
    #             continue
    #     return min_obj, best_trees, best_ids



    def score_candidates(self, candidates, transmat=None, mode=None, nproc=1):
            if mode is None:
                mode = self.mode 
            total_likelihood = 0
      
            best_scores = []
            best_tree_ids = {}
            # observed_data = {}
            all_best_scores = []
            refined_cands = {}
            for i,c in enumerate(self.clonotypes):
                    # print(f'clonotype {c}')
                    ts = TribalSub( isotype_weights=transmat,alpha=self.alpha, nworkers=nproc)
                    if mode=="refine_ilp":
                        all_scores = ts.forest_mode(candidates[c], mode=mode)
                    else:
                        all_scores = ts.forest_mode_loop(candidates[c], mode=mode)
                    sl = ScoreList(all_scores)
                    best_obj, best_scores = sl.find_best_scores()
                    total_likelihood += best_obj
                    # best_score, best_trees, best_ids = self.find_best_scores(all_scores)
                    all_best_scores.append(best_scores)
           
                    # total_likelihood += best_score
                    best_tree_ids[c] = best_scores[0].tree.id
                    # for j,t in enumerate(best_trees):
                    #     tree_name = f"{c}_{j}"
                    #     observed_data[tree_name] = self.obs_states[c]
                    #     t.set_name(tree_name)
                    bs = best_scores[0]
                    lin_tree = bs.tree
                    leaf_isotypes = {l: bs.isotypes[l] for l in lin_tree.get_leafs() }
                    msa = {l: bs.labels[l] for l in lin_tree.get_leafs() }
                    leaf_isotypes[lin_tree.root] = bs.isotypes[lin_tree.root]
                    msa[lin_tree.root] = bs.labels[lin_tree.root]
                    refined_cands[c] = LineageForest(msa, leaf_isotypes, [lin_tree])
               
                    # all_best_trees += best_trees
                    
                    # print(f"{i}: {c} Best tree: {best_tree.id} Tree Score: {best_score} Total Score: {total_likelihood}")
            # lin_forest = LineageForest()
            # lin_forest.generate_from_list(all_best_trees)
            return total_likelihood,all_best_scores, best_tree_ids, refined_cands 
    
    @staticmethod
    def check_convergence(old, new, threshold):
        return np.abs(old -new) < threshold





    # def run(self, nproc=1):

    #     best_log_like_score = np.inf
    
    #     best_trans =None
    #     best_states = None
    #     best_fit_scores = {}
    #     log_like_scores = {}
    #     init_mat_lists = [self.init_transmat.copy()]
    #     # for i in range(self.restarts-1):
    #     #     init_mat_lists.append(tm.add_noise(self.init_transmat.copy(), self.rng, mu=self.mu, sigma=self.sigma))

    #     # results = []
    #     # for t in init_mat_lists:
    #     #     results.append(self.fit(t, nproc))


    #     # with Pool(nproc) as p:
    #     #     results = p.map(self.fit, init_mat_lists)
    #     restart=0
    #     for fit_score, exp_log_like, tmat, state_probs in results:     
    #         print(f"\nFit Phase Complete for restart {restart}!\nFit Score: {fit_score} Exp Log Like {exp_log_like}")
    #         if exp_log_like < best_log_like_score:
    #             best_log_like_score = exp_log_like
    #             best_trans = tmat
    #             best_states = state_probs
    #             best_iter = restart
    #         best_fit_scores[restart] = fit_score
    #         log_like_scores[restart] = exp_log_like
    #         restart+= 1


    #     return best_fit_scores[best_iter], best_trans,best_states, log_like_scores

    
    def fit(self, transmat, nproc=1):
  
     
        ''' resolve the polytomys using a basic sankoff cost function
            
            then generate an estimate of the transition matrix from EM algorithm
            do until convergence of transition matrix and best trees:
                then for each clonotype,
                    run tribal polytomy with the updated transition matrix
                    select all best trees
                refit the transition matrix   
            return a single best tree for each clonotype

        '''
        best_tree_ids = None 
        old_score = np.Inf
        cand_tmat = []
        cand_scores = []
        best_trees = []
        # transmat = self.transmat
        stay_probs = np.linspace(0.55,0.95, self.restarts)
        # for i in range(self.restarts):
        # init_mat_lists = [  tm.gen_trans_mat(stay_probs[i], self.n_isotypes) for i in range(self.restarts)]

        # for i in range(self.niterations):
        for i in range(self.restarts):
            print(f"\nStarting Cycle {i}...")

            transmat =  tm.gen_trans_mat(stay_probs[i], self.n_isotypes)
            
            # candidates = self.intialize_candidates(best_tree_ids)
            candidates = self.intialize_candidates()
            init_score, best_scores, best_tree_ids, refined_cands = self.score_candidates(candidates, transmat=transmat, nproc=nproc)
            for _ in range(self.niterations):
                
                print("\nFitting transition matrix...")
                for j in range(self.niterations):

                   
                    # cur_log_like, state_probs, transmat= EMProbs(lin_forest, transmat, self.states).fit(obs_data)
                    transmat, state_probs = MaxLike(self.n_isotypes).infer(best_scores) 
                    updated_score = 0
                    for scores in best_scores:
                        updated_score += scores[0].compute_score(-1*np.log(transmat))
                    
                    current_score, best_scores, best_tree_ids, refined_cands = self.score_candidates(refined_cands, transmat=transmat, nproc=nproc, mode="score")
                    print(f"iteration: {j} init score: {updated_score} current score: {current_score}")
                    if self.check_convergence(current_score, updated_score, self.threshold):
                        cand_tmat.append((transmat, state_probs))
                        cand_scores.append(current_score)
                        best_trees.append(best_scores)
                        transmat = tm.add_noise(transmat.copy(), self.rng, mu=self.mu, sigma=self.sigma)
                        current_score, best_scores, best_tree_ids, refined_cands = self.score_candidates(refined_cands, transmat=transmat, nproc=nproc, mode="score")

                        break
          
            # DrawStateDiag(transmat).heatmap(f"/scratch/projects/tribal/experimental_data/hoehn_paper/Mouse_2/tribal/refine_ilp/temp/transmat_{i}.png")
            # ts = TribalSub( isotype_weights=transmat,alpha=self.alpha, nworkers=nproc)
            # all_scores = ts.forest_mode(candidates[c], mode=self.mode)
 
            


            # if i > self.min_iterations and self.check_convergence( old_score, updated_score, self.threshold):
            #     break
                   
    
        min_value = min(cand_scores)
        min_index = cand_scores.index(min_value)
        transmat, state_probs = cand_tmat[min_index]
        best_scores = best_trees[min_index]
        best_scores = ScoreList([b[0] for b in best_scores])
        return min_value, transmat, state_probs, best_scores
        # return best_fit_scores[best_iter], best_trans,best_states, log_like_scores
        # return  current_score, updated_score, transmat, state_probs
               

        

    # def search(self, lin_forest):

    #     ''' resolve the polytomys using a basic sankoff cost function
    #         then generate a max likelihood estimate of the transition matrix from EM algorithm
    #         do until convergence of transition matrix and best trees:
    #             then for each clonotype,
    #                 run tribal polytomy search with the updated transition matrix
    #                 select all best trees
    #             refit the transition matrix   
    #         return a single best tree for each clonotype

    #     '''

    #     total_likelihood = 0
    #     best_lin_trees = {}
    #     for i,tree in enumerate(lin_forest.get_trees()):
    #             clono = tree.name
    #             alignment= self.clonotypes[clono].alignment
    #             isotypes = self.clonotypes[clono].isotypes
        
    
    #             best_score, best_lin_tree, best_labels, best_iso = TribalPoly( n_isotypes=self.n_isotypes,transmat=transmat,
    #                                                 alpha=self.alpha).greedy_hill_climbing(tree, alignment, isotypes)
    #             best_lin_trees[clono] = {'tree': best_lin_tree, 'labels': best_labels, 'isotypes':best_iso}
    #             total_likelihood += best_score 
              
                
    #             print(f"{i}: {tree.name}:  Tree Score: {best_score} Total Score: {total_likelihood}")
  
    #     return total_likelihood, best_lin_trees
  



  
      
  

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
    # alignment = ut.read_fasta(align_fname)
    # alignment = {key: list(value.strip()) for key,value in alignment.items()}


    #only treat isotype file as a fasta file if .fasta is in the name, otherwise, we assume it is a csv file dictionary
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
    parser.add_argument("--mu", type=float, default=0.07, help="mean of gaussian white noise to add for distortion")
    parser.add_argument("--sigma", type=float, default=0.05, help="std of gaussian white noise to add for distortion")
    parser.add_argument("--nworkers", type=int, default=2, help="number of workers to use in the event in multiple restarts")
    parser.add_argument("--max_cand", type=int, default = 20,  help="max candidate tree size per clonotype")
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("--restarts",  type=int, default=1, help="number of restarts")
    parser.add_argument("--mode", choices=["score", "refine", "refine_ilp", "search"], default="score")
    # parser.add_argument( "-o", "--outpath", type=str, help="path to directory where output files should be saved")
    parser.add_argument("--score", type=str, help="filename where the score file should be saved")
    parser.add_argument("--transmat_infer", type=str, help="filename where the inferred transition matrix should be saved")
    parser.add_argument("--state_probs", type=str, help="filename where the inferred state probabilities should be saved")
    parser.add_argument("--heatmap", type=str, help="filename where the {png,pdf} of transition matrix should be saved")
    parser.add_argument("--propmap", type=str, help="filename where the {pdf,png} of isotype proportions should be saved")

    # parser.add_argument("--save_all_restarts", type=str, help="path where all restarts should be saved")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
  
    # fpath = "/scratch/projects/tribal/benchmark_pipeline/sim_data/recomb/direct/cells35/size75/rep1/2.0/0.365"
    # args = parser.parse_args(["--clonotypes","/scratch/projects/tribal/benchmark_pipeline/tmats_exp/clonotypes25.txt", 
    #     "--encoding", "/scratch/projects/tribal/benchmark_pipeline/sim_encoding.txt",
    #     "--alpha", "0.8",
    #     "-p", fpath,
    #     "--max_cand", "25",
    #     "--niter" , "20",
    #     "--restarts", "15",
    #     "--root", "naive",
    #     "--tree_path", fpath,
    #     "--nworkers", "7",
    #     "--fasta", "GCsim_dedup.fasta",
    #     "--isotypes", "GCsim.isotypes",
    #     "--candidates", "dnapars/outtree",
    #     "--mode", "refine",
    #     "--heatmap", "test/heatmap_ml_refine.pdf",
    #     "--transmat_inf", "test/transmat_ml_refine.txt"
    #     ])

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
                alpha = args.alpha, 
                seed = args.seed,
                isotype_encoding= iso_encoding,
                max_cand= args.max_cand,
                niter = args.niter,
                threshold=args.thresh,
                not_trans_prob= 1-args.jump_prob,
                restarts=args.restarts,
                mu = args.mu,
                sigma=args.sigma,
                mode = args.mode
                )
    

  
    obj_score, transmat, state_probs,  best_trees= tr.fit(args.nworkers)


    print("\nTRIBAL Complete!, saving results...")

 

    # best_trees.pickle_scores("test/trees.pickle")
    # for i,score in enumerate(best_trees):
    #     score.tree.save_png(f"test/trees/clonotype{i+1}.png", score.isotypes, show_legend=True)


    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write(str(obj_score))

    if args.transmat_infer is not None:
        np.savetxt(args.transmat_infer, transmat)
    if args.state_probs is not None:
        np.savetxt(args.state_probs, state_probs)
    
    if args.heatmap is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).heatmap(args.heatmap)
   
    # if args.propmap is not None:
    #     DrawStateDiag(transmat, state_probs, rev_encoding).state_heatmap(args.propmap)


    # if args.outpath is not None:
    #     lin_forest.save_trees(args.outpath)
 
  

 
    # if args.save_all_restarts is not None:
    #     with open(f"{args.save_all_restarts}/expectation_log_like.csv", "w+") as file:
    #         file.write("restart,max_exp_log_like\n")
    #         for i in all_transmats:
    #             np.savetxt( f"{args.save_all_restarts}/transmat_restart{i}.txt",all_transmats[i])
    #             np.savetxt(f"{args.save_all_restarts}/state_probs_restart{i}.txt",all_state_probs[i])
    #             DrawStateDiag(all_transmats[i], all_state_probs[i], rev_encoding).save(f"{args.save_all_restarts}/state_diagram_restart{i}.png")
    #             file.write(f"{i},{all_log_like[i]}\n")



        
    


  









