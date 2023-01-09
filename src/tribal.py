
import networkx as nx
import numpy as np
import sys, os, re
import argparse 
import pickle
from multiprocessing import Pool
from copy import deepcopy
import utils as ut
from tree_utils import TribalTree
from ete3 import Tree
from trans_matrix import TransMat
from em_weight_matrix import EMProbs
from lineage_tree import LineageTree, LineageForest
import init_transmat as tm

from tribal_sub import TribalSub
from draw_state_diagram import DrawStateDiag

class Tribal:
    def __init__(self, clonotypes, transmat=None, alpha=0.9, 
                alphabet= ("A", "C", "G", "T","N", "-"), 
                isotype_encoding=None, seed= 1026, 
                max_cand=50, niter=10,
                threshold=0.1, restarts=5,
                n_isotypes = 7, not_trans_prob=0.65, mu=0.07, sigma=0.05 ):
        
        self.clonotypes = clonotypes
        self.isotype_encoding = isotype_encoding
        self.alphabet = alphabet
        self.alpha = alpha
        self.not_trans_prob = not_trans_prob
        print(self.not_trans_prob)
        if transmat is None:
            if isotype_encoding is not None:
                self.n_isotypes = len(isotype_encoding)
            else:
                self.n_isotypes = n_isotypes 
            print(f"generating transitiom matrix with not jump prob {self.not_trans_prob}")
            self.transmat = tm.gen_trans_mat(self.not_trans_prob, self.n_isotypes)
      
        else:
            self.transmat = transmat
        
        self.init_transmat = self.transmat.copy()

    
        self.states = np.arange(self.transmat.shape[0])

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
        candidates = {}
        for c in self.clonotypes:
            if self.clonotypes[c].size() > self.max_cand:
                cand = self.rng.choice(self.clonotypes[c].size(), self.max_cand, replace=False).tolist()
            else:
                cand = [i for i in range(self.clonotypes[c].size())]
            
            candidates[c]= LineageForest( alignment=self.clonotypes[c].alignment, isotypes=self.clonotypes[c].isotypes)
            for index in cand:
                candidates[c].add(deepcopy(self.clonotypes[c][index]))
            if best_tree_ids is not None:
                if best_tree_ids[c] not in cand:
                    candidates[c].add(deepcopy(self.clonotypes[c][best_tree_ids[c]]))

        return candidates 
    

    def score_candidates(self, candidates, transmat=None):
          
            total_likelihood = 0
            best_trees = []
            best_tree_ids = {}
            for i,c in enumerate(self.clonotypes):
                    ts = TribalSub( isotype_weights=transmat,alpha=self.alpha)
                    best_score,_ = ts.forest_mode_loop(candidates[c], mode="refine")
                    best_score = best_score[0]
                    best_tree = best_score.tree
                    total_likelihood += best_score.objective 
                    best_tree_ids[c] = best_tree.id
                    best_tree.set_name(c)
                    best_trees.append(best_tree)
                    
                    # print(f"{i}: {c} Best tree: {best_tree.id} Tree Score: {best_score} Total Score: {total_likelihood}")
            lin_forest = LineageForest()
            lin_forest.generate_from_list(best_trees)
            return total_likelihood, best_tree_ids, lin_forest
    
    @staticmethod
    def check_convergence(old, new, threshold):
        return np.abs(old -new) < threshold


    # def run(self):
    #     best_tree_id_pointer = None
    #     best_log_like_score = np.NINF
    #     best_lin_forest = None
    #     all_transmats ={}
    #     all_state_probs = {}
    #     best_trans =None
    #     best_states = None
    #     best_fit_scores = {}
    #     log_like_scores = {}
    #     for i in range(self.restarts):

    #         print(f"\nStarting fitting phase for restart {i}.....")
    #         fit_score, best_tree_ids, lin_forest, exp_log_like = self.fit()
    #         print(f"\nFit Phase Complete for restart {i}!\nFit Score: {fit_score} Exp Log Like {exp_log_like}\n")
    #         if exp_log_like > best_log_like_score:
    #             best_log_like_score = exp_log_like

    #             best_trans = self.transmat.copy()
    #             best_states = self.state_probs.copy()
    #             best_lin_forest = lin_forest
    #             best_tree_id_pointer = best_tree_ids 
    #             best_iter = i
    #         best_fit_scores[i] = fit_score
    #         log_like_scores[i] = exp_log_like

    #         all_transmats[i] =self.transmat.copy()
    #         all_state_probs[i] = self.state_probs.copy()
    #         self.transmat = self.init_transmat
    #         self.transmat =tm.add_noise(self.transmat, self.rng, mu=0.05, sigma=0.01, min_prob=0.01)


            


        # search_score, best_lin_trees =self.search(lin_forest)
        # print(f"\nSearch Phase Complete!\nSearch Score: {search_score}")

        # return best_fit_scores[best_iter], best_lin_forest, best_trans,best_states, all_transmats, all_state_probs, log_like_scores


    def run(self, nproc=1):

        best_log_like_score = np.NINF
        best_lin_forest = None
        all_transmats ={}
        all_state_probs = {}
        best_trans =None
        best_states = None
        best_fit_scores = {}
        log_like_scores = {}
        init_mat_lists = [self.init_transmat.copy()]
        for i in range(self.restarts-1):
            init_mat_lists.append(tm.add_noise(self.init_transmat.copy(), self.rng, mu=self.mu, sigma=self.sigma, min_prob=0.01))

        with Pool(nproc) as p:
            results = p.map(self.fit, init_mat_lists)
        restart=0
        for fit_score, best_tree_ids, lin_forest, exp_log_like, tmat, state_probs in results:     
            print(f"\nFit Phase Complete for restart {restart}!\nFit Score: {fit_score} Exp Log Like {exp_log_like}")
            if exp_log_like > best_log_like_score:
                best_log_like_score = exp_log_like

                best_trans = tmat
                best_states = state_probs
                best_lin_forest = lin_forest
                best_tree_id_pointer = best_tree_ids 
                best_iter = restart
            best_fit_scores[restart] = fit_score
            log_like_scores[restart] = exp_log_like

            all_transmats[restart] = tmat
            all_state_probs[restart] = state_probs
            restart+= 1



            


        # search_score, best_lin_trees =self.search(lin_forest)
        # print(f"\nSearch Phase Complete!\nSearch Score: {search_score}")

        return best_fit_scores[best_iter], best_lin_forest, best_trans,best_states, all_transmats, all_state_probs, log_like_scores

    
    def fit(self, transmat):
  
     
        ''' resolve the polytomys using a basic sankoff cost function
            
            then generate a max likelihood estimate of the transition matrix from EM algorithm
            do until convergence of transition matrix and best trees:
                then for each clonotype,
                    run tribal polytomy with the updated transition matrix
                    select all best trees
                refit the transition matrix   
            return a single best tree for each clonotype

        '''
        best_tree_ids = None 
        old_score = np.Inf
        for i in range(self.niterations):
            print(f"\nStarting Cycle {i}...")
            # if i== 0:
            #     transmat = None
            # else:
            #     transmat = self.transmat

            
            candidates = self.intialize_candidates(best_tree_ids)
            current_score, best_tree_ids, lin_forest = self.score_candidates(candidates, transmat=transmat)

            print("\nFitting transition matrix...")
            cur_log_like, state_probs, transmat= EMProbs(lin_forest, transmat, self.states).fit(self.obs_states)
            # self.transmat = new_tmat

         
            print(f"\nCycle {i} Complete: Current Score: {current_score} Previous Score: {old_score} Curr EM Log Like: {cur_log_like}")
            
            if i > min(7, self.niterations):

                if self.check_convergence( old_score, current_score, self.threshold):
                  break
                   
            old_score = current_score 
        return  current_score, best_tree_ids, lin_forest, cur_log_like, transmat, state_probs
               

        

    def search(self, lin_forest):

        ''' resolve the polytomys using a basic sankoff cost function
            then generate a max likelihood estimate of the transition matrix from EM algorithm
            do until convergence of transition matrix and best trees:
                then for each clonotype,
                    run tribal polytomy search with the updated transition matrix
                    select all best trees
                refit the transition matrix   
            return a single best tree for each clonotype

        '''

        total_likelihood = 0
        best_lin_trees = {}
        for i,tree in enumerate(lin_forest.get_trees()):
                clono = tree.name
                alignment= self.clonotypes[clono].alignment
                isotypes = self.clonotypes[clono].isotypes
        
    
                best_score, best_lin_tree, best_labels, best_iso = TribalPoly( n_isotypes=self.n_isotypes,transmat=transmat,
                                                    alpha=self.alpha).greedy_hill_climbing(tree, alignment, isotypes)
                best_lin_trees[clono] = {'tree': best_lin_tree, 'labels': best_labels, 'isotypes':best_iso}
                total_likelihood += best_score 
              
                
                print(f"{i}: {tree.name}:  Tree Score: {best_score} Total Score: {total_likelihood}")
  
        return total_likelihood, best_lin_trees
  



  
      
  

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


def create_input( path,  tree_path, clonotype, root, seq_fasta_fname, trees_fname, iso_fasta_fname, iso_encoding=None, start_iso=None):

    tree_fname =f"{tree_path}/{clonotype}/{trees_fname}"
    align_fname = f"{path}/{clonotype}/{seq_fasta_fname}"
    iso_fname =f"{path}/{clonotype}/{iso_fasta_fname}"
    tree_list = create_trees(tree_fname)

     

    alignment = ut.read_fasta(align_fname)
    alignment = {key: list(value.strip()) for key,value in alignment.items()}


    #only treat isotype file as a fasta file if .fasta is in the name, otherwise, we assme it is a csv file dictionary
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
            # ete_tree.name = args.root
            # print(ete_tree)
            nx_tree= convert_to_nx(ete_tree, args.root)
            # print(list(nx_tree.edges))
            # ttree = TribalTree(nx_tree, root=args.root, is_rooted=True)
            cand_trees.append(nx_tree)
        return cand_trees


def pickle_save(obj, fname):
        with open(fname, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--clonotypes", required=False, type=str,
        help="filename with list of clonotype ids")
    parser.add_argument("-p", "--path", type=str, required=True  )
    parser.add_argument("-i", "--isotypes",  type=str, default= "isotype.fasta",
        help="filename of isotype fasta file within each clonotype directory")
    parser.add_argument("-t", "--transmat", required=False, type=str,
        help="filename of input transition matrix")
    parser.add_argument("-r", "--root", required=False, default="naive",
        help="the common id of the root in all clonotypes")
    parser.add_argument(  "--tree_path", type=str, required=True, help="path to directory where candidate trees are saved")

    parser.add_argument("--candidates", type=str, default="outtree", help="filename containing newick strings for candidate trees")
    parser.add_argument("--niter", type=int, help="max number of iterations in the fitting phase", default=10)
    parser.add_argument("--thresh", type=float, help="theshold for convergence in fitting phase" ,default=0.1)
    parser.add_argument("--mu", type=float, default=0.07, help="mean of gaussian white noise to add for distortion")
    parser.add_argument("--sigma", type=float, default=0.05, help="std of gaussian white noise to add for distortion")
    parser.add_argument("--nworkers", type=int, default=1, help="number of workers to use in the event in multiple restarts")

    parser.add_argument("--max_cand", type=int, default = 10,  help="max candidate tree size per clonotype")
    parser.add_argument("-s", "--seed", type=int, default=1026)
    parser.add_argument("-e", "--encoding", type=str, help="text file isotype states listed in germline order")
    parser.add_argument("--n_isotypes", type=int, default=7, help="the number of isotypes states to use if isotype encoding file is not providied")

    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("-j", "--jump_prob", type=float, default=0.25, help="for inititalization of transition matrix if not provided")
    parser.add_argument("--restarts",  type=int, default=1, help="number of restarts")
    parser.add_argument("--trees_fname", type=str, )
    
    
    
    parser.add_argument( "--fasta", type=str, default= "concat.aln.fasta", help="filename of input MSA in fasta file")
    parser
    # parser.add_argument("--mode", type=str, choices=["fit", "search"],
    #             help="fit only resolves polytomyies and does not perform tree moves, search will search tree space for the candidate set")
    
    
    parser.add_argument( "-o", "--outpath", type=str, help="path to directory where output files should be saved")
    parser.add_argument("--score", type=str, help="filename where the score file should be saved")
    parser.add_argument("--transmat_infer", type=str, help="filename where the inferred transition matrix should be saved")
    parser.add_argument("--state_probs", type=str, help="filename where the inferred state probabilities should be saved")
    parser.add_argument("--diagram", type=str, help="filename where the png of transition matrix should be saved")
    parser.add_argument("--diagram_pdf", type=str, help="filename where the pdf of transition matrix should be saved")
    parser.add_argument("--heatmap", type=str, help="filename where the heatmap pdf of transition matrix should be saved")
    parser.add_argument("--save_all_restarts", type=str, help="path where all restarts should be saved")

    # parser.add_argument( "--sequences", type=bool, action="store_true", help="if ancestral sequences should be saved")
    # parser.add_argument("-n", "--newick", type=bool, action="store_true",  help="if newick string should be saved")
    # parser.add_argument("--score",  type=bool, action="store_true")
    # parser.add_argument("--iso_infer", type=bool, action="store_true")

    # parser.add_argument("--save_candidates", type=bool, action="store_true")
    args= parser.parse_args()

    print(f"jump prob: {args.jump_prob}")
    # path = "/scratch/projects/tribal/real_data"
    # dataset = "GCB_NP_2"
 

    # args =parser.parse_args([
    #     "-c", f"{path}/{dataset}/clonotypes.txt",
    #     "-p", f"{path}/{dataset}/input",
    #     "-r", "naive",
    #     "-j", "0.35",
    #     "-e", f"{path}/mouse_isotype_encoding.txt",
    #     "-s", "3",
    #     # "--score", f"{path}/{dataset}/tribal/score.txt",
    #     "-o", f"{path}/test",
    #     "--alpha", "0.75",
    #     "--max_cand", "5",
    #     "--niter", "3",
    #     "--thresh", "0.1",
    #     "--tree_path", f"{path}/{dataset}/dnapars",
    #     "--save_all_restarts", f"{path}/test",
    #     "--diagram", f"{path}/test/state_diagram.png",
    #     "--restarts", "3",
    #     "--nworkers", "3"

    # ])

    if args.encoding is not None:
        iso_encoding, start_iso, n_isotypes = create_isotype_encoding(args.encoding)
        rev_encoding = {val: key for key, val in iso_encoding.items()}
    else:
        n_isotypes = args.n_isotypes
        iso_encoding = None
        start_iso= None 
        rev_encoding = None
    
    
    if args.transmat is not None:
    #     transmat= TransMat(n_isotypes=n_isotypes).fit(args.jump_prob)
    # else:
        transmat = np.loadtxt(args.transmat)
    else:
        transmat= None

    # if n_isotypes != transmat.shape[0]:
    #     raise ValueError("Isotype states in transition matrix does not match the number of states in encoding file")
    
    
    if args.clonotypes is not None:
        clonotypes = []
        with open(args.clonotypes, 'r+') as file:
            for line in file:
                clonotypes.append(line.strip())

    else:
         clonotypes = [it.name for it in os.scandir(args.path) if it.is_dir()]

    clonodict = {}
    for c in clonotypes:
        print(args.isotypes)
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
                sigma=args.sigma
                )
    

  
    # obj_score, lin_forest, transmat, state_probs, all_transmats, all_state_probs, all_log_like = tr.run()
    obj_score, lin_forest, transmat, state_probs, all_transmats, all_state_probs, all_log_like = tr.run(args.nworkers)

    # print(f"Tribal Object Score: {obj_score}")
    print("\nTRIBAL Complete!, saving results...")

    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write(str(obj_score))

    if args.transmat_infer is not None:
        np.savetxt(args.transmat_infer, transmat)
    if args.state_probs is not None:
        np.savetxt(args.state_probs, state_probs)
    
    if args.diagram is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).save(args.diagram)
    if args.diagram_pdf is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).save_pdf(args.diagram_pdf)
    if args.heatmap is not None:
            DrawStateDiag(transmat, state_probs, rev_encoding).heatmap(args.heatmap)


    if args.outpath is not None:
        lin_forest.save_trees(args.outpath)
        # save_results(args.outpath, lin_tree_dict, pngs=True, isotype_mapping=rev_encoding)
  

 
    if args.save_all_restarts is not None:
        with open(f"{args.save_all_restarts}/expectation_log_like.csv", "w+") as file:
            file.write("restart,max_exp_log_like\n")
            for i in all_transmats:
                np.savetxt( f"{args.save_all_restarts}/transmat_restart{i}.txt",all_transmats[i])
                np.savetxt(f"{args.save_all_restarts}/state_probs_restart{i}.txt",all_state_probs[i])
                DrawStateDiag(all_transmats[i], all_state_probs[i], rev_encoding).save(f"{args.save_all_restarts}/state_diagram_restart{i}.png")
                file.write(f"{i},{all_log_like[i]}\n")



        
    


  









