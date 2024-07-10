
import numpy as np
from utils import pickle_load
from draw_state_diagram import DrawStateDiag
import argparse, sys
class MaxLike:
    def __init__(self, n_isotypes=7, pseudo_count =1) -> None:
        
        self.n_iso = n_isotypes
        self.trans_probs = np.full((self.n_iso, self.n_iso), pseudo_count)
        self.trans_probs =np.triu(self.trans_probs)
        self.state_probs = np.ones(self.n_iso)

        
    def infer(self, score_list):
        '''
        given score list object (inluding trees and isotype labels)
        infer the ML isotype transitions and proportions.
        '''
        for scores in score_list:
            for s in scores:
                tree = s.tree
                iso_labs = s.isotypes
                self.update_counts(tree, iso_labs)

        
        #add pseudo counts of unobserved transitions in the data

                
            # if self.state_probs[i] ==0:
            #     self.state_probs[i] = self.pseudo_counts 
            # for j in range(self.trans_probs.shape[1]):
            #     if j >= i and self.trans_probs[i,j] ==0:

            #         self.trans_probs[i,j]= self.pseudo_counts
        row_sums = self.trans_probs.sum(axis=1)
    
        self.trans_probs = self.trans_probs/ row_sums[:, np.newaxis]
        self.state_probs = self.state_probs/self.state_probs.sum()

        return self.trans_probs, self.state_probs

    

    def update_counts(self, tree, isotypes):

        nodes  = tree.preorder_traversal()
        for n in nodes:
     
            from_state = isotypes[n]
            self.state_probs[from_state] = self.state_probs[from_state] + 1
            children = tree.children(n)
            if len(children) > 0:
                for c in children:
                    to_state = isotypes[c]
                    self.trans_probs[from_state, to_state] =self.trans_probs[from_state, to_state] + 1
        
      

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--path", type=str, required=True, help="path to the directory containing input files")
    parser.add_argument("-k", "--clonotypes", required=True, type=int, help="number of clonotypes in the experiment")
    parser.add_argument("-e", "--encoding", type=str, help="text file isotype states listed in germline order")
    parser.add_argument("--pseudo", type=float, help="pseudo count value", default=0.001)

    parser.add_argument("--n_isotypes", type=int, default=7, help="the number of isotypes states to use if isotype encoding file is not provided and input isotypes are encoded numerically")
    parser.add_argument("--transmat_infer", type=str, help="filename where the inferred transition matrix should be saved")
    parser.add_argument("--state_probs", type=str, help="filename where the inferred state probabilities should be saved")
    parser.add_argument("--heatmap", type=str, help="filename where the {png,pdf} of transition matrix should be saved")
    parser.add_argument("--propmap", type=str, help="filename where the {pdf,png} of isotype proportions should be saved")
    parser.add_argument("--method", choices = ["dnapars", "dnaml", "igphyml"])

    # parser.add_argument("--save_all_restarts", type=str, help="path where all restarts should be saved")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    # args = parser.parse_args([
    #     "-p", "/scratch/projects/tribal/benchmark_pipeline/sim_data/tmat_inf/direct/cells35/size25/rep4/2.0/0.365",
    #     "-e", "/scratch/projects/tribal/benchmark_pipeline/sim_encoding.txt",
    #     "-k", "25",
    #     "--transmat_infer", "/scratch/projects/tribal/test/transmat.txt",
    #     "--state_probs", "/scratch/projects/tribal/test/state_probs.txt",
    #     "--heatmap", "/scratch/projects/tribal/test/heatmap.png",
    # ])

    if args.encoding is not None:
        iso_encoding, start_iso, n_isotypes = create_isotype_encoding(args.encoding)
        rev_encoding = {val: key for key, val in iso_encoding.items()}
    else:
        n_isotypes = args.n_isotypes
        iso_encoding = None
        start_iso= None 
        rev_encoding = None
    

    path = args.path

    clonotypes = [i+1 for i in range(args.clonotypes)]
    # score_list = []
    # for c in clonotypes:
    #     scores = pickle_load(f"{path}/{c}/dnapars/best_results.pickle")
    #     score_list.append(scores[0])

    score_list = [pickle_load(f"{path}/{c}/{args.method}/best_results.pickle") for c in clonotypes]
    for i,s in enumerate(score_list):
        print(f"clonotype {i+1} : {len(s)}")

    ml = MaxLike(n_isotypes=n_isotypes, pseudo_count=args.pseudo)
    transmat, state_probs = ml.infer(score_list)



    

    if args.transmat_infer is not None:
        np.savetxt(args.transmat_infer, transmat)
    if args.state_probs is not None:
        np.savetxt(args.state_probs, state_probs)
    
    if args.heatmap is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).heatmap(args.heatmap)
   
    if args.propmap is not None:
        DrawStateDiag(transmat, state_probs, rev_encoding).state_heatmap(args.propmap)


    print("done")





