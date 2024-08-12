'''
script to run command line interface of tribal
'''

import sys
import argparse
import pandas as pd 
import numpy as np
from tribal import preprocess, Tribal







def parse_tuple(value):
    values = value.split(',')
    return float(values[0]), float(values[1])

def main():
    parser = argparse.ArgumentParser(description="Tribal CLI Tool")
    # parser.add_argument("-f", "--file", type=str)
    # args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    # print(args.file)
    
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands")

    # Preprocess command
    parser_preprocess = subparsers.add_parser("preprocess", help="Preprocess data")
    parser_preprocess.set_defaults(func=preprocessor)
  
    parser_preprocess.add_argument("-d", "--data", type=str,
            help="filename of csv file with the sequencing data", required=True)
    parser_preprocess.add_argument("-r", "--roots", type=str,
            help="filename of csv file with the root sequences", required=True)
    parser_preprocess.add_argument("-e", "--encoding", type=str,
            help="filename isotype encodings", required=True)
    parser_preprocess.add_argument( "--min-size", type=int, default=4,
            help="minimum clonotype size (default 4)")
    parser_preprocess.add_argument("--dataframe",  type=str,
            help="path to where the filtered dataframe with additional sequences and isotype encodings should be saved.")
    parser_preprocess.add_argument("-o", "--out", type=str, required=True,
            help="path to where pickled clonotype dictionary input should be saved")
    parser_preprocess.add_argument("-j", "--cores", type=int, default=1,
            help="number of cores to use (default 1)")
    parser_preprocess.add_argument("--heavy", action="store_true", 
            help= "only use the heavy chain and ignore the light chain")
    parser_preprocess.add_argument("-v", "--verbose", action="store_true",
            help="print additional messages")
 

    

     # tribal command
    parser_tribal = subparsers.add_parser("fit", help="B cell lineage tree inference")
    parser_tribal.set_defaults(func=fit)
    parser_tribal.add_argument("-c", "--clonotypes", type=str, required=True,
                               help="path to pickled clonotypes dictionary of parsimony forests, alignments, and isotypes" )
    parser_tribal.add_argument("--n_isotypes", type=int, required=True, help="the number of isotypes states to use")
    parser_tribal.add_argument("--stay-prob", type=parse_tuple, default=(0.55,0.95), help="the lower and upper bound of not class switching, example: 0.55,0.95")
    parser_tribal.add_argument("-t", "--transmat", required=False, type=str,
        help="optional filename of isotype transition probabilities")
    # parser_tribal.add_argument("-r", "--root", required=False, default="naive",
    #     help="the common id of the root in all clonotypes")
    parser_tribal.add_argument("--niter", type=int, help="max number of iterations during fitting", default=10)
    parser_tribal.add_argument("--thresh", type=float, help="theshold for convergence in during fitting" ,default=0.1)
    parser_tribal.add_argument("-j", "--cores", type=int, default=1, help="number of cores to use")
    parser_tribal.add_argument("--max-cand", type=int, default = 50,  help="max candidate tree size per clonotype")
    parser_tribal.add_argument("-s", "--seed", type=int, default=1026, help="random number seed")
    parser_tribal.add_argument("--restarts",  type=int, default=1, help="number of restarts")
    parser_tribal.add_argument("--mode", choices=["score", "refinement"], default="refinement", 
                               help="mode for fitting B cell lineage trees, one of 'refinment' or 'score' ")
    parser_tribal.add_argument("--score", type=str, help="filename where the objective values file should be saved")
    parser_tribal.add_argument("--transmat-infer", type=str, help="filename where the inferred transition matrix should be saved")
    parser_tribal.add_argument("--verbose", action="store_true",
            help="print additional messages.")
    parser_tribal.add_argument("--pickle",  help="path where the output dictionary of LineageTree lists should be pickled"  )
    parser_tribal.add_argument("--write-results",  help="path where all optimal solution results are saved"  )

 

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

def preprocessor(args):

    with open(args.encoding, "r+") as file:
        isotypes = [ line.strip() for line in file]
           
  
    input_df = pd.read_csv(args.data)
    roots = pd.read_csv(args.roots)
    clonotypes, df = preprocess(input_df, 
                                roots,
                                isotypes, 
                                min_size=args.min_size, 
                                use_light_chain=not args.heavy,
                                cores=args.cores, 
                                verbose=args.verbose)
    
    if args.out is not None:
        pd.to_pickle(clonotypes, args.out)
    
    if args.dataframe is not None:
        df.to_csv(args.dataframe, index=False)


def fit(args):


    clonodict = pd.read_pickle(args.clonotypes)
    

    tr= Tribal( n_isotypes=args.n_isotypes,
                seed = args.seed,
                max_cand= args.max_cand,
                niter = args.niter,
                threshold=args.thresh,
                restarts=args.restarts,
                stay_probs= args.stay_prob,
                verbose = args.verbose
                )
    

    if args.transmat is not None:

        transmat = np.loadtxt(args.transmat)
    else:
        transmat= None


  
    shm_score, csr_likelihood, best_scores, transmat= tr.fit(clonodict, 
                                                             transmat=transmat, 
                                                             mode = args.mode,
                                                             cores=args.cores)


    print("\nTRIBAL Complete!, saving results...")


    if args.score is not None:
        with open(args.score, 'w+') as file:
            file.write(f"shm,csr\n")
            file.write(f"{shm_score},{csr_likelihood}\n")

    if args.transmat_infer is not None:
        np.savetxt(args.transmat_infer, transmat)
   
    if args.pickle is not None:
        pd.to_pickle(best_scores, args.pickle)
    
    if args.write_results is not None:
            _write_results(args.write_results, best_scores, clonodict)

    # if args.all_optimal_solutions is not None:
    #     save_results(args.all_optimal_solutions, best_scores, pngs=True )





def _write_results(outpath, best_trees, clonotypes):

   
    for clono, tree_list in best_trees.items():
        clonotype = clonotypes[clono]
        tree_list.write_all(outpath, clonotype.isotype_encoding)

