import argparse
import numpy as np 
import utils as ut



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", required=False, type=str)
    parser.add_argument("-s", "--sequences", type=str, required=False  )
    parser.add_argument("-i", "--isotypes",  type=str, default= "isotype.fasta",
        help="filename of isotype fasta file within each clonotype directory")
    parser.add_argument( "--transmat", required=False, type=str,
        help="filename of input transition matrix")
    parser.add_argument( "-o", "--output" ,required=False, type=str,
        help="filename of input transition matrix")

    

    args= parser.parse_args()
    # pth = "/scratch/projects/tribal/benchmark_pipeline/sim_data/replicates/cells35/size25/rep1/2.0/0.365/4"
    # tmat_fname ="/scratch/projects/tribal/benchmark_pipeline/sim_data/replicates/cells35/size25/rep1/2.0/0.365/tribal/transmat.txt"
    # parents_fname = f"{pth}/GCsim_collapsed_tree.parents"
    # iso_fname = f"{pth}/GCsim_collapsed_tree.isotypes"
    # seq_fname =  f"{pth}/GCsim_collapsed_tree.sequences"


    transmat = np.loadtxt(args.transmat)
    weights, states = ut.convert_transmat_to_weights(transmat)
    parents = ut.read_dict(args.tree)
    isotypes = ut.read_dict(args.isotypes)
    seqs = ut.read_dict(args.sequences)

    seqs = {key: np.array(list(val)) for key, val in seqs.items()}
    isotypes = {key: int(val) for key, val in isotypes.items()}

    edge_list = []
    ancestral_nodes = []
    for child, val in parents.items():
        if val != "":
            edge_list.append((val, child))
        if "seq" in val and val not in ancestral_nodes:
            ancestral_nodes.append(val)
    
    # for a in ancestral_nodes:
    #     edge_list.append((a,a))
        

    

    pars_score = 0
    iso_score = 0
    for u,v in edge_list:
        pars_score += ut.hamming_distance(seqs[u],seqs[v])
        
        iso_score += weights[isotypes[u], isotypes[v]]
    
    # print(f"Parsimony: {pars_score} Isotype Score: {iso_score}")
    with open(args.output, "w+" ) as file:
        file.write(f"{pars_score},{iso_score}\n")
    # print("Done!")
    

    
    
    
