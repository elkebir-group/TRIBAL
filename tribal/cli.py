######### PREPROCESS#########

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
    return iso_encoding

def main(args):
    df  =pd.read_csv(args.file)
    roots =pd.read_csv(args.roots)
    iso_encoding = create_isotype_encoding(args.encoding)
    clonodict, df = preprocess(df, roots, isotype_encoding = iso_encoding,
                                 min_size=args.min_size,
                                 verbose = args.verbose, heavy=args.heavy)
    if args.pickle is not None:
        pd.to_pickle(clonodict, args.pickle)
    
    if args.dataframe is not None:
        df.to_csv(args.dataframe, index=False)

    if args.clonotypes is not None:
        with open(args.clonotypes, "w+") as file:
            for j in df["Clonotype"].unique():
                file.write(j + "\n")





if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file", type=str,
        help="filename of csv file with the clone data")
    parser.add_argument("-r", "--roots", type=str,
        help="filename of csv file with the root sequences")
    parser.add_argument("-e", "--encoding", type=str,
        help="filename isotype encodings")
    parser.add_argument( "--min-size", type=int, default=5,
        help="minimum clonotype size")
    parser.add_argument("--dataframe",  type=str,
        help="path to where the filtered dataframe with additional sequences and isotype encodings should be saved.")
    parser.add_argument("--clonotypes",  type=str,
        help="path to where a list of the included clontoypes should be saved.")
    parser.add_argument("--verbose", action="store_true",
        help="print additional messages.")
    parser.add_argument("-P", "--pickle", type=str,
        help="path to where pickled clonotype dictionary input should be saved")
    parser.add_argument("-j", "--cores", type=int, default=1,
        help="number of cores to use")
    parser.add_argument("--heavy", action="store_true", 
        help= "only use the heavy chain and ignore the light chain")
 


    args= parser.parse_args()
    # args = parser.parse_args([
    #     "-f", "Human2/human_data.csv",
    #     "-r", "Human2/human_data_root_seq.csv",
    #     "-e", "Human2/human_encoding.txt",
    #     # "-P", "Human/clonotypes.pkl",
    #     "--verbose"

    # ])

    main(args)




########## TRIBAL cli
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





def expand_tree(T, leafset):
   

    remap = {n: f"{n}_int" for n in T if n in leafset and T.out_degree[n] > 0}
    T = nx.relabel_nodes(T, remap)
    for child, parent, in remap.items():
        T.add_edge(parent, child) 
    return T     
  
def create_trees(cand_fname, format=0):
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
            
            try:
                ete_tree = Tree(nw, format=format)
            except:
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
    parser.add_argument("--edge-list", type=str, help="filename containing edge list of input tree")
    parser.add_argument("--mode", choices=["score", "refine", "refine_ilp", "search"], default="score")
    parser.add_argument("-e", "--encoding", type=str, required=True)
    parser.add_argument("--tree",  type=str, help="outputfile of best tree")
    parser.add_argument( "--fasta", type=str, help="filename where reconstructed ancestral sequences should be saved as fasta file")
    parser.add_argument( "--png", type=str, help="filename where to save a png of the optimal tree")
    parser.add_argument( "--sequences", type=str, help="filename where reconstructed ancestral sequences should be saved as csv file")
    parser.add_argument("--score",  type=str, help="filename of the objective function value objective function value")
    parser.add_argument("--reversible",  action="store_true", 
                        help="a flag to indicate the standard 0/1 cost function is used (the number of isotype changes is minimized and irreversiblility is ignored)")
    parser.add_argument("--compute-seq",  action="store_true", 
                        help="if the small parsimony problem should be solved for sequences")
    parser.add_argument("--iso_infer",  type=str, help="filename of the inferred isotypes for the internal nodes")
    parser.add_argument("--all_optimal_sol",  help="path where all optimal solution results are saved"  )
    parser.add_argument("--nworkers", type=int, default=1, help="number of workers to use in the event of multiple input candidate trees")
    parser.add_argument("--pickle_best", type=str, help="filename to pickle the best results")
    parser.add_argument("--pickle_all", type=str, help="filename to pickle the best results")
    parser.add_argument("--nwk-format", type=int, default=0, help="ete3 newick format for input newick tree file, default=0")


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    # ttype = "direct"
    # n = 65
    # r = 1
    # clonotype = 63
    # dist = "clonaltree_a0"
    # inpath = f"benchmark_pipeline/sim_data/recomb/{ttype}/cells{n}/size75/rep{r}/2.0/0.365/{clonotype}"

    # args = parser.parse_args([
    #     "-r", "naive",
    #     "--nwk-format" , "1",
    #     "-i", f"{inpath}/GCsim.isotypes",
    #     "-a", f"{inpath}/GCsim_dedup.fasta",
    #     "--edge-list", f"{inpath}/{dist}/tree.nwk.csv",
    #     "-e", "benchmark_pipeline/sim_encoding.txt",
    #     "--fasta", f"{inpath}/{dist}/seq.fasta",
    #     "--iso_infer", f"{inpath}/{dist}/isotype.csv",
    #     "--pickle_all", f"{inpath}/{dist}/all_results.pkl",
    #     "--compute-seq"
    # ]

    # )
  
    ###CLI for TRibalSub
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
        cand = create_trees(args.candidates, args.nwk_format)
        lin_forest = LineageForest(alignment=alignment, isotypes=isotypes_filt)
        lin_forest.generate_from_list(cand, args.root)
      
    elif args.edge_list is not None:
        edge_list = []
        with open(args.edge_list, "r+" ) as file:
            for line in file:
                edge = line.strip().split(",")
                edge_list.append((edge[0],edge[1]))
        T = expand_tree(nx.DiGraph(edge_list), set(alignment.keys()) - set([args.root]))
        lin_forest = LineageForest(alignment=alignment, isotypes=isotypes_filt)
        lin_forest.generate_from_list([T], args.root)


        

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

    tr = TribalSub(isotype_weights, 0.9, timeout=args.timeout, 
                   nworkers=args.nworkers, root_id=args.root,
                     reversible=args.reversible, compute_seq = args.compute_seq)

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
    