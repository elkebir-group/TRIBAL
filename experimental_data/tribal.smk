configfile: "tribal.yml"
import pandas as pd 
import sys
sys.append("../src")

def get_files(dirn, modes, fname):
    targets = []
    for d in config["dataset"]:
        df =pd.read_csv(f"{d}/igphyml/name_isotype_mapping.csv")
        clonotypes = df['clone_id'].unique()
        for c in clonotypes:
            for m in modes:
                for s in config["script"]:
                    targets.append(f"{d}/{dirn}/{s}/{m}/{c}/{fname}")
    return targets 


rule all:
    input: 
        expand( "{dataset}/tribal_recomb/input_forest.pickle",
            dataset = config["dataset"]
        ),
        expand("{dataset}/tribal_recomb/{script}/transmat.txt",
            dataset = config["dataset"],
            script = config["script"]
        ),
        get_files("tribal_recomb", ["score", "refine_ilp"], "tree.txt")

rule prep_dnapars:
    input: 
        clonotypes="{dataset}/clonotypes_igphyml.txt",
        encoding = "mouse_isotype_encoding.txt",
    params:
        fasta = "heavy.aln.fasta",
        tree_path = "./{dataset}",
        inpath = "./{dataset}/recomb_input",
        root = "naive",
    output: "{dataset}/tribal_recomb/input_forest.pickle"
    shell:
        "python ../src/newicks_to_lineageforest.py --clonotypes {input.clonotypes} "
        "--encoding {input.encoding} -p {params.inpath} "
        "--tree_path {params.tree_path} -r {params.root} "
        "--fasta {params.fasta} -o {output} "



        
        # expand(  "{dataset}/{mode}/{clonotype}/tree.png",
        #     dataset=config["dataset"],
        #     mode =  config["refine_modes"],
        #     clonotype = get_clonotypes(),
        # ),
              

rule get_clonotypes:
    input: "{dataset}/igphyml/name_isotype_mapping.csv"
    output: "{dataset}/clonotypes_igphyml.txt"
    run:
        df =pd.read_csv(input[0])
        clonotypes = df['clone_id'].unique()
        with open(output[0], "w+") as file:
            for c in clonotypes:
                file.write(f"{c}\n")

rule tribal:
    input: 
        forest = "{dataset}/tribal_recomb/input_forest.pickle",
        encoding = "mouse_isotype_encoding.txt",
    params:
        max_cand = config["max_cand"],
        niter = config["niter"],
        thresh = config["threshold"],
        seed =  config["seed"],
        root = "naive",
        restarts = config["nrestarts"],
        mode = "refine_ilp",
        script = lambda wildcards: config['script'][wildcards.script]
        # inpath = "./{dataset}/recomb_input",
        # tree_path = "./{dataset}/dnapars",
        # fasta= "heavy.aln.fasta"
    threads: config["nworkers"]
    output:
        transmat = "{dataset}/tribal_recomb/{script}/transmat.txt",
        stateprobs = "{dataset}/tribal_recomb/{script}/state_probs.txt",
        heatmap = "{dataset}/tribal_recomb/{script}/transmat.pdf",
        propmap = "{dataset}/tribal_recomb/{script}/state_probs.pdf",
        score = "{dataset}/tribal_recomb/{script}/fit_score.csv"
    log:
        std = "{dataset}/tribal_recomb/{script}/fit.log",
        err = "{dataset}/tribal_recomb/{script}/fit.err.log"
    benchmark: "{dataset}/tribal_recomb/{script}/benchmark.log"
    shell:
        "nice -n 5 python {params.script} --forest {input.forest} "
        "-r {params.root} -e {input.encoding} -s {params.seed}  "
        " --niter {params.niter} --mode {params.mode} "
        "--thresh {params.thresh} --transmat_infer {output.transmat}  "
        "--state_probs {output.stateprobs} --score {output.score}  "  
        "--nworkers {threads} --restarts {params.restarts} "
        "--heatmap {output.heatmap} --propmap {output.propmap} > {log.std} 2> {log.err} "    


rule tribal_refine:
    input: 
        forest = "{dataset}/tribal_recomb/input_forest.pickle",
        encoding = "mouse_isotype_encoding.txt",
        transmat =  "{dataset}/tribal_recomb/{script}/transmat.txt",
    output: 
        forest="{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/forest.pickle",
        tree= "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/tree.txt",
        seq= "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/inferred_seq.csv",
        fasta= "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/inferred_seq.fasta",
        isotypes= "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/isotypes.csv",
        score= "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/scores.csv",
        png = "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/tree.png",
        opts = directory("{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/opt_trees"),
        rf_dist="{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/best_rf_dist.csv",
        node_degree ="{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/parsimony_node_degree.csv",
        refine_degree ="{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/refine_node_degree.csv",
    params:
        root = "naive",
        # mode = "refine_ilp"
    threads: config['nworkers']
    log: 
        run= "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/refine.log",
        err ="{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/refine.err.log"  
    benchmark: "{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/refine.benchmark.log" 
    shell:
        "nice -n 5 python ../src/tribal_tree.py "
        " -r {params.root} "
        " -f {input.forest} "
        "-c {wildcards.clonotype} "
        "-t {input.transmat} "
        "-e {input.encoding} "
        "--sequences {output.seq} "
        "--fasta {output.fasta} "
        "--score {output.score} "
        "--iso_infer {output.isotypes} "
        "--png {output.png} "
        "--mode {wildcards.mode} "
        "--all_optimal_sol {output.opts} "
        "--best_tree_diff {output.rf_dist} "
        "--nworkers {threads} "
        "--tree {output.tree} "
        "--pars_degree {output.node_degree} "
        "--refine_degree {output.refine_degree} "
        "-o {output.forest} > {log.run} 2> {log.err} "  

   
rule newick_strings:
    input:    
        scores="{dataset}/tribal_recomb/{script}/{mode}/{clonotype}/forest.pickle",
        mapping =  "{dataset}/igphyml/seq_mappings/{clonotype}.mapping.csv",
     output:
        newicks = "{dataset}/tribal_recomb/{script}/{mode}/newick/{clonotype}.nwk.csv"
      run:
          import utils as ut 
          import lineage_tree as lt 
          seq_mapping = ut.read_dict(input.mapping)
          rev_mapping = {val: key for key,val in seq_mapping.items()}
          scores  = lt.load(input.forest)
          for score in scores:
            id = score.lin_tree.id 



           
    






