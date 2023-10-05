configfile: "igphyml_config.yml"
import pandas as pd 


def get_files():
    targets = []
    for d in config["datasets"]:
        df =pd.read_csv(f"{d}/igphyml/name_isotype_mapping.csv")
        clonotypes = df['clone_id'].unique()
        for clone in clonotypes:
                    targets.append(f"{d}/igphyml/pngs/{clone}.png")
    return targets 

rule all:
   input:  
        expand("{dataset}/igphyml/data_lineages_gy.tsv_igphyml_stats.txt",
            dataset = config['datasets']
        ),
        expand("{dataset}/igphyml/likelihoods.csv",
            dataset = config['datasets']
        ),
        get_files()


rule clean_airr:
    input:
        metadata = "{dataset}/{dataset}_dandelion_metadata.tsv",
        table = "{dataset}/{dataset}_dandelion_table.tsv",
        clonotypes = "{dataset}/clonotypes.txt"

    output:
        indata= "{dataset}/igphyml/input_data.tsv",
        mapping= "{dataset}/igphyml/name_isotype_mapping.csv",
    shell:
        "python preprocess_airr.py "
        "-m {input.metadeta} -t {input.table} "
        " -c {input.clonotypes} -o {output.indata} "
        "-p {output.mapping} "


# snakemake --use-conda  -s igphyml_pipeline.smk -j 1
rule run_igphyml:
    input: 
        airr = "{dataset}/igphyml/input_data.tsv"
    output:
       clones= "{dataset}/igphyml/clones.rds",
       trees = "{dataset}/igphyml/trees.rds",
    conda: "r_dowser"
    params:
        nproc = 7
    script:
        'igphyml_run.R'

rule write_trees:
    input: 
        trees = "{dataset}/igphyml/trees.rds",
    output:
       all_files = "{dataset}/igphyml/all_files.rds"
    params:     
        outdir = "./{dataset}/igphyml/results"
    conda: "r_dowser"
    script:
        'igphyml_write_trees.R'

rule extract_likelihood:
    input:     
        trees = "{dataset}/igphyml/trees.rds"
    output:
        likelihood = "{dataset}/igphyml/likelihoods.csv"
    conda: "r_dowser"
    script:
        "extract_likelihood.R"

rule convert_to_linforests:
    input: 
        mapping= "{dataset}/igphyml/name_isotype_mapping.csv",
        tree = "{dataset}/igphyml/results/{clone}.tree.csv",
        sequences =  "{dataset}/igphyml/results/{clone}.sequence.csv",
        isotypes = "{dataset}/recomb_input/{clone}/isotype.fasta",
        encoding = "mouse_isotype_encoding.txt"
    output:
        linforest = "{dataset}/igphyml/results/{clone}.pickle",
        png = "{dataset}/igphyml/pngs/{clone}.png",
        seq_mapping = "{dataset}/igphyml/seq_mappings/{clone}.mapping.csv",
    params:
        root = "Germline"
    run:
        import sys 
        sys.path.append("../src")
        import lineage_tree as lt 
        import utils as ut 
        import networkx as nx 
        encoding_dict = {}
        with open(input.encoding, "r+" ) as file:
            for i,line in enumerate(file):
                if "m" and "d" in line.lower():
                    encoding_dict["Ighm/d"] =i
    
                else:
                    encoding_dict[line.strip()] =i 
        
        iso_dict = ut.read_fasta(input.isotypes)
        iso_encodings = {key : encoding_dict[val] for key, val in iso_dict.items()}
        iso_encodings["naive"] =0       
        tree =nx.DiGraph()
        root = params.root
        edges = ut.read_edge_list(input.tree)
        
        seq = ut.read_dict(input.sequences)
        seq  = {key: list(value.strip()) for key, value in seq.items()}

        tree.add_edges_from(edges)
        root_name = [n for n in tree if tree.in_degree[n]==0][0]
        def reroot_tree(tree):
            if len(list(tree.neighbors(root))) == 0:
                G = tree.to_undirected()
                H = nx.dfs_tree(G,source=root)
                if H.out_degree[root_name]==0:
                    H.remove(root_name)
                return H 
            else:
                return tree 
        tree = reroot_tree(tree)
        df = pd.read_csv(input.mapping)
        df['seq_name'] = df['seq_name'].str.replace("_","")
        df  = df[df["clone_id"]==wildcards.clone]
        seq_dict = dict(zip(df["sequence_id"], df["seq_name"]))
        seq_dict[root] = "naive"
        ut.save_dict(seq_dict, output.seq_mapping)

    
        alignment = {}
        for n, s in seq.items():
            if n in seq_dict:
                alignment[seq_dict[n]] = s 
            else:
                alignment[n] = s


        lin_tree = lt.LineageTree(tree, "naive", name=wildcards.clone)
        lin_tree.relabel(seq_dict)
        lin_tree.save_png(output.png, iso_encodings, hide_underscore=False)
        lin_forest = lt.LineageForest(alignment, iso_encodings, [lin_tree])
        lin_forest.save_forest(output.linforest)

      













# rule prepdata:
#   input: "{dataset}/igphyml/input_data.tsv"
#   output: 
#      lineages = "{dataset}/igphyml/data_lineages.tsv"
#   params:
#     dirname = "data"
#   log: "{dataset}/igphyml/data_prep.log"
#   shell:
#      "BuildTrees.py -d {input} --outname {params.dirname} --log {log} --format airr"


# rule generate_topologies:
#     input:  "{dataset}/igphyml/data_lineages.tsv"
#     output: "{dataset}/igphyml/data_lineages_gy.tsv"
#     params: 
#         app = config["ig_path"]
#     threads: 10
#     log:
#         std=  "{dataset}/igphyml/generate_topologies.log",
#         err=  "{dataset}/igphyml/generate_topologies.err.log"
#     benchmark: "{dataset}/igphyml/generate_topologies.benchmark"
#     shell:
#         "nice -n 5 {params.app} --repfile {input} -m GY --outrep {output} --run_id gy --threads {threads} > {log.std} 2> {log.err}"

# rule igphyml:
#     input:  "{dataset}/igphyml/data_lineages_gy.tsv"
#     output: "{dataset}/igphyml/data_lineages_gy.tsv_igphyml_stats.txt"
#     threads: 10
#     log:
#         std=  "{dataset}/igphyml/igphyml.log",
#         err=  "{dataset}/igphyml/igphyml.err.log"
#     benchmark: "{dataset}/igphyml/igphyml.benchmark"
#     params: 
#         app = config["ig_path"]
#     shell:
#         "nice -n 5 {params.app} --repfile {input} -m HLP --threads {threads}  > {log.std} 2> {log.err}"