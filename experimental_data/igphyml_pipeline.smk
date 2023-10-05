configfile: "igphyml_config.yml"


rule all:
   input:  
        expand("{dataset}/igphyml/data_lineages_gy.tsv_igphyml_stats.txt",
            dataset = config['datasets']
        ),
        expand("{dataset}/igphyml/trees.rds",
            dataset = config['datasets']
        )


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

#snakemake --use-conda --conda-frontend mamba --conda-prefix /scratch/data/leah/anaconda3/envs/r_dowser -s igphyml_pipeline.smk -j 1
rule run_igphyml:
    input: 
        airr = "{dataset}/igphyml/input_data.tsv"
    output:
       clones= "{dataset}/igphyml/clones.rds",
       trees = "{dataset}/igphyml/trees.rds",
       results = directory("{dataset}/igphyml/reults")
    conda: "/scratch/data/leah/anaconda3/envs/r_dowser"
    params:
        nproc = 7
    script:
        'igphyml_run.R'




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