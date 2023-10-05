configfile: "dnapars.yml"
import pandas as pd 



def get_files(dirn, fname):
    targets = []
    for d in config["datasets"]:
        df =pd.read_csv(f"{d}/igphyml/name_isotype_mapping.csv")
        clonotypes = df['clone_id'].unique()
        for c in clonotypes:
            targets.append(f"{d}/{dirn}/{c}/{fname}")
    return targets 




rule all: 
    input: 
        get_files("recomb_input", "heavy.aln.fasta"),
        get_files("dnapars", "dnapars.log")




rule recode_sequences:
    input: 
        fasta= "{dataset}/igphyml/data/{c}.fasta",
        mapping= "{dataset}/igphyml/name_isotype_mapping.csv"
    output:
        alignment = "{dataset}/recomb_input/{c}/heavy.aln.fasta",
        isotype ="{dataset}/recomb_input/{c}/isotype.fasta"
    shell:
        "python recode_seqs.py -f {input.fasta} -m {input.mapping} "
        "-i {output.isotype} -o {output.alignment} -c {wildcards.c} "


rule convert_phylip:
    input: "{dataset}/recomb_input/{c}/heavy.aln.fasta"
    output: "{dataset}/recomb_input/{c}/heavy.algn.phylip"
    shell:
        "seqmagick convert {input} {output}"


rule dnapars_config:
    input: "{dataset}/recomb_input/{c}/heavy.algn.phylip"
    output: "{dataset}/dnapars/{c}/dnapars.cfg"
    shell:
        "mkconfig {input} dnapars > {output} "


rule dnapars_run:
    input: "{dataset}/dnapars/{c}/dnapars.cfg",
    params: 
        direct="{dataset}/dnapars/{c}"
    output:  "{dataset}/dnapars/{c}/dnapars.log"
    shell:       
        """
        (
        cd {params.direct}
        dnapars < dnapars.cfg > dnapars.log
        )
        """
   



# rule preprocess:
#     input: 
#         data_fname = "{run_dir}/{run_dir}_dandelion_table.tsv",
#         root_fname = "{run_dir}/{run_dir}_root_sequences.csv"
#     params: 
#         min_size = 5,
#         pth = "/scratch/projects/tribal/real_data/{run_dir}/input"
#     output: 
#         summary = "{run_dir}/clonotype_summary.csv",
#         id_mapping = "{run_dir}/barcode_id_mapping.csv"
#     script:
#         "scripts/preprocess.R"


# rule align:
#     input: 
#         fasta = "{run_dir}/input/{clonotype}/{seq}.fasta",
#         summary = "{run_dir}/clonotype_summary.csv"
#     output: "{run_dir}/input/{clonotype}/{seq}.aln.fasta"
#     log:
#         err = "{run_dir}/input/{clonotype}/{seq}.muscle.err.log",
#         std ="{run_dir}/input/{clonotype}/{seq}.muscle.log"
#     shell:
#         "./muscle5.1.linux_intel64 -align {input.fasta} -output {output} > {log.std} 2> {log.err}"


# rule concat_fasta:
#     input:
#         light= "{run_dir}/input/{clonotype}/light.aln.fasta",
#         heavy= "{run_dir}/input/{clonotype}/heavy.aln.fasta"
#     output:  "{run_dir}/input/{clonotype}/concat.aln.fasta"
#     shell:
#         "python ../scripts/concat_fasta.py -l {input.light} -y {input.heavy} -o {output}"


     
       


# rule tribal_infer:
#     input: 
#         clonotypes = "{run_dir}/clonotypes.txt",
#         encoding = "mouse_isotype_encoding.txt",
#     params:
#         inpath = "/scratch/projects/tribal/real_data/{run_dir}/input",
#         alpha = 0.75,
#         max_cand = config["max_cand"],
#         niter = config["niter"],
#         thresh = config["threshold"],
#         seed =  config["seed"],
#         root = "naive",
#         outpath = "/scratch/projects/tribal/real_data/{run_dir}/tribal/{jump_prob}",
#         tree_path = "/scratch/projects/tribal/real_data/{run_dir}/dnapars",
#         nworkers = config['nworkers'],
#         restarts = config["nrestarts"],
#         mu = config["mu"],
#         sigma = config["sigma"]
#     output:
#         transmat = "{run_dir}/tribal/{jump_prob}/transmat.txt",
#         stateprobs = "{run_dir}/tribal/{jump_prob}/state_probs.txt",
#         score = "{run_dir}/tribal/{jump_prob}/fit_score.txt",
#         diagram = "{run_dir}/tribal/{jump_prob}/transmat.png",
#         diagram_pdf = "{run_dir}/tribal/{jump_prob}/transmat.pdf",
#     log:
#         std = "{run_dir}/tribal/{jump_prob}/fit.log",
#         err = "{run_dir}/tribal/{jump_prob}/fit.err.log"
#     shell:
#         "python ../src/tribal.py -c {input.clonotypes} -p {params.inpath} "
#         "-r {params.root} -e {input.encoding} -s {params.seed} "
#         "-o {params.outpath} --alpha {params.alpha} --niter {params.niter} "
#         "--thresh {params.thresh} --transmat_infer {output.transmat}  --tree_path {params.tree_path}  "
#         "--state_probs {output.stateprobs} --score {output.score} --diagram_pdf {output.diagram_pdf}  "
#         "--nworkers {params.nworkers} --restarts {params.restarts} --mu {params.mu} --sigma {params.sigma} "
#         "--save_all_restarts {params.outpath} -j {wildcards.jump_prob} "
#         "--diagram {output.diagram} > {log.std} 2> {log.err} "


# rule state_pdf:
#     input: 
#         transmat =  "{run_dir}/tribal/{jump_prob}/transmat.txt",
#         encoding = "mouse_labels.txt",
#     output: 
#         heatmap ="heatmaps/{run_dir}_{jump_prob}.pdf",
#         state_heatmap="heatmaps/{run_dir}_{jump_prob}_isotype_prop.pdf"
#     shell:
#         "python ../src/draw_state_diagram.py -t {input.transmat}  -e {input.encoding} "
#         "--heatmap {output.heatmap} --statemap {output.state_heatmap} "

# rule save_heatmap:
#     input: 
#         infer =  "{run_dir}/tribal/0.25/transmat.txt",
#         encoding = "mouse_isotype_encoding.txt",
#     output: 
#         heatmap ="heatmaps/{run_dir}.pdf",
#         png  ="heatmaps/{run_dir}.png"
#     run:
#         import numpy as np 
#         import sys 
#         sys.path.insert(0, '/scratch/projects/tribal/src')
#         from draw_state_diagram import DrawStateDiag

#         infer = np.loadtxt(input.infer)
#         ds = DrawStateDiag(infer)
#         ds.heatmap(output.heatmap)
#         ds.heatmap(output.png)






