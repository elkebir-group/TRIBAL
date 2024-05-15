import pandas as pd
import argparse
import utils as ut  
from io import StringIO
import subprocess

def convert_to_string(seqs, root):
    mystr = f">naive\n{root}\n"
    for key, val in seqs.items():
        mystr += f">{key}\n{val}\n"
    return mystr

def align(seqs, root):

    #!cat seqs.fasta mafft --quiet - > aligned.fasta
    seq_str = convert_to_string(seqs,root)
    child = subprocess.Popen(["mafft", "--quiet", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    child.stdin.write(seq_str.encode())
    child_out = child.communicate()[0].decode("utf8")

    child.stdin.close()
    aligned_seqs = ut.read_fasta(child_out)
 
    return aligned_seqs
  

def run_dnapars(cfg_str, phylip_str):
    # Prepare input and output handles
    input_handle = StringIO(phylip_str)
    output_handle = StringIO()
    treefile_handle = StringIO()

    # Create the configuration file buffer
    cfg_handle = StringIO(cfg_str)
    
    # Run dnapars
    dnapars_cline = ["dnapars"]
    process = subprocess.Popen(dnapars_cline, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=phylip_str + "\n" + cfg_str)
    
    if process.returncode != 0:
        raise RuntimeError(f"dnapars failed with error code {process.returncode}: {stderr}")
    
    return stdout, stderr

def align_clonotypes(df, roots):
    clonotype_alignments= {}
    for j in df["Clonotype"].unique():
        clono = df[df["Clonotype"]== j]
        heavy_seqs = dict(zip(clono["seq"],clono["Heavy Chain Variable Seq"]))
    
        heavy_root = roots[roots["Clonotype"]==j]["Heavy Chain Root"].values[0]
  
        heavy_chain_align = align(heavy_seqs, heavy_root)
        light_seqs = dict(zip(clono["seq"],clono["Light Chain Variable Seq"]))
    
        light_root = roots[roots["Clonotype"]==j]["Light Chain Root"].values[0]
        light_chain_align = align(light_seqs, light_root)
        concat_seqs = {}
        for key in heavy_seqs:
            concat_seqs[key] = heavy_chain_align[key].upper() + light_chain_align[key].upper()
        clonotype_alignments[j] = concat_seqs
    return clonotype_alignments

def main(args):
    df  =pd.read_csv(args.file)
   
    print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
    df = df.groupby("Clonotype").filter(lambda x: len(x) >= args.min_size)
    df['seq'] = df.groupby('Clonotype').cumcount() + 1
    df['seq'] = 'seq' + df['seq'].astype(str)
    print(f"After filtering, the number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
    print(df.head())
    df.columns.values[0] = "cellid"
    roots =pd.read_csv(args.roots)
    roots.columns.values[0] = "Clonotype"
    alignments = align_clonotypes(df, roots)
    print(len(alignments))
    print("done")
   
    
        











if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file", type=str,
        help="filename of csv file with the clone data")
    parser.add_argument("-r", "--roots", type=str,
        help="filename of csv file with the root sequences")
    parser.add_argument( "--min-size", type=int, default=4,
        help="minimum clonotype size")
    # parser.add_argument("-i", "--iso", type=str,
    #     help="filename of output fasta file with isotypes")
    # parser.add_argument("-m", "--mapping", type=str,
    #     help="mapping of short sequence using for dnapars and original sequence id")

    args= parser.parse_args()
    args = parser.parse_args([
        "-f", "Human/human_data.csv",
        "-r", "Human/human_data_root_seq.csv",

    ])

    main(args)


  
    # seqs = dict(zip(df['seq_short'], df['sequence']))
    # isos = dict(zip(df['seq_short'], df['c_call_short']))

    # id_map = df[['seq_short', 'sequence_id']]
    # id_map = id_map[id_map["seq_short"] != "naive"]
    # for key in isos:
    #     if isos[key] =="D":
    #         isos[key] = "M"

    # if args.seq is not None:
    #     ut.write_fasta(args.seq, seqs)
    # if args.iso is not None:

    #     ut.write_fasta(args.iso, isos)
    # if args.mapping is not None:
    #     id_map.to_csv(args.mapping, index=False)
    

    

import subprocess
import io
from Bio import AlignIO

def run_mafft(input_fasta):
    # Run MAFFT and capture the output
    mafft_cline = ["mafft", "--auto", "--quiet", "-"]
    process = subprocess.Popen(mafft_cline, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=input_fasta)
    
    if process.returncode != 0:
        raise RuntimeError(f"MAFFT failed with error code {process.returncode}: {stderr}")
    
    return stdout

def convert_alignment_to_phylip(alignment_str, input_format="fasta"):
    # Create a file-like object from the string
    input_handle = io.StringIO(alignment_str)
    
    # Read the alignment
    alignment = AlignIO.read(input_handle, input_format)
    
    # Write the alignment in PHYLIP format to a string buffer
    output_handle = io.StringIO()
    AlignIO.write(alignment, output_handle, "phylip")
    
    # Get the content of the PHYLIP formatted alignment
    phylip_str = output_handle.getvalue()
    
    return phylip_str



# Example usage
input_fasta = """>seq1
AGCTAGCTAG
>seq2
CGTAGCTAGC
>seq3
TTAGCTAGCT"""

# Step 1: Run MAFFT to get the aligned sequences
aligned_sequences = run_mafft(input_fasta)

# Step 2: Convert the aligned sequences to PHYLIP format
phylip_str = convert_alignment_to_phylip(aligned_sequences)

# Step 3: Generate the configuration file as a string
cfg_str = "O\n1\nS\nY\n4\n5\n.\nY"

# Step 4: Run dnapars using the PHYLIP formatted data and the configuration file
dnapars_stdout, dnapars_stderr = run_dnapars(phylip_str, cfg_str)

print("dnapars output:", dnapars_stdout)
print("dnapars error (if any):", dnapars_stderr)
