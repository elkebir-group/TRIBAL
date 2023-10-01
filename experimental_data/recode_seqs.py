import pandas as pd
from collections import Counter
import argparse 
import sys 

sys.path.append("/scratch/projects/tribal/src")
from utils import read_fasta, write_fasta

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta", required=True, type=str,
        help="fasta file")
    parser.add_argument("-m", "--mapping", required=True, type=str,
        help="sequence_id to seq_name mapping")
    parser.add_argument("-i", "--isotypes", required=False, type=str,
        help="dandelion table file")
    parser.add_argument("-o", "--output", required=False, type=str,
        help="output alignment")
    parser.add_argument("-c", "--clonotype", required=True, type=str)
 
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    alignment= read_fasta(args.fasta)
  
    mapping = pd.read_csv(args.mapping)

    df = mapping[mapping['clone_id']==args.clonotype]

    id_to_name = {}
    iso_labs = {}
   
    for index,row in df.iterrows():
        id = row["sequence_id"]
        name = row["seq_name"]
        iso = row["isotype"]
        id_to_name[id] = name
        if iso in ["Ighm", "Ighd"]:
            iso = "Ighm/d"
        iso_labs[id] = iso
    new_alignment = {}
    isotypes = {"naive": "Ighm/d"}
    for a, seq in alignment.items():
        if a not in id_to_name:
            new_alignment["naive"] = seq
    
        else:
            nm =id_to_name[a].split("_")
            label = "".join(nm)
            new_alignment[label] = seq
            isotypes[label] = iso_labs[a]
    
    if args.output is not None:
        write_fasta(args.output, new_alignment)
    
    if args.isotypes is not None:
        write_fasta(args.isotypes, isotypes)
    
    







    