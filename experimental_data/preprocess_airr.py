import pandas as pd
from collections import Counter
import argparse 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--clonotypes", required=True, type=str,
        help="clonotypes file")
    parser.add_argument("-m", "--metadata", required=True, type=str,
        help="dandelion metadata file")
    parser.add_argument("-t", "--table", required=True, type=str,
        help="dandelion table file")
    parser.add_argument("-o", "--output", required=False, type=str,
        help="output lineages file")
    parser.add_argument("-p", "--mapping", required=False, type=str,
        help="output lineages file")
    
    parser.add_argument("-d", "--dataset", type=str, required=True, help="name of dataset")

    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    # dataset = "GCB_NP_2"
    # dpath = f"/scratch/projects/tribal/experimental_data/{dataset}"
    dpath = args.path
    dataset = args.dataset


    def read_tsv(fname, first_col="index"):

        firstline  =True
        lines = []
        line_lengths = []
        with open(fname, "r+") as file:
            for line in file:
                if firstline:
                    header = line.strip().split("\t")
                    ncols = len(header)
                
                    firstline = False 
                else:
                    line = line.strip().split("\t")
                    line_lengths.append(len(line))
                    lines.append(line)
        header[0] = first_col
        df = pd.DataFrame(lines, columns=header)
        return(df)

    metadat = read_tsv(args.metadata)

    clonotypes = pd.read_csv(args.clonotypes, names=["clone_id"])

    meta_filt = pd.merge(metadat, clonotypes, on='clone_id', how='inner')

    df = meta_filt
    result = df[df['v_call'].str.startswith('IGH')]

    # Group by 'clone_id' and 'germline_alignment_d_mask', and calculate the count within each group
    grouped = result.groupby(['clone_id', 'germline_alignment_d_mask']).size().reset_index(name='count')

    # Find the germline sequence with the maximum count within each 'clone_id' group
    max_seq = grouped.groupby('clone_id').apply(lambda x: x[x['count'] == x['count'].max()])

    # Reset the index to clean up the result DataFrame
    max_seq.reset_index(drop=True, inplace=True)

    # Display the result with the germline sequence having the maximum count
    # max_seq

    num_germlines = max_seq.groupby("clone_id").size()
    bad_germlines = num_germlines[num_germlines > 1]
    bad_clonotypes =bad_germlines.index.values

    germline_to_keep = max_seq[~(max_seq['clone_id'].isin(bad_clonotypes))]
    germline_to_keep = germline_to_keep.drop('count', axis=1)

    meta_filt2 = pd.merge(meta_filt, germline_to_keep, on=['clone_id', 'germline_alignment_d_mask'], how="inner")
    meta_filt2= meta_filt2.set_index("index")
    
    if args.output is not None:
        meta_filt2.to_csv(args.output, index=False, sep="\t")

    tabdat = pd.read_csv(args.table, sep="\t")

    tabdat= tabdat.rename(columns={'Unnamed: 0': 'cellId', 'Clonotype': 'clone_id','Heavy Chain Isotype (Expression)': 'isotype' })



    tab_filt = tabdat[['cellId', 'clone_id',  'isotype']]

    meta_iso =tab_filt.merge(meta_filt2, on=['cellId', 'clone_id'], how="inner")

    meta_iso =meta_iso[['sequence_id', 'isotype', 'clone_id', 'cellId']]


    meta_iso['seq_name'] = 'seq_' + (meta_iso.groupby('clone_id').cumcount() + 1).astype(str)
    meta_iso = meta_iso.sort_values(['clone_id', 'seq_name'])
    meta_iso = meta_iso.sort_values(['clone_id', 'seq_name'])
    if args.mapping is not None:
        meta_iso.to_csv(args.mapping, index=False)
