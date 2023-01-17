import argparse
import networkx as nx 
import pandas as pd


def read_dict(fname):
    mydict= {}
    with open(fname, 'r') as file:
        for line in file:
            line_list =line.strip().split(",")
            mydict[line_list[0]] = line_list[1]
    return mydict

def get_levels(row):
    clono = row['clonotype']
    mut_nodes = row['nodes'].split(",")
    clono_path = f"{tree_path}/{clono}/{lamb}" 
    edges = read_dict(f"{clono_path}/tree.txt")
    tree = nx.DiGraph()
    root = "naive"
    tree.add_edges_from([(val, key) for key, val in edges.items()])
    
    path_lengths =[nx.shortest_path_length(tree, root, n.strip()) for n in mut_nodes]
    
    return(min(path_lengths))




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", required=False, type=str)
    parser.add_argument("-o", "--outfile", required=False, type=str)
    parser.add_argument("-t", "--treepath", required=False, type=str)
    parser.add_argument("-l", "--lamb", required=False, type=float, default=0.75)



    args= parser.parse_args()

    lamb =args.lamb
    # pth = "/scratch/projects/tribal/real_data/day_14"
    # args= parser.parse_args([
    #     "-i", f"{pth}/day_14_W33L.csv",
    #     "-o", f"{pth}/day_14_W33L_level_0.75.csv",
    #     "--treepath", f"{pth}/tribal_search"
    # ])

    tree_path = args.treepath
    

    df = pd.read_csv(args.infile, skiprows=1, header=None, names=['clonotype', 'nodes'])
    df['min_mut_level'] = df.apply(get_levels, axis=1)
    df.to_csv(args.outfile, index=False)


    

    

