import argparse
import networkx as nx 
import pandas as pd
import itertools 

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

def get_clonos(dfx, dfy=None):
    if dfy is None:
        return dfx['clonotype'].tolist()
    else:
        setx = set(dfx['clonotype'].tolist())
        sety =set(dfy['clonotype'].tolist())
        return list(setx.intersection(sety))

def helper(tree, root, mut_nodes_temp, desc):
    children = list(tree.successors(root))
    mut_sucessors = set(children).intersection(set(mut_nodes_temp))
    for u in mut_sucessors:
        to_remove = []
        for v in mut_nodes_temp:
            if v not in mut_sucessors:
                
                if nx.has_path(tree, u,v ):
                    desc[u].append(v)
                    to_remove.append(v)
        for v in to_remove:
            mut_nodes_temp.remove(v)
    if len(mut_nodes_temp) ==0:
        return
    for c in children:
        helper(tree, c, mut_nodes_temp, desc)
    
    
    


def analyze_single_mut(clono, mut, tpath, df):
    tree = get_tree(clono, tpath)
    #traverese the tree in breadth first search
    dist=1
    mut_nodes = df[df['clonotype']==clono]['nodes'].values[0]
    mut_nodes = mut_nodes.split(",")
    mut_nodes =[m.strip() for m in mut_nodes]
  
    
    # [0].strip.split(",")
    mut_nodes_temp = [m for m in mut_nodes]
    leaf_nodes = [n for n in mut_nodes if "seq" in mut_nodes]
    # internal_nodes = [n for n in mut_nodes if n not in leaf_nodes]
    # incomp =0
    root = "naive"
    desc ={m: [] for m in mut_nodes}
    helper(tree, root,mut_nodes_temp, desc)
    return desc, tree

def first_muts(desc):
    total = 0
    fmuts = []
    for key, val in desc.items():
        if len(val) > 0:
            fmuts.append(key)
            total += (1+len(val))
    if total != len(desc):

        for key in desc:
            if 'seq' in key:
                in_list = False 
                for key2,val2 in desc.items():
                    if key in val2:
                        in_list = True
                        break
                if not in_list:

                    fmuts.append(key)
  
        
    return fmuts

            
def compare_dicts(desc1, desc2, tree):
    num_branches = 0
    num_ancestral_1 = 0
    num_ancestral_2 = 0
    fmuts1= first_muts(desc1)
    fmuts2 = first_muts(desc2)

    for u in fmuts1:
        for v in fmuts2:
            if nx.has_path(tree,u,v):
                num_ancestral_1+=1
            elif nx.has_path(tree,v,u):
                num_ancestral_2+=1
            else:
                num_branches +=1
    
    return num_branches, num_ancestral_1, num_ancestral_2



       
            


    # for u,v in itertools.product(mut_nodes):
    #     if len(nx.shortest_path(tree,u,v)):

    # for u,v in itertools.combinations(internal_nodes,2):
    #     pt1 = nx.shortest_path(tree, u,v)
    #     pt2 = nx.shortest_path(tree, v,u)
    #     if len(pt1)==0 and len(pt2)==0:
    #         incomp +=1
    #     if len(pt1)



    # while True:

    #     level_desc = nx.descendants_at_distance(tree, "naive",dist )
    #     incomp_nodes = mut_nodes.intersection(level_desc)


        
def get_tree(clono, tpath, lamb=0.75):
    clono_path = f"{tpath}/{clono}/{lamb}" 
    edges = read_dict(f"{clono_path}/tree.txt")
    tree = nx.DiGraph()

    tree.add_edges_from([(val, key) for key, val in edges.items()])
    return tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", required=False, type=str)
    parser.add_argument("-o", "--outfile", required=False, type=str)
    parser.add_argument("-t", "--treepath", required=False, type=str)
    parser.add_argument("-l", "--lamb", required=False, type=float, default=0.75)



    args= parser.parse_args()
    # mut = ""
    lamb =args.lamb
    pth = "/scratch/projects/tribal/real_data"
    # args= parser.parse_args([
    #     "-i", f"{pth}/day_14_W33L.csv",
    #     "-o", f"{pth}/day_14_W33L_level_0.75.csv",
    #     "--treepath", f"{pth}/tribal_search"
    # ])

   
    



    # df['min_mut_level'] = df.apply(get_levels, axis=1)

    # all_trees = {c: get_tree for c in df['clonotype']}
 
    muts = ["W33L", "S66N","K59R" ]
    datasets = ["day_14", "GCB_NP_1", "GCB_NP_2"]
    df_dict = {d: {} for d in datasets }
    for d in datasets:
        for m in muts:
                infile = f"{pth}/mut_levels_in/{d}_{m}.csv"
                try:
                    df_dict[d][m] = pd.read_csv(infile, skiprows=1, header=None, names=['clonotype', 'nodes'])
                except:
                    continue

  
    num_clonos = {d: {} for d in datasets}
   
    # incomp
    for d in datasets:
        dpath = f"{pth}/{d}/tribal_search"
        with open(f"{pth}/paired_analysis/{d}_paired_muts.csv", "w+") as file:
            file.write("clonotype,mutx,muty,num_branches,num_x_anc_y,num_y_anc_x\n")

            for x,y in itertools.combinations_with_replacement(muts, 2):
                if x not in df_dict[d] or y not in df_dict[d]:
                        continue
                if x ==y:
                    clonos = get_clonos(df_dict[d][x])
                    num_clonos[d][x,y] = len(clonos)
                 
                    for c in clonos:
                        desc,_ = analyze_single_mut(c,x, dpath, df_dict[d][x])
                        branches = first_muts(desc)
                        num_branches= len(branches)
                        
                            
                        # if num_branches == 2:
                        file.write(f"{c},{x},{y},{num_branches},NA,NA\n")
                        # print(desc)
                else:
                        clonos = get_clonos(df_dict[d][x], df_dict[d][y])
                        num_clonos[d][x,y] = len(clonos)
                        for c in clonos:
                            # if d=="day_14" and c=="B_147_6_100_148_1_41" and x =="W33L" and y=="S66N":
                            #     print("here")
                            # print(f"{d},{c},{x},{y}") 
                            desc1, t1 = analyze_single_mut(c,x, dpath, df_dict[d][x])
                            desc2, t2= analyze_single_mut(c,y, dpath, df_dict[d][y])
                            num_branches, num_ancestralx, num_ancestraly= compare_dicts(desc1, desc2, t1)
                            file.write(f"{c},{x},{y},{num_branches},{num_ancestralx},{num_ancestraly}\n")


           
               
            # else:
            #     
            
            #  
    with open(f"{pth}/paired_analysis/nclonotypes_by_paired_muts.csv", "w+") as file:
        file.write("dataset,mutx,mutx,nclonotypes\n")
        for d in num_clonos:
                for x,y in num_clonos[d]:
                    file.write(f"{d},{x},{y},{num_clonos[d][x,y]}\n")

    # print(num_clonos)
            #     incomp[x,y] = {}
            #     losses[x] ={}
            #     for c in clonos:
            #         incomp[x,y][c], losses[x][c] = analyze_single_mut(c,x)
            # else:
            #     clonos = get_clonos(x)
            #     #get count of incomparables 
            #     #get count of loss 
            #     pass 
                





    # df.to_csv(args.outfile, index=False)


    

    

