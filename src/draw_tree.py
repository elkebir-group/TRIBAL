import networkx as nx
import pydot 
from utils import read_dict 
import argparse 

class DrawTree:
    def __init__(self, parents, isotypes, color_encoding=None, show_legend=False, isotype_encoding=None) -> None:
        
        self.parents = parents
        self.isotypes = isotypes

        if color_encoding is None:
            self.color_encoding =  {
                0 : "antiquewhite",
                1 : "dimgrey",
                2 : "darkcyan",
                3 : "cornflowerblue",
                4 : "darksalmon",
                5 : "mediumseagreen",
                6 : "orangered",
                7 : "orchid",
                8 : "darkgoldenrod",
                9 : "thistle2"
            }
        else:
            self.color_encoding = color_encoding

        self.graph = pydot.Dot("tribal_tree", graph_type="digraph", bgcolor="white")

        # Add nodes
        iso_trans = {}
        for u in self.color_encoding:
            for v in self.color_encoding:
                iso_trans[u,v] =0 
       
        total_trans = 0
        used_isotypes = []
        for p in self.parents:
            iso = int(self.isotypes[p])
            if iso not in used_isotypes:
                used_isotypes.append(iso)
            col = self.color_encoding[iso]
            self.graph.add_node(pydot.Node(str(p), shape="circle", color=col, style='filled'))
        
        for key, val in self.parents.items():
            if val != "":
                new_edge = pydot.Edge(dst=key, src=val, color="black")
                iso_u = int(self.isotypes[val])
                iso_v =int(self.isotypes[key])
                iso_trans[iso_u,iso_v] += 1
                total_trans += 1
              
                self.graph.add_edge(new_edge)
        
        if show_legend:
            used_isotypes.sort()
            if isotype_encoding is not None:
                for u in used_isotypes:
                    iso_name = isotype_encoding[u]
                    self.graph.add_node(pydot.Node(iso_name, shape="circle", color=self.color_encoding[u], style='filled'))
                
                for u,v in iso_trans:
                    if iso_trans[u,v] > 0:
                            u_name = isotype_encoding[u]
                            v_name = isotype_encoding[v]
                            lab = str(round(iso_trans[u,v]/total_trans,2))
                            new_edge = pydot.Edge(src=u_name, dst=v_name, color="black", label=lab)
                            self.graph.add_edge(new_edge)

    
    def save(self, fname):
        self.graph.write_png(fname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", required=True, type=str,
        help="filename of tree file in parent edge list form")
    parser.add_argument("-i", "--isotypes", required=True, type=str,
        help="filename of input file containing inferred isotypes")
    parser.add_argument("-e", "--encoding", required=False, type=str,
        help="filename of input transition matrix")
    parser.add_argument("-l", "--legend", type=bool, action="store_true")
    parser.add_argument("-o", "--outfile", type=str, default="TRIBAL_Tree.png")
    args= parser.parse_args()


    iso_labels = read_dict(args.isotypes)
    parents = read_dict(args.tree)
    if args.encoding is not None:
        with open(encoding_fname, 'r+') as file:
            isotype_encoding = {}
            counter = 0
            for line in file:
                isotype_encoding[counter] = line.strip()
                counter += 1
    else:
        isotype_encoding = None

    dt = DrawTree(parents, iso_labels, show_legend=args.legend,isotype_encoding=isotype_encoding )
    
    dt.save(args.outfile)


# clono = "B_75_1_1_148_1_40"
# encoding_fname = "/scratch/projects/tribal/real_data/mouse_isotype_encoding.txt"

# for i in range(6):
    # parents_fname = f"/scratch/projects/tribal/real_data/day_14/tribal_fit/{clono}/candidate{i}.tree.txt"
    # iso_fname = f"/scratch/projects/tribal/real_data/day_14/tribal_fit/{clono}/candidate{i}.isotypes"
    # out_fname = f"/scratch/projects/tribal/real_data/day_14/tribal_fit/{clono}/candidate{i}.png"



  







    

        
