import networkx as nx
import pygraphviz as pgv
from utils import read_dict 
import argparse 

class DrawTree:
    def __init__(self, parents, isotypes=None, color_encoding=None, root="naive", 
                 show_legend=False, isotype_encoding=None, 
                 show_labels=True, hide_underscore=True, use_dashed=False) -> None:
        self.parents = parents
        if isotypes is not None:
            self.isotypes = isotypes
        else:
            nodes = [key for key in parents] + [val for key, val in parents.items()]
            nodes =set(nodes)
            self.isotypes = {n: 0 for n in nodes }
        
    
     
        

        self.root= root

        if color_encoding is None:
            self.color_encoding =  {
                -1: "#FFFFFF",
                0 : "#808080",
                1 : "#FFEDA0",
                2 : "#FD8D3C",
                3 : "#E31A1C",
                4 : "#800026",
                5 : "#6A51A3",
                6 : "#74C476",
                7 : "mediumseagreen",
                8 : "darkgoldenrod",
                9 : "thistle1"
            }
          
        else:
            self.color_encoding = color_encoding

        self.graph = pgv.AGraph(directed=True, bgcolor="white")
        def add_node(p):
            if p in self.isotypes:
                iso = int(self.isotypes[p])
            else:
                iso = -1
            if iso not in used_isotypes:
                used_isotypes.append(iso)
            fill_color = self.color_encoding[iso]
            outline_color = "black"  # Set outline color here
            lab = str(p)
            if p == self.root:
                lab = "r"
       

            if "_" in lab and hide_underscore:
                    lab =lab.split("_")[0]
       
            if not show_labels:
                lab = ""
            self.graph.add_node(str(p), label=lab, shape="circle", style='filled', penwidth="1", fillcolor=fill_color, color=outline_color)


        used_isotypes = []
        for p in self.parents:
            add_node(p)
        
        for key, val in self.parents.items():
            if val != "":
                if val == self.root:
                    add_node(val)
                new_edge = (val, key)
                try:
                    if self.isotypes[val] == self.isotypes[key]:

                        self.graph.add_edge(*new_edge, color="black")
                    else:
                        if self.use_dashed:
                            self.graph.add_edge(*new_edge, color="black", style="dashed")
                        else:
                            self.graph.add_edge(*new_edge, color="black")
                except:
                    self.graph.add_edge(*new_edge, color="black")
        
        if show_legend:
            used_isotypes.sort()
            # if isotype_encoding is not None:
            for u in used_isotypes:
                    if u < 0:
                        continue
                    # iso_name = self.color_encoding[u]
                    self.graph.add_node(f"i_{u}", label=str(u), shape="circle", fillcolor=self.color_encoding[u], style='filled', color='black')
                
        #         for u,v in iso_trans:
        #             if iso_trans[u,v] > 0:
        #                     u_name = isotype_encoding[u]
        #                     v_name = isotype_encoding[v]
        #                     lab = str(round(iso_trans[u,v]/total_trans,2))
        #                     new_edge = pydot.Edge(src=u_name, dst=v_name, color="black", label=lab)
        #                     self.graph.add_edge(new_edge)

    
    def save(self, fname):
        self.graph.draw(fname, prog="dot")
    

    def save_pdf(self, fname):


         self.graph.draw(fname, prog="dot") 
# Or, save it as a DOT-file:
# graph.write_raw("output_raw.dot")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", required=True, type=str,
        help="filename of tree file in parent edge list form")
    parser.add_argument("-i", "--isotypes", required=True, type=str,
        help="filename of input file containing inferred isotypes")
    parser.add_argument("-e", "--encoding", required=False, type=str,
        help="filename of input transition matrix")
    parser.add_argument("-l", "--legend", action="store_true")
    parser.add_argument("-o", "--outfile", type=str)
    parser.add_argument("--pdf", type=str)
    args= parser.parse_args()

    # clono = "B_75_1_1_148_1_40"
    # encoding_fname = "/scratch/projects/tribal/real_data/mouse_isotype_encoding.txt"

    # i = 0
    # parents_fname = f"/scratch/projects/tribal/real_data/day_14/tribal_fit/{clono}/candidate{i}.tree.txt"
    # iso_fname = f"/scratch/projects/tribal/real_data/day_14/tribal_fit/{clono}/candidate{i}.isotypes"
    # out_fname = f"/scratch/projects/tribal/real_data/day_14/tribal_fit/{clono}/candidate{i}.png"
    
    # args= parser.parse_args([
    #     "-t", "bcr-phylo-benchmark/sim_data/replicates/cells35/size25/rep1/2.0/0.365/1/GCsim_collapsed_tree.parents",
    #     "-i", "bcr-phylo-benchmark/sim_data/replicates/cells35/size25/rep1/2.0/0.365/1/GCsim_collapsed_tree.isotypes",
    #     "-o", "bcr-phylo-benchmark/sim_data/replicates/cells35/size25/rep1/2.0/0.365/1/GCsim_collapsed_tree.png"


    # ])

    
    # args= parser.parse_args([
    #     "-t", "bcr-phylo-benchmark/sim_data/replicates/cells35/size25/rep1/2.0/0.365/tribal_refine/1/0.75/inferred_collapsed_tree.parents",
    #     "-i", "bcr-phylo-benchmark/sim_data/replicates/cells35/size25/rep1/2.0/0.365/tribal_refine/1/0.75/inferred_collapsed_tree.isotype",
    #     "-o", "bcr-phylo-benchmark/sim_data/replicates/cells35/size25/rep1/2.0/0.365/tribal_refine/1/0.75/inferred_collapsed_tree.png"


    # ])

    iso_labels = read_dict(args.isotypes)
    parents = read_dict(args.tree)
    if args.encoding is not None:
        with open(args.encoding, 'r+') as file:
            isotype_encoding = {}
            counter = 0
            for line in file:
                isotype_encoding[counter] = line.strip()
                counter += 1
    else:
        isotype_encoding = None

    dt = DrawTree(parents, iso_labels, show_legend=args.legend,isotype_encoding=isotype_encoding )
    
    dt.save(args.outfile)

    if args.pdf is not None:
        dt.save_pdf(args.pdf)










  







    

        

