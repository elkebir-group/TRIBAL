import numpy as np
import pydot 
import argparse
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt

class DrawStateDiag:
    def __init__(self, transmat, state_probs=None, isotype_mapping=None, color_encoding=None, cmap_name="rocket_r") -> None:
    
        if isotype_mapping is None:
            isotype_mapping = {i: i for i in range(transmat.shape[0])}
        
        self.isotype_mapping = isotype_mapping
        
        self.cmap = sns.color_palette(cmap_name, as_cmap=True)
        self.states = [i for i in range(transmat.shape[0])]
        self.transmat = transmat
  
        if color_encoding is None:
            self.color_encoding =  {
                0 : "antiquewhite",
                1 : "turquoise1",
                2 : "darkcyan",
                3 : "cornflowerblue",
                4 : "darksalmon",
                5 : "mediumseagreen",
                6 : "orangered",
                7 : "orchid",
                8 : "darkgoldenrod",
                9 : "thistle1"
            }

        self.graph = pydot.Dot("state_diagram", graph_type="digraph", bgcolor="white")
        added_nodes = []
        for s in self.states:
            lab =isotype_mapping[s]
            if state_probs is not None:
                if state_probs[s] >= 0.00:
          
                    added_nodes.append(s)
                    lab +=  "\n" + str(round(state_probs[s],3))

            col = self.color_encoding[s]
            self.graph.add_node(pydot.Node(isotype_mapping[s], shape="circle", color=col, style='filled',label=lab))
        

        for s in self.states:
            # if s not in added_nodes:
            #     continue
            for t in self.states:
                if self.transmat[s,t] > 0.01:# and t in added_nodes:
                    # if t not in added_nodes and s == t:
                    #     continue
         
    
            
                    lab = str(round(self.transmat[s,t],3))
                    new_edge = pydot.Edge(dst=isotype_mapping[t], src=isotype_mapping[s], color="black", label=lab)
                
                    self.graph.add_edge(new_edge)
    
     
    def heatmap(self,fname):
        plt.figure()
        labs= [self.isotype_mapping[s] for s in self.states]
        df = pd.DataFrame(self.transmat, index=labs, columns=labs)
        fig =sns.heatmap(df, annot=True, fmt=".03f", cmap=self.cmap)

        fig.set(xlabel="Isotype", ylabel="Isotype")
        plt.savefig(fname)
      

        
    def save(self, fname):
        self.graph.write_png(fname)
    
    def save_pdf(self, fname):


        self.graph.write_pdf(fname) 



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--transmat", required=True, type=str,
        help="filename of tree file in parent edge list form")
    parser.add_argument("-s", "--state_probs", required=False, type=str,
        help="filename of tree file in parent edge list form")

    parser.add_argument("-e", "--encoding", required=False, type=str,
        help="filename of input transition matrix")

    parser.add_argument("-o", "--outfile", type=str)
    parser.add_argument("--pdf", type=str)
    parser.add_argument("--heatmap", type=str)
    args= parser.parse_args()

    # args = parser.parse_args([
    #     "-t", "/scratch/projects/tribal/benchmark_pipeline/sim_data/transmats/transmat1.txt",
    #     "-e", "/scratch/projects/tribal/benchmark_pipeline/sim_encoding.txt",
    #     "--heatmap", "src/test/heatmap.png"

    # ])


    tmat = np.loadtxt(args.transmat)
    if args.state_probs is not None:
        state_probs = np.loadtxt(args.state_probs)
    else:
        state_probs = None
    if args.encoding is not None:

        with open(args.encoding, 'r+') as file:
                    isotype_encoding = {}
                    counter = 0
                    for line in file:
                        isotype_encoding[counter] = line.strip()
                        counter += 1
    else:
        isotype_encoding = None


    ds = DrawStateDiag(tmat, state_probs, isotype_encoding)
    if args.pdf is not None:
        ds.save_pdf(args.pdf)
    
    if args.outfile is not None:
        ds.save(args.outfile)

    if args.heatmap is not None:
        ds.heatmap(args.heatmap)