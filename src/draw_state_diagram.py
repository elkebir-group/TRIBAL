import numpy as np
import pydot 
import argparse
class DrawStateDiag:
    def __init__(self, transmat, state_probs, isotype_mapping=None, color_encoding=None) -> None:
        
        if isotype_mapping is None:
            isotype_mapping = {i: i for i in range(transmat.shape)}
        

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
            if state_probs[s] >= 0.00:
          
                added_nodes.append(s)
                lab = isotype_mapping[s] + "\n" + str(round(state_probs[s],3))

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
        
    def save(self, fname):
        self.graph.write_png(fname)
    
    def save_pdf(self, fname):


        self.graph.write_pdf(fname) 



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--transmat", required=True, type=str,
        help="filename of tree file in parent edge list form")
    parser.add_argument("-s", "--state_probs", required=True, type=str,
        help="filename of tree file in parent edge list form")

    parser.add_argument("-e", "--encoding", required=False, type=str,
        help="filename of input transition matrix")

    parser.add_argument("-o", "--outfile", type=str, default="transmat.png")
    parser.add_argument("--pdf", type=str, default="transmat.pdf")
    args= parser.parse_args()


    tmat = np.loadtxt(args.transmat)
    state_probs = np.loadtxt(args.state_probs)
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