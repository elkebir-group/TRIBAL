import numpy as np
import pydot 
class DrawStateDiag:
    def __init__(self, transmat, state_probs, isotype_mapping=None, color_encoding=None) -> None:
        
        if isotype_mapping is None:
            isotype_mapping = {i: i for i in range(transmat.shape)}
        

        self.states = [i for i in range(transmat.shape[0])]
        self.transmat = transmat
        
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

        self.graph = pydot.Dot("state_diagram", graph_type="digraph", bgcolor="white")
        added_nodes = []
        for s in self.states:
            if state_probs[s] >= 0.00:
          
                added_nodes.append(s)
                lab = isotype_mapping[s] + "\n" + str(round(state_probs[s],2))

                col = self.color_encoding[s]
                self.graph.add_node(pydot.Node(isotype_mapping[s], shape="circle", color=col, style='filled',label=lab))
        

        for s in self.states:
            # if s not in added_nodes:
            #     continue
            for t in self.states:
                if self.transmat[s,t] > 0.01:# and t in added_nodes:
                    # if t not in added_nodes and s == t:
                    #     continue
         
    
            
                    lab = str(round(self.transmat[s,t],2))
                    new_edge = pydot.Edge(dst=isotype_mapping[t], src=isotype_mapping[s], color="black", label=lab)
                
                    self.graph.add_edge(new_edge)
        
    def save(self, fname):
        self.graph.write_png(fname)
    

# tmat = np.loadtxt("/scratch/projects/tribal/real_data/test/transmat.txt")
# state_probs= np.loadtxt("/scratch/projects/tribal/real_data/test/state_probs.txt")
# encoding=  "/scratch/projects/tribal/real_data/mouse_isotype_encoding.txt"
# with open(encoding, 'r+') as file:
#             isotype_encoding = {}
#             counter = 0
#             for line in file:
#                 isotype_encoding[counter] = line.strip()
#                 counter += 1

# ds = DrawStateDiag(tmat, state_probs, isotype_encoding)
# ds.save("/scratch/projects/tribal/real_data/test/transmat.png")