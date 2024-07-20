import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from dataclasses import dataclass 
import numpy as np
from itertools import combinations

class SteinerTree:
    '''
    Given an instance of the Steiner Minimal Tree Problem (G,S,c ) where G=(V,E), S is a subset of V and 
    c is an edge weight mapping, finds a Steiner Tree in G where the sum of the edge weights is minimum
    c = lambda* seq_weights + (1-lambda)*iso_weights 
    '''
    def __init__(self, G,T, S,  seq_weights, iso_weights,node_mapping, 
                 tree_to_graph, tree_out_degree, pars_score=None, root=0, lamb=0.9, threads=3) -> None:
        
        #set solve parameters 
        self.m = gp.Model('steiner')
        self.m.Params.LogToConsole = 0
        self.m.Params.Threads = threads


        self.terminals = S
        self.root= root
        self.nodes = list(G.nodes)
        self.edges = list(G.edges)
        self.seq_weights = seq_weights
        self.internal_nodes = [n for n in self.nodes if n not in self.terminals and n!= self.root ]
        self.in_nodes = {j : [] for j in self.nodes}
        self.out_nodes = {j : [] for j in self.nodes}
        for i,j in self.edges:
                self.in_nodes[j].append(i)
            
                self.out_nodes[i].append(j)
        
        #compute edge weights from the sequence and isotype weights 
    
        self.c = {e:iso_weights[e]for e in self.edges}


        #create a tuple of terminals and edges 
        self.flow_dest = [(t,i,j) for i,j in self.edges for t in self.terminals]

       
             
        
        #add a continuous flow variable 
        self.f = self.m.addVars(self.flow_dest, name='flow', lb=0.0)

        #add a binary variable indicating if edge i,j is included in the Steiner tree 
        self.x = self.m.addVars(self.edges, vtype=GRB.BINARY, name="edges")

        
        for u,v in T.edges:
             self.m.addConstr(sum(self.x[i,j] for i in tree_to_graph[u] for j in tree_to_graph[v] if (i,j) in self.edges ) ==1)
        #we need this variable if we want to prevent unifurcations 
        # self.z = self.m.addVars(self.internal_nodes, vtype=GRB.BINARY, name="node use")  

         #minimize the sum of the edge weights of the Steiner Tree in G
        self.m.setObjective( sum(self.c[i,j]*self.x[i,j] for [i,j] in self.edges) , GRB.MINIMIZE)

             

        # need to send 1 unit of flow from the root to each terminal
        for t in self.terminals:

            #enforce that a node only has flow if edge i,j is included in the minimal Steiner tree
            self.m.addConstrs(
                (self.f[t, i, j] <= self.x[i,j] for i,j in self.edges), "flow upper bound")
            
            #flow conversvation on the internal nodes
            for v in self.internal_nodes:
                self.m.addConstr(sum(self.f[t,i,v] for i in self.in_nodes[v])==  \
                                 sum(self.f[t,v,j] for j in self.out_nodes[v]), "flow conservation") 
            
            #ensures 1 unit of flow reaches each termminal
            self.m.addConstr(sum(self.f[t,i,t] for i in self.in_nodes[t])==1)
            
            self.m.addConstr(sum(self.f[t,t,j] for j in self.out_nodes[t])==0)
            
            #ensure 1 unit of flow designated for each terminal leaves the root
            self.m.addConstr(sum(self.f[t,self.root,j] for j in self.out_nodes[self.root])==1)

        

    def post_process(self, T):
    #  return T
        '''
        Remove any unifurcations that have 0 branch length
        '''
        unifurcations = [n for n in T if T.out_degree[n]==1 and n != self.root]

        to_remove = []
        for u in unifurcations:
            
            child = list(T.neighbors(u))[0]

            if self.c[u,child] ==0: # and self.seq_weights[u,child]==0:
                to_remove.append((u,child))
        
        for u,v in to_remove:
              parent = list(T.predecessors(u))
              if len(parent) > 0:
                if u.split("_")[0] ==parent[0].split("_")[0]:
                    parent = parent[0]
                    T.remove_node(u)
                    T.add_edge(parent, v)


        return T
                   
 
    
    
    def run(self):
            T = nx.DiGraph()
            self.m.optimize()
            if self.m.Status == GRB.OPTIMAL:
                solution = self.m.getAttr('X', self.x)
                flow = self.m.getAttr('X', self.f)
                score = self.m.objVal
   
                for i,j in self.edges:
                    total_flow = sum(flow[t,i,j] for t in self.terminals)
                    if solution[i,j] > 0.5 and total_flow >0:
                        T.add_edge(i,j)
            

            else:
                 print("Warning model infeasible")
                 score = np.Inf
       
     
            T = self.post_process(T)
            return score, T
    


def name_node(tree, node, label,  is_poly=False, is_leaf=False): 

        lab = str(node) + "_" + str(label)

        if is_poly:
             lab += "_p"
        return  lab

def decode_node(node, is_leaf=False):
        codes =node.split("_")

        if is_leaf:
             name = codes[0]
        else:
             name = node

        if codes[-1] == "p":
             isotype = codes[-2]
        else:
             isotype = codes[-1]
        return int(isotype), name  


#takes in a tree and constructs a network 
class ConstructGraph:
    def __init__(self, iso_costs, isotypes, root_identifier="root") -> None:
        
        self.Graphs = []
        self.iso_costs = iso_costs

        self.iso_labs = isotypes
        self.root_identifier = root_identifier
        self.node_isotypes = {}

        # for i,j in self.iso_costs:
        #      if j >= i:
        #           if np.isinf(self.iso_costs[i,j]):
        #                print(self.iso_costs[i,j])
        #                print(self.iso_costs)
        #           assert(not np.isinf(self.iso_costs[i,j]))
                  

    def build(self, LinTree ):

        T = LinTree.T
        root = LinTree.root
        id = LinTree.id
      
        G = nx.DiGraph()
        seq_weights = {}
        iso_weights = {}

        nodes_to_states = {}
        postorder =  list(nx.dfs_postorder_nodes(T, source=root))
    
        out_degree_max = {}
        leafset = []
        node_mapping ={}
        tree_to_graph = {n: [] for n in T.nodes}
        node_out_degree = {n: 0 for n in T.nodes}
        for n in postorder:
            node_out_degree[n] =T.out_degree[n]
           
            node_extensions= []

            if T.out_degree[n] ==0:
                iso = self.iso_labs[n]

                new_node = name_node(id,n,iso, is_leaf=T)
                G.add_node(new_node)
                node_mapping[new_node] = n
                tree_to_graph[n].append(new_node)
              

                out_degree_max[new_node] = T.out_degree[n]
                node_out_degree
                leafset.append(new_node)
                self.node_isotypes[new_node] = iso
   
                nodes_to_states[n] = [self.iso_labs[n]]
                  
            elif T.out_degree[n] > 0: 
                desc = LinTree.find_leaf_descendants(n,T)
                upper_bound = max(self.iso_labs[d] for d in desc)
                nodes_to_states[n] =[]
                for i in range(0, upper_bound+1):
                    if i > 0 and n in self.iso_labs:
                         break
                    if self.root_identifier == n:
                         new_node = self.root_identifier
                    else:
                        new_node = name_node(id,n,i)
                    nodes_to_states[n].append(i)

            
                    G.add_node( new_node)
                    node_mapping[new_node] = n
                    tree_to_graph[n].append(new_node)
                    
                    out_degree_max[new_node] = T.out_degree[n]
                    node_extensions.append(new_node)

                    self.node_isotypes[new_node] = i
                    for v in T.neighbors(n):
                        is_leaf=T.out_degree[v]==0
                        for j in nodes_to_states[v]:
                            if i <= j: #and T.out_degree[n] > 2 and  self.iso_costs[int(i),int(j)] > 0:
                                G.add_edge(new_node, name_node(id,v,j, is_leaf=is_leaf))
                                seq_weights[new_node,name_node(id,v,j,is_leaf=is_leaf)] = 0

                                # seq_weights[new_node,name_node(id,v,j,is_leaf=is_leaf)] = hamming_distance(seq_labs[n], seq_labs[v])
                                if self.iso_costs[i,j] == np.Inf:
                                        
                                        print(f"warning, impossible edge i: {i} j: {j} cost: {self.iso_costs[i,j]}")
                         
                #sequence weights are always 0 because it's the name sequence
            # if n != self.root_identifier and list(T.predecessors(n))[0] != self.root_identifier:
                for u,v in combinations(node_extensions, 2):
                        if self.node_isotypes[u] < self.node_isotypes[v]:
                            G.add_edge(u,v)
                            seq_weights[u,v] = 0
                        elif self.node_isotypes[v] < self.node_isotypes[u]:
                            G.add_edge(v,u)
                            seq_weights[v,u] = 0
            

        for u,v in G.edges:
            s = self.node_isotypes[u]
            t =  self.node_isotypes[v]
            iso_weights[u,v] = self.iso_costs[int(s),int(t)]

        
        fg = FlowGraph(id, G, seq_weights, iso_weights, self.node_isotypes, node_mapping, tree_to_graph, node_out_degree)
   
        self.Graphs.append(fg)
        return fg
    
        
    def decodeTree(self, tree, root_id="naive"):
         isotypes = {}
         #old as keys, new names as values
         relabeling = {}
         for n in tree.nodes:
            if n == root_id:
                name = root_id 
                iso = 0
        
            else:
                is_leaf = tree.out_degree[n]==0
                iso, name = decode_node(n, is_leaf)
            relabeling[n] = name
            isotypes[name] = iso
         tree= nx.relabel_nodes(tree, relabeling)
         return tree, isotypes
              
         


@dataclass
class FlowGraph:
     id: int 
     G: nx.DiGraph
     seq_weights: dict 
     iso_weights: dict
     isotypes: dict
     node_mapping: dict 
     tree_to_graph: dict 
     node_out_degree: dict 



     def find_terminals(self):
          return  [v for v in self.G if self.G.out_degree[v]==0]
     
     def save_graph(self, fname):
        color_encoding =  {
                0 : "#f0f0f0",
                1 : "#FFEDA0",
                2 : "#FD8D3C",
                3 : "#E31A1C",
                4 : "#800026",
                5 : "mediumseagreen",
                6 : "#74C476",
                7 : "#6A51A3",
                8 : "darkgoldenrod",
                9 : "thistle1"
            }
        # Draw the graph
        if '.pdf' in fname:
             ext = 'pdf'
        else:
            ext ='png'
        pgv_graph = nx.nx_agraph.to_agraph(self.G)


        node_colors = {n: color_encoding[val] for n,val in self.isotypes.items()}

        for node in pgv_graph.nodes():
            
            pgv_graph.get_node(node).attr["fillcolor"] = node_colors[node]
            pgv_graph.get_node(node).attr["style"] = "filled"
            if node != "naive":
                pgv_graph.get_node(node).attr["label"] =  node.split("_")[0]
            else:
                pgv_graph.get_node(node).attr["label"] = "r"

        pgv_graph.draw(fname, prog="dot", format=ext)

 
                         
                     
                    
             


                
        
    

                     
                 
                 
    
                  

          


