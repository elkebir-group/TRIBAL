import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from utils import hamming_distance
from lineage_tree import LineageTree
from dataclasses import dataclass 
import numpy as np
from itertools import combinations

class SteinerTree:
    '''
    Given an instance of the Steiner Minimal Tree Problem (G,S,c ) where G=(V,E), S is a subset of V and 
    c is an edge weight mapping, finds a Steiner Tree in G where the sum of the edge weights is minimum
    c = lambda* seq_weights + (1-lambda)*iso_weights 
    '''
    def __init__(self, G, S, seq_weights, iso_weights,node_mapping, tree_to_graph, tree_out_degree, root=0, lamb=0.9, threads=3) -> None:
        
        #set solve parameters 
        self.m = gp.Model('steiner')
        self.m.Params.LogToConsole = 0
        self.m.Params.Threads = threads


        self.terminals = S
        self.root= root
        self.nodes = list(G.nodes)
        self.edges = list(G.edges)
        
        self.internal_nodes = [n for n in self.nodes if n not in self.terminals and n!= self.root ]
        self.in_nodes = {j : [] for j in self.nodes}
        self.out_nodes = {j : [] for j in self.nodes}
        for i,j in self.edges:
                self.in_nodes[j].append(i)
            
                self.out_nodes[i].append(j)
        
        #compute edge weights from the sequence and isotype weights 
        self.c = {e:lamb* seq_weights[e] + (1-lamb) *iso_weights[e]for e in self.edges}


        #create a tuple of terminals and edges 
        self.flow_dest = [(t,i,j) for i,j in self.edges for t in self.terminals]
        
        #add a continuous flow variable 
        self.f = self.m.addVars(self.flow_dest, name='flow', lb=0.0)

        #add a binary variable indicating if edge i,j is included in the Steiner tree 
        self.x = self.m.addVars(self.edges, vtype=GRB.BINARY, name="edges")

        #we need this variable if we want to prevent unifurcations 
        # self.z = self.m.addVars(self.internal_nodes, vtype=GRB.BINARY, name="node use")  

         #minimize the sum of the edge weights of the Steiner Tree in G
        self.m.setObjective( sum(self.c[i,j]*self.x[i,j] for [i,j] in self.edges) , GRB.MINIMIZE)


        #enforce that there is a single MRCA child of the root node 
        self.m.addConstr(sum(self.x[self.root,j] for j in self.out_nodes[self.root] ) <= 1)

        for n, deg in tree_out_degree.items():
             nodes = tree_to_graph[n]
             if deg ==1:
                  continue
             self.m.addConstr(sum(self.x[i,j] for i in nodes for j in self.out_nodes[i]) \
                               - sum(self.x[i,j] for j in nodes for i in self.in_nodes[j]) <= deg -1)
             

        
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
            
            # self.m.addConstr(sum(self.f[t,t,j] for j in self.out_nodes[t])==0)
            
            #ensure 1 unit of flow designated for each terminal leaves the root
            self.m.addConstr(sum(self.f[t,self.root,j] for j in self.out_nodes[self.root])==1)

        ## constraint to prevent unifurcations in internal nodes 
        # for i in self.internal_nodes:
        # self.m.addConstr(sum(sum(self.f[t,i,j] for j in self.out_nodes[i]) for t in self.terminals) >= 2*self.z[i])
        # self.m.addConstr(1e5*self.z[i] >= sum(sum(self.f[t,i,j] for j in self.out_nodes[i]) for t in self.terminals))

        # for n,val in degree_max.items():
        #         if len(self.out_nodes[n]) > 0:
        #             self.m.addConstr(sum(self.x[n,j] for j in self.out_nodes[n]) <= val)

        #ensure only 1 node is selected when multiple candidate isotypes exist
        # for node_list in mut_exc_list:
        #     cand_edges_in = []
        #     for n in node_list:
        #         predeccessors = G.predecessors(n)
        #         for p in predeccessors:
        #              cand_edges_in.append((p,n))     
            # self.m.addConstr(sum(self.x[i,j] for i,j in cand_edges_in )<=1)
    

        

        
 
    
    
    def run(self):
            T = nx.DiGraph()
            self.m.optimize()
            if self.m.Status == GRB.OPTIMAL:
                solution = self.m.getAttr('X', self.x)
                flow = self.m.getAttr('X', self.f)
                score = self.m.objVal
   
                for i,j in self.edges:
               
                    if solution[i,j] > 0.5:
                        T.add_edge(i,j)
                    # else:
                         
                    #     for t in self.terminals:
                    #         if flow[t,i,j] > 0.0:
                    #            print(f"{t}:  edge {i} -> {j}")
            
                    #    break

            else:
                 print("Warning model infeasible")
                 score = np.Inf
    
            
            return score, T
    

#treeid_nodeid_isotype_polytomy
def name_node(tree, node, label,  is_poly=False, is_leaf=False): 

        lab = str(node) + "_" + str(label)
        # if not is_leaf:
        #     lab = str(tree) + "_" + lab
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

        for i,j in self.iso_costs:
             if j >= i:
                  if np.isinf(self.iso_costs[i,j]):
                       print(self.iso_costs[i,j])
                       print(self.iso_costs)
                  assert(not np.isinf(self.iso_costs[i,j]))
                  

    def build(self, LinTree, seq_labs ):

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
                            if i <= j:
                                G.add_edge(new_node, name_node(id,v,j, is_leaf=is_leaf))
                                seq_weights[new_node,name_node(id,v,j,is_leaf=is_leaf)] = hamming_distance(seq_labs[n], seq_labs[v])
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
    
    # def combineGraphs(self):
    #      graphs = [fg.G for fg in self.Graphs]
    #      combined_graph = nx.compose_all(graphs)
    #      seq_weights = {}
    #      iso_weights = {}
    #      for fg in self.Graphs:
    #           seq_weights.update(fg.seq_weights)
    #           iso_weights.update(fg.iso_weights)
        
    #      return FlowGraph(0, combined_graph, seq_weights, iso_weights)
        
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
        # Generate layout
        # pos = nx.spring_layout(self.G)

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

 
                         
                     
# T = nx.DiGraph()
# T.add_edges_from([(0, 8), (8,7), (7,1), (7,2), (7,3), (7,10), (4,5), (4,6), (8,4),(8,9)])
# isotypes  = {0: 0,  1:3, 2:0, 3:1, 5:2, 6:0, 9:3, 10:3}

# # T.add_edges_from([('naive', 'w'), ('w','x'), ('w','y'), ('w','z')]) # (0,5)
# # T= nx.relabel_nodes(T, {1:"f", 2:"b", 3:"c", 4:"d", 5:"e"})
# # isotypes  = {'naive': 0,  'x':3, 'y':1, 'z':1, 5:3}


# # T.add_edges_from([('naive', 0), (0,1), (1,2), (1,3), (1,4),(1,6), (0,5)])  # (0,5)
# # isotypes  = {'naive': 0, "b":1, "c":3, "d":1, "e":2}  #5:2
# iso_costs = {}
# for u in range(4):
#      for v in range(4):
#           if u > v:
#                continue
#           elif u == v:
#                 iso_costs[(u,v)] = 1
#           elif v == u + 1:
#                iso_costs[(u,v)] = 2
#                 # iso_costs[(u,v)] = 0
#           else:
#                iso_costs[(u,v)] = 3
#             #    iso_costs[(u,v)] = 10
#             #    iso_costs[(u,v)] = 0   

# # for key, val in iso_costs.items():
# #      iso_costs[key] = -1*np.log(val)
# # for key, val in iso_costs.items():
# #      print(f"{key}: {val}")
# ext = "png"
# sequences = {n: ['a,c'] for n in T.nodes}
# sequences = {5: ['T', 'T', 'A'], 6: ['T', 'T', 'T'], 1: ['C','C','G'], 10: ['C','C','G'], 2: ['C', 'C', 'G'], 3:['G', 'C', 'C'], 0: ['A', 'T', 'T'], 9:['G','T','T']}

# lt =LineageTree(T,0)
# lt.save_png(f"test/start_tree.{ext}", isotypes)
# pars_score, seq_labels =lt.sequence_parismony(sequences)

# cg = ConstructGraph(iso_costs, isotypes, root_identifier="naive")
# fg = cg.build(lt, seq_labels) 
# fg.save_graph(f"test/flow_graph.{ext}")

# score, T2= SteinerTree( fg.G, fg.seq_weights, fg.iso_weights, root="0_0", degree_max=fg.degree_max).run()
# lt2 = LineageTree(T2, "naive", 0)
# print(f"score: {score}")
# lt2.save_png(f"test/out_tree.{ext}", fg.isotypes, show_legend=False)
# print("DONE!")






# lt2 = LineageTree(T2, "naive", 0)
# print(f"score: {score}")
# lt2.save_png("test/out_tree3.png", fg.isotypes)


# T.add_edges_from([ ('naive', 'r') ,('r', 'a') , ('r','b'), ('r', 'c'), ('r', 'd'), ('r','e')])
# #iso_costs = {(0,0): 0.91, (0,1): 1.6, (0,2): 1.6, (0,3): 1.6, (1,1): 0.11, (1,2):  3, (1,3): 3, (2,2): 0.51, (2,3): 0.43, (3,3): 0  }

# iso_costs = {(0,0): 0.55, (0,1): 0.44, (0,2): 0.005, (0,3): 0.005,
#               (1,1): 0.55, (1,2):  0.44, (1,3): 0.05, 
#               (2,2): 0.55, (2,3): 0.44, 
#               (3,3): 1  }
# iso_costs_log = {}
# for key in iso_costs:
#      iso_costs[key] = -1*np.log(iso_costs[key])


# isotypes  = {'naive': 0, 'a': 1, 'b':1, 'c': 2, 'd':3, 'e':3}


# # # lt2 =LineageTree(T,'root', 1)




# score, T2= SteinerTree( fg.G, fg.seq_weights, fg.iso_weights, fg.degree_bound, fg.mut_exclusive_list, root="naive").run()



# T.add_edges_from([ ('root', 'r') , ('r','a'), ('r', 'b'), ('r', 'c'), ('r', 'd')])

# T.add_edges_from([ ('root', 'r') , ('r','a'), ('r', 'b'), ('r', 'c'), ('r', 'd'), ('r', 'e'), ('r', 'f'), ('r', 'g')])
# T.add_edges_from([ ('naive', 'r') ,('r', 'a') , ('r','b'), ('r', 'c'), ('c', 'd'), ('c', 'e'), ('c', 'f')])



# iso_costs = {}
# for i in range(10):
#      for j in range(10):
#           if j >= i:
#                if i ==j:
#                     iso_costs[(i,j)] =0
#                else:
#                     iso_costs[(i,j)] =1


# isotypes  = {'naive': 0, 'a': 2, 'b':2, 'd': 3, 'f':2, 'e': 3}








# _ = cg.build(lt1, sequences)
# fg = cg.combineGraphs()

# print(list(fg.G.nodes))
# print(list(fg.G.edges))



# score, T = SteinerTree(fg.G, fg.seq_weights, fg.iso_weights, root='root', lamb=0).run()
# print(score)
# print(list(T.edges))

                    
             


                
        
    

                     
                 
                 
    
                  

          


