


#!/usr/bin/env python3.7

# Copyright 2022, Gurobi Optimization, LLC

# Solve a multi-commodity flow problem.  Two products ('Pencils' and 'Pens')
# are produced in 2 cities ('Detroit' and 'Denver') and must be sent to
# warehouses in 3 cities ('Boston', 'New York', and 'Seattle') to
# satisfy supply/demand ('inflow[h,i]').
#
# Flows on the transportation network must respect arc capacity constraints
# ('capacity[i,j]'). The objective is to minimize the sum of the arc
# transportation costs ('cost[i,j]').

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from lineage_tree import LineageTree
from utils import hamming_distance
from collections import Counter
from dataclasses import dataclass 
import numpy as np
import matplotlib.pyplot as plt

class SteinerTree:
    def __init__(self, G, seq_weights, iso_weights, outdegree_bound, mut_exc_list=[], root=0, lamb=0.9, threads=3) -> None:
        self.m = gp.Model('steiner')
        self.m.Params.LogToConsole = 0
        self.m.Params.Threads = 3
        
        self.terminals = [v for v in G if G.out_degree[v]==0]
        self.nodes = list(G.nodes)
        self.root= root
        self.internal_nodes = [n for n in self.nodes if n not in self.terminals and n!= self.root ]
        self.edges = list(G.edges)
        self.flow_dest = [(t,i,j) for i,j in self.edges for t in self.terminals]
        #create a tuple of terminals and edges 
        self.f = self.m.addVars(self.flow_dest, name='flow', lb=0.0)
        self.x = self.m.addVars(self.edges, vtype=GRB.BINARY, name="edges")
        self.z = self.m.addVars(self.internal_nodes, vtype=GRB.BINARY, name="node use")
        self.iso_weights = iso_weights
        self.seq_weights = seq_weights

        self.in_nodes = {j : [] for j in self.nodes}
        self.out_nodes = {j : [] for j in self.nodes}
        for i,j in self.edges:
                self.in_nodes[j].append(i)
            
                self.out_nodes[i].append(j)

        # self.m.setObjective(lamb* sum(self.seq_weights[i,j]*self.x[i,j] for i,j in self.edges) + 
        #                     (1-lamb)* sum(self.iso_weights[i,j]*self.x[i,j] for [i,j] in self.edges) , GRB.MINIMIZE)

        self.m.setObjective( sum(self.iso_weights[i,j]*self.x[i,j] for [i,j] in self.edges) , GRB.MINIMIZE)


        for i in self.internal_nodes:
            self.m.addConstr(sum(self.x[i,j] for j in self.out_nodes[i]) >= 2*self.z[i])
        
            self.m.addConstr(1e5*self.z[i] >= sum(self.x[i,j] for j in self.out_nodes[i]))

        #ensure only 1 node is selected when multiple candidate isotypes exist
        for node_list in mut_exc_list:
            cand_edges_in = []
            for n in node_list:
                predeccessors = G.predecessors(n)
                for p in predeccessors:
                     cand_edges_in.append((p,n))

                     
            self.m.addConstr(sum(self.x[i,j] for i,j in cand_edges_in )<=1)
    

        for n,val in outdegree_bound.items():
                self.m.addConstr(sum(self.x[n,j] for j in self.out_nodes[n]) <= val -1)
        
        # for i in self.internal_nodes:
        #      self.m.addConstr(sum(self.x[i,j] for j in self.out_nodes[i]) - sum(self.x[j,i] for j in self.in_nodes[i]) >=1)

        #if edge is selected, then there must be flow to at least on terminal through that edge
        for i,j in self.edges:
            self.m.addConstr(sum(self.f[t,i,j] for t in self.terminals) >= self.x[i,j])

        for t in self.terminals:

            self.m.addConstrs(
                (self.f[t, i, j] <= self.x[i,j] for i,j in self.edges), "flow upper bound")
            for v in self.internal_nodes:
                self.m.addConstr(sum(self.f[t,i,v] for i in self.in_nodes[v])== sum(self.f[t,v,j] for j in self.out_nodes[v]), "flow conservation") #flow conservation
            self.m.addConstr(sum(self.f[t,i,t] for i in self.in_nodes[t])==1)
            self.m.addConstr(sum(self.f[t,t,j] for j in self.out_nodes[t])==0)
            self.m.addConstr(sum(self.f[t,self.root,j] for j in self.out_nodes[self.root])==1)
 
    
    
    def run(self):
            T = nx.DiGraph()
            self.m.optimize()
            if self.m.Status == GRB.OPTIMAL:
                solution = self.m.getAttr('X', self.x)
                flow = self.m.getAttr('X', self.f)
            
                score = self.m.objVal
        
          

                # for i,j in self.special_arcs:
                #     print('%s -> %s: %g' % (i, j, polytomies[ i, j]))
                
                for i,j in self.edges:
                    if solution[i,j] > 0:
                        T.add_edge(i,j)
                        # for t in self.terminals:
                        #     if flow[t,i,j] > 0:
                    #    print(f"edge {i} -> {j}")
            
                    #    break

            else:
                 print("Warning model infeasible")
                 score = np.Inf
    
            # for t in self.terminals:
            #      for i,j in self.edges:
            #           if flow[t,i,j] > 0:
            #                print(f"terminal: {t} flow {i} -> {j}: {flow[t,i,j]}")
                      
            
            return score, T
    

                # print(states)
       
                # print("used edges in network")
                # for i, j in self.all_arcs:
                #         if solution[ i, j] > 0:
                #             print('%s -> %s: %g' % (i, j, solution[ i, j]))
                
                # print("polytomies used")
            #     for i,j in self.special_arcs:
            #         if polytomies[i,j] > 0:
            #             # print('%s -> %s' % (i, j))
            #             t= int(i.split("_")[1])
                        
            #             for k in self.in_nodes[i]:
            #                 if solution[k,i] > 0:
                              
            #                     child =k.split("_")[0]
            #                     if t in poly_map:
            #                         poly_map[t].append(child)
            #                     else:
            #                         poly_map[t] = [child]
                                         
            # return score, states, poly_map

        

        
# G = nx.DiGraph()       
# G.add_nodes_from(['r', 'd', 'e', 'f', 'g', 'a', 'b', 'c'])
# G.add_edges_from([('r', 'd'), ('r', 'e'), ('d', 'a'), ('d','b'), ('d', 'c'), ('d','g'), ('g', 'a'), ('g', 'b'), ('e','a'), ('e', 'b'), ('e', 'c'), ('r','f') ])
# labels = {'r': 0, 'f': 2, 'a': 1, 'b': 1, 'c': 2, 'd': 0, 'e': 1, 'g': 1 }
# iso_weights = {}
# seq_weights = {}

# costs = {(0,0): 0.2, (0,1): 1, (0,2): 50, (1,1): 0.01, (1,2):  0.1 }
# for i,j in G.edges:
#     if labels[i] ==labels[j]:
#         iso_weights[i,j] =1
#     else:
#         iso_weights[i,j]=2
# # for i,j in G.edges:
# #     iso_weights[i,j] =costs[(labels[i],labels[j])]
#     seq_weights[i,j]=0
# print(iso_weights)
# score = SteinerTree(G,seq_weights, iso_weights, root='r', lamb=0).run()
# print(score)


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
        # print(f"Building graph for tree {id}....")
        # if id ==8:
        #      print('here')
        G = nx.DiGraph()
        seq_weights = {}
        iso_weights = {}
        max_outdegree = {}
        min_state = {}
        nodes_to_states = {}
        postorder =  list(nx.dfs_postorder_nodes(T, source=root))
        mut_exc_list =[]
        for n in postorder:
            # if n == self.root_identifier:
            #      print('here')
            node_exc_list = []
            if T.out_degree[n] ==0:
                iso = self.iso_labs[n]

                new_node = name_node(id,n,iso, is_leaf=T)
                G.add_node(new_node)
                self.node_isotypes[new_node] = iso
                min_state[n] = self.iso_labs[n]
                nodes_to_states[n] = [self.iso_labs[n]]
                  
            elif T.out_degree[n] > 0: 
                poly_nodes_added = []
                degree = T.out_degree[n]
                min_state[n] = min([min_state[v] for v in T.neighbors(n)])
                nodes_to_states[n] =[]

                for i in range(0, min_state[n]+1):
              
                    if i > 0 and n in self.iso_labs:
                         break
                    if self.root_identifier == n:
                         new_node = self.root_identifier
                    else:
                        new_node = name_node(id,n,i)
                    nodes_to_states[n].append(i)

            
                    G.add_node( new_node)
                    node_exc_list.append(new_node)

                    self.node_isotypes[new_node] = i
                    for v in T.neighbors(n):
                        is_leaf=T.out_degree[v]==0
                        for j in nodes_to_states[v]:
                            if i <= j:
                                G.add_edge(new_node, name_node(id,v,j, is_leaf=is_leaf))
                                seq_weights[new_node,name_node(id,v,j,is_leaf=is_leaf)] = hamming_distance(seq_labs[n], seq_labs[v])
                                if self.iso_costs[i,j] == np.Inf:
                                        
                                        print(f"warning, impossible edge i: {i} j: {j} cost: {self.iso_costs[i,j]}")
                                # iso_weights[new_node, name_node(id,v,j,is_leaf=is_leaf)] = self.iso_costs[i,j]
                
                #insert polytomy nodes and add edges to Graph
                if degree > 2:
                   
                    poly_node_states = []
                    for v in T.neighbors(n):
                         poly_node_states = poly_node_states + nodes_to_states[v]
                    max_node = max(poly_node_states)
                    node_counts = Counter(poly_node_states)

                    
                    for j in range(1, max_node+1):
                            # if node_counts[j] > 1:
                                p_node = name_node(id,n,j,True)
                                max_outdegree[p_node] = degree
                                self.node_isotypes[p_node] = j
                                G.add_node(p_node)
                                poly_nodes_added.append(p_node)
                                for i in nodes_to_states[n]:
                                    if i < j:
                                        G.add_edge(name_node(id,n,i), name_node(id,n,j,True))
                                        seq_weights[name_node(id,n,i), name_node(id,n,j,True)] = 0
                                        if self.iso_costs[i,j] == np.Inf:
                                                print(f"warning, impossible edge i: {i} j: {j} cost: {self.iso_costs[i,j]}")

                                        # iso_weights[name_node(id,n,i), name_node(id,n,j,True)] = self.iso_costs[i,j]
                                
                                for k in range(1, j):
                                     if name_node(id,n,k,True) in G.nodes:
                                          G.add_edge(name_node(id,n,k,True), p_node)
                                          seq_weights[name_node(id,n,k,True),p_node] = 0

                                
                                for v in T.neighbors(n):
                                        is_leaf=T.out_degree[v]==0
                                        for s in nodes_to_states[v]:
                                             if j <= s:
                                                G.add_edge(name_node(id,n,j,True), name_node(id,v,s, is_leaf=is_leaf))
                                                seq_weights[name_node(id,n,j,True),name_node(id,v,s,is_leaf=is_leaf)] = hamming_distance(seq_labs[n], seq_labs[v])
                                                iso_weights[name_node(id,n,j,True), name_node(id,v,s,is_leaf=is_leaf)] = self.iso_costs[j,j]

                                        # if j in nodes_to_states[v]:                           
                                        #     G.add_edge(name_node(id,n,j,True), name_node(id,v,j, is_leaf=is_leaf))
                                        #     seq_weights[name_node(id,n,j,True),name_node(id,v,j,is_leaf=is_leaf)] = hamming_distance(seq_labs[n], seq_labs[v])
                                        #     # iso_weights[name_node(id,n,j,True), name_node(id,v,j,is_leaf=is_leaf)] = self.iso_costs[j,j]

            if len(node_exc_list) >0:
                 mut_exc_list.append(node_exc_list)   
        
        #turn back on
        #post-process graph to remove potential unifurcations
        bad_nodes = []
        while True:
            for n in G.nodes:
            
                if n == self.root_identifier or "Germline" in n:
                    continue
                if len(list(G.neighbors(n))) == 1:
                    bad_nodes.append(n)
            if len(bad_nodes) ==0:
                 break
            for b in bad_nodes:
                G.remove_node(b)
                del self.node_isotypes[b]
                del max_outdegree[b]
            bad_nodes =[]
          
        for u,v in G.edges:
            s = self.node_isotypes[u]
            t =  self.node_isotypes[v]
            iso_weights[u,v] = self.iso_costs[int(s),int(t)]
        
        fg = FlowGraph(id, G, seq_weights, iso_weights, max_outdegree, self.node_isotypes, mut_exc_list)
        self.Graphs.append(fg)
        return fg
    
    def combineGraphs(self):
         graphs = [fg.G for fg in self.Graphs]
         combined_graph = nx.compose_all(graphs)
         seq_weights = {}
         iso_weights = {}
         for fg in self.Graphs:
              seq_weights.update(fg.seq_weights)
              iso_weights.update(fg.iso_weights)
        
         return FlowGraph(0, combined_graph, seq_weights, iso_weights)
        
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
     degree_bound: dict
     isotypes: dict
     mut_exclusive_list: list


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
    # Draw pygraphviz graph and save as PDF
        pgv_graph.draw(fname, prog="dot", format=ext)

        # Save as PDF
        # plt.savefig(fname)

                         
                     
# T = nx.DiGraph()

# T.add_edges_from([ ('naive', 'r') ,('r', 'a') , ('r','b'), ('r', 'c'), ('r', 'd'), ('r','e')])
# #iso_costs = {(0,0): 0.91, (0,1): 1.6, (0,2): 1.6, (0,3): 1.6, (1,1): 0.11, (1,2):  3, (1,3): 3, (2,2): 0.51, (2,3): 0.43, (3,3): 0  }
# sequences = {n: ['a,a'] for n in T.nodes}
# iso_costs = {(0,0): 0.55, (0,1): 0.44, (0,2): 0.005, (0,3): 0.005,
#               (1,1): 0.55, (1,2):  0.44, (1,3): 0.05, 
#               (2,2): 0.55, (2,3): 0.44, 
#               (3,3): 1  }
# iso_costs_log = {}
# for key in iso_costs:
#      iso_costs[key] = -1*np.log(iso_costs[key])


# isotypes  = {'naive': 0, 'a': 1, 'b':1, 'c': 2, 'd':3, 'e':3}

# lt1 =LineageTree(T,'naive', 0)
# # # lt2 =LineageTree(T,'root', 1)
# lt1.save_png("test/start_tree.png", isotypes)
# cg = ConstructGraph(iso_costs, isotypes, root_identifier="naive")

# fg = cg.build(lt1, sequences) 
# fg.save_graph("test/G1.png")
# score, T2= SteinerTree( fg.G, fg.seq_weights, fg.iso_weights, fg.degree_bound, fg.mut_exclusive_list, root="naive").run()

# lt2 = LineageTree(T2, "naive", 0)
# print(f"score: {score}")
# lt2.save_png("test/out_tree3.png", fg.isotypes)

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

                    
             


                
        
    

                     
                 
                 
    
                  

          


# class PolytomyResolver:
#     def __init__(self, tree_nodes, weights, scores, states, candidate_state):
#         # Create optimization model
#         self.m = gp.Model('netflow')
#         self.tree_nodes = tree_nodes
#         self.nodes = ['source', 'sink'] + [n for n in tree_nodes]
#         self.states = [str(s) for s in states] 
#         self.score_nodes = []
#         self.capacity_upper = {}
#         self.capacity_lower = {}
#         self.fixed_costs = {}
#         self.costs = {}
#         for n in tree_nodes:
#             self.capacity_upper[("source", n )] = 1
#             self.capacity_lower[("source", n )] = 0
#             self.costs[("source", n)] = 0
#             for state in scores[n]:
#                 node_name = f"{n}_{state}"
#                 self.score_nodes.append(node_name)
#                 self.nodes.append(node_name)
       
        
#                 self.capacity_upper[(n,node_name)] =1
#                 self.capacity_lower[(n,node_name)] =0
#                 self.costs[(n,node_name)] = scores[n][state]

        
#                 self.capacity_upper[(node_name,"sink")] = 1
#                 self.capacity_lower[(node_name,"sink")] = 0
#                 self.costs[(node_name, "sink")] = weights[candidate_state, int(state)]
#                 for s in self.states:
#                     if int(s) <= state and int(s) != candidate_state:
                 
#                         self.capacity_upper[(node_name,f"state_{s}")] = 1
#                         self.capacity_lower[(node_name,f"state_{s}")] = 0
             
#                         self.costs[(node_name,f"state_{s}")] = weights[int(s), int(state)]
                
#         self.special_arcs = []
#         #TODO: need to add n-1 layers of states
#         for s in self.states:
#             if int(s) >= candidate_state:
#                 node_name = f"state_{s}"
#                 self.special_arcs.append((node_name, "sink"))
       
#                 self.nodes.append(node_name)
#                 self.capacity_upper[(node_name,"sink")] = len(tree_nodes)-1
#                 self.capacity_lower[(node_name,"sink")] =  2 #0
#                 self.fixed_costs[(node_name, "sink")] = weights[candidate_state, int(s)]
        
   
#         self.all_arcs, self.cap_upper = gp.multidict(self.capacity_upper)

#         # self.inflow = {a: 0 for a in self.all_arcs}
#         # self.inflow[('edges','sink')] = len(tree_nodes)
#         # self.inflow['edge','source'] = -1*len(tree_nodes)

#         self.normal_arcs = [i for i in self.all_arcs if i not in self.special_arcs]

#         # self.flow = self.m.addVars(self.all_arcs, lb=self.capacity_lower, ub=self.capacity_upper, vtype=GRB.INTEGER, name='flow')
#         self.flow = self.m.addVars(self.all_arcs, vtype=GRB.INTEGER, name='flow')

#         self.polytomy = self.m.addVars(self.special_arcs, vtype=GRB.BINARY, name='polytomy')



#         # Arc-capacity constraints
#         self.m.addConstrs(
#             (self.flow[ i, j] <= self.capacity_upper[i, j] for i, j in self.normal_arcs), "cap_upper")

#         self.m.addConstrs(
#             (self.flow[ i, j] >= self.capacity_lower[i, j] for i, j in self.normal_arcs), "cap_lower")
        
#         #polytomy is opened only if flow is greater than lower bound of capacity
#         self.m.addConstrs(
#             (self.flow[ i, j]>= self.capacity_lower[i, j]*self.polytomy[i,j] for i, j in self.special_arcs), "cap_special_lower")
        
#         self.m.addConstrs(
#             (self.flow[ i, j] <= self.capacity_upper[i, j]*self.polytomy[i,j] for i, j in self.special_arcs), "cap_special_upper")
        
#         #Flow conservation 
#         # m.addConstr(quicksum(x[i,j] for j in (set(V) - set(S))) >= 2)
#         self.in_nodes = {j : [] for j in self.nodes}
#         self.out_nodes = {j : [] for j in self.nodes}
#         for j in self.nodes:
#             for (i,k) in self.all_arcs:
#                 if k==j:
#                     self.in_nodes[j].append(i)
#                 if i == j:
#                     self.out_nodes[j].append(k)
        
#         self.m.addConstr(sum(self.flow['source',j] for j in self.out_nodes['source']) == len(tree_nodes))


#         #flow leaving the source cannot exceed the number of tree nodes 
#         self.m.addConstr(sum(self.flow[j, 'sink'] for j in self.in_nodes['sink']) == len(tree_nodes))
#         for  j in self.nodes:
#             if j not in ['source', 'sink']:
#                 self.m.addConstr(sum(self.flow[i, j]  for i in self.in_nodes[j]) == sum(self.flow.sum(j, i) for i in self.out_nodes[j]))
        
#         #flow leaving the source cannot exceed the number of tree nodes 
   

#         self.m.setObjective(sum(self.costs[i,j]*self.flow[i,j] for i,j in self.normal_arcs) + 
#                         sum(self.fixed_costs[i,j]*self.polytomy[i,j] for i,j in self.special_arcs)  , GRB.MINIMIZE)
        

#     def run(self):
#             states = {}
#             poly_map = {}
#             score= None
#             self.m.optimize()
#             if self.m.Status == GRB.OPTIMAL:
#                 solution = self.m.getAttr('X', self.flow)
#                 polytomies = self.m.getAttr('X', self.polytomy)

#                 # for i,j in self.special_arcs:
#                 #     print('%s -> %s: %g' % (i, j, polytomies[ i, j]))
                
#                 for i in self.tree_nodes:
#                     for j in self.out_nodes[i]:
#                         if solution[i,j] > 0:
#                             states[i] = int(j.split("_")[1])

#                 # print(states)
#                 score = self.m.objVal
              
#                 # print("used edges in network")
#                 # for i, j in self.all_arcs:
#                 #         if solution[ i, j] > 0:
#                 #             print('%s -> %s: %g' % (i, j, solution[ i, j]))
                
#                 # print("polytomies used")
#                 for i,j in self.special_arcs:
#                     if polytomies[i,j] > 0:
#                         # print('%s -> %s' % (i, j))
#                         t= int(i.split("_")[1])
                        
#                         for k in self.in_nodes[i]:
#                             if solution[k,i] > 0:
                              
#                                 child =k.split("_")[0]
#                                 if t in poly_map:
#                                     poly_map[t].append(child)
#                                 else:
#                                     poly_map[t] = [child]
                                         
#             return score, states, poly_map




        
    
            
            


# tree_nodes = ['a', 'b', 'c']
# states = [0,1]
# weights = {(0, 0): 0, (0, 1): 1, (1,0): 0, (1,1): 0}
# scores= {'a': {0: 3, 1: 4}, 'b': {0: 2, 1: 5}, 'c': {0: 2, 1: 11}}
# candidate_state = states[0]


# tree = nx.DiGraph()

# tree.add_edges_from([(0,1), (0,2), (0,3), (1,4), (1,5), (2,6), (2,7), (3,8), (3,9)])

# tree_nodes = [1,2,3]
# candidate_state = 0
# scores = {1: {0: 2, 1: 1}, 2: {0:2, 1:0}, 3: {0: 2, 1: 2, 2: 0}}

# states = [0,1,2]
# # isotypes = {4:1, 5:2, 6:1, 7:2, 8:1, 9:2}

# weights = {}
# for s in states:
#     for t in states:
#         if s > t:
#             weights[s,t] = 1000000
#         elif s==t:
#             weights[s,t] =0
#         else:
#             weights[s,t] =1




# pr = PolytomyResolver(tree_nodes, weights, scores, states, candidate_state)
# pr.run()

       
        
