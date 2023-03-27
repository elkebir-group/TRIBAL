


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


class SteinerTree:
    def __init__(self, G, seq_weights, iso_weights, root=0, lamb=0.5) -> None:
        self.m = gp.Model('steiner')
        self.terminals = [v for v in G if G.out_degree[v]==0]
        self.nodes = list(G.nodes)
        self.root= root
        self.internal_nodes = [n for n in self.nodes if n not in self.terminals and n!= self.root ]
        self.edges = list(G.edges)
        self.flow_dest = [(t,i,j) for i,j in self.edges for t in self.terminals]
        #create a tuple of terminals and edges 
        self.f = self.m.addVars(self.flow_dest, name='flow', lb=0.0)
        self.x = self.m.addVars(self.edges, vtype=GRB.BINARY, name="edges")
        self.iso_weights = iso_weights
        self.seq_weights = seq_weights

        self.in_nodes = {j : [] for j in self.nodes}
        self.out_nodes = {j : [] for j in self.nodes}
        for i,j in self.edges:
                self.in_nodes[j].append(i)
            
                self.out_nodes[i].append(j)

        self.m.setObjective(lamb* sum(self.seq_weights[i,j]*self.x[i,j] for i,j in self.edges) + 
                            (1-lamb)* sum(self.iso_weights[i,j]*self.x[i,j] for [i,j] in self.edges) , GRB.MINIMIZE)

        for t in self.terminals:
            self.m.addConstrs(
                (self.f[t, i, j] <= self.x[i,j] for i,j in self.edges), "flow upper bound")
            for v in self.internal_nodes:
                self.m.addConstr(sum(self.f[t,i,v] for i in self.in_nodes[v])== sum(self.f[t,v,j] for j in self.out_nodes[v]), "flow conservation") #flow conservation
            self.m.addConstr(sum(self.f[t,i,t] for i in self.in_nodes[t])==1)
            self.m.addConstr(sum(self.f[t,t,j] for j in self.out_nodes[t])==0)
            self.m.addConstr(sum(self.f[t,self.root,j] for j in self.out_nodes[self.root])==1)
        print("model initialized..")
    
    
    def run(self):
            T = nx.DiGraph()
            self.m.optimize()
            if self.m.Status == GRB.OPTIMAL:
                solution = self.m.getAttr('X', self.x)
                flow = self.m.getAttr('X', self.f)
            
            score = self.m.objVal
            print(score)
          

                # for i,j in self.special_arcs:
                #     print('%s -> %s: %g' % (i, j, polytomies[ i, j]))
                
            for i,j in self.edges:
                if solution[i,j] > 0:
                        #    print(f"edge {i} -> {j}")
                           T.add_edge(i,j)

            
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



def name_node(tree, node, label,  is_poly=False): 
        # lab = str(tree) + "_" + str(node) + "_" + str(label)
        lab = str(node) + "_" + str(label)
        if is_poly:
             lab += "_p"
        return  lab




#takes in a tree and constructs a network 
class ConstructGraph:
    def __init__(self, LinTree, iso_costs, sequence, isotypes) -> None:
          
        self.G = nx.DiGraph()
        self.T = LinTree.T
        self.root = LinTree.root
        self.id = LinTree.id
   
     
        self.iso_costs = iso_costs
        self.seq_weights = {}
        self.iso_weights= {}
        self.seq_labs = sequence 
        self.iso_labs = isotypes


    def build(self):
        min_state = {}
        nodes_to_states = {}
        postorder =  list(nx.dfs_postorder_nodes(self.T, source=self.root))
        for n in postorder:
            if self.T.out_degree[n] ==0:
                iso = self.iso_labs[n]
                new_node = name_node(self.id,n,iso)
                self.G.add_node(new_node)
                min_state[n] = self.iso_labs[n]
                nodes_to_states[n] = [self.iso_labs[n]]
                  
            elif self.T.out_degree[n] > 0: 
                poly_nodes_added = []
                degree = self.T.out_degree[n]
                min_state[n] = min([min_state[v] for v in self.T.neighbors(n)])
                nodes_to_states[n] =[]
                for i in range(0, min_state[n]+1):
                    if i > 0 and n in self.iso_labs:
                         break
                    new_node = name_node(self.id,n,i)
                    nodes_to_states[n].append(i)

            
                    self.G.add_node( new_node)
                    for v in self.T.neighbors(n):
                        for j in nodes_to_states[v]:
                        
                            self.G.add_edge(new_node, name_node(self.id,v,j))
                            self.seq_weights[new_node,name_node(self.id,v,j)] = hamming_distance(self.seq_labs[n], self.seq_labs[v])
                            self.iso_weights[new_node, name_node(self.id,v,j)] = self.iso_costs[i,j]
                
                #insert polytomy nodes and add edges to Graph
                if degree > 2:
                   
                    poly_node_states = []
                    for v in self.T.neighbors(n):
                         poly_node_states = poly_node_states + nodes_to_states[v]
                    max_node = max(poly_node_states)
                    node_counts = Counter(poly_node_states)

                    
                    for j in range(1, max_node+1):
                            if node_counts[j] > 1 and node_counts[j] < degree:
                                poly_nodes_added.append(name_node(self.id,n,j,True))
                                for i in nodes_to_states[n]:
                                    if i < j:
                                        self.G.add_edge(name_node(self.id,n,i), name_node(self.id,n,j,True))
                                        self.seq_weights[name_node(self.id,n,i), name_node(self.id,n,j,True)] = 0
                                        self.iso_weights[name_node(self.id,n,i), name_node(self.id,n,j,True)] = self.iso_costs[i,j]
                                
                                for v in self.T.neighbors(n):
                                        if j in nodes_to_states[v]:
                            
                                            self.G.add_edge(name_node(self.id,n,j,True), name_node(self.id,v,j))
                                            self.seq_weights[name_node(self.id,n,j,True),name_node(self.id,v,j)] = hamming_distance(self.seq_labs[n], self.seq_labs[v])
                                            self.iso_weights[name_node(self.id,n,j,True), name_node(self.id,v,j)] = self.iso_costs[j,j]
        return self.G, self.seq_weights, self.iso_weights
                                 

                      
                     
T = nx.DiGraph()
# T.add_edges_from([ ('root', 'r') , ('r','a'), ('r', 'b'), ('r', 'c'), ('r', 'd')])

T.add_edges_from([ ('root', 'r') , ('r','a'), ('r', 'b'), ('r', 'c'), ('r', 'd'), ('r', 'e'), ('r', 'f'), ('r', 'g')])

iso_costs = {(0,0): 0.91, (0,1): 1.6, (0,2): 1.6, (0,3): 1.6, (1,1): 0.11, (1,2):  3, (1,3): 3, (2,2): 0.51, (2,3): 0.43, (3,3): 0  }
sequences = {'root': ['a', 'a'], 'r': ['a', 'a'], 'a': ['a', 'a'], 'b': ['a', 'a'], 'c': ['a', 'a'] , 'd': ['a','a'], 'e': ['a','a'], 'f': ['a','a'], 'g': ['a','a']}


isotypes  = {'root': 0, 'a': 2, 'b':2, 'c': 2, 'd': 3, 'f':3, 'e': 3, 'g': 1}
lt =LineageTree(T,'root', 0)
 
G, s_weights, i_weights = ConstructGraph( lt, iso_costs, sequences, isotypes).build() 

print(list(G.nodes))
print(list(G.edges))

score, T = SteinerTree(G,s_weights, i_weights, root='root_0', lamb=0).run()
print(score)
print(list(T.edges))

                    
            # else:
            #      min_state[n] = min(min_state[v] for v in self.T.neighbors[n])
            #      for i in range(0, min_state[n]+1):
            #          new_node = name_node(n,i)
            #          nodes_to_states.append(i)
            #        for i in range(1, min_state[n] + 1):
            #             new_poly = name_poly(n,i)
            #             for v in self.T.neighbors[n]:
            #                  if i in node_states[v]:
            #                     self.G.add_edge(new_poly, name_node(v,i))
            #                     self.seq_weights[new_node,name_node[v,i]] = hamming_distance(self.seq_labs[n], self.seq_labs[v])
            #                     self.iso_weights[new_node, name_node(v,i)] = self.iso_costs[i,i]

                        


                
        
    

                     
                 
                 
    
                  

          


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

       
        
