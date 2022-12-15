


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
import numpy as np 

class PolytomyResolver:
    def __init__(self, tree_nodes, weights, scores, states, candidate_state):
        # Create optimization model
        self.m = gp.Model('netflow')
        self.tree_nodes = tree_nodes
        self.nodes = ['source', 'sink'] + [n for n in tree_nodes]
        self.states = [str(s) for s in states] 
        self.score_nodes = []
        self.capacity_upper = {}
        self.capacity_lower = {}
        self.fixed_costs = {}
        self.costs = {}
        for n in tree_nodes:
            self.capacity_upper[("source", n )] = 1
            self.capacity_lower[("source", n )] = 0
            self.costs[("source", n)] = 0
            for state in scores[n]:
                node_name = f"{n}_{state}"
                self.score_nodes.append(node_name)
                self.nodes.append(node_name)
       
        
                self.capacity_upper[(n,node_name)] =1
                self.capacity_lower[(n,node_name)] =0
                self.costs[(n,node_name)] = scores[n][state]

        
                self.capacity_upper[(node_name,"sink")] = 1
                self.capacity_lower[(node_name,"sink")] = 0
                self.costs[(node_name, "sink")] = weights[candidate_state, int(state)]
                for s in self.states:
                    if int(s) != candidate_state:
                        self.capacity_upper[(node_name,f"state_{s}")] = 1
                        self.capacity_lower[(node_name,f"state_{s}")] = 0
             
                        self.costs[(node_name,f"state_{s}")] = weights[int(s), int(state)]
                
        self.special_arcs = []
        #TODO: need to add n-1 layers of states
        for s in self.states:
            if int(s) != candidate_state:
                node_name = f"state_{s}"
                self.special_arcs.append((node_name, "sink"))
       
                self.nodes.append(node_name)
                self.capacity_upper[(node_name,"sink")] = len(tree_nodes)-1
                self.capacity_lower[(node_name,"sink")] =  2 #0
                self.fixed_costs[(node_name, "sink")] = weights[candidate_state, int(s)]
        
   
        self.all_arcs, self.cap_upper = gp.multidict(self.capacity_upper)

        # self.inflow = {a: 0 for a in self.all_arcs}
        # self.inflow[('edges','sink')] = len(tree_nodes)
        # self.inflow['edge','source'] = -1*len(tree_nodes)

        self.normal_arcs = [i for i in self.all_arcs if i not in self.special_arcs]

        # self.flow = self.m.addVars(self.all_arcs, lb=self.capacity_lower, ub=self.capacity_upper, vtype=GRB.INTEGER, name='flow')
        self.flow = self.m.addVars(self.all_arcs, vtype=GRB.INTEGER, name='flow')

        self.polytomy = self.m.addVars(self.special_arcs, vtype=GRB.BINARY, name='polytomy')



        # Arc-capacity constraints
        self.m.addConstrs(
            (self.flow[ i, j] <= self.capacity_upper[i, j] for i, j in self.normal_arcs), "cap_upper")

        self.m.addConstrs(
            (self.flow[ i, j] >= self.capacity_lower[i, j] for i, j in self.normal_arcs), "cap_lower")
        
        #polytomy is opened only if flow is greater than lower bound of capacity
        self.m.addConstrs(
            (self.flow[ i, j]>= self.capacity_lower[i, j]*self.polytomy[i,j] for i, j in self.special_arcs), "cap_special_lower")
        
        self.m.addConstrs(
            (self.flow[ i, j] <= self.capacity_upper[i, j]*self.polytomy[i,j] for i, j in self.special_arcs), "cap_special_upper")
        
        #Flow conservation 
        # m.addConstr(quicksum(x[i,j] for j in (set(V) - set(S))) >= 2)
        self.in_nodes = {j : [] for j in self.nodes}
        self.out_nodes = {j : [] for j in self.nodes}
        for j in self.nodes:
            for (i,k) in self.all_arcs:
                if k==j:
                    self.in_nodes[j].append(i)
                if i == j:
                    self.out_nodes[j].append(k)
        
        self.m.addConstr(sum(self.flow['source',j] for j in self.out_nodes['source']) == len(tree_nodes))


        #flow leaving the source cannot exceed the number of tree nodes 
        self.m.addConstr(sum(self.flow[j, 'sink'] for j in self.in_nodes['sink']) == len(tree_nodes))
        for  j in self.nodes:
            if j not in ['source', 'sink']:
                self.m.addConstr(sum(self.flow[i, j]  for i in self.in_nodes[j]) == sum(self.flow.sum(j, i) for i in self.out_nodes[j]))
        
        #flow leaving the source cannot exceed the number of tree nodes 
   

        self.m.setObjective(sum(self.costs[i,j]*self.flow[i,j] for i,j in self.normal_arcs) + 
                        sum(self.fixed_costs[i,j]*self.polytomy[i,j] for i,j in self.special_arcs)  , GRB.MINIMIZE)
        

    def run(self):
            states = {}
            poly_map = {}
            score= None
            self.m.optimize()
            if self.m.Status == GRB.OPTIMAL:
                solution = self.m.getAttr('X', self.flow)
                polytomies = self.m.getAttr('X', self.polytomy)

                # for i,j in self.special_arcs:
                #     print('%s -> %s: %g' % (i, j, polytomies[ i, j]))
                
                for i in self.tree_nodes:
                    for j in self.out_nodes[i]:
                        if solution[i,j] > 0:
                            states[i] = int(j.split("_")[1])

                # print(states)
                score = self.m.objVal
              
                # print("normal used edges in network")
                # # for i, j in self.all_arcs:
                # #         if solution[ i, j] > 0:
                # #             print('%s -> %s: %g' % (i, j, solution[ i, j]))
                
                # print("polytomies used")
                for i,j in self.special_arcs:
                    if polytomies[i,j] > 0:
                        # print('%s -> %s' % (i, j))
                        t= int(i.split("_")[1])
                        
                        for k in self.in_nodes[i]:
                            if solution[k,i] > 0:
                              
                                child =int(k.split("_")[0])
                                if t in poly_map:
                                    poly_map[t].append(child)
                                else:
                                    poly_map[t] = [child]
                                         
            return score, states, poly_map




        
    
            
            


# tree_nodes = ['a', 'b', 'c']
# states = [0,1]
# weights = {(0, 0): 0, (0, 1): 1, (1,0): 0, (1,1): 0}
# scores= {'a': {0: 3, 1: 4}, 'b': {0: 2, 1: 5}, 'c': {0: 2, 1: 11}}
# candidate_state = states[0]


# tree = nx.DiGraph()

# tree.add_edges_from([(0,1), (0,2), (0,3), (1,4), (1,5), (2,6), (2,7), (3,8), (3,9)])

tree_nodes = [1,2,3]
candidate_state = 0
scores = {1: {0: 2, 1: 1}, 2: {0:2, 1:0}, 3: {0: 2, 1: 2, 2: 0}}

states = [0,1,2]
# isotypes = {4:1, 5:2, 6:1, 7:2, 8:1, 9:2}

weights = {}
for s in states:
    for t in states:
        if s > t:
            weights[s,t] = 1000000
        elif s==t:
            weights[s,t] =0
        else:
            weights[s,t] =1




pr = PolytomyResolver(tree_nodes, weights, scores, states, candidate_state)
pr.run()

       
        
