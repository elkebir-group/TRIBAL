import pyomo
import pyomo.environ as pyo
import pyomo.opt
import networkx as nx 
import numpy as np 
import logging

class MPTR:
    def __init__(self, G, T, S, edge_weight, tree_to_graph, root=0, threads=3) -> None:
        
        self.G= G
        self.T = T
        self.terminals = S 
        self.tree_to_graph = tree_to_graph
        self.root= root 
        self.nodes = set(G.nodes)
        self.edges = set(G.edges)
        self.c = {e: edge_weight[e]for e in self.edges}
        self.internal_nodes = [n for n in self.nodes if n not in self.terminals and n!= self.root ]
        self.flow_dest = {(t,e) for e in self.edges for t in self.terminals}


    def createModel(self):
 
             
        self.m = pyo.ConcreteModel()


        self.m.f = pyo.Var(self.terminals, self.edges, domain=pyo.NonNegativeReals)
        self.m.x = pyo.Var(self.edges, domain=pyo.Binary)

        def obj_rule(m):
            return sum(self.c[e] *m.x[e] for e in self.edges)
        
        self.m.OBJ = pyo.Objective(rule=obj_rule, sense=pyo.minimize)

        def flow_upper_bound(m, t,u,v):

            return m.f[t,u,v] <= m.x[u,v]
        
        self.m.UpperBound = pyo.Constraint(self.terminals, self.edges, rule=flow_upper_bound)

        def transitory_rule(m, u, v):
            return sum(m.x[i, j] for i in self.tree_to_graph[u] for j in self.tree_to_graph[v] if (i, j) in self.edges) == 1

        self.m.TransitoryConstraint = pyo.Constraint(self.T.edges, rule=transitory_rule)


        
        
        def flow_conservation_rule(m,t,v):
            incoming= list(self.G.predecessors(v))
            outgoing = list(self.G.successors(v))
            if v == self.root:
                return sum(m.f[t,(v,o)] for o in outgoing) ==1
            elif v in self.terminals:  
                return pyo.Constraint.Skip
            
            else:  #internal node 
                return sum(m.f[t,(i,v)] for i in incoming)== sum(m.f[t,(v,o)] for o in outgoing)
            

        self.m.flow_conservation = pyo.Constraint(self.terminals, self.nodes, 
                                                  rule=flow_conservation_rule, name="flow conservation")
 
        def terminal_incoming(m, t):
            incoming= self.G.predecessors(t)
            return sum(m.f[t,(i,t)] for i in incoming)==1
        self.m.term_incoming = pyo.Constraint(self.terminals, rule=terminal_incoming, 
                                           name="terminal incoming flow")

        def terminal_outgoing(m,t):
            outgoing = list(self.G.successors(t))
            if len(outgoing) == 0:
                return pyo.Constraint.Skip
            return sum(m.f[t,(t,o)]for o in outgoing) ==0
        
        self.m.term_outgoing = pyo.Constraint(self.terminals, rule=terminal_outgoing, 
                                           name="terminal outgoing flow")
        

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
     

            logging.getLogger('pyomo.core').setLevel(logging.ERROR)
            self.createModel()
            solver = pyomo.opt.SolverFactory("glpk")
    
            results = solver.solve(self.m, tee=False, keepfiles=False)

            if (results.solver.status != pyomo.opt.SolverStatus.ok):
                print('Check solver not ok?')
            if (results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal):  
                print('Check solver optimality?') 

            if results.solver.termination_condition == pyo.TerminationCondition.optimal:
                solution = {e: pyo.value(self.m.x[e]) for e in self.edges}
                flow = {(t, e): pyo.value(self.m.f[t, e]) for t, e in self.flow_dest}
                score = pyo.value(self.m.OBJ)  # Assuming you have an objective defined
                T = nx.DiGraph()
                for e in self.edges:
                  
                    total_flow = sum(flow[t, e] for t in self.terminals)
               
                    if solution[e] > 0.5 and total_flow > 0:
                        # print(f"{'glpk'}: {e} {total_flow}")
                        T.add_edge(*e)
        
        
            return score,T
        