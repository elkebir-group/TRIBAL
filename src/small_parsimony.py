import networkx as nx
import numpy as np
from scipy.special import logsumexp
from trans_matrix import TransMat
from itertools import product
from copy import deepcopy
# from polytomy_resolver import PolytomyResolver as pr
np.seterr(all="ignore")

class SmallParsimony:
    def __init__(self, T,  root,  alphabet=["A", "C", "G", "T", "N", "-"], cost=None):
        
        self.T = T
        self.root = root
        self.cost = cost
        self.alphabet = alphabet

        if self.cost is None:
            self.cost= {}

            for i in self.alphabet:
                for j in self.alphabet:
                    if (i,j) not in self.cost:
                        if i == "-" and j == "-":
                            self.cost[i,j] =np.Inf
                        if i == "-" or j == "-":
                            self.cost[i,j] = 1
                        elif  i ==j:
                            self.cost[i,j] = 0
                        
                        else:
                            self.cost[i,j] = 1

        self.postorder = self.postorder_traversal()
        self.nodes = list(self.T.nodes())

      
        # self.att_name = attribute
        # if attribute == "SEQ":
        #     self.att=sequences
        # else:
        #     self.att = isotypes
         #nx.get_node_attributes(self.T, attribute)
        
        
    

    def postorder_traversal(self):
        return list(nx.dfs_postorder_nodes(self.T, source=self.root))
    

    def preorder_traversal(self):
        return list(nx.dfs_preorder_nodes(self.T, source=self.root))
    
    def is_leaf(self, node):
        return self.T.out_degree(node) ==0

 

    # def sankoff_polytomy2(self, isotypes, weights, states):
    #     self.att= isotypes 
    #     dp_matrix = {n : {}  for n in self.nodes}
    #     dp_bt = {n : {} for n in self.nodes}

    #     #save the state of the polytomy and the proposed descendents
    #     poly_bt = {n : {} for  n in self.nodes}


    #     for n in self.postorder:
        
    #         if self.is_leaf(n):
    #             dp_matrix[n][isotypes[n]] =0
              
            
    #         else:
    #              #if node is binary, we perform regular Sankoff
    #             if self.T.out_degree[n] == 2:
                    
    #                 for s in states:
    #                     dp_bt[n][s] = {}
    #                     total_score = 0
                        
    #                     for c in self.T.successors(n):
                        
    #                         best_score = np.Inf
    #                         for t in states:
    #                             if t in dp_matrix[c]:
    #                                 score = weights[s,t] + dp_matrix[c][t]
    #                                 if score < best_score:
    #                                     best_score = score
    #                                     best_state = t
    #                         dp_bt[n][s][c] = best_state
    #                         total_score += best_score
    #                     if total_score < np.Inf:
    #                         dp_matrix[n][s] = total_score
    #                     else:
    #                         dp_bt[n][s][c] = {}
    #             else:
    #                 #we simulatenously resolve the polytomy and complete the dp table
    #                 best_score = np.Inf
    #                 child_order = list(self.T.successors(n))
    #                 child_scores = {c: dp_matrix[c] for c in child_order }
                    
    #                 #we only need to consider internal states that are less than or equal to the min possible state 
    #                 #in the children
                    
                
    #                 min_state = max(states)
    #                 for c in child_scores:
    #                     for t in child_scores[c]:
    #                         if t < min_state:
    #                             min_state = t
    #                 inner_node_states = [s for s in states if s <= min_state]
                    
    #                 for s in inner_node_states:
                        
    #                     score, opt_states, poly_map = pr(child_order, weights, child_scores, states, candidate_state=s ).run()

    #                     if score < best_score:
    #                         best_score = score 
    #                         dp_matrix[n][s] = score
    #                         dp_bt[n][s] = opt_states
    #                         poly_bt[n][s] = poly_map
                            
    #                         #fill out dp table for state s
    #                 # dp_matrix[n][s] = best_score 
    #                 # dp_bt[n][s] = {}
    #                 # if best_score < np.Inf:
    #                 #     for c,t in zip(child_order, best_state):
    #                 #         dp_bt[n][s][c] = t
    #                 #     poly_bt[n][s] = deepcopy(best_poly_dict)

        
    #     return dp_matrix, dp_bt, poly_bt 
    
    def sankoff_poly(self, isotypes, weights, states, ighd=None):
        total_nodes_added = 0
        self.att= isotypes 
        dp_matrix = {n : 0 for n in self.nodes}
        dp_bt = {}



    


        for n in self.postorder:
            
            if self.is_leaf(n):
                dp_matrix[n] =0
                dp_bt[n] = self.att[n]
            
            
            
            else:
                children = [c  for c in self.T.successors(n)]
                if n == self.root:
                    dp_bt[n] = 0
                else:
                    child_states = [dp_bt[c] for c in children]
                    if ighd is not None:
                 
                        if any([c ==ighd for c in child_states]):
                            if  all([c ==ighd for c in child_states]):
                                dp_bt[n] = ighd 
                            else:
                                dp_bt[n] = min(states)
                        else:
                            dp_bt[n] = min(child_states)
                    else:
                        dp_bt[n] = min(child_states)
             

                s = dp_bt[n]
                state_counts = {t:0 for t in states}
                state_lists = {t: [] for t in states}
                for c in children:
                    state_counts[dp_bt[c]] += 1
                    state_lists[dp_bt[c]].append(c)
            
                if len(children) > 2:
                    for t in state_counts:
                        
                        if state_counts[t] > 1 and t != s:
                            
                            #add in new node to partially resolve polytomy
                            poly_node = len(list(self.T.nodes)) + 1
                            while poly_node in self.T:
                                poly_node += 1
                            total_nodes_added += 1
                            dp_bt[poly_node] = t
                            dp_matrix[n] += weights[s,t]
                            for kid in state_lists[t]:
                                    self.T.remove_edge(n, kid)
                                
                                    self.T.add_edge(n, poly_node)
                                    self.T.add_edge(poly_node, kid)
                                    dp_matrix[n] += dp_matrix[kid] + weights[t,t] 
                        else:
                            for kid in state_lists[t]:
                                dp_matrix[n] += dp_matrix[kid] + weights[s,t]
                
                else:
                    
                    for c in children:
                        t = dp_bt[c]
                        dp_matrix[n] += dp_matrix[c] + weights[s,t]
        # print(total_nodes_added)
        return dp_matrix, dp_bt

       


               

        

        
    
    # def sankoff_polytomy(self, isotypes, weights, states):
    #     self.att= isotypes 
    #     dp_matrix = {n : {}  for n in self.nodes}
    #     dp_bt = {n : {} for n in self.nodes}

    #     #save the state of the polytomy and the proposed descendents
    #     poly_bt = {n : {} for  n in self.nodes}


    #     for n in self.postorder:
        
    #         if self.is_leaf(n):
    #             dp_matrix[n][isotypes[n]] =0
              
            
    #         else:
    #              #if node is binary, we perform regular Sankoff
    #             if self.T.out_degree[n] == 2:
                    
    #                 for s in states:
    #                     dp_bt[n][s] = {}
    #                     total_score = 0
                        
    #                     for c in self.T.successors(n):
                        
    #                         best_score = np.Inf
    #                         for t in states:
    #                             if t in dp_matrix[c]:
    #                                 score = weights[s,t] + dp_matrix[c][t]
    #                                 if score < best_score:
    #                                     best_score = score
    #                                     best_state = t
    #                         dp_bt[n][s][c] = best_state
    #                         total_score += best_score
    #                     if total_score < np.Inf:
    #                         dp_matrix[n][s] = total_score
    #                     else:
    #                         dp_bt[n][s][c] = {}
    #             else:
    #                 #we simulatenously resolve the polytomy and complete the dp table
    #                 net_best_score = np.Inf
    #                 child_order = list(self.T.successors(n))
    #                 child_scores = {c: dp_matrix[c] for c in child_order }
                    
    #                 #we only need to consider internal states that are less than or equal to the min possible state 
    #                 #in the children
                    
                
    #                 min_state = max(states)
    #                 for c in child_scores:
    #                     for t in child_scores[c]:
    #                         if t < min_state:
    #                             min_state = t
    #                 inner_node_states = [s for s in states if s <= min_state]
                    
    #                 for s in inner_node_states:
                        
    #                     net_score, opt_states, poly_map = pr(child_order, weights, child_scores, states, candidate_state=s ).run()

    #                     if net_score < net_best_score:
    #                         net_best_score = net_score 
    #                         best_opt_states = opt_states
    #                         bst_poly_map = poly_map
    #                         dp_matrix[n][s] = net_score
    #                         dp_bt[n][s] = opt_states
    #                         poly_bt[n][s] = poly_map
    #                     # print(f"state{s}: {net_best_score}")  
    #                 child_order = list(self.T.successors(n))
    #                 state_cand =[]
    #                 for c in child_order:
    #                     state_cand.append(list(dp_matrix[c].keys()))
    #                 for s in inner_node_states:
    #                         best_score = np.Inf
    #                         best_poly_dict = {}
                        
    #                         for cand in product(*state_cand):
    #                             cand_poly= {}
    #                             cand_score = 0
    #                             state_counts = {t:0 for t in states}
    #                             state_lists = {t: [] for t in states}
                        
    #                             #initialize the score with normal sankoff
    #                             for c,t in zip(child_order, cand):
    #                                 cand_score += weights[s,t] + dp_matrix[c][t]
                                  
    #                                 state_counts[t] += 1
    #                                 state_lists[t].append(c)
                                
    #                             if cand_score < best_score:
    #                                 best_score = cand_score 
    #                                 best_state = cand
    #                                 best_poly_dict = {}

                             
    #                             if any([state_counts[t] ==len(child_order) for t in state_counts]):
    #                                 continue
                             
    #                             #now attempt to introduce poltyomies to see if the score can be improved
                            
    #                             cand_score  = 0
    #                             poly_cand = {}
    #                             for t in state_counts:
    #                                 if state_counts[t] ==0:
    #                                     continue
                                
    #                                 subscore = sum(dp_matrix[kid][t] for kid in state_lists[t])
                                   
                        
    #                                 if state_counts[t] > 1:
                                       
    #                                 #try every possible state to resolve the polytomy, tracking the polytomy that is best
    #                                     best_poly_score = np.Inf
    #                                     for j in states:
                                            
    #                                         poly_score = weights[s,j] + weights[j,t]* state_counts[t] + subscore
    #                                         if poly_score < best_poly_score:
                                             
    #                                             poly_cand[t] = j
    #                                             best_poly_score =poly_score

    #                                     cand_score += best_poly_score

    #                                     if t in poly_cand:
    #                                         if poly_cand[t] in cand_poly:
    #                                             cand_poly[poly_cand[t]] += state_lists[t]

    #                                         else:
    #                                             cand_poly[poly_cand[t]] = state_lists[t]
    #                                     # for kid in state_lists[t]:
    #                                     #     if t in poly_cand:

    #                                     #         cand_poly[kid] = poly_cand[t]
                                     
                                    
    #                                 elif state_counts[t] == 1:
    #                                     cand_score += weights[s,t] + subscore

    #                                     #do normal sankoff 
                                
    #                             if cand_score < best_score:
    #                                 best_score = cand_score
    #                                 best_state = cand 
    #                                 best_poly_dict = deepcopy(cand_poly)
                            
    #                         #fill out dp table for state s
    #                         dp_matrix[n][s] = best_score 
    #                         dp_bt[n][s] = {}
    #                         if best_score < np.Inf:
    #                             for c,t in zip(child_order, best_state):
    #                                 dp_bt[n][s][c] = t
    #                             poly_bt[n][s] = deepcopy(best_poly_dict)
        
        
    #     return dp_matrix, dp_bt, poly_bt 

    # def polytomy_backtrace(self,dp_mat, dp_bt, poly_bt, start_state=0):
    #     tree = self.T.copy()
    #     labels = {}
    #     best_score = np.Inf
    #     if start_state is None:
    #         for s in dp_mat[self.root]:
    #             if dp_mat[self.root][s] < best_score:
    #                 best_score = dp_mat[self.root][s]
    #                 labels[self.root] = s
    #     else:
    #         labels[self.root] = start_state
    #         best_score = dp_mat[self.root][start_state]

    #     total_nodes_added =0 
    #     for n in self.preorder_traversal():
    #             if self.is_leaf(n):
    #                 labels[n] = self.att[n]
    #                 continue 
    #             state = labels[n]
 
    #             #update the state of children nodes
    #             for kid in dp_bt[n][state]:

    #                 labels[kid] = dp_bt[n][state][kid]
                
    #             #add in any inferred polytomies and assign their label
    #             if len(poly_bt[n]) >0:
    #                 if len(poly_bt[n][state]) > 0:
                        
    #                         for j in poly_bt[n][state]:
    #                             poly_node = str(len(list(tree.nodes)) + 1)
    #                             labels[poly_node] = j
    #                             # print(f"adding node {poly_node}")
    #                             total_nodes_added += 1
    #                             for kid in poly_bt[n][state][j]:
    #                                 tree.remove_edge(n, kid)
                                
    #                                 tree.add_edge(n, poly_node)
    #                                 tree.add_edge(poly_node, kid)
            
                         
    #     # print(f"total nodes added: {total_nodes_added}")
    #     return best_score, labels, tree

    def polytomy_resolver(self, leaf_labels,  weights, states, ighd=None):
            dp_mat, dp_bt = self.sankoff_poly(leaf_labels, weights, states, ighd=ighd)

            # # dp_mat, dp_bt, dp_poly = self.sankoff_polytomy2(leaf_labels, weights, states)
            # score, labels, tree = self.polytomy_backtrace(dp_mat, dp_bt,dp_poly, start_state)

 
            return dp_mat[self.root], dp_bt, self.T




    def sankoff_dp(self, pos):
        dp_matrix = {n : {}  for n in self.nodes}
     
        for n in self.postorder:

            if self.is_leaf(n):
                for a in self.alphabet:

                    if self.att[n][pos] == a:
                        dp_matrix[n][a] = 0
                    else:
                        dp_matrix[n][a] = np.Inf
            
  

            elif n == self.root:
                base = self.att[n][pos]
                for a in self.alphabet:
                    if a == base:
                        dp_matrix[n][a] = np.sum([[self.min_cost(c,a, dp_matrix) for c in self.T.successors(n)]])
                    else:
                        dp_matrix[n][a] = np.Inf
            
            else:
                
          
                for a in self.alphabet:
                    score_array = np.array([self.min_cost(c,a, dp_matrix) for c in self.T.successors(n)])
                    
                
                    dp_matrix[n][a] = score_array.sum()

        opt_score = np.array([dp_matrix[self.root][a] for a in self.alphabet]).min()
        
        return opt_score, dp_matrix

    @staticmethod
    def argmin( my_dict):
        minval = np.Inf
        
        for key,val in my_dict.items():
            if val <= minval:
                minval = val
                minarg = key
            
        return minarg


    def traceback(self, dp_matrix):
        labels = {}
        labels[self.root] = self.argmin(dp_matrix[self.root])
        for c in self.T.successors(self.root):
            self.traceback_helper(labels[self.root], c, dp_matrix, labels)
        
        return labels
        
     
            
    def traceback_helper(self,state, root, dp_matrix, labels):
        self.find_state(state, root, dp_matrix, labels)
        for c in self.T.successors(root):
            self.traceback_helper(labels[root], c,dp_matrix, labels)
    


    def find_state(self, state, parent,  dp_matrix, labels):

        min_cost = np.Inf
        for a in self.alphabet:
            if parent == self.root:
                trans_cost = dp_matrix[parent][state]
            else:
                trans_cost = self.cost[a,state] + dp_matrix[parent][a]
            if trans_cost < min_cost:
                min_cost = trans_cost
                labels[parent] = a

    @staticmethod 
    def concat_labels(all_labs, nodes):
        seq_assign = {}
        for k in nodes:
            seq_assign[k] = [all_labs[i][k] for i in all_labs if k in all_labs[i]]
          
            
        return seq_assign


    def sankoff(self, sequences):
        self.att = sequences
        if type(self.att[self.root]) != int:
            seq_length = len(self.att[self.root])
    
        else:
            seq_length = 1
            self.att = {key: [val] for key,val in self.att.items()}

        all_labels = {}
 
        min_scores =np.zeros(seq_length, dtype=float)
        for i in range(seq_length):
            min_scores[i], dp_mat = self.sankoff_dp(i)
            label_dict = self.traceback(dp_mat)
            all_labels[i] = label_dict
        
        node_labels= self.concat_labels(all_labels, self.nodes)
        if seq_length == 1:
            for key in node_labels:
                if key in self.att:
                    node_labels[key] = self.att[key]
            node_labels = {key: val[0] for key,val in node_labels.items()}
        # nx.set_node_attributes(self.T, node_labels, self.att_name)

        return min_scores.sum(),node_labels

            

    def min_cost(self, child, to_state, dp_matrix):
   

        scores = np.array([self.cost[to_state,a] + dp_matrix[child][a]for a in self.alphabet])
        
        return scores.min()
    

    def fitch_dp(self):
   
        opt_states = {}
        for n in self.postorder:
            if len(list(self.T.successors(n))) ==0:
                opt_states[n] = int(self.att[n])
            
            else:
                cand_states = []
                for s in self.T.successors(n):
                    cand_states.append( opt_states[s])
                
                opt_states[n]= np.min(np.unique(cand_states))
        return opt_states
                
    
    def fitch_score(self, dp_mat):
        score = 0
        for n in self.postorder:
            for u in self.T.predecessors(n):
                score += (dp_mat[n] != dp_mat[u])
        
        return score
    
    def likelihood_score(self, isotypes, transMat):
        likelihood = {}
        #likelihood n, s is the likelihood of the subtree rooted at node n when taking character state s
  
        isotype_states = np.arange(6)
  
   


        for n in self.postorder:

            for s in isotype_states:

                #initialize the base case
                if len(list(self.T.neighbors(n))) ==0:
                        likelihood[n,s] = np.log(1*(s==isotypes[n]))
                       
                else:
                    partial_total = 0
                    for u in self.T.neighbors(n):
                        state_val = []
                        for t in isotype_states:
                            state_val.append( np.log(transMat[s,t])  + likelihood[u,t])
                        partial_total += logsumexp(state_val)
                    likelihood[n,s] = partial_total
        
        return likelihood[0,0], likelihood
                


    
    def fastml_dp(self, transMat):
        Ldict = {n : {} for n in self.postorder}
        Cdict = {n : {} for n in self.postorder}
        isotypes = np.arange(6)
      
        for n in self.postorder:
            if len(list(self.T.neighbors(n))) ==0:
                for u in isotypes:
                    Cdict[n][u] = [self.att[n]]
                
                for v in isotypes:
                        Ldict[n][v] = np.log(transMat[v,Cdict[n][v][0]])
            else:
                for i in isotypes:
                    best_like = np.NINF
                
                    for j in isotypes:
                        val_array = [np.log(transMat[i,j])]
                        for child in self.T.successors(n):
                            val_array.append(Ldict[child][j])
                        
                        like_j = np.sum((val_array))
                        if like_j >  best_like:
                            Cdict[n][i] = [j]
                            best_like = like_j 
                        elif like_j == best_like and like_j > np.NINF:
                            Cdict[n][i].append(j)
                        elif like_j == best_like:
                            Cdict[n][i] =[]
                    
                    Ldict[n][i]  = best_like 
        return Ldict, Cdict

                        

    def fastml_traceback(self, Ldict, Cdict):
        score = Ldict[self.root][0]

        for n in self.preorder_traversal():
            if n == self.root:
                opt_state = {self.root: [0]}
            else:
                parent = list(self.T.predecessors(n))[0]
                opt_state_parent = opt_state[parent]
                opt_state[n] = []
                for o in opt_state_parent:
                    for c in Cdict[n][o]:
                        opt_state[n].append(c)
        
        return score, opt_state

   



    def fastml(self, isotypes, transMat):
        self.att = isotypes 
        Ldict, Cdict = self.fastml_dp(transMat)
        score, opt_states = self.fastml_traceback(Ldict, Cdict)

        return score, opt_states
    

    def fitch(self, isotypes, transMat=None):
        self.att = isotypes 
        opt_states = self.fitch_dp()
        score = self.fitch_score(opt_states)
        # score, likelihood = self.likelihood_score(isotypes, transMat)
        isotype_states = np.arange(6)
        poss_states = {}
        for n in self.postorder:
          

            poss_states[n] = [s for s in isotype_states if likelihood[n,s] > np.NINF ]
  
        



        return score, poss_states
        


# tree = nx.DiGraph()

# tree.add_edges_from([  (0,1), (0,2), (0,3), (1,4), (1,5), (2,6), (2,7), (3,8), (3,9)])

# states = [0,1,2]
# isotypes = {4:1, 5:2, 6:1, 7:2, 8:2, 9:2}

# weights = {}
# for s in states:
#     for t in states:
#         if s > t:
#             weights[s,t] = np.Inf
#         elif s==t:
#             weights[s,t] =0
#         else:
#             weights[s,t] =1

# sp = SmallParsimony(tree, 0, states, weights)
# score, labels, tree = sp.sankoff_poly(isotypes, weights, states)
# # score,labels, tree = sp.polytomy_backtrace(dp_mat, dp_bt, dp_poly)
# print(score)
# for key, value in labels.items():
#     print(f"node: {key} label: {value}")
        



# if __name__ == "__main__":

#     tree = nx.DiGraph()
#     tree.add_edges_from([(0,1), (0,4), (1,2), (1,3)])

#     cost_mat = {}
#     alphabet = ("a", "g", "c", "t")
#     seq = {2: "cc", 3: "gg", 4: "tt", 0:"cc"}
#     for a in alphabet:
#         for b in alphabet:
#             if a == b:
#                 cost_mat[a,b] =0
#             elif a in ["a", "g"] and b in ["a", "g"]:
#                 cost_mat[a,  b] = 1
#             elif a in ["c", "t"] and b in ["c", "t"]:
#                 cost_mat[a,b] =1
#             else:
#                 cost_mat[a,b] =3
#     nx.set_node_attributes(tree, seq, "sequence")

#     sk = SmallParsimony(tree,  alphabet, 0, "sequence", cost_mat)
#     opt_score, labels = sk.sankoff()
#     print(f"Parisomony score: {opt_score}")

#     tree2 = nx.DiGraph()
#     tree2.add_edges_from([(0,1), (0,4), (1,2), (1,3)])

#     cost_mat = {}
#     for i in range(7):
#         for j in range(7):
#             if i < j:
#                 cost_mat[str(i),str(j)] = 1
#             elif i ==j:
#                 cost_mat[str(i),str(j)] = 0
#             else:
#                 cost_mat[str(i),str(j)] = np.Inf
#     alphabet = [str(i) for i in range(7)]
#     seq = {2: "2", 3: "1", 4: "5", 0:"0"}


    # nx.set_node_attributes(tree2, seq, "state")
    # sk2 = SmallParsimony(tree2, alphabet, 0, "state", cost_mat)
    # score, labels = sk2.fitch()
    # print(labels)
    # sk_score, sk_labels = sk2.sankoff()
    # print(sk_labels)

    # print(f"Fitch parisomony score: {score}\nSankoff parsimony score: {sk_score}")

   



    

    
