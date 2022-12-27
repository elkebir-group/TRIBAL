import numpy as np
import networkx as nx 
import itertools
from scipy.special import logsumexp
class EMProbs:
    def __init__(self, tree, root, states,transmat=None,  maxiter=100, epsilon=0.1) -> None:

        self.T = tree
        self.root =root 
        self.states = states 
        self.tmat = transmat 
     
        self.maxiter = maxiter
        
        self.epsilon = epsilon 
        self.postorder = self.postorder_traversal()
        self.preorder = self.preorder_traversal()
    
 
    
        self.beta = {n: {} for n in self.preorder}
        self.alpha = {n: {} for n in self.preorder}
        self.beta_parent = {n: {self.parent(n): {}} for n in self.preorder}

        self.beta_subtree_minus_child = {n: {c: {} for c in self.children(n) }for n in self.preorder}
        self.log_tmat = np.log(self.tmat)

        self.Amat = {n: {} for n in self.preorder}
        self.Bmat = {n: {s: {} for s in self.states} for n in self.preorder}

       



    def postorder_traversal(self):
        return list(nx.dfs_postorder_nodes(self.T, source=self.root))
    

    def preorder_traversal(self):
        return list(nx.dfs_preorder_nodes(self.T, source=self.root))

    def parent(self,n):
        if n == self.root:
            return None
        return list(self.T.predecessors(n))[0]
    
    def children(self,n):
        return list(self.T.neighbors(n))
    

    def is_leaf(self, node):
        return self.T.out_degree(node) ==0
    

    def upward_recursion(self):
        for n in self.postorder:
 
            if self.is_leaf(n):
                for j in self.states:

                    self.beta[n][j] = np.log(int(self.obs_states[n]==j))

            else:
                for c in self.children(n):
                    for m in self.states:
                        arr= []
                        for j in self.states:
                            arr.append(self.log_tmat[m,j] + self.beta[c][j])
                        
                        self.beta_parent[c][n][m] = logsumexp(arr)
                
                for m in self.states:
                    self.beta[n][m] = np.sum([self.beta_parent[c][n][m] for c in self.children(n)])
                
                
                for c in self.children(n):
                    for m in self.states:
                        self.beta_subtree_minus_child[n][c][m] = self.beta[n][m] - self.beta_parent[c][n][m]
            
    def downward_recursion(self):
        for n in self.preorder:
            if n ==self.root:
                for m in self.states:
                    self.alpha[n][m] = np.log(int(m==0))

            else:
                parent = self.parent(n)
          
                for m in self.states:
                    arr = []
                    for j in self.states:
                        arr.append(self.alpha[parent][j] + self.log_tmat[j,m] + self.beta_subtree_minus_child[parent][n][j])
                    arr =[a for a in arr if not np.isnan(a)]
                    if len(arr) > 0:
                        self.alpha[n][m] = logsumexp(arr)
                    else:
                        self.alpha[n][m] = np.NINF

    def compute_conditionals(self):
        for n in self.preorder:
            if n ==1:
                print("here")
           
            for m in self.states:
                arr = []
                for j in self.states:
                    arr.append(self.alpha[n][j] + self.beta[n][j])
                

                self.Amat[n][m] =self.alpha[n][m] + self.beta[n][m] - logsumexp(arr)
                if n != self.root:
                    parent= self.parent(n)
                    for j in self.states:
                        if self.log_tmat[j,m] != np.NINF:
    
                            self.Bmat[n][m][j] = self.alpha[parent][j] + self.beta_subtree_minus_child[parent][n][j] + \
                                                self.beta[n][m] + self.log_tmat[j,m] - logsumexp(arr)
                        else:
                            self.Bmat[n][m][j] = np.NINF 
    

        

    def e_step(self):
        self.upward_recursion()

        self.downward_recursion()
        self.compute_conditionals()

        
    def check_convergence(self, curr, old ):
        return  np.abs(curr - old ) <= self.epsilon
      

    
    def  m_step(self):
    
       
        #log like expectation of for root being in state 0
        exp_log_like = np.exp(self.Amat[self.root][0]*np.log(1))
        new_tmat = np.zeros((len(self.states), len(self.states)))
        state_probs = np.zeros(len(self.states))

        for node in self.Amat:
        # if n == self.root:
        #     continue
            for state,val in self.Amat[node].items():
                state_probs[state] += np.exp(val)
            
        state_probs[state_probs==0] = 0.01
        state_probs = state_probs/state_probs.sum()
  
        for node in self.Bmat:
            print(node)
            for child_state in self.Bmat[node]:
                for parent_state, val in self.Bmat[node][child_state].items():
                
                    if not np.isnan(val):
                        print(val)
                        new_tmat[parent_state, child_state] += np.exp(val)
        #  for n in self.postorder:
        #     if n == self.root:
        #         continue
        #     parent = self.parent(n)

        #     if self.is_leaf(n):
        #         for v in self.states:
        #             coeff = np.exp(self.Amat[parent][v])
        #             new_tmat[v,self.obs_states[n]]+=coeff
        #             exp_log_like += coeff*self.log_tmat[v, self.obs_states[n]]
        #     else:
        #         for v in self.states:
        #             for u in self.states:
        #                 coeff = np.exp(self.Bmat[n][u][v])
        #                 exp_log_like += coeff*self.log_tmat[v, u]
        #                 new_tmat[v,u]+= coeff
            
           
      
            #add pseudo counts for unobserved possible transitions 
        for s in self.states:
            for t in self.states:
                    if s <=t and new_tmat[s,t]==0:
                        new_tmat[s,t] = 0.01
        norm_term = new_tmat.sum( axis=1).reshape(-1,1)

        new_tmat = new_tmat/norm_term

        
        return exp_log_like, state_probs, new_tmat





    def fit(self, obs_states):
        self.obs_states = obs_states
        old_log_like = np.Inf
        for i in range(10):
            self.e_step()
            cur_log_like, state_probs, new_tmat = self.m_step()
            self.tmat = new_tmat
    
            self.log_tmat = np.log(self.tmat)

            if self.check_convergence(cur_log_like, old_log_like):
                break
            else:
                old_log_like = cur_log_like
        return cur_log_like, state_probs, new_tmat

    
# tree = nx.DiGraph()

# tree.add_edges_from([  (0,1), (0,2), (0,3), (1,4), (1,5), (2,6), (2,7), (3,8), (3,9)])

# states = [0,1,2,3,4,5,6]
# isotypes = {4:1, 5:2, 6:1, 7:2, 8:2, 9:2}
# init_probs = {s:0 for s in states}
# init_probs[0] =1
# state_probs = {s: 0 for s in states}

# for s in states:
#     state_probs[s] = 1/len(states)


# transmat = np.loadtxt("/scratch/projects/tribal/real_data/mouse_transmat2.txt")


# em = EMProbs(tree, 0, states, transmat, state_probs, init_probs )
# em.fit(isotypes)








    