import numpy as np
import networkx as nx 
import itertools
from scipy.special import logsumexp
class EMProbs:
    def __init__(self, lineage_forest, transmat, states=None,  maxiter=100, epsilon=0.1) -> None:

        self.lin_trees = lineage_forest.forest
        self.tmat = transmat 
        if states is None:
            self.states = np.arange(self.tmat.shape[0])
        else:
            self.states = states 
       
        self.log_tmat = np.log(self.tmat)
     
        self.maxiter = maxiter
        self.epsilon = epsilon
  
    
        #intialize dynamic programming dictionaries
    
        self.beta = {t.name: {n: {} for n in t.nodes()} for t in self.lin_trees}
        self.alpha = {t.name: {n: {} for n in t.nodes()} for t in self.lin_trees}
        self.beta_parent = {t.name: {n: {t.parent(n): {}} for n in t.nodes()} for t in self.lin_trees}

        self.beta_subtree_minus_child = {t.name: {n: {c: {} for c in t.children(n) }for n in t.nodes()}for t in self.lin_trees}


        self.Amat = {t.name: {n: {} for n in t.nodes()} for t in self.lin_trees}
        self.Bmat =  {t.name: {n: {s: {} for s in self.states} for n in t.nodes()} for t in self.lin_trees}



        

   
    

    def upward_recursion(self, tree):
        id = tree.name 
        for n in tree.postorder_traversal():
 
            if tree.is_leaf(n):
                for j in self.states:

                    self.beta[id][n][j] = np.log(int(self.obs_states[n]==j))

            else:
                for c in tree.children(n):
                    for m in self.states:
                        arr= []
                        for j in self.states:
                            arr.append(self.log_tmat[m,j] + self.beta[id][c][j])
                        
                        self.beta_parent[id][c][n][m] = logsumexp(arr)
                
                for m in self.states:
                    self.beta[id][n][m] = np.sum([self.beta_parent[id][c][n][m] for c in tree.children(n)])
                
                
                for c in tree.children(n):
                    for m in self.states:
                        self.beta_subtree_minus_child[id][n][c][m] = self.beta[id][n][m] - self.beta_parent[id][c][n][m]
            
    def downward_recursion(self, tree):
        id = tree.name 
        for n in tree.preorder_traversal():
            if tree.is_root(n):
                for m in self.states:
                    self.alpha[id][n][m] = np.log(int(m==0))

            else:
                parent = tree.parent(n)
          
                for m in self.states:
                    arr = []
                    for j in self.states:
                        arr.append(self.alpha[id][parent][j] + self.log_tmat[j,m] + self.beta_subtree_minus_child[id][parent][n][j])
                    arr =[a for a in arr if not np.isnan(a)]
                    if len(arr) > 0:
                        self.alpha[id][n][m] = logsumexp(arr)
                    else:
                        self.alpha[id][n][m] = np.NINF

    def compute_conditionals(self, tree):
        id = tree.name 
        for n in tree.preorder_traversal():
  
           
            for m in self.states:
                arr = []
                for j in self.states:
                    arr.append(self.alpha[id][n][j] + self.beta[id][n][j])
                

                self.Amat[id][n][m] =self.alpha[id][n][m] + self.beta[id][n][m] - logsumexp(arr)
                if not tree.is_root(n):
                    parent= tree.parent(n)
                    for j in self.states:
                        if self.log_tmat[j,m] != np.NINF:
    
                            self.Bmat[id][n][m][j] = self.alpha[id][parent][j] + self.beta_subtree_minus_child[id][parent][n][j] + \
                                                self.beta[id][n][m] + self.log_tmat[j,m] - logsumexp(arr)
                        else:
                            self.Bmat[id][n][m][j] = np.NINF 
    

        

    def e_step(self, obs_states):
        for tree in self.lin_trees:
            self.obs_states = obs_states[tree.name]
            self.upward_recursion(tree)

            self.downward_recursion(tree)
            self.compute_conditionals(tree)

        
    def check_convergence(self, curr, old ):
        return  np.abs(curr - old ) <= self.epsilon
      

    
    def  m_step(self):
    
       
        #log like expectation of for root being in state 0

        
        exp_log_like = 0
        new_tmat = np.zeros((len(self.states), len(self.states)))
        state_probs = np.zeros(len(self.states))

        for t in self.Amat:
            for node in self.Amat[t]:
                for state,val in self.Amat[t][node].items():
                    state_probs[state] += np.exp(val)
            
        state_probs[state_probs==0] = 0.01
        state_probs = state_probs/state_probs.sum()

        for t in self.Bmat:
            
            for node in self.Bmat[t]:
        
                for child_state in self.Bmat[t][node]:
                    for parent_state, val in self.Bmat[t][node][child_state].items():
                        if parent_state > child_state:
                            continue
                        if not np.isnan(val):
                            coeff = np.exp(val)
                            new_tmat[parent_state, child_state] += coeff
                            new_term = coeff* self.log_tmat[parent_state, child_state]
                            if np.isnan(new_term):
                                print("here")
                            exp_log_like += coeff* self.log_tmat[parent_state, child_state]

  
        #add pseudo counts for possible unobserved transitions 
        for s in self.states:
            for t in self.states:
                    if s <=t and new_tmat[s,t]==0:
                        new_tmat[s,t] = 1
        norm_term = new_tmat.sum( axis=1).reshape(-1,1)

        new_tmat = new_tmat/norm_term

        
        return exp_log_like, state_probs, new_tmat





    def fit(self, obs_states):
    
        old_log_like = np.Inf
        for i in range(self.maxiter):
            self.e_step(obs_states)
            cur_log_like, state_probs, new_tmat = self.m_step()
            self.tmat = new_tmat
    
            self.log_tmat = np.log(self.tmat)
            print(f"Iteration: {i} Current LogLike: {cur_log_like} Prior LogLike: {old_log_like}")
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








    