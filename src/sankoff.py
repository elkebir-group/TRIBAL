import networkx as nx
import numpy as np

class Sankoff:
    def __init__(self, T, cost, alphabet, root):
        
        self.T = T
        self.root = root
        self.cost = cost
        self.alphabet = alphabet

        self.postorder = self.postorder_traversal()
        self.nodes = list(self.T.nodes())

      

        self.seq= nx.get_node_attributes(self.T, "sequence")
    

    def postorder_traversal(self):
        return list(nx.dfs_postorder_nodes(self.T, source=self.root))
    
    def is_leaf(self, node):
        return self.T.out_degree(node) ==0

    def sankoff_dp(self, pos):
        dp_matrix = {n : {}  for n in self.nodes}
        #TODO: adjust root for specified sequence
        for n in self.postorder:

            if self.is_leaf(n):
                for a in self.alphabet:

                    if self.seq[n][pos] == a:
                        dp_matrix[n][a] = 0
                    else:
                        dp_matrix[n][a] = np.Inf
            
            else:
                
          
                for a in self.alphabet:
                    score_array = np.array([self.min_cost(c,a, dp_matrix) for c in self.T.successors(n)])
                    dp_matrix[n][a] = score_array.sum()

        opt_score = np.array([dp_matrix[0][a] for a in self.alphabet]).min()
        
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
        labels[0] = self.argmin(dp_matrix[self.root])
        for c in self.T.successors(self.root):
            self.traceback_helper(labels[0], c, dp_matrix, labels)
        
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
        seq_assign = {k : '' for k in nodes}
        for k in nodes:
            seq_assign[k] += "".join([all_labs[i][k] for i in all_labs])
        return seq_assign


    def run(self):
        seq_length = len(self.seq[self.root])
        all_labels = {}
 
        min_scores =np.zeros(seq_length, dtype=int)
        for i in range(seq_length):
            min_scores[i], dp_mat = self.sankoff_dp(i)
            label_dict = self.traceback(dp_mat)
            all_labels[i] = label_dict
        
        
        
        node_labels= self.concat_labels(all_labels, self.nodes)

        return min_scores.sum(),node_labels

            

    def min_cost(self, child, to_state, dp_matrix):
   

        scores = np.array([self.cost[to_state,a] + dp_matrix[child][a]for a in self.alphabet])
        
        return scores.min()


if __name__ == "__main__":

    tree = nx.DiGraph()
    tree.add_edges_from([(0,1), (0,4), (1,2), (1,3)])

    cost_mat = {}
    alphabet = ("a", "g", "c", "t")
    seq = {2: "cc", 3: "gg", 4: "tt", 0:"cc"}
    for a in alphabet:
        for b in alphabet:
            if a == b:
                cost_mat[a,b] =0
            elif a in ["a", "g"] and b in ["a", "g"]:
                cost_mat[a,b] = 1
            elif a in ["c", "t"] and b in ["c", "t"]:
                cost_mat[a,b] =1
            else:
                cost_mat[a,b] =3
    nx.set_node_attributes(tree, seq, "sequence")

    sk = Sankoff(tree, cost_mat, alphabet, 0)
    opt_score, labels = sk.run()

    

    
