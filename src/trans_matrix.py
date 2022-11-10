
import numpy as np
class TransMat:
    def __init__(self, rng=None, n_isotypes=6):
        
        self.n_isotypes = n_isotypes
        self.isotypes = np.arange(n_isotypes)
        
        if rng is None:
            self.rng =np.random.default_rng(1026)
        else:
            self.rng = rng
        self.trans_mat =np.zeros((n_isotypes, n_isotypes))
    


    def fit(self, jump_prob=0.5):     
        for i in range(self.n_isotypes):
            if i < self.trans_mat.shape[1]-1:
                self.trans_mat[i,i] = 1 - jump_prob
            else:
                self.trans_mat[i,i] =1

            number_of_paths = self.trans_mat.shape[1]- i -1
            for j in range(i+1, self.n_isotypes):
                self.trans_mat[i,j:] =jump_prob/number_of_paths
        
        return self.trans_mat

    def fit_local(self, jump_prob=0.5):
        for i in range(self.n_isotypes):
            if i < self.trans_mat.shape[1]-1:
                self.trans_mat[i,i] = 1 - jump_prob
                self.trans_mat[i, i+1] = jump_prob
            else:
                self.trans_mat[i,i] =1

        
        return self.trans_mat
    
  
    def fit_dirichlet(self,conc=2):
     
    
        k = self.trans_mat.shape[1]
        for i in range(self.n_isotypes):
            alpha =np.full(shape=(k-i), fill_value=conc)
            dir_sample = self.rng.dirichlet(alpha, 1).reshape(-1)
            self.trans_mat[i, i:] = dir_sample
        return self.trans_mat


