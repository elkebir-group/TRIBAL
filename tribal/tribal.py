import numpy as np
from copy import deepcopy
from max_likelihood_trans_probs import MaxLike
import init_transmat as tm
from lineage_tree import LineageTree, LineageTreeList
from base_tree import ParsimonyForest
from multiprocessing import Pool


class Tribal:
    """
    Infer a B cell lineage tree for each clonotype and shared isotype transition probabilities. 

    Attributes
    ----------

    n_isotypes: int
       the number of isotype states
    
    alphabet: tuple
        the valid alphabet for BCR sequences (default:  ("A", "C", "G", "T","N", "-"))

    sankoff_cost_function: dict
        the cost function to use for running the Sankoff algorithm during ancestral BCR sequence reconstruction.
        the keys must be all pairs from the provided alphabet. defaults to using the standard cost
        function (1 for a mismatch, 0 for match) is no cost function provided.
    
    seed: int
        random number seed used to randomly downsample candidate trees within each iterations when size of a parsimony forests
        is greater than max_cand
        
    max_cand: int
        the maximum allowable size of a parsimony forest to consider within each coordinate descent iterations.
        if the size of each maximum parsimony forest is less than max_cand then downsampling does not occur.
   
     niter: int
        the maximum number of coordinate descent iterations to perform if convergence criteria is not met (default 10)
    
    threshold: float
        The tolerance for convergence of the CSR objective (default 0.5)
    
    restarts: int
        the number of restarts, i.e, different initialization of the isotype transition probabilities (default 5)

    stay_probs: tuple
        the lower bound and upper bound for the initialization of the probability of not class switching.
        the initalization values of this parameter are determined by usined np.linspace on this range with the
        restarts parameter. 

    mode: str
        One of 'refinement' or 'score' (default: 'refinement'). The mode determines how lineage trees 
        and ancestral isotype stands are inferred when optimzing the class-switch recombination objective.
        In 'refinment' mode, TRIBAL solves the Most Parsimonious Tree Refinement (MPTR) problem for each candidate
        tree in the maximum parsimony forest given fixed isotype transition probabilities. In 'score' mode, TRIBAL
        infers the ancestral isotype states using the Sankoff algorithm using a cost function with weights determined
        by the isotype transition probability matrix. 

    Notes
    ----
    TRIBAL takes as input a maximum parsimony forest (a set of trees that minimizes the SHM score for a 
    multiple sequenced alignment of the concatentated heavy and light chain (optional) variable
    region sequences) and the encoded isotype of each sequenced B cell for k clonotypes.   See the Preprocessor class for help preparing the input data. 
    It then infers a B cell lineage tree(s) for each clonotype that minimizes the somatic hypermutation (SHM) parsimony score and then maximizes the class swtich recombination (CSR)
    likelihood score. It also infers the optimal isotype transition probabilities that jointly maximizes the CSR likelihood.

    TRIBAL uses a coordinate ascent algorithm to alternately infer the optimal isotype transition probabilities and then
    infer a representative B cell lineage tree for each clonotype. This proceeds until the CSR likelihood objective
    convergences within a tolerance defined by threshold. Multiple restarts are performed with different initial
    isotype transition probabilities. The initializations are defined by the stay_probs, a tuple of floats the defines
    an upper and lower bound on the diagonal of the isotype transition probability matrix. Given a number of restarts and stay probs, 
    np.linspace(lower bound, upper bound, restarts) is used to define a stay probability, i.e., the probability of 
    not class switching, for each restart. The remain entries in each row are either 0 if the transition violates
    class switching constraints, e.g., IgG -> IgM, or initialized uniformly (1 - stay_prob)/# of valid transitions from start 
    state. 

    The default mode for TRIBAL if 'refinement', during the coordinate descent step TRIBAL will find the 
    most parsimonious tree refinement of each tree  in the maximum parsimony forest and ancestral isotypes given the current isotype transition probabilites.
    In score mode, TRIBAL will not modify the input maximum parsimony trees and will only infer the ancestral
    isotypes using the Sankoff algorithm with a cost function given by the negative log of the isotype
    transition probabilities. 


    If the size the of the maximum parsimony forest is very large, the parsimony forest can optionally
    be downsampled to a size of max_cand during each iteration of the algorithm. The best tree found so far is always included
    in the next iteration to ensure convergence. 



      
    Examples
    --------
    Here is an example of how to use the TRIBAL class::

        import tribal
        clonotypes = tribal.load(np_klh_data)

        tr = tribal.Tribal(7, max_cand=50, niter=1, restarts=3, theshold=1, mode='refinement')
        shm_obj, csr_obj, lineage_trees, transmat = tr.fit(clonotypes, cores=10, threads=1)

    """




    def __init__(self, 
                n_isotypes = 7, 
                alphabet= ("A", "C", "G", "T","N", "-"), 
                sankoff_cost_function = None, 
                seed= 1026, 
                max_cand=50, 
                niter=10,
                threshold=0.5, 
                restarts=5,
                stay_probs=(0.55,0.95), 
                mode="refinement",
                verbose =False ):

        self.mode = mode
        self.n_isotypes = n_isotypes 
        self.alphabet = alphabet
        self.sankoff_cost_function = sankoff_cost_function
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.stay_probs = stay_probs 
        self.min_iterations = min(4, niter)
        self.max_cand = max_cand
        self.states = np.arange(self.n_isotypes)
        self.threshold = threshold
        self.niterations = niter
        self.restarts = restarts 
        self.verbsose = verbose 



    

        #don't remember if I need these or not
        self.candidates = {}
       
        self.obs_states = {key: val.isotypes for key,val in self.clonotypes.items()}




     
    

    def intialize_candidates(self, best_tree_ids=None):
        '''
        randomly initialize a set of candidate trees up to max_cands 
        for each clonotype, including the best trees found so far is dictionary
        is given.'''
        candidates = {}
        for c in self.clonotypes:
            # print(f"clontoype {c}: size: {self.clonotypes[c].size()}")
            if self.clonotypes[c].size() > self.max_cand:
                cand = self.rng.choice(self.clonotypes[c].size(), self.max_cand, replace=False).tolist()
            else:
                cand = [i for i in range(self.clonotypes[c].size())]
            
            candidates[c]= ParsimonyForest( alignment=self.clonotypes[c].alignment, isotypes=self.clonotypes[c].isotypes)
            for index in cand:
                candidates[c].add(deepcopy(self.clonotypes[c][index]))
            if best_tree_ids is not None:
                for i in best_tree_ids[c]:
                    if i not in cand:
                        candidates[c].add(deepcopy(self.clonotypes[c][i]))

        return candidates 
    

    def csr_optimize(self, candidates, transmat=None, mode=None):
            if mode is None:
                mode = self.mode 
      
            #forest returns dict with clonotype as key and LineageTreeList as value
            all_scores = self.forest_mode(candidates, transmat, mode=mode)

            top_scores = {}
            total_likelihood = 0
            for c, sl in all_scores.items():
                best_obj, best_scores = sl.find_all_best_trees()
                top_scores[c] = best_scores
                total_likelihood += best_obj
            return total_likelihood, top_scores
          
    
    @staticmethod
    def check_convergence(old, new, threshold):
        return np.abs(old -new) < threshold


    def fit(self, clonotypes, mode="refinement", nproc=1, threads=1):
        """
        Run TRIBAL on an input data


        """
       
        self.nproc = nproc
        transmat = self.infer_probabilities(clonotypes, mode)
        csr_likelihood, best_scores = self.csr_optimize(clonotypes, transmat, mode)
        shm_score = 0
        for c, bst_lst in best_scores.items():
            alignment = clonotypes[c].alignment
            for lt in bst_lst:
                lt.ancestral_sequence_reconstruction(alignment, self.alphabet)
            shm_score += min(lt.shm_obj for lt in bst_lst)
        return  shm_score, csr_likelihood, best_scores, transmat


    


    def infer_probabilities(self):
  
     
        ''' 
            jointly infer isotype transition probabilities for a set of clonotypes
            and a B cell lineage tree that maximizes the CSR likelihood
            
        '''
        
        cand_tmat = []
        cand_scores = []
        best_trees = []
        stay_probs = np.linspace(*self.stay_probs, self.restarts)
        for i in range(self.restarts):
            if self.verbose:
                print(f"\nStarting Cycle {i}...")

            transmat =  tm.gen_trans_mat(stay_probs[i], self.n_isotypes)
            best_tree_ids = None 
            old_score = np.Inf
            
            for j in range(self.niterations):
 
                candidates = self.intialize_candidates(best_tree_ids)
                current_score, best_scores, best_tree_ids, _ = self.score_candidates(candidates, transmat=transmat)
                all_best_scores.append(best_scores)
            

                best_tree_ids[c] = [score.tree.id for score in best_scores]
        
                bs = best_scores[0]
                lin_tree = bs.tree
                leaf_isotypes = {l: bs.isotypes[l] for l in lin_tree.get_leafs() }
                msa = {l: bs.labels[l] for l in lin_tree.get_leafs() }
                leaf_isotypes[lin_tree.root] = bs.isotypes[lin_tree.root]
                msa[lin_tree.root] = bs.labels[lin_tree.root]
                refined_cands[c] = ParsimonyForest(msa, leaf_isotypes, [lin_tree])
               

                #return total_likelihood,all_best_scores, best_tree_ids, refined_cands 
                if self.verbose:
                    print(f"iteration: {j} old score: {old_score} current score: {current_score}")

                if self.check_convergence(current_score, old_score, self.threshold):
                        cand_tmat.append((transmat, state_probs))
                        cand_scores.append(current_score)
                        best_trees.append(best_scores)
                        break
                else:
                    old_score = current_score
                    transmat, state_probs = MaxLike(self.n_isotypes).infer(best_scores) 
        
        min_value = min(cand_scores)
        min_index = cand_scores.index(min_value)
        transmat, state_probs = cand_tmat[min_index]
        best_scores = best_trees[min_index]
        best_scores = LineageTreeList([b[0] for b in best_scores])
        
        return min_value, transmat, state_probs, best_scores

         

    #candidates should be a dict of lineage forests with clonotypes as key
    def forest_mode(self, transmat, candidates, mode="refinement"):
            
            arg_vals = []
            for c in candidates:
                lin_forest = candidates[c]
                isotype_labels = lin_forest.isotypes 
                for t in lin_forest.get_trees():
                    lt = LineageTree(clonotype=c, tree=t)
                    arg_vals.append((transmat, lt, isotype_labels ))


            if mode == "refinement":
                mode_func = self.refine
            else:
                mode_func = self.score 
     

            if self.nproc <= 1:
                results = [mode_func(*args) for args in arg_vals]
            else:
                with Pool(self.nproc) as pool:
                    results = pool.starmap(mode_func, arg_vals) 
            
            all_scores = {c: LineageTreeList() for c in candidates}
            for lt in results:
                all_scores[lt.clonotype].append(lt)
            
            return all_scores
      



    def score(self, transmat, lt, isotype_labels):
        lt.isotype_parsimony(isotype_labels, cost=transmat, states=self.states)
        return lt 




    def refine(self, transmat, lt, isotype_labels):
        lt.refinement(isotype_labels, cost=transmat)
        return lt
     

                      

      
        





        
    


  









