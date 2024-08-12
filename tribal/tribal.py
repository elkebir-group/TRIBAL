"""The TRIBAL algorithm."""

from multiprocessing import Pool
import numpy as np
from .max_likelihood_trans_probs import MaxLike
from .init_transmat import gen_trans_mat
from .lineage_tree import LineageTree, LineageTreeList
from .clonotype import Clonotype



class Tribal:
    """
    A class to infer a B cell lineage tree for each clonotype and shared isotype transition probabilities. 

    Attributes
    ----------
    n_isotypes : int
       the number of isotype states
    alphabet : tuple
        the valid alphabet for BCR sequences, defaults to  ("A", "C", "G", "T","N", "-")
    sankoff_cost_function : dict
        the cost function to use for running the Sankoff algorithm during ancestral BCR sequence reconstruction.
        the keys must be all pairs from the provided alphabet. defaults to using the standard cost
        function (1 for a mismatch, 0 for match)
    seed: int
        random number seed used to randomly downsample candidate trees within each iterations when size of a parsimony forests
        is greater than max_cand
    max_cand: int
        the maximum allowable size of a parsimony forest to consider within each coordinate descent iterations.
        if the size of each maximum parsimony forest is less than max_cand then downsampling does not occur.
     niter: int
        the maximum number of coordinate descent iterations to perform if convergence criteria is not met, defaults to 10)
    threshold : float
        The tolerance for convergence of the CSR objective, defaults to 0.5
    restarts : int
        the number of restarts, i.e, different initialization of the isotype transition probabilities, defaults to 10
    stay_probs : tuple
        the lower bound and upper bound for the initialization of the probability of not class switching.
        the initalization values of this parameter are determined by using np.linspace on this range with the
        restarts parameter, defaults to (0.55, 0.95)

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
        
        from tribal import Tribal, clonotypes

        #clonotypes dictionary includes the following clonotypes
        #["Clonotype_1036", "Clonotype_1050", "Clonotype_10884","Clonotype_1457", "Clonotype_755", "Clonotype_322"]
        
        #the clonotype data contains the following isotypes encoded from 0 to 7
        isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']
        tr = Tribal(n_isotypes=len(isotypes), verbose=True, restarts=2, niter=15)
        
        #run in refinement mode
        shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="refinement", cores=6)

    
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
                verbose =False ):


        self.n_isotypes = n_isotypes 
        self.alphabet = alphabet
        self.sankoff_cost_function = sankoff_cost_function
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.stay_probs = stay_probs 
        self.max_cand = max_cand
        self.threshold = threshold
        self.niterations = niter
        self.restarts = restarts 
        self.verbose  = verbose 

        if self.verbose:
            param_str  = "\nTRIBAL parameters:\n"
            param_str += f"number of isotypes: {self.n_isotypes}\n"
            param_str += f"sequence alphabet: {self.alphabet}\n"
            param_str += f"maximum candidates: {self.max_cand}\n"
            param_str += f"stay probabilities: {self.stay_probs}\n"
            param_str += f"random number seed: {self.seed}\n"
            param_str += f"restarts: {self.restarts}\n"
            param_str += f"iterations: {self.niterations}\n"
            param_str += f"convergence threshold: {self.threshold}\n"
            param_str += f"verbose: {self.verbose}\n"

            print(param_str)



    

     
    

    def _intialize_candidates(self, clonotypes, best_tree_ids=None):
        # '''
        # randomly initialize a set of candidate trees up to max_cands 
        # for each clonotype, including the best trees found so far is dictionary
        # is given.'''
       
        candidates = {}
        for c in clonotypes:
 
            if clonotypes[c].size() > self.max_cand:
                cand = self.rng.choice(clonotypes[c].size(), self.max_cand, replace=False).tolist()
            else:
                cand = [i for i in range(clonotypes[c].size())]
            
            candidates[c]= Clonotype(id = clonotypes[c].id, alignment=clonotypes[c].alignment, isotypes=clonotypes[c].isotypes)
            
            #TODO: remove deep copy
            for index in cand:
                candidates[c].add(clonotypes[c][index])
            if best_tree_ids is not None:
                for i in best_tree_ids[c]:
                    if i not in cand:
                        candidates[c].add(clonotypes[c][i])

        return candidates 
    

    def _csr_optimize(self, candidates, transmat=None, mode=None):

      
            #forest returns dict with clonotype as key and LineageTreeList as value
            all_scores = self._forest_mode(candidates, transmat)

            top_scores = {}
            total_likelihood = 0
            for c, sl in all_scores.items():
                best_obj, best_scores = sl.find_all_best_trees()
                top_scores[c] = best_scores
                total_likelihood += best_obj
            return total_likelihood, top_scores
          
    
    @staticmethod
    def _check_convergence(old, new, threshold):
        return np.abs(old -new) < threshold


    def fit(self, clonotypes, mode="refinement", transmat=None, cores=1):
        """
        Run TRIBAL on a dictionary of clonotypes and infer B cell lineage tree(s)
        for each clonotype and a shared istoype transition probability matrix. 

        Parameters
        ----------
        clonotypes: dict
            a dictionary of Clonotypes each containing a parsimony forest, isotypes and multiple sequence alignment    
        mode: str
            the mode for optimizing the class switch recombination (CSR) likelihood, one of ["refinement", "score"]. In
            'refinement' mode, TRIBAL solves the most parsiminious tree refinement (MPTR) problem for each 
            candidate tree in the parsimony forest. In 'score' mode, TRIBAL  infers the ancestral isotypes using
            the Sankoff algorithm with the weights coming from the isotype transition probabilities. 
        transmat: list
            a optional isotype transition probabilty matrix to infer a B cell lineage
            tree(s) for each clonotype. If not provided, the isotype transition probabilites
            are inferred from the data. 
            
        cores: int
            The number of cores to use (default 1)



        Examples
        --------
        Here are examples of how to run the fit function::

            from tribal import Tribal, clonotypes

            #clonotypes dictionary includes the following clonotypes
            #["Clonotype_1036", "Clonotype_1050", "Clonotype_10884","Clonotype_1457", "Clonotype_755", "Clonotype_322"]
            
            isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']
            tr = Tribal(n_isotypes=len(isotypes), verbose=True, restarts=2, niter=15)
            
            #run in refinement mode
            shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="refinement", cores=6)

            #run in scoring mode
            shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="score", cores=6)


            #given a user-specified isotype transition probability matrix 
            from tribal import probabilites
            shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes =clonotypes,
                                                                        transmat= probabilites,
                                                                        mode="refinement", cores=6)


                
        Returns
        -------

        shm_score: float
            a float with the total somatic hypermutation (SHM) parsimony score for all clonotypes
        
        csr_likelihood: float
            a float with the total class switch recombination (CSR) likelihood score for all clonotypes

        best_scores: dict 
            a dictionary of LineageTreeLists containing all optimal LineageTrees per clonotype

        transmat: np.array
            a numpy array of isotype transition probabilities


        """



       
        self.mode = mode 
        self.nproc = cores
        if len(clonotypes)==0:
            raise ValueError("The number of clonotypes is 0. Recheck data and try again.")
        if self.verbose:
            print(f"\nStarting TRIBAL with {len(clonotypes)} clonotypes...")
        
        if transmat is None:
            if self.verbose:
                print("Inferring isotype transition probabilities...")
            transmat = self._infer_probabilities(clonotypes)
        
        else:
            if self.verbose:
                print("Using provided isotype transition probabilities for CSR optimization...")
            if transmat.shape[0] != self.n_isotypes or transmat.shape[1] != self.n_isotypes:
                raise ValueError("User provided isotype transition probability matrix does not match the number of isotype states")
        if self.verbose:
            print("\nOptimizing the CSR Likelihood...")
        csr_likelihood, best_scores = self._csr_optimize(clonotypes, transmat=transmat)

        if self.verbose:
            
            print("Reconstructing the ancestral sequences...")
        
        shm_score = 0
        for c, bst_lst in best_scores.items():
            alignment = clonotypes[c].alignment
            min_score = np.Inf
            for lt in bst_lst:
                lt.ancestral_sequence_reconstruction(alignment, self.alphabet)
                if lt.shm_obj < min_score:
                    min_score = lt.shm_obj
            shm_score += min_score
        
        if self.verbose:
            print(f"SHM Score: {shm_score} CSR Likelihood: {csr_likelihood}")
            print("The TRIBE has spoken!")
        return  shm_score, csr_likelihood, best_scores, transmat


    


    def _infer_probabilities(self, clonotypes):
  
     
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
                print(f"\nStarting restart {i}...")

            transmat =  gen_trans_mat(stay_probs[i], self.n_isotypes)
            best_tree_ids = None 
            old_score = np.Inf
            
            for j in range(self.niterations):
 
                candidates = self._intialize_candidates(clonotypes, best_tree_ids)
                current_score, all_best_scores = self._csr_optimize(candidates, transmat=transmat)
                # all_best_scores.append(best_scores)
            

                best_tree_ids =  {c: all_best_scores[c]._get_ids() for c in all_best_scores}
        
                if self.verbose:
                    print(f"iteration: {j} old score: {old_score} current score: {current_score}")

                if self._check_convergence(current_score, old_score, self.threshold) or j == self.niterations -1:
                        cand_tmat.append(transmat)
                        cand_scores.append(current_score)
                        best_trees.append(all_best_scores)
                        break
                else:
                    old_score = current_score
                    transmat, _ = MaxLike(self.n_isotypes).infer(all_best_scores) 
        
        min_value = min(cand_scores)
        min_index = cand_scores.index(min_value)
        transmat = cand_tmat[min_index]
        
        return  transmat

         

    #candidates should be a dict of lineage forests with clonotypes as key
    def _forest_mode(self,  candidates, transmat):
            
            arg_vals = []
            for c in candidates:
                lin_forest = candidates[c]
                isotype_labels = lin_forest.isotypes 
                for t in lin_forest.get_forest():
                    lt = LineageTree(clonotype=c, tree=t)
                    arg_vals.append((transmat, lt, isotype_labels ))


            if self.mode == "refinement":
                mode_func = self._refine
            else:
                mode_func = self._score 
     

            if self.nproc <= 1:
                results = [mode_func(*args) for args in arg_vals]
            else:
                with Pool(self.nproc) as pool:
                    results = pool.starmap(mode_func, arg_vals) 
            
            all_scores = {c: LineageTreeList() for c in candidates}
            for lt in results:
                all_scores[lt.clonotype].append(lt)
            
            return all_scores
      



    def _score(self, transmat, lt, isotype_labels):
        lt.isotype_parsimony(isotype_labels, transmat=transmat)
        return lt 




    def _refine(self, transmat, lt, isotype_labels):
        lt.refinement(isotype_labels, transmat=transmat)
        return lt
     

                      

      
        





        
    


  









