
import networkx as nx
import numpy as np


from copy import deepcopy
import utils as ut
from itertools import repeat
from lineage_tree import LineageTree, LineageTreeList
from base_tree import Clonotype, BaseTree
from multiprocessing import Pool
from steiner_tree import ConstructGraph, SteinerTree



class TreeSelection:
    def __init__(self, 
                isotype_weights=None, 
                alpha=0.9, 
                n_isotypes=7,
                cost_function=None, 
                alphabet= ("A", "C", "G", "T","N", "-"), 
                nworkers=1, 
                root_id="naive",
                reversible=False,
                compute_seq= False ):

        

        


        if isotype_weights is None:

            if not reversible:
                rev_val = np.Inf
            else:
                rev_val = 0

         
            print("isotype weights is none, using standard cost function")
                #use unweighted sankoff cost function
            self.states = [i for i in range(n_isotypes)]
            self.iso_weights = {}

            for s in range(n_isotypes):
                for t in range(n_isotypes):
                    if s > t:
                        self.iso_weights[s,t] = rev_val
                    elif s < t:
                        self.iso_weights[s,t] = 1
                    else:
                        self.iso_weights[s,t] =0



        else:
            if type(isotype_weights) is np.ndarray:
                self.iso_weights, self.states = ut.convert_transmat_to_weights(isotype_weights)
            else:
                self.iso_weights = isotype_weights
                self.states = list(set([s for s,t in self.iso_weights]))


        self.alpha=alpha
        self.alphabet = alphabet

        self.cost_function = cost_function
        self.nworkers = nworkers
        self.root_id = root_id
        self.compute_seq= compute_seq
      

    
    def score(self, lin_tree, alignment, isotype_labels):
        
        seq_score, seq_labels = lin_tree.sequence_parismony(alignment, 
                                                    alphabet=self.alphabet, 
                                                    cost_function=self.cost_function, compute=self.compute_seq)
        iso_score, iso_labels = lin_tree.isotype_parsimony(isotype_labels,  weights=self.iso_weights, states=self.states)
        obj = self.compute_score(seq_score, iso_score)



        return   LineageTree(obj, seq_score, iso_score, seq_labels, iso_labels, lin_tree)


    def refinement_heuristic(self, lin_tree, alignment, isotype_labels):
        iso_score, iso_labels = lin_tree.isotype_parsimony_polytomy(isotype_labels, weights=self.iso_weights, states=self.states )

        seq_score, seq_labels = lin_tree.sequence_parismony(alignment, 
                                                    alphabet=self.alphabet, 
                                                    cost_function=self.cost_function)
        obj = self.compute_score(seq_score, iso_score)
        return LineageTree(obj, seq_score, iso_score, seq_labels, iso_labels, lin_tree)

    def refinement(self, lin_tree, alignment, isotype_labels):

  
            cg = ConstructGraph(self.iso_weights, isotype_labels, root_identifier=self.root_id)
           

            seq_score_prior, seq_labels = lin_tree.sequence_parismony(alignment)
            fg = cg.build(lin_tree, seq_labels)

            st = SteinerTree(fg.G, lin_tree.T, fg.find_terminals(), fg.seq_weights, 
                             fg.iso_weights,fg.node_mapping, fg.tree_to_graph,
                               fg.node_out_degree, pars_score = seq_score_prior, root=self.root_id, lamb=self.alpha, threads=1 )
            obj, tree = st.run()

            out_tree, out_iso = cg.decodeTree(tree)
            

            out_lt = BaseTree(out_tree, "naive", lin_tree.id, lin_tree.name)
    
            if len(alignment) > 0:
                seq_score, seq_labels = out_lt.sequence_parismony(alignment)
                assert seq_score_prior ==seq_score
            else:
                seq_labels = {}
                seq_score = 0



            sc = LineageTree(obj, seq_score, obj, seq_labels, out_iso, out_lt)
            sc.check_score(self.iso_weights)

            return sc 

      


    
    def compute_score(self, pars, iso):
        return self.alpha*pars + (1-self.alpha)*iso




        
            


    def forest_mode_loop(self, lin_forest, alignment=None, isotypes=None, mode="score"):
            
    

            all_results = {}
            
            if alignment is None:
                alignment = lin_forest.alignment
            if isotypes is None:
                isotype_labels = lin_forest.isotypes 
            else:
                isotype_labels = isotypes

            if mode =="refine":
                mode_func = self.refinement_heuristic
            elif mode == "refine_ilp":
                mode_func = self.refinement
            elif mode == "search":
                mode_func = self.search 
            else:
                mode_func = self.score 
     
            all_results = [] 
            for lin_tree in lin_forest.get_trees():
         
                result = mode_func(lin_tree,alignment, isotype_labels)
                all_results.append(  result)
                       
                      

      
            return  all_results



    def forest_mode(self, lin_forest, alignment=None, isotypes=None, mode="score", ntrees=1):
            
         

            all_results = {}
            
            if alignment is None:
                alignment = lin_forest.alignment
            if isotypes is None:
                isotype_labels = lin_forest.isotypes 
            else:
                isotype_labels = isotypes

            if mode =="refine":
                mode_func = self.refinement_heuristic
            elif mode == "refine_ilp":
                mode_func = self.refinement
            elif mode == "search":
                mode_func = self.search 
            else:
                mode_func = self.score 

        
            with Pool(self.nworkers) as pool:
                all_results = pool.starmap(mode_func,zip(lin_forest.get_trees(), repeat(alignment), repeat(isotype_labels)))     

    
            return all_results



    