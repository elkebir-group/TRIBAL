"""This module provides functions and classes to handle lineage trees."""
import os
from dataclasses import dataclass, field
from functools import total_ordering
import pickle
from typing import Dict, Tuple
import numpy as np
from .base_tree import BaseTree
from .small_parsimony import SmallParsimony
from .expansion_graph import ConstructGraph
from .mptr import MPTR
from .draw_tree import DrawTree
from .utils import save_dict, write_fasta

@dataclass
@total_ordering
class LineageTree:
    """
    A class to model B cell lineage trees.

    Attributes
    ----------
    clonotype : str
        The name of the clonotype.
    tree : BaseTree
        The rooted tree topology.
    csr_obj : float 
        The current CSR likelihood of the Lineage tree (default: 0).
    isotype : dict
        The isotype labels of the Lineage tree nodes.
    shm_obj : float 
        The current SHM Parsimony score of the Lineage tree (default: 0).
    sequences : dict
        The BCR sequences of the Lineage tree nodes.

    Notes
    -----
    A B cell lineage tree is a rooted tree with nodes labeled by BCR sequences (concatenated heavy and light chain) and
    isotypes. The somatic hypermutation (SHM) parsimony score is the total number of base substitutions within the sequences, and 
    the class switch (CSR) likelihood is the likelihood of the isotype labels given an isotype transition probability matrix.
    """
    
    clonotype: str
    tree: BaseTree= None
    csr_obj: float = 0
    isotypes: dict = field(default_factory=dict)
    shm_obj: float = 0
    sequences: dict = field(default_factory=dict)

    def __post_init__(self):
        """Initialize the objective tuple and root."""
        self.objective = (self.shm_obj,self.csr_obj)
        self.root = self.tree.root

    def _validate_item(self, item):
        """Check comparison items is an instance of LineageTree."""
        if isinstance(item, (LineageTree)):
            return item
        raise TypeError(
            f"Score expected, item got {type(item).__name__}"
        )
    def __eq__(self, __value: object) -> bool:
        """Check if two lineage tree items have the same objective score."""
        item = self._validate_item(__value)
        return self.objective ==item.objective
   
    def __lt__(self, __value: object) -> bool:
        """Check if a LineageTree is less than another."""
        item = self._validate_item(__value)
        return self.objective <item.objective

    def __str__(self):
        """To string method."""
        mystr = f"B cell lineage tree for {len(self.tree.get_leafs())} cells"
        mystr += f"\nRoot id: {self.root}"
        mystr += f"# of nodes: {len(self.tree.T.nodes)}"
        mystr +=  f"\nObjective\n----------\nSHM: {self.shm_obj}\nCSR: {self.csr_obj}"
        return mystr

    def to_pickle(self, fname):
        """Pickle the B cell lineage tree.

        Parameters
        ----------
        fname : str
            the path to where the pickled object should be saved,
        """
        with open(fname, 'wb') as file:
            pickle.dump(self, file)

    def _seq_len(self):
        for _, val in self.sequences.items():
            return len(val)
    
    def compute_csr_likelihood(self, transmat):
        """Compute CSR likelihood of a lineage tree for a given isotype transition probability matrix.

        Parameters
        ----------
        transmat : numpy.array
            The isotype transition probability used to compute the CSR likelihood.

        Returns
        -------
            The class switch recombination (CSR) likelihood.
        """
        transmat = -np.log(transmat)
       
        iso = self.isotypes
        score = 0
        try:
            nodes = self.tree.preorder_traversal()
            for n in nodes:
                t = iso[n]
                for c in self.tree.children(n):
                    s = iso[c]
                    score += transmat[t, s]
        except:
            raise ValueError("Invalid isotypes or tree. \
                             First run isotype parsimony or refinement functions.")
        
        self.csr_obj = score
        return self.csr_obj
    

    def get_id(self):
        """Get the internal id of the tree topology. Useful for mapping refined trees back the to the unrefined tree in the parsimony forest.

        Returns
        -------
        : int
            the internal id of the tree topology
        """
        return self.tree.id

    def ancestral_sequence_reconstruction(self, alignment,
                           alphabet=("A", "C", "G", "T","N", "-"), 
                           cost_function=None):
        """Infer the ancestral BCR sequences of the internal nodes given an alignment.

        Parameters
        ----------
        alignment : dict
            a dictionary with leaf labels and root id as keys and the BCR sequence as value.
            
        alphabet : tuple
            the valid alphabet for BCR sequences, default:  ("A", "C", "G", "T","N", "-")
        
        cost_function : dict|None
            the cost function for substitution of a single nucleotide base. If None, the 
            standard 0/1 cost function is used for matches and mismatches. If dictionary
            all pairs of the elements in the alpabet should be be includes in the keys.
        
        Examples
        --------
        Here is an example of how to reconstruct ancestral sequences::
        
        ```python  
            from tribal import clonotypes, LineageTree
            id = "Clonotype_1036"
            clonotype = clonotypes[id]
            forest = clonotype.get_forest()
            lt = LineageTree(id=id, tree = forest[0])
            shm_score, sequences = lt.ancestral_sequence_reconstruction(clonotype.alignment)
            print(lt)
        ```

        Returns
        -------
        : float
            a float with the somatic hypermutation (SHM) parsimony score for the lineage tree
        
        : dict
            a dictionary containing the BCR sequence labels of the lineage tree

        """     
        alignment = {k: list(alignment[k]) for k in alignment}
        sp = SmallParsimony(self.tree, 
                            alphabet= alphabet,
                            cost = cost_function)
        self.shm_obj, sequences = sp.sankoff(alignment)
        self.sequences = {key : "".join(value) for key, value in sequences.items()}

        return self.shm_obj, self.sequences
    


    def isotype_parsimony(self, isotype_labels:dict, transmat:np.array):
        """Infer the isotype of the B cell lineage tree using weighted parsimony.

        Parameters
        ----------
        isotype_labels : dict
            a dictionary with leaf labels and root id as keys and isotypes as values.  
        transmat : numpy.array
            the isotype transition probability used to compute the CSR likelihood.  

        Examples
        --------
        Here is an example of how to infer isotypes::

        ```python
            from tribal import clonotypes, probabilities, LineageTree

            id = "Clonotype_1036"
            clonotype = clonotypes[id]
            forest = clonotype.get_forest()
            lt = LineageTree(id=id, tree = forest[0] )
            csr_likelihood, isotypes = lt.isotype_parsimony(isotype_labels= clonotype.isotypes,
                                                                        transmat=probabilities )
            print(lt)
        ```

        Returns
        -------
        csr_obj : float
            a float with the class switch recombination likelihood score for the lineage tree.  
        isotypes : dict
            a dictionary containing the isotypes of the lineage tree.  
        """
        transmat = -np.log(transmat)
        states = list(range(transmat.shape[0]))
        sp = SmallParsimony(self.tree, alphabet=states,cost=transmat)
        self.csr_obj, self.isotypes = sp.sankoff(isotype_labels)


    def refinement(self, isotype_labels: Dict[str, str], transmat: np.ndarray) -> Tuple[float, Dict[str, str]]:
        """Solves the most parsimonious tree refinement problem (MPTR).

        Parameters
        ----------
        isotype_labels : dict
            A dictionary with leaf labels and root id as keys and isotypes as values.
        transmat : numpy.array
            The isotype transition probability used to compute the CSR likelihood.

        Examples
        --------
        Here is an example of how to refine a lineage tree:

        ```python
        from tribal import clonotypes, probabilities, LineageTree

        id = "Clonotype_1036"
        clonotype = clonotypes[id]
        forest = clonotype.get_forest()
        lt = LineageTree(id=id, tree=forest[0])
        csr_likelihood, isotypes = lt.refinement(isotype_labels=clonotype.isotypes, transmat=probabilities)
        print(lt)
        ```

        Returns
        -------
        float
            A float with the class switch recombination likelihood score for the lineage tree.  
            
        dict
            A dictionary containing the isotypes of the lineage tree.  
        """
        cost = -np.log(transmat)
        cg = ConstructGraph(cost, isotype_labels, root_identifier=self.root)
        fg = cg.build(self.tree)

        st = MPTR(fg.G,
                  self.tree.T,
                  fg.find_terminals(),
                  fg.iso_weights,
                  fg.tree_to_graph,
                  root=self.root)
        
        self.csr_obj, tree = st.run()

        tree, self.isotypes = cg.decodeTree(tree)
        self.tree = BaseTree(tree, self.root, self.tree.id, self.tree.name)

        return self.csr_obj, self.isotypes

    def draw(self,
            fname,
            isotype_encoding=None,
            show_legend=False,
            show_labels=True,
            hide_underscore=True,
            color_encoding = None,
            dot = False):
        """Visualization of the current B cell lineage tree saves as a png or pdf.

        Parameters
        ----------
        fname : str
            The filename where the visualization should be saved.  
        isotype_encoding : list
            The list of the isotype labels to use.  
        show_legend : bool
            Optionally display the legend of the isotype colors (default=True).  
        show_labels : bool
            label the nodes by the sequence label
        hide_underscore : bool
            internal nodes that undergo refinement will have an underscore and copy number appended
            to the label. Setting this to true hides the underscore during visualization and retains 
            only the original label.
        color_encoding : dict, optional
            optional dictionary that maps isotype encoding to a color, if None, the default
            color palette is used.
        dot : bool
            if the file should be saved as a dot file, otherwise it will be saved
            as a png or pdf, depending on the file exentsion of fname
        """
        parents = self.get_parents()
        dt = DrawTree(parents, 
                    self.isotypes,
                    show_legend=show_legend,
                    root=self.root,
                    isotype_encoding=isotype_encoding,
                    show_labels=show_labels,
                    hide_underscore=hide_underscore,
                    color_encoding=color_encoding)
        if not dot:
            dt.save(fname)
        else:
            dt.save_dot(fname)


    def postorder_traversal(self) -> list:
        """Perform a postorder traversal of the lineage tree."""
        return self.tree.postorder_traversal()


    def preorder_traversal(self) -> list:
        """Perform a preorder traversal of the lineage tree."""
        return self.tree.preorder_traversal()

    def parent(self,n):
        """Identify the parent of a specified node.

        Parameters
        ----------
        n: str | int
            id of query node.

        Returns
        -------
           the parent of query node n.
        """
        return self.tree.parent(n)
    
    def children(self, n):
        """
        Identify the set of children of a specified node.

        Parameters
        ----------
        n : str
            ID of the query node.  

        Returns
        -------
            A list of children of node `n`.
        """
        return self.tree.children(n)
    
    def is_leaf(self,n):
        """Check if node is a leaf."""
        return self.tree.is_leaf(n)
    
        
    def get_leafs(self):
        """
        Identify the leafset of the lineage tree.

        Returns
        -------
            the leafset of the lineage tree.

        """
        return self.tree.get_leafs()
    
    def get_parents(self):
        """Identify the part of each node in the lineage tree.

        Returns
        -------
            a mapping of each node in the lineage tree to its parent node.

        """
        return self.tree.get_parents()
    
   
    def save_tree(self,fname):
        """Write the parent dictionary of the lineage tree to a file.

        Parameters
        ----------
        fname : str
            filename where the file should be saved 
        """
        parents = self.get_parents()
        save_dict( parents, fname)
    
    def save_edges(self, fname):
        """Write the edge list of a lineage tree to a file.

        Parameters
        ----------
        fname : str
            filename to where edge list should be saved 
        """
        self.tree.save_edges(fname)

    
    def get_edge_dataframe(self):
        """Obtain the edge list of the lineage tree as a pandas.DataFrame."""
        return self.tree.get_edge_df()
    
    def write(self, outpath:str, isotype_encoding=None, tree_label=None):
        """Write the lineage tree data to files.

        Parameters
        ----------
        outpath : str
            the path to file the files should be written.  
        isotype_encoding : list, optional
            the ordered isotype labels  
        """
        if tree_label is None:
            tree_label = ""
        else:
            tree_label = f"{tree_label}"
        clono_name = self.clonotype
        if isotype_encoding is not None:
            isotypes = {}
            for key, iso in self.isotypes.items():
                if iso >=0 and iso < len(isotype_encoding):
                    isotypes[key] = isotype_encoding[iso]
                else:
                    isotypes[key] = iso
            isotype_encoding = isotype_encoding
  
        else:
            isotypes = self.isotypes
            isotype_encoding = None
    

        if not os.path.exists(outpath):
            os.makedirs(outpath)

        write_fasta(f"{outpath}/{clono_name}_sequences{tree_label}.fasta", self.sequences)
        save_dict(f"{outpath}/{clono_name}_isotypes{tree_label}.csv",isotypes)
        self.save_edges(f"{outpath}/{clono_name}_edge_list{tree_label}.txt")
        self.draw(f"{outpath}/{clono_name}_tree{tree_label}.png",
                    isotype_encoding=isotype_encoding,
                    hide_underscore=False,
                    show_legend=True)
    
class LineageTreeList(list):
    """Extends class list in order to store and manipulate list of LineageTrees."""

    def __str__(self) -> str:
        """To string method."""
        return f"{type(self)} with {len(self)} LineageTree(s)."
    
    def append(self, item):
        """Append a LineageTree to the LineageTreeList."""
        super().append(self._validate_item(item))
    
    #only Score objects are allowd in
    def _validate_item(self, item):
        if isinstance(item, (LineageTree)):
            return item
        raise TypeError(
            f"Score expected, item got {type(item).__name__}"
        )

    def write(self, fname:str, sep=","):
        """Write the objective scores of all LineageTrees in the list to a file.

        Parameters
        ----------
        fname : str
            path to where file should be written.
        sep : str
            the seperator (default: ",").
        """
        with open(fname,'w+') as file:
            file.write(f"id{sep}shm_score{sep}csr_likelihood\n")

            for score in self:
            
                file.write(f"{score.tree.id}{sep}{score.objective}{sep}{score.seq_score}{sep}{score.iso_score}{sep}\n")
    
    def find_best_tree(self):
        """Find the LineageTree with optimal score.
         
        If there are multiple optimal solutions, it returns the first one it finds.  
        See `find_all_best_trees` to get all optimal solutions in the LineageTreeList or `sample_best_trees` to randomly 
        sample from among the optimal LineageTrees.

        Returns
        ----------
        float
            the optimal CSR likelihood of the best LineageTree
        
        LineageTree
            a LineageTree with optimal CSR likelihood
        """
        min_score = min(self, key=lambda x: x.objective).csr_obj

        # Filter objects with the minimum score
        min_score_object = [obj for obj in self if round(obj.csr_obj, 5) == round(min_score,5)][0]
        return min_score, [min_score_object]
    
    def find_all_best_trees(self):
        """Find the LineageTree(s) with optimal scores.

        Returns
        -------
        float
            the optimal CSR likelihood of the best LineageTree.
        
        LineageTreeList
            a LineageTreeList of LineageTree with optimal CSR likelihood.
        """
        min_score = min(self, key=lambda x: x.csr_obj).csr_obj


        min_score_object = [obj for obj in self if round(obj.csr_obj, 5) == round(min_score,5)]
        min_score_object = LineageTreeList(min_score_object)
        return min_score, min_score_object
    
    def _get_all_trees(self):

        return [x.tree for x in self]
    
    def _get_ids(self):
        return [x.get_id() for x in self]
    
    def to_pickle(self, fname):
        """Pickle the LineageTreeList.

        Parameters
        ----------
        fname : str
            the path to where the pickled object should be saved.
        """
        with open(fname, 'wb') as file:
            pickle.dump(self, file)
    
    def sample_best_tree(self, rng=None, seed=1016):
        """Find a LineageTree with optimal score.
        
        If there are multiple optimal solutions, it randomly samples among the lineage trees with optimal solutions.  
        See `find_all_best_trees` to get all optimal solutions in the LineageTreeList or `find_best_tree` to randomly 
        sample from among the optimal LineageTrees.

        Parameters
        ----------
        rng : numpy.random.Generator
            a numpy random number generator to use for sampling. (default: None)
        seed : int
            a random number seed to use to initialize a numpy.random.Generator (default: 1016)
        
        Returns
        ----------
        LineageTree
            a randomly sampled LineageTree with optimal CSR likelihood
        """
        if rng is None:
            rng = np.random.default_rng(seed)
        _, best_scores = self.find_all_best_trees()
        sampled_index = np.random.choice(len(best_scores))
        return best_scores[sampled_index]
    

   
    def write_all(self, outpath, isotype_encoding=None):
        """Write the data for all LineageTrees in the list to files.

        Parameters
        ----------
        outpath : str
            the path to where the files should be written.
        isotype_encoding : list
            the ordered isotype labeles
        """
        for i,lt in enumerate(self):
            lt.write(outpath,isotype_encoding=isotype_encoding, tree_label=i)
    