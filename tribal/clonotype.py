"""A module for modeling clonotype data."""

import pickle
from dataclasses import dataclass, field
from .base_tree import BaseTree

@dataclass
class Clonotype:
    """A class to model the clonotype input data. 

    Attributes
    ----------
    id : str
        the label of the clonotype
    alignment : dict
        the multiple sequence alignment for the BCR sequences
    isotypes : dict
        the encoded isotypes of the sequenced B cells
    forest : list
        a list of networkx.Digraphs containing the maximum parsimony forest for a 
        the multiple sequence alignment
    isotype_encoding : list
        an ordered list of the labels of the isotypes. This is important because naming
        conventions of isotype states, i.e., IgM, M, IghM, ighm, vary across datasets.

    """
    id: str   #id label of the clonotype
    alignment: dict = None
    isotypes: dict = None
    forest: list = field(default_factory=list)
    isotype_encoding: list = None

    def generate_from_list(self, tree_list, root=None):
        """Populate the parsimony forst from a list of networkx digraphs.
        
        Parameters
        ----------
        tree_list : list
            a list of nx.DiGraphs containing the trees to populate the parsimony forest
        root : str, optional
            the root id of each B cell lineage tree
        
        """
        for i,t in enumerate(tree_list):
            if isinstance(t, BaseTree):
                t.set_id(i)
                self.add(t)
            else:
                self.add(BaseTree(t,root,i))

    def add(self, tree):
        """Add a tree to the parsimony forest.
        Parameters
        ----------
        tree : nx.DiGraph
            a tree to add to the parsimony forest
        """
        self.forest.append(tree)

    def __getitem__(self, key):
        """Slice a tree in the parsimony forest."""
        return self.forest[key]
    
    def size(self):
        """Return the size of the parsimony forest."""
        return len(self.forest)
  
    def get_forest(self):
        """Return the parsimony forest."""
        return self.forest
    
    def save(self, fname):
        """Pickle the clonotype.
        
        Parameters
        ----------
        fname : str
            the filename where the clonotype should be pickled.
        """
        with open(fname, 'wb') as file:
            pickle.dump(self, file)


