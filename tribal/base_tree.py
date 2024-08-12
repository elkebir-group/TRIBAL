"""A module to model tree topology."""

from dataclasses import dataclass
import networkx as nx
from pandas import DataFrame


@dataclass
class BaseTree:
    """Dataclass for the tree toplogy of a B cell lineage tree.

    Attributes
    ----------
    T : nx.Digraph
        a rooted tree 
    root : str
        the unique identifier of the root
    id : int, optional
        the identifier of the tree
    name : str, optional
        the name of the tree
    """

    T: nx.DiGraph
    root: str
    id: int = 0
    name: str = None
    

    def postorder_traversal(self) -> list:
        """
        Perform a postorder traversal of the tree.

        Postorder traversal means visiting the children nodes first, and then the parent node.

        Returns
        -------
        list
            List of nodes in postorder.
        """
        return list(nx.dfs_postorder_nodes(self.T, source=self.root))

    def preorder_traversal(self) -> list:
        """
        Perform a preorder traversal of the tree.

        Preorder traversal means visiting the parent node first, and then the children nodes.

        Returns
        -------
        list
            List of nodes in preorder.
        """
        return list(nx.dfs_preorder_nodes(self.T, source=self.root))

    def nodes(self):
        """Return the list of nodes of the tree."""
        return list(self.T.nodes)
    
    def edges(self):
        """Return the list of edges of the tree."""
        return list(self.T.edges)
    
    def parent(self, n):
        """Get the parent node of the given node.

        Parameters
        ----------
        n : node
            The node for which to find the parent.

        Returns
        -------
        str | None
            The parent node if it exists, otherwise None.
        """
        preds = list(self.T.predecessors(n))
        if len(preds) == 0:
            return None
        
        return preds[0]

    def relabel(self, labels):
        """Relabels a tree in place.
        
        Parameters
        ----------
        labels : dict
            a dictionary with current nodes as keys and new labels as values. May include only 
            a subset of nodes.
        """
        self.T = nx.relabel_nodes(self.T, labels)
        if self.root in labels:
            self.root = labels[self.root]
    
    def children(self,n):
        """Return the set of children of a specified node `n`."""
        return list(self.T.neighbors(n))
    
    def is_leaf(self, node):
        """Check if node is a leaf of the tree.

        Parameters
        ----------
        node : str | int
            a node in the lineage tree
        
        Returns
        -------
        bool
            leaf status of the specified node 
        """
        return self.T.out_degree(node) ==0
    
    def get_leafs(self):
        """Return the leafset of the lineage tree."""
        return [n for n in self.T if self.is_leaf(n)]
    
    def is_root(self, node):
        """Check if node is the root.

        Parameters
        ----------
        node : str | int
            a node in the lineage tree
        
        Returns
        -------
        bool
            indicator if the node is the root

        """
        if type(node) == int:
            node =str(node)
        return node == self.root
    
    def get_parents(self):
        """Obtain the parent of each node.
        
        Returns
        -------
        dict
            node as key and parent node as child. "" indicates no parent, i.e., the root.
        """
        parents = {}
        for n in self.T:
            parent = self.parent(n)
            if parent is not None:
                parents[n] = parent
        return parents
    
    def set_id(self,id):
        """Update the id of the tree.
        
        Parameters
        -------
        id : str,int
            the new id of the tree
        """
        self.id =id
    
    def set_name(self, name):
        """Update the name of the tree.
        
        Parameters
        -------
        name : str
            the new id of the tree
        """
        self.name= name
    
    @staticmethod
    def find_leaf_descendants(node, graph):
        """Obtain all descendants of a node in a graph that is a leaf.
        
        Parameters
        ----------
        node : str, int
        graph : nx.DiGraph

        Returns
        -------
            the set of leaf descendants

        """
        leaf_descendants = set()

        # Helper function to traverse the graph
        def dfs(current_node):
            nonlocal leaf_descendants
            # If the current node is a leaf, add it to the set
            if graph.out_degree(current_node) == 0:
                leaf_descendants.add(current_node)
            else:
                # Traverse all child nodes recursively
                for child_node in graph.successors(current_node):
                    dfs(child_node)

        # Start the depth-first search from the specified node
        dfs(node)
        return leaf_descendants

    def get_clade_set(self, tree):
        """Get the clades of the tree."""
        clade_set = []
        for node in tree:
            clade_set.append(self.find_leaf_descendants(node, tree))
    
        return(set(map(frozenset, clade_set)))
    

    def save_edges(self, fname):
        """Write the edge list to a file.

        Parameters
        ----------
        fname : str
            the file where the edge list should be saved.
        """
        with open(fname, "w+") as file:
            for u,v in self.T.edges:
                file.write(f"{u},{v}\n")
    
    def get_edge_df(self):
        """Get the edge list as a pandas.DataFrame."""
        u_list =[u for u,v in self.T.edges]
        v_list = [v for u,v in self.T.edges]
        return DataFrame({'parent': u_list, 'child': v_list})
