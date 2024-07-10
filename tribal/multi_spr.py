import networkx as nx

class MultSPR:
    def __init__(self, ttree, isotypes, states=range(9)):

        self.T = ttree.T
        self.root = ttree.root
        self.nodes = ttree.nodes()
        self.isotypes = isotypes 



        self.pseudo_root = list(self.T.neighbors(self.root))[0]

        self.cand_nodes = [n for n in self.nodes if n != self.root and n != self.pseudo_root]

        # self.generate_valid_cand()


        self.pre_order = list(nx.dfs_preorder_nodes(self.T, source=self.root))
        self.has_state = {s: [] for s in states}
        self.has_child_with_state = {s: [] for s in states}
        self.has_desc_with_state = {s: [] for s in states}
        for key, val in isotypes.items():
            self.has_state[val].append(key)
        
        for n in self.pre_order:
            if n != self.root:
                for c in self.T.neighbors(n):
                    child_state = isotypes[c]
                    self.has_child_with_state[child_state].append(n)

        
        for n in self.pre_order:
            if n != self.root:
                for d in nx.dfs_preorder_nodes(self.T, source=n):
                    desc_state = isotypes[d]
                    self.has_desc_with_state[desc_state].append(n)
        self.cand_moves = self.generate_cand()

 

        # print("done")


    def check_valid(self, node, state):
        return (node in self.has_state[state] or node in self.has_desc_with_state[state]) and not self.is_leaf(node) 
    
    def get_parent(self, n):
        return list(self.T.predecessors(n))[0]
    
    def is_leaf(self, n):
        return self.T.out_degree(n) ==0


    def not_neighbors(self, n, u):
        return self.get_parent(n) != self.get_parent(u) 

      
    
    def generate_cand(self):
        cand_moves = []
        for n in self.cand_nodes:
        
            node_state = self.isotypes[n]
            for u in self.cand_nodes:
                if self.get_parent(n) != u and u != n:
                    if self.check_valid(u,node_state) and self.not_neighbors(n,u) and  self.get_parent(n) != u:
                        cand_moves.append((n, u))
            #also allow every node to attempt to be a direct descendent of the the pseudo root
            if self.get_parent(n) != self.pseudo_root:
                cand_moves.append((n, self.pseudo_root))
        return cand_moves



    
    def spr(self, cut_node, regraft_node):
        
        tree = self.T.copy()

        #TODO: if regraft node is a descendent of cut node: 
        
        descendents = list(nx.dfs_preorder_nodes(self.T, source=cut_node))  
        if regraft_node in descendents:
            grandparent = self.get_parent(cut_node)
            new_parent = nx.shortest_path(self.T, source=cut_node, target=regraft_node)[1]
            tree.remove_edge(cut_node, new_parent)
            tree.remove_edge(grandparent, cut_node)
            tree.add_edge(grandparent, new_parent)
            tree.add_edge(regraft_node, cut_node)

            neighbors = list(tree.neighbors(cut_node))
            if len(neighbors) ==1:
                parent = list(tree.predecessors(cut_node))[0]
                tree.remove_node(cut_node)
                tree.add_edge(parent, neighbors[0])



        # if self.get_parent(regraft_node) == cut_node:
        #     grandparent = self.get_parent(cut_node)
        #     tree.remove_edge(cut_node, regraft_node)
        #     tree.remove_edge(grandparent, cut_node)
        #     tree.add_edge(grandparent, regraft_node)
        #     tree.add_edge(regraft_node, cut_node)

        
        else:

            cut_node_parent = self.get_parent(cut_node)
            cut_node_grandparent = self.get_parent(cut_node_parent)
            tree.remove_edge(cut_node_parent, cut_node)
            tree.add_edge(regraft_node, cut_node)
            
            #check for an remove any degree 2 nodes
            cut_node_parent_children = list(tree.neighbors(cut_node_parent))        
            if len(cut_node_parent_children) == 1:
        
                tree.remove_node(cut_node_parent)
                tree.add_edge(cut_node_grandparent, cut_node_parent_children[0])
        
        # for n in tree:
        #     if (n,n) in list(tree.edges):
        #         print(f"self loop: cut node: {cut_node} regraft node: {regraft_node}")

     
        # if not nx.is_tree(tcopy):
        #     print(f"not a tree: cut node: {cut_node} regraft node: {regraft_node}")

        # if not nx.is_connected(tcopy):
        #     print(f"not connected: cut node: {cut_node} regraft node: {regraft_node}")
        
        return tree
        



        
    def __iter__(self):
       ''' Returns the Iterator object '''
       return MultSPRIterator(self)

class MultSPRIterator:

    def __init__(self, spr):
        self.spr = spr

    def __next__(self):



        if len(self.spr.cand_moves) > 0:
            cut_node, regraft_node  = self.spr.cand_moves.pop()
            # print(f"cut node {cut_node} regraft node {regraft_node}")
            spr_tree = self.spr.spr(cut_node, regraft_node)
            return spr_tree

       
        else:
            raise StopIteration
        
    


  
        


    # def spr(self):
    #     cand_trees = []
    #     cand_nodes = [n for n in self.nodes if n != self.root ]

    #     for n in list(self.T.successors(self.root)):
    #         case2_trees = self.case2_spr(n)
    #         cand_trees += case2_trees

    #     for n in cand_nodes:

    #         if n not in list(self.T.successors(self.root)):
                
           
    #             case1_trees = self.case1_spr(n)
    #             case3_trees = self.case3_spr(n)
    #             cand_trees += (case1_trees + case3_trees)

    
    #     to_delete = []
    #     for i, c1 in enumerate(cand_trees):
    #         for j, c2 in enumerate(cand_trees):
    #             if i < j:
    #                 if set(list(c1.edges)) == set(list(c2.edges)):
    #                     to_delete.append(j)
                      
    #     cand_tree = [cand_trees[i] for i in range(len(cand_trees)) if i not in to_delete]       
    #     return cand_tree



# test_tree = nx.DiGraph()
# test_tree.add_edges_from([(0,6), (0,8), (6,7), (6,1), (6,7), (7,2), (7,3), (8,4), (8,5)])
# #states = {n :0 for n in test_tree.nodes}
# #test_tree.add_edges_from([(0,1), (0,2), (1,3), (1,4), (1,5), (2,6), (2,7), (0,8)])
# states = {0 : 0, 1 : 1, 2: 2, 3: 2, 4:0, 5:2, 6:1, 7:2, 8:0}



# nx.set_node_attributes(test_tree, states, "isotype")

# tm = TreeMoves(test_tree, 0)
# one_spr_list = tm.spr()

# tm.root_tree(test_tree,0)
# tm.all_nni()
# tm.nni_moves(test_tree, (2,3))




        

    # def check_nni_improvement(unrooted_tree, edge):
    #     return True 


    # def all_nni(self):
    #     nni_trees = []
    #     unrooted_tree  = self.T.to_undirected()
     
    #     cand_edges = self.find_internal_edges(unrooted_tree)
    #     for c in cand_edges:
    #         if check_improvement(unrooted_tree, c):
    #             nni_trees.append(nni_move(unrooted_tree,c ))
            

    # @staticmethod
    # def find_internal_edges(unrooted_tree):
    #     internal_vertices = [n for n in unrooted_tree.nodes if unrooted_tree.degree[n] >= 2]
    #     internal_edges = []
    #     for v in internal_vertices:
    #         for u in unrooted_tree.neighbors(v):
    #             if unrooted_tree.degree[u] >= 2:
    #                 if v < u:
    #                     if (v,u) not in internal_edges:
    #                         internal_edges.append((v,u))
    #                 else:
    #                     if not (u,v) in internal_edges:
    #                         internal_edges.append((u,v))



    #     return internal_edges
    

    # @staticmethod
    # def nni_moves(unrooted_tree, edge):
    #     out_trees = [] 
    #     left, right = edge

    #     unrooted_tree.remove_edge(left, right)

    #     tree1 = unrooted_tree.copy()
    #     tree2 = unrooted_tree.copy()
    #     e = find_internal_edges(tree1)
    #     left_subtree_roots = list(unrooted_tree.neighbors(left))
    #     right_subtree_roots = list(unrooted_tree.neighbors(right))

    #     if len(left_subtree_roots) >= 1 and len(right_subtree_roots) >= 1:


    #         tree1.remove_edge(right,right_subtree_roots[0])
    #         tree1.remove_edge(left, left_subtree_roots[0])


    #         tree1.add_edge(left,right_subtree_roots[0])
    #         tree1.add_edge(right,left_subtree_roots[0])
    #         tree1.add_edge(left,right)
    #         out_trees.append(tree1)

    #     if len(left_subtree_roots) ==2  and len(right_subtree_roots)==2:
    #         tree2.remove_edge(right,right_subtree_roots[0])
    #         tree2.remove_edge(left, left_subtree_roots[1])

    #         tree2.add_edge(right, left_subtree_roots[1])
    #         tree2.add_edge(left, right_subtree_roots[0])
    #         tree2.add_edge(left,right)
    #         out_trees.append(tree2)
        
    #     return out_trees
        





    














    # @staticmethod
    # def root_tree(unrooted_tree, root):
    #     #do a bfs on the graph 
    #     rooted_tree = nx.bfs_tree(unrooted_tree, source=root)

    #     return rooted_tree



    # def find_neighborhood(self):
    #     pass 
    

    # def rooted_spr(self):
    #     pass


