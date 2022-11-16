import networkx as nx
import numpy as np
from itertools import combinations
class USPR:
    def __init__(self, tree, rng=None, min_radius=2, max_radius=2):

        self.T = tree
  
        self.nodes = list(self.T.nodes())
        if rng is None:
            self.rng = np.random.default_rng()
        else:
            self.rng = rng
        self.min_radius = min_radius
        self.max_radius = max_radius
        #compute the number of spr moves as a sanity check
        self.taxa = [n for n in self.nodes if self.T.degree[n]==1]
        self.ntaxa = len(self.taxa)

        self.num_spr = 2*(self.ntaxa -3)*(2*self.ntaxa -7)
        # print(f"Expected number of spr moves: {self.num_spr}")

        # self.generate_valid_cand()
        self.generate_cand()
        # print(f"Total number of 2 edge spr moves: {len(self.spr_cand)}")
        self.generate_cand_nni()
        self.shuffle(self.nni_cand)
        self.shuffle(self.spr_cand)

  
   
        num_nni_moves = len(self.nni_cand)
        num_2edge_spr = len(self.spr_cand)
        assert self.num_spr >= (num_nni_moves + num_2edge_spr)

        # print(f"Total number of NNI moves: {len(self.nni_cand)}")
        # print(f"Computed number of spr moves: {num_nni_moves + num_2edge_spr}")


    def shuffle(self, mylist):
        self.rng.shuffle(mylist)
    @staticmethod  
    def is_adjacent(e1, e2):
        common_vertices = set(e1).intersection(set(e2))
        return len(common_vertices) > 0
    
    def edge_radius(self, e1,e2):
        shortest_path = np.Inf
        for u in e1:
            for v in e2:
                path = nx.shortest_path(self.T, u,v)
                if len(path) -1 < shortest_path:
                    shortest_path = len(path)-1
    
        return shortest_path

                
                

    def generate_cand(self):
        edges = list(self.T.edges)
        self.spr_cand = []
        #get all ordered pairs of edges
        for i in range(len(edges)):
            for j in range( len(edges)):
                if i == j:
                    continue
                pruned = edges[i]
                regraft = edges[j]

                #check if the pruned edge is adjacent to the regraft edge
                if not self.is_adjacent(pruned,regraft):
                    #this excludes 6(n-2) edge pairs with SPR identity moves
                    e_rad = self.edge_radius(pruned, regraft)
                    if e_rad > 1 and e_rad <= self.max_radius and e_rad >= self.min_radius:

                        self.spr_cand.append((pruned, regraft))
               

       

    def generate_cand_nni(self):
        self.nni_cand = []
        internal_edges = []
        for n in self.T.nodes:
            if self.T.degree[n] ==3:
                for u in self.T.neighbors(n):
                    if self.T.degree[u] ==3:
                        if (u,n) not in internal_edges:
                            internal_edges.append((n,u))
        for left, right in internal_edges:
   
      
            left_subtrees =[l for l in self.T.neighbors(left) if not l == right]
            left_swap = left_subtrees[0]
            right_subtrees = [r for r in self.T.neighbors(right) if r != left]
        
            for r in right_subtrees:
                self.nni_cand.append((left, left_swap, right, r))

  



    def spr(self, pruned, regraft):
        T =self.T.copy()
        #prune the cut edge resulting in T with two components
        shortest_path = np.Inf
        for u in pruned:
            path = nx.shortest_path(T, u, regraft[0])
            if len(path) < shortest_path:
                shortest_path = len(path) 
                cut_node = u

        for v in pruned:
            if v != cut_node:
                reattach_node = v
               
        T.remove_edge(cut_node, reattach_node)

        #remove any newly created internal nodes with degree 2
        if T.degree[cut_node] == 2:
            neighbors = list(T.neighbors(cut_node))
            T.add_edge(neighbors[0], neighbors[1])
            T.remove_node(cut_node)

        new_node = len(self.nodes) + 1
        while True:
            if new_node in T:
                new_node +=1
            else:
                break


        #split regraft edge in half
        T.remove_edge(regraft[0], regraft[1])
        for u in regraft:
            T.add_edge(u, new_node)
        
        #regraft the pruned edge 
        T.add_edge(new_node, reattach_node)
        
        return T

        


    def nni(self, left, left_swap, right, right_swap):
        T = self.T.copy()
        left_neighbors = list(T.neighbors(left_swap))
        for l in left_neighbors:
            T.remove_edge(l, left_swap)
    
        right_neighbors = list(T.neighbors(right_swap))
        for r in right_neighbors:
            T.remove_edge(r, right_swap)
        for l in left_neighbors:

            T.add_edge(l, right_swap)
        for r in right_neighbors:
            T.add_edge(r, left_swap)
        return T





        


        
    def __iter__(self):
       ''' Returns the Iterator object '''
       return USPRIterator(self)

class USPRIterator:

    def __init__(self, uspr):
        self.uspr = uspr

    def __next__(self):

        if len(self.uspr.spr_cand) > 0 and len(self.uspr.nni_cand) > 0:
            if self.uspr.rng.random() > 0:
                return self.spr_next()
            else:
                return self.nni_next()


        elif len(self.uspr.spr_cand) > 0:
            return self.spr_next()
      

        # elif len(self.uspr.nni_cand) > 0:
        #     return self.nni_next()  


        else:
            raise StopIteration
        
    
    def spr_next(self):
        pruned, regraft = self.uspr.spr_cand.pop()
        spr_tree = self.uspr.spr(pruned,regraft)
        return spr_tree 

    def nni_next(self): 
        left, left_swap, right, right_swap = self.uspr.nni_cand.pop()
        nni_tree = self.uspr.nni(left, left_swap, right, right_swap)
        return nni_tree


# tree = nx.Graph()
# tree.add_edges_from([('a','c'), ('b','c'), ('c', 'd'), ('d', 'e'), 
#                     ('d', 'f'), ('f', 'g'), ('f','j'), ('g', 'h'), ('g', 'i') ])

# us = iter(USPR(tree) )
# count = 0
# while True:
   
#     try:
#         foo = next(us)
#     except:
#         print(count)
#         break
#     count += 1




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


