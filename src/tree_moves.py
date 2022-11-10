import networkx as nx
import numpy as np
class TreeMoves:
    def __init__(self, T, root, seq_att_name= "sequence", isotypes):
        self.T = T
        self.root = root
        self.nodes = list(self.T.nodes())
        self.rng = np.random.default_rng(5)
        self.iso_att_name = iso_att_name
        self.seq_att_name = seq_att_name
        self.isotypes = isotypes
    


  
    def regrafting(self, spr_tree, cut_node, internal_node_label, grafting_node):
          parent_node = list(spr_tree.predecessors(grafting_node))[0]
          spr_tree.remove_edge(parent_node, grafting_node)
          
          internal_node_state = min(self.isotypes[cut_node], self.isotypes[grafting_node])
          spr_tree.add_node(internal_node_label, isotype=internal_node_state)
          spr_tree.add_edge(parent_node, internal_node_label)
          spr_tree.add_edge(internal_node_label, cut_node)
          spr_tree.add_edge(internal_node_label, grafting_node)

          

          return spr_tree

    def case1_spr(self, cut_node):
        #TODO: fix case1 to consider multifurcating trees
        #case 1 is that the parent of of cut node is not the root 
        # and it will be grafted to a preexisting edge
        if cut_node ==4:
            print(cut_node)
        tree = self.T.copy()
        one_spr_trees = []
        parent = list(tree.predecessors(cut_node))[0]
        children = list(self.T.successors(parent))
        other_child = [c for c in children if c != cut_node][0]
        grandparent = list(tree.predecessors(parent))[0]


        tree.remove_edge(parent,cut_node)
        tree.add_edge(grandparent,other_child)
        tree.remove_node(parent)

        cut_state = self.isotypes[cut_node]
        sucessors = list(nx.dfs_preorder_nodes(tree, source=cut_node))
        move_cand = [n for n in tree.nodes if (n not in sucessors)
                        and self.isotypes[n] <= cut_state and n != self.root and n != other_child]
        
        for m in move_cand:
                 
            spr_tree = tree.copy()
            if m ==3:
                print(cut_node)
            self.regrafting(spr_tree, cut_node, parent,m)
            # print(spr_tree.edges)
            one_spr_trees.append(spr_tree)
        
        return one_spr_trees

        

    def case2_spr(self,cut_node):
        #case two is when a subtree of the root is the cut edge
        tree = self.T.copy()
        one_spr_trees = []
        #parent is the root
        children = list(self.T.successors(self.root))
        other_child = [c for c in children if c != cut_node][0]
        grandchildren = tree.successors(other_child)
        tree.remove_node(other_child)
        tree.remove_edge(self.root, cut_node)
        
        for g in grandchildren:
            tree.add_edge(self.root, g)
        
        cut_state = self.isotypes[cut_node]

        # root_children = self.T.successors(self.root)
        sucessors = list(nx.dfs_preorder_nodes(tree, source=cut_node))
        move_cand = [n for n in tree.nodes if (n not in sucessors)
                        and self.isotypes[n] <= cut_state and n != self.root]
        
        for m in move_cand:
                 
            spr_tree = tree.copy()
            self.regrafting(spr_tree, cut_node,other_child,m)
            # print(spr_tree.edges)
            one_spr_trees.append(spr_tree)
        
        return one_spr_trees
        

        
    def case3_spr(self, cut_node):
        #the cut node that is not a subtree of the root is regrafted
        # to a new root
        tree = self.T.copy()

        parent = list(tree.predecessors(cut_node))[0]
        children = list(self.T.successors(parent))
        other_child = [c for c in children if c != cut_node][0]
        grandparent = list(tree.predecessors(parent))[0]


        tree.remove_edge(parent,cut_node)
        tree.add_edge(grandparent,other_child)
        tree.remove_node(parent)
        root_successors = list(tree.successors(self.root))
        for s in root_successors:
            tree.remove_edge(self.root, s)
            tree.add_edge(parent, s)
        
        tree.add_edge(self.root, cut_node)
        tree.add_edge(self.root, parent)
        # print(list(tree.edges()))

        return [tree]
        


    def spr(self):
        cand_trees = []
        cand_nodes = [n for n in self.nodes if n != self.root ]

        for n in list(self.T.successors(self.root)):
            case2_trees = self.case2_spr(n)
            cand_trees += case2_trees

        for n in cand_nodes:

            if n not in list(self.T.successors(self.root)):
                
           
                case1_trees = self.case1_spr(n)
                case3_trees = self.case3_spr(n)
                cand_trees += (case1_trees + case3_trees)

    
        to_delete = []
        for i, c1 in enumerate(cand_trees):
            for j, c2 in enumerate(cand_trees):
                if i < j:
                    if set(list(c1.edges)) == set(list(c2.edges)):
                        to_delete.append(j)
                      
        cand_tree = [cand_trees[i] for i in range(len(cand_trees)) if i not in to_delete]       
        return cand_tree



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


