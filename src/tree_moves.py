import networkx as nx
class TreeMoves:
    def __init__(self, T, root):
        self.T = T
        self.root = root
    


    def check_nni_improvement(unrooted_tree, edge):
        return True 


    def all_nni(self):
        nni_trees = []
        unrooted_tree  = self.T.to_undirected()
     
        cand_edges = self.find_internal_edges(unrooted_tree)
        for c in cand_edges:
            if check_improvement(unrooted_tree, c):
                nni_trees.append(nni_move(unrooted_tree,c ))
            

    @staticmethod
    def find_internal_edges(unrooted_tree):
        internal_vertices = [n for n in unrooted_tree.nodes if unrooted_tree.degree[n] >= 2]
        internal_edges = []
        for v in internal_vertices:
            for u in unrooted_tree.neighbors(v):
                if unrooted_tree.degree[u] >= 2:
                    if v < u:
                        if (v,u) not in internal_edges:
                            internal_edges.append((v,u))
                    else:
                        if not (u,v) in internal_edges:
                            internal_edges.append((u,v))



        return internal_edges
    

    @staticmethod
    def nni_moves(unrooted_tree, edge):
        out_trees = [] 
        left, right = edge

        unrooted_tree.remove_edge(left, right)

        tree1 = unrooted_tree.copy()
        tree2 = unrooted_tree.copy()
        e = find_internal_edges(tree1)
        left_subtree_roots = list(unrooted_tree.neighbors(left))
        right_subtree_roots = list(unrooted_tree.neighbors(right))

        if len(left_subtree_roots) >= 1 and len(right_subtree_roots) >= 1:


            tree1.remove_edge(right,right_subtree_roots[0])
            tree1.remove_edge(left, left_subtree_roots[0])


            tree1.add_edge(left,right_subtree_roots[0])
            tree1.add_edge(right,left_subtree_roots[0])
            tree1.add_edge(left,right)
            out_trees.append(tree1)

        if len(left_subtree_roots) ==2  and len(right_subtree_roots)==2:
            tree2.remove_edge(right,right_subtree_roots[0])
            tree2.remove_edge(left, left_subtree_roots[1])

            tree2.add_edge(right, left_subtree_roots[1])
            tree2.add_edge(left, right_subtree_roots[0])
            tree2.add_edge(left,right)
            out_trees.append(tree2)
        
        return out_trees
        





    














    @staticmethod
    def root_tree(unrooted_tree, root):
        #do a bfs on the graph 
        rooted_tree = nx.bfs_tree(unrooted_tree, source=root)

        return rooted_tree



    def find_neighborhood(self):
        pass 
    

    def rooted_spr(self):
        pass

test_tree = nx.Graph()
test_tree.add_edges_from([(0,2), (1,2), (2,3), (3,7), (3,6), (6,4), (6,5)])

tm = TreeMoves(test_tree, 0)

tm.root_tree(test_tree,0)
tm.all_nni()
tm.nni_moves(test_tree, (2,3))
