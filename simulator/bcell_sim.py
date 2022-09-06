
import numpy as np
import networkx as nx
class BCellSim:
    def __init__(self,n, m_light, m_heavy, n_isotypes=7, beta=1/3, prop=1/3, seed=1026) -> None:
        
        self.rng = np.random.default_rng(seed)
        self.n = n+1
        self.root = 0
        self.alphabet = ["A", "C", "G", "T"]
        self.alpha_to_index = {a : i for i,a in enumerate(self.alphabet)}
        self.ml = m_light
        self.mh = m_heavy
        self.m = self.ml+self.mh
        self.hotspots = self.sample_sites(self.m, prop)

        #construct rate matrices
        self.R_cold = beta*np.ones((4,4),dtype=float)
        for i in range(self.R_cold.shape[0]):
            self.R_cold[i,i] = -3
        self.R_hot = self.R_cold 
        
        #initialize the root
        self.labels = {}
        self.labels[0] =self.generate_seq(m_light) + self.generate_seq(m_heavy)

        self.isotypes = np.arange(n_isotypes)
        self.isotype_labels = {}
        self.isotype_labels[0] = 0




    def generate_seq(self, n):
        return self.rng.choice(self.alphabet, n).tolist()

    def generate(self):
        self.tree = self.random_tree()

        self.draw_branch_lengths(5)

        self.evolve()
        # nx.set_node_attributes(self.tree, self.seq_labels)

        return self.tree, self.labels, self.isotype_labels

    @staticmethod
    def root_tree(tree, root):
        tree = nx.dfs_tree(tree,root)
     
        root_children = list(tree.successors(root))
        for c in root_children:
            grandchildren = tree.successors(c)
            for g in grandchildren:
                tree.add_edge(root, g)

        tree.remove_node(c)

        return tree
    
    


    def random_tree(self):
        tree = nx.Graph()
        ids = [i for i in range(self.n+1)]
        n = len(ids)
        center = 2*n -3 
        for i in range(n):
            tree.add_edge(i, center)
      
        next_node = n
       
        for i in range(n-3):
            pair = self.rng.choice(ids,2, replace=False)

            for p in pair:
                tree.add_edge(p, next_node)
                tree.remove_edge(p, center)
                tree.add_edge(p, next_node)
                tree.add_edge(center, next_node)
                ids.remove(p)
            ids.append(next_node)
            next_node += 1
        
        tree = self.root_tree(tree, 0)

        return tree

   

    
    def sample_sites(self, m, prop):
        npos = int(prop*m)
        return self.rng.choice(m,npos, replace=False)


    def class_switch(self, iso):
        
        return  self.rng.choice(self.isotypes[self.isotypes >= iso],1)[0]


    

    def draw_branch_lengths(self, scale=1):
        self.branch_lengths = {}
        for e in self.tree.edges:
            self.branch_lengths[e] = self.rng.exponential(scale)


    
    def simulate_base(self, Q,br_length, anc_base):
        trans_matrix = np.exp(Q*br_length)   
        probs = trans_matrix[self.alpha_to_index[anc_base],]
        probs = probs/probs.sum()
        new_base = self.rng.choice(self.alphabet, size=1, p=probs )  
        return new_base[0]

    def evolve(self):
        nodes = list(nx.dfs_preorder_nodes(self.tree,0))
        nodes = [n for n in nodes if n != self.root]
        
        for n in nodes:
            parent = list(self.tree.predecessors(n))[0] 
            seq = self.labels[parent].copy()  
            iso = self.isotype_labels[parent]    
            #TODO: vectorize this for hot and cold spots
            for i in range(self.m):
     
                br_length =  self.branch_lengths[(parent,n)]
                if i in self.hotspots:
                    seq[i] = self.simulate_base(self.R_hot, br_length, seq[i])  
                else:
                    seq[i] = self.simulate_base(self.R_cold, br_length, seq[i])
            
            self.isotype_labels[n] = self.class_switch(iso)
            self.labels[n] = seq


sm = BCellSim(7, 300, 250 )
tree, seqs, isos = sm.generate() 
print("done")              
    # def evolve(self, prop):
    #     lc_root = self.labels[0]["L"]
    #     hc_root = self.labels[0]['H']
    #     light_pos  = self.sample_sites(len(lc_root)), prop)
    #     heavy_pos = self.sample_sites(len(hc_root), prop)
    #     events = ['class_switch', 'shm', 'both']
   
    #     node_list = [0]
    #     while len(node_list) > 0:
    #         curr_node = node_list.pop()
    #         cur_iso = self.labels[curr_node]["I"]

    #         children = self.tree.successors(curr_node)

    #         for c in children:
    #             #events occurring on incoming edges to the node
    #             parent_lc = self.labels[curr_node]["L"].copy()
    #             parent_hc = self.labels[curr_node]['H'].copy()
    #             parent_iso = self.labels[0]["I"]

    #             event = self.rng.random_choice(events, 1)
             
    #             if event == "class_switch":
    #                 new_iso = self.class_switch(cur_iso)
    #                 self.labels[c] = {"L" : parent_lc, "H" : parent_hc, "I":new_iso}
                 
                  
    #             elif event == "shm":
    #                 ncycles = self.rng.random_choice(4, 1)
    #                 self.self.shm(child_lc, child_hc, light_pos, heavy_pos, ncycles)
    #                 self.labels["I"]
    #             else:
    #                 ncycles = self.rng.random_choice(4, 1)
    #                 new_iso, new_seq = self.both(self.isotype[curr_node], parent_seq, light_pos, heavy_pos, ncycles)
    #                 self.isotype[c]  = new_iso
    #                 self.seq_labels[c] = new_seq
    #             node_list.append(c)






    # def shm(self, seq, light_pos, heavy_pos, cycles=1):
    #     #mutate either heavy, light chain

    #     new_seq = parent_seq.copy()

    #     for i in range(cycles):

    #         if self.rng.random() < 0.5:
    #             positions = light_pos
    #             chain = "L"
    #         else:
    #             poisitions = heavy_pos
    #             chain = "H"
    #             parent_base = seq[i]
            
    #         mod_pos = self.rng.random_choice(positions, 1)

    #         choice_base = [a for a in self.alphabet if a != parent_base]
    #         base = self.rng.random_choice(choice_base)
       
    #         new_seq[pos][chain] = choice_base
    #     return new_seq

    
    # def both(self,iso, seq, light_pos, heav_pos, cycle):
    #     new_iso, = self.class_switch(iso)
    #     new_seq = self.shm(seq, light_pos, heav_pos, cycle)

    #     return new_iso, new_seq

    
    # def generate_tree(self):
    #    sequence = rng.random_choice(self.n, self.n-2)
    #    self.tree = nx.from_prufer_sequence(sequence)
    #    self.
    

