
import numpy as np

class BCellSim:
    def __init__(self, seed,n, ml, mh) -> None:
        
        self.rng = np.random.default_rng(seed)
        self.n = n
        self.alphabet = ["A", "C", "G", "T"]
        self.ml = ml
        self.mh = mh
        self.generate_tree()

        self.isotypes = np.arange(6)
        self.seq_labels = {}
        self.seq_labels[0] = {'L': self.generate_seq(ml), 'H' : self.generate_seq(mh)}
   

    
    def generate_tree(self):
       sequence = rng.random_choice(self.n, self.n-2)
       self.tree = nx.from_prufer_sequence(sequence)
       self.
    

    def generate_seq(self, n):
        return self.rng.random_choice(self.alphabet, n)


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
        
        self.tree = self.root_tree(tree, self.root)
   

    
    def sample_sites(self, m, prop):
        npos = int(prop*m)
        return self.rng.choice(ml,npos, replace=False)

    @staticmethod 
    def class_switch(iso):
        new_iso = self.rng.random_choice([i for i in range(iso, 6)],1)
        return new_iso, seq
    
    def shm(self,iso,  seq, positions, cycles=1):
        mod_pos = self.rng.random_choice(positions, cycles, replace=F)
        new_seq = parent_seq.copy()
        for i in mod_pos:
            parent_base = seq[i]

            choice_base = [a for a in self.alphabet if a != parent_base]
            base = self.rng.random_choice(choice_base)
       
            new_seq[pos] = choice_base
        return iso, new_seq

    
    def both(self,iso, seq, positions, cycle):
        new_iso, seq = self.class_switch(iso, seq, positions)
        new_iso, new_seq = self.shm(new_iso, seq, positions)

        return new_iso, new_seq


    def evolve(self, root_seq, prop):

        light_pos self.sample_sites(len(root_seq), prop)
        events = ['class_switch', 'shm', 'both']
        curr_node = 0
        while True:
            children = self.tree.successors(curr_node)
            for c in children:
                event = self.rng.random_choice(events, 1)
                if event == "class_switch":
                    self.isotype[c] = self.class_switch(self.isotype[curr_node])
                elif event == "shm":
                    ncycles = self.rng.random_choice(4, 1)
                    self.seq_labels






    