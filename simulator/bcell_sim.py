
import numpy as np
import networkx as nx
from scipy.linalg import expm
import argparse
import pickle 
import copy 
import itertools

class BCellSim:
    def __init__(self,n, m_light=320, m_heavy=355, root_lc=None, root_hc =None, 
                 pos_lc=None, pos_hc=None, prop_hot = 0.3,
                 n_isotypes=7, cold_rate=0.002, isotype_rate_matrix = None, 
                 hot_rate=1, expected_branch_length =5/(320+355), seed=1026) -> None:
        

        #branch lengths are the expected number of subsitutions per site
        expected_branch_length = 2.5/(m_light+m_heavy)
        self.changes = 0
        self.rng = np.random.default_rng(seed)
        self.n = n

        self.root = 0
        self.alphabet = ["A", "C", "G", "T"]
        
        self.alpha_to_index = {a : i for i,a in enumerate(self.alphabet)}
        

        self.labels = {}
        self.isotypes = np.arange(n_isotypes)
        if isotype_rate_matrix is None:
            self.Q_isotype =np.zeros((n_isotypes, n_isotypes))
            
            for i in range(n_isotypes):
                self.Q_isotype[i,:] = (self.isotypes >= i)/ np.sum((self.isotypes >= i))
          
        else:
            self.Q_isotype = isotype_rate_matrix

 
        self.char = ['L', 'H', 'I']
        self.chains = ['L', 'H']
        self.labels = {i : {} for i in self.char}
        self.m = {'L' : m_light, 'H': m_heavy}
        self.labels['I'][0]=0
        self.position_type = {}
        self.tree = None
        self.branch_lengths = None
        
        for root, chain in zip([root_lc, root_hc],self.chains):
            if root is None:
        
                self.labels[chain][0] =self.generate_seq(self.m[chain])
           
            else:
                self.labels[chain][0] =  root.split("")
        
               
        for pos_type, chain in zip([pos_lc, pos_hc],self.chains):
            
            if pos_type is None:
                self.position_type[chain]= self.rng.choice(['c', 'h'], size=self.m[chain],p=[1-prop_hot, prop_hot])
           
            else:
                self.position_type[chain] =  pos_type.tolower().split("")
        

        self.Q_cold = self.construct_rate_matrix(cold_rate)
        self.Q_hot = self.construct_rate_matrix(hot_rate)
    
        self.exp_br_len =expected_branch_length
        



    @staticmethod
    def construct_rate_matrix(rate):
        Q = rate * np.ones((4, 4))
     
   
        for i in range(Q.shape[0]):
            Q[i,i] = -3*rate
        return Q
    

    def generate_seq(self, n):
        return self.rng.choice(self.alphabet, n).tolist()

    def generate(self):
        self.tree = self.random_tree()

        self.branch_lengths =self.draw_branch_lengths(self.exp_br_len)

        self.evolve()
       
        self.pairwise_distance()
        print(self.pw_dist)



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
        n = self.n + 1
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

    def simulate_isotype(self, parent_iso):
        trans_prob = self.Q_isotype[parent_iso,:]
        new_iso = self.rng.choice(self.isotypes, 1, p=trans_prob)
        return new_iso[0]
    
    def draw_branch_lengths(self, scale=0.001):
        branch_lengths = {}
        for e in self.tree.edges:
            branch_lengths[e] = self.rng.exponential(scale)
        return branch_lengths

    def simulate_base(self, Q,br_length, anc_base):
        trans_matrix = expm(Q*br_length)   
        probs = trans_matrix[self.alpha_to_index[anc_base],]
        
        new_base = self.rng.choice(self.alphabet, size=1, p=probs )  
        if anc_base != new_base:
            self.changes += 1
        return new_base[0]

    def evolve(self):
        nodes = list(nx.dfs_preorder_nodes(self.tree,0))
        nodes = [n for n in nodes if n != self.root]
        
        for n in nodes:
            parent = list(self.tree.predecessors(n))[0] 
            iso = self.labels['I'][parent] 
            self.labels['I'][n] = self.simulate_isotype(iso)
            br_length =  self.branch_lengths[(parent,n)]
         
            for chain in self.chains:
              
                seq = copy.copy(self.labels[chain][parent])
                  
                for i in range(len(seq)):
    
                    if self.position_type[chain][i]=='h':
                        seq[i] = self.simulate_base(self.Q_hot, br_length, seq[i])  
                    else:
                        seq[i] = self.simulate_base(self.Q_cold, br_length, seq[i])
            
            
                self.labels[chain][n] = seq


    def save_labels(self, fname):
        with open(fname, 'w+') as file:
            file.write("node,isotype,light_chain,heavy_chain")
            for n in self.tree.nodes:
                lc = "".join(self.labels['L'][n])
                hc = "".join(self.labels['H'][n])
                file.write(f"\n{n},{self.labels['I'][n]},{lc},{hc}")

    def pairwise_distance(self):
        leaf_nodes = [n for n in self.tree.nodes if self.tree.out_degree(n)==0]
        self.pw_dist = 0
        self.dist= {}
        for n1, n2 in itertools.combinations(leaf_nodes, 2):
            chain_total = 0
            for c in self.chains:
                seq1 = np.array(self.labels[c][n1])
                seq2 = np.array(self.labels[c][n2])
                val =(seq1 != seq2).sum()

                chain_total += val
            self.dist[n1,n2] = chain_total 
            self.pw_dist += chain_total
        print(self.pw_dist)




    
    def save_tree(self, path):

        leafs = [n for n in self.tree.nodes if len(list(self.tree.successors(n))) ==0]
                          
        with open(path, "w+") as file:
            file.write(f"{len(list(self.tree.edges))} #edges\n")
            for u,v in list(self.tree.edges):
                file.write(f"{u} {v}\n")
            file.write(f"{len(leafs)} #leaves\n")
            for l in leafs:
                file.write(f"{l}\n")
            
    
    def alignment_save(self, fname):
        with open(fname, "w+") as file:
            file.write("id,sequence,isotype")
            for n in self.tree.nodes:
                if len(list(self.tree.successors(n))) ==0 or n ==self.root:
                    lc = self.labels["L"][n]
                    lc = "".join(lc)
                    hc = self.labels["H"][n]
                    hc = "".join(hc)
                    iso = self.labels["I"][n]
                    file.write(f"\n{n},{lc+hc},{iso}")


    def pickle_save(self, fname):
        with open(fname, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__=="__main__":

    parser= argparse.ArgumentParser()

    parser.add_argument("-n", "--taxa", type=int, default=7,
                        help="number of taxa to simulate")
    parser.add_argument("-l", "--light", type=int, default=350,
                        help="number of sites in light chain to simulate")
    parser.add_argument( "--heavy", type=int, default=320,
                        help="number of sites in light chain to simulate")
    parser.add_argument("-r", "--root", type=str, 
                        help="filename containing the root sequence for the light and heavy chain")
    parser.add_argument( "--config", type=str,
                        help="name of yaml config file")
    parser.add_argument("-p", "--prop_hot", type=float, default=0.3,
                        help="proportion of sites that should be hot sites")
    parser.add_argument("-i", "--nisotypes", type=int,
                        help="the number of isotypes to include")
    parser.add_argument("-c", "--cold_rate", type=float, default=0.03,
                        help="mutation rate for coldspots ")
    parser.add_argument( "--hot_rate", type=float, default=1.0,
                        help="mutation rate for hotspots")
    parser.add_argument("--iso_trans", type=str,
                        help="file containing the desired isotype transition rate")
    parser.add_argument("-s", "--seed", type=int, default=1026,
                        help="random number seed")
    parser.add_argument("-b", "--br_length", type=float, default=0.01,
                        help="expected number of substitions per site")   

    parser.add_argument("-o", "--output", type=str,
                        help="name of output file to pickle the simulator")  
    parser.add_argument( "--labels", type=str, 
                        help="name of output file to write the labels of the nodes")
    parser.add_argument( "--tree", type=str, 
                        help="name of output file for tree")
    parser.add_argument( "--alignment", type=str, 
                        help="input alignment for methods")
    
    
    pth = "simulator"
    args= parser.parse_args([
        "-l", "328",
        "--heavy", "376",
        "--tree", f"{pth}/tree.txt",
        "--labels", f"{pth}/labels.csv",
        "--output", f"{pth}/sim.pickle",
        "--alignment",f"{pth}/alignment.csv"
    
    ])

    sm = BCellSim(n = args.taxa,
        m_light = args.light,
        m_heavy =  args.heavy,   
        prop_hot = args.prop_hot,
        cold_rate = args.cold_rate,
        hot_rate = args.hot_rate,
        seed = args.seed,
    )

    sm.generate()
    if args.tree is not None:
        sm.save_tree(args.tree)
    if args.labels is not None:
        sm.save_labels(args.labels)
    if args.output is not None:
        sm.pickle_save(args.output)
    if args.alignment is not None:
        sm.alignment_save(args.alignment)




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
    

