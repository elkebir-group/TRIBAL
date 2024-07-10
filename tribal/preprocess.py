import pandas as pd
import utils as ut  
from io import StringIO
import subprocess
from Bio import AlignIO
import tempfile
import os, re
from lineage_tree import ParimonyForest
import networkx as nx
from ete3 import Tree
import multiprocessing as mp

CONFIG = """
J
1
10
O
1
4
5
.
Y
"""

class Preprocessor:
    def __init__(self, isotypes,  min_size=4, use_light_chain=True, verbose=False, ) -> None:
        """
        A class used to preprocess data for TRIBAL.

        Attributes
        ----------
            isotypes: list
                Dictionary that maps the isotype labels, i.e., IghM, to a linear encoding.
            min_size: int
                The mininum number of B cells needed to form a clonotype (default 4).
            use_light_chain: bool
                 Should the light chain be included in the BCR sequences (default True)
            verbose: bool  
                Print more verbose output (default False)
        
        Notes
        -----
        The isotype_encoding is an ordered list of the isotypes. Since there is no standard
        snytax for isotype labels, e.g., IgM, IghM, M, IGM), this list must be provided. must
        include all isotype labels expected in the input data and be ordered to reflect the linear
        ordering of the heavy chain constant genes in the genome.
      


        
        Examples
        --------
        Here is an example of how to use the Preprocessor class::

            isotypes = ['IgM', 'IgG', 'IgE', 'IgA']
            preprocessor = Preprocessor(isotypes, min_size=5, use_light_chain=True)
            clonotypes, df_filtered = preprocessor.process(df, roots, cores=5)
            print(len(clonotypes))
        """

        self.verbose = verbose
        self.isotype = {iso: i for i, iso in enumerate(isotypes)}
        self.min_size = min_size
        self.use_light_chain = True


    def preprocess(self, df: pd.DataFrame, roots: pd.DataFrame, cores:int =1, dnapars_path=None):
        """
        Preprocess the input data to prepare data for TRIBAL. The clonotypes will first be filtered to ensure each 
        clontoype has at least min_size cells. Each retained clonotype is aligned to the root sequence and then maximum 
        parsimony forest is enumerated for the B cells with that clonotype.

        Parameters
        ----------

        df: pd.Dataframe
            A datframe with columns including 'Clonotype', 'Heavy Chain Variable Seq', 'Light Chain Variable Seq',
            'Heavy Chain V Allele', 'Light Chain V Allele', and 'Heavy Chain Isotype'.  The 'Light Chain' columns are optional.
            
        roots: pd.Dataframe
            A datframe with columns including 'Clonotype', 'Heavy Chain Root', 'Light Chain Root'
             The 'Light Chain Root' column is only required if use_light_chains is True.
        
        cores: int
            The number of cores to use (default 1)
        
        dnapars_path: str
            The path to the dnapars executable if dnapars is not on system path.

        
        Returns
        -------
        dict
            a dictionary containing the input data for TRIBAL

        """

        self.check_df(df, roots)

        #first filter out clonotypes smaller than min size
        if self.verbose:
            print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
        
        df = df.groupby("Clonotype").filter(lambda x: len(x) >= self.min_size)

        if self.verbose:
            print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
    
        df = self.filter_alleles(df, "Heavy Chain V Allele")
        if self.use_light_chain:
            df = self.filter_alleles(df, "Light Chain V Allele")
        df = df[ df['Heavy Chain Isotype'].notna()]

        df = df.groupby("Clonotype").filter(lambda x: len(x) >= self.min_size)

        if self.verbose:
            print(f"The number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
    
        if self.verbose:
            print(f"After filtering, the number of cells is {df.shape[0]} and the number of clonotypes is {df['Clonotype'].nunique()}.")
        
        #prep dnapars sequence ids
        df['seq'] = df.groupby('Clonotype').cumcount() + 1
        df['seq'] = 'seq' + df['seq'].astype(str)

        df.columns.values[0] = "cellid"
        roots.columns.values[0] = "Clonotype"
        df["isotype"] = df['Heavy Chain Isotype'].map(self.isotypes)


    
        clonodict = {}
        tree_size_dict = {}

        
        instances  = [ (j, df.copy(), roots.copy(), dnapars_path) for j in df["Clonotype"].unique()]
        with mp.Pool(cores) as pool:
    
            results = pool.starmap(self.run, instances)
        
        for j, linforest in results:
            clonodict[j] = linforest
            tree_size_dict[j] = linforest.size()
        
            
        df["ntrees"] = df["Clonotype"].map(tree_size_dict)
        return clonodict, df
    
    
    @staticmethod
    def check_df(df, roots):
 
        return True
    
    
    def run(self, j, df, roots, dnapars_path):
        clono = df[df["Clonotype"]== j]
        isotypes = dict(zip(clono["seq"], clono["isotype"]))
        isotypes["naive"] = 0
        heavy_seqs = dict(zip(clono["seq"],clono["Heavy Chain Variable Seq"]))

        heavy_root = roots[roots["Clonotype"]==j]["Heavy Chain Root"].values[0]

        heavy_chain_align = self.align(heavy_seqs, heavy_root)
        if self.use_light_chain:
            light_seqs = dict(zip(clono["seq"],clono["Light Chain Variable Seq"]))

            light_root = roots[roots["Clonotype"]==j]["Light Chain Root"].values[0]



            light_chain_align = self.align(light_seqs, light_root)
        alignment = {}
        for key in heavy_chain_align:
            if self.use_light_chain:
                alignment[key] = heavy_chain_align[key].upper() + light_chain_align[key].upper()
            else:
                alignment[key] = heavy_chain_align[key].upper()

        

        mapping = dict(zip(clono["seq"].values,clono['cellid'].values ))
        if self.verbose:
            print(f"Running dnapars for clonotype {j} with {len(alignment)} sequences...")
        align_str = self.convert_to_string(alignment)
        phylip_str = self.convert_alignment_to_phylip(alignment_str=align_str, input_format="fasta")

        outtrees = self.run_dnapars(phylip_str, dnapars_path)
        tree_list = self.create_trees(outtrees)


    
        linforest = ParimonyForest(alignment=alignment, isotypes=isotypes, mapping=mapping)
        if self.verbose:
            print(f"Clonotype {j} has {len(tree_list)} max parsimony trees.")
        linforest.generate_from_list(tree_list, root="naive")
        return j,linforest

    @staticmethod
    def run_dnapars(phylip_content, dnapars_path=None):
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Write PHYLIP content to a temporary "infile"
            infile_path = os.path.join(temp_dir, "infile")
            with open(infile_path, "w") as infile:
                infile.write(phylip_content)
            
            # Write config content to a temporary "config" file
            # config_path = os.path.join(temp_dir, "dnapars.cfg")
            # with open(config_path, "w") as config:
            #     config.write(config_content)
            # outfile = os.path.join(temp_dir, "outfile")
            # Run dnapars from the temporary directory
            if dnapars_path is None:
                dnapars_cline = ["dnapars"]
            else:
                dnapars_cline = [dnapars_path]
            # process = subprocess.Popen(dnapars_cline, cwd=temp_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            process = subprocess.Popen(dnapars_cline, cwd=temp_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            stdout = process.communicate(input=CONFIG)
            
            # Read the contents of outtree file
            outtree_path = os.path.join(temp_dir, "outtree")
            with open(outtree_path, "r") as outtree:
                outtree_content = outtree.read()
        
        return outtree_content

    @staticmethod
    def convert_to_string(seqs, root=None):
        if root is not None:
            mystr = f">naive\n{root}\n"
        else:
            mystr = ""
        for key, val in seqs.items():
            mystr += f">{key}\n{val}\n"
        return mystr

    @staticmethod
    def align(self, seqs, root):

        #!cat seqs.fasta mafft --quiet - > aligned.fasta
        seq_str = self.convert_to_string(seqs,root)
        child = subprocess.Popen(["mafft", "--quiet", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        child.stdin.write(seq_str.encode())
        child_out = child.communicate()[0].decode("utf8")

        child.stdin.close()
        aligned_seqs = ut.read_fasta(child_out)
    
        return aligned_seqs
    
    @staticmethod
    def convert_alignment_to_phylip(alignment_str, input_format="fasta"):
        # Create a file-like object from the string
        input_handle = StringIO(alignment_str)
        
        # Read the alignment
        alignment = AlignIO.read(input_handle, input_format)
        
        # Write the alignment in PHYLIP format to a string buffer
        output_handle = StringIO()
        AlignIO.write(alignment, output_handle, "phylip")
        
        # Get the content of the PHYLIP formatted alignment
        phylip_str = output_handle.getvalue()
        
        return phylip_str



    @staticmethod
    def align_clonotypes(self, df, roots):
        clonotype_alignments= {}
        for j in df["Clonotype"].unique():
            clono = df[df["Clonotype"]== j]
            heavy_seqs = dict(zip(clono["seq"],clono["Heavy Chain Variable Seq"]))
        
            heavy_root = roots[roots["Clonotype"]==j]["Heavy Chain Root"].values[0]
    
            heavy_chain_align = self.align(heavy_seqs, heavy_root)
            light_seqs = dict(zip(clono["seq"],clono["Light Chain Variable Seq"]))
        
            light_root = roots[roots["Clonotype"]==j]["Light Chain Root"].values[0]
            light_chain_align = self.align(light_seqs, light_root)
            concat_seqs = {}
            for key in heavy_chain_align:
                concat_seqs[key] = heavy_chain_align[key].upper() + light_chain_align[key].upper()
            clonotype_alignments[j] = concat_seqs
        return clonotype_alignments


    @staticmethod
    def convert_to_nx(self, ete_tree, root):
        nx_tree = nx.DiGraph()
        internal_node = 1
        internal_node_count = 0
        for node in ete_tree.traverse("preorder"):

            if node.name == "":
                node.name = internal_node
                internal_node_count += 1
                internal_node += 1
            if node.is_root():
                root_name =node.name


            for c in node.children:
                if c.name == "":
                    c.name = str(internal_node)
                    internal_node += 1
                    internal_node_count += 1
        

                nx_tree.add_edge(node.name, c.name)
    
        if len(list(nx_tree.neighbors(root))) == 0:
    
            nx_tree.remove_edge(root_name, root)
            nx_tree.add_edge(root, root_name)
        

        return nx_tree
    
    @staticmethod
    def create_trees(self, outtrees):
        cand_trees = []
        exp = '\[.*\]'
        file = StringIO(outtrees)
        # with open(cand_fname, 'r') as file:
        nw_strings = []
        nw_string = ""
        for nw in file:
                line = nw.strip()
                nw_string += line
                if ";" in line:
                    
                    nw_strings.append(nw_string)
                    nw_string = ""

        for nw in nw_strings:

            nw = re.sub(exp, '', nw)
            

            ete_tree = Tree(nw, format=0)

            nx_tree= self.convert_to_nx(ete_tree, "naive")
            
            cand_trees.append(nx_tree)
        return cand_trees


        


    @staticmethod
    def filter_alleles(self, df, col):
        grouped = df.groupby(['Clonotype', col]).size().reset_index(name='count')

        # # Identify the allele with the highest count for each Clonotype
        idx = grouped.groupby('Clonotype')['count'].idxmax()
        max_alleles = grouped.loc[idx]

        # # Merge with the original dataframe to filter the desired rows
        filtered_df = pd.merge(df, max_alleles[['Clonotype', col]], on=['Clonotype', col])

        return filtered_df


  






  





