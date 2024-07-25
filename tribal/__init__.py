from .tribal import Tribal
from .preprocess import preprocess
from .lineage_tree import LineageTree, LineageTreeList
from .base_tree import BaseTree, Clonotype
from .utils import write_fasta, save_dict



import pkg_resources
import pandas as pd 
import numpy as np

dat_file = pkg_resources.resource_filename('tribal', 'data/df.pkl')
roots_file = pkg_resources.resource_filename('tribal', 'data/roots.pkl')
clono_file = pkg_resources.resource_filename('tribal', 'data/clonotypes.pkl')
transmat_file = pkg_resources.resource_filename('tribal', 'data/transmat.txt')

df= pd.read_pickle(dat_file)
roots = pd.read_pickle(roots_file)
clonotypes = pd.read_pickle(clono_file)
probabilites = np.loadtxt(transmat_file)


