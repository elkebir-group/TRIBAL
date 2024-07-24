from .tribal import Tribal
from .preprocess import preprocess
from .lineage_tree import LineageTree, LineageTreeList
from .base_tree import BaseTree, Clonotype


import pkg_resources
import pandas as pd 
import numpy as np

dat_file = pkg_resources.resource_filename('tribal', 'data/data.csv')
roots_file = pkg_resources.resource_filename('tribal', 'data/roots.csv')
clono_file = pkg_resources.resource_filename('tribal', 'data/clonotypes.pkl')
transmat_file = pkg_resources.resource_filename('tribal', 'data/transmat.txt')

df= pd.read_csv(dat_file)
roots = pd.read_csv(roots_file)
clonotypes = pd.read_pickle(clono_file)
probabilites = np.loadtxt(transmat_file)


