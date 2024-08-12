"""Initialization of tribal module."""
import pkg_resources
import pandas as pd
import numpy as np
from .tribal import Tribal
from .preprocess import preprocess
from .lineage_tree import LineageTree, LineageTreeList
from .base_tree import BaseTree
from .clonotype import Clonotype

dat_file = pkg_resources.resource_filename('tribal', 'data/df.pkl')
roots_file = pkg_resources.resource_filename('tribal', 'data/roots.pkl')
clono_file = pkg_resources.resource_filename('tribal', 'data/clonotypes.pkl')
transmat_file = pkg_resources.resource_filename('tribal', 'data/transmat.txt')
lin_tree_file = pkg_resources.resource_filename('tribal', 'data/lineage_tree.pkl')
lin_tree_list_file = pkg_resources.resource_filename('tribal', 'data/lineage_tree_list.pkl')

df= pd.read_pickle(dat_file)
roots = pd.read_pickle(roots_file)
clonotypes = pd.read_pickle(clono_file)
probabilites = np.loadtxt(transmat_file)
lineage_tree = pd.read_pickle(lin_tree_file)
lineage_tree_list = pd.read_pickle(lin_tree_list_file)
