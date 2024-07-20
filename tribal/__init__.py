from .tribal import Tribal
from .preprocess import preprocess


import pkg_resources
import pandas as pd 

dat_file = pkg_resources.resource_filename('tribal', 'data/data.csv')
roots_file = pkg_resources.resource_filename('tribal', 'data/roots.csv')
clono_file = pkg_resources.resource_filename('tribal', 'data/clonotypes.pkl')

df= pd.read_csv(dat_file)
roots = pd.read_csv(roots_file)
clonotypes = pd.read_pickle(clono_file)


