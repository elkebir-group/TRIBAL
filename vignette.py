"""
Vignette to demonstrate the capabilities of the tribal package.
"""

from tribal.preprocess import preprocess
from tribal import probabilites,df, roots
from tribal import Tribal

from tribal import lineage_tree, lineage_tree_list


if __name__ == "__main__":


    print(lineage_tree)
    print(lineage_tree_list)
    print(len(lineage_tree_list))
    print(probabilites.shape)
  

            
    isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']

    clonotypes, df_filt = preprocess(df, roots,isotypes, cores=3, verbose=True )


    tr = Tribal(n_isotypes=len(isotypes), verbose=True, restarts=1, niter=15)
            
    #run in refinement mode
    shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="refinement", cores=6)

