"""
Vignette to demonstrate the capabilities of the tribal package.
"""

from tribal.preprocess import preprocess
from tribal import probabilites,df, roots
from tribal import Tribal
import pandas as pd 



if __name__ == "__main__":


    from tribal import lineage_tree, lineage_tree_list
    print(lineage_tree)
    print(lineage_tree_list)
    print(len(lineage_tree_list))
    print(probabilites.shape)
  
    # clonotype_list = ["Clonotype_1036", "Clonotype_1050", "Clonotype_10884", 
    #                   "Clonotype_1457", "Clonotype_755", "Clonotype_322"]

            
    isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']



    # df = df[df['clonotype'].isin(clonotype_list)]
    clonotypes, df_filt = preprocess(df, roots,isotypes, cores=3, verbose=True )

    pd.to_pickle(clonotypes, "tribal/data/clonotypes.pkl")
    

    tr = Tribal(n_isotypes=len(isotypes), verbose=True, restarts=1, niter=15)
            
    #run in refinement mode
    shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="refinement", cores=6)

    # lineage_tree_list = best_scores['Clonotype_1036']
    # pd.to_pickle(lineage_tree_list, "tribal/data/lineage_tree_list.pkl")

    # pd.to_pickle(lineage_tree_list[0], "tribal/data/lineage_tree.pkl")


    # #run in scoring mode
    # shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="score", cores=6)


    # #given a user-specified isotype transition probability matrix 
    # shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes =clonotypes,
    #                                                             transmat= probabilites,
    #                                                             mode="refinement", cores=6)
