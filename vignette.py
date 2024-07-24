
from tribal.preprocess import preprocess
from tribal import  clonotypes, probabilites,df, roots
from tribal import Tribal




if __name__ == "__main__":



    clonotype_list = ["Clonotype_1036", "Clonotype_1050", "Clonotype_10884", 
                      "Clonotype_1457", "Clonotype_755", "Clonotype_322"]

            
    isotypes = ['IGHM', 'IGHG3', 'IGHG1', 'IGHA1','IGHG2','IGHG4','IGHE','IGHA2']



    # df = df[df['clonotype'].isin(clonotype_list)]
    # clonotypes, df_filt = preprocess(df, roots,isotypes, cores=3, verbose=True )
    # import pandas as pd 
    # pd.to_pickle(clonotypes, "tribal/data/clonotypes.pkl")
    
    print(len(clonotypes))
    tr = Tribal(n_isotypes=len(isotypes), verbose=True, restarts=1, niter=15)
            
    #run in refinement mode
    shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="refinement", cores=6)

    # #run in scoring mode
    # shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="score", cores=6)


    # #given a user-specified isotype transition probability matrix 
    # shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes =clonotypes,
    #                                                             transmat= probabilites,
    #                                                             mode="refinement", cores=6)
