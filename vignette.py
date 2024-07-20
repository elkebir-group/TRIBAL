
from tribal.preprocess import preprocess
from tribal import df, roots, clonotypes
from tribal import Tribal


if __name__ == "__main__":



    # clonotype_list = ["Clonotype_1036", "Clonotype_1050", "Clonotype_10884"]
    # use_light_chain = True
    # verbose = True
    isotypes = ['IGHM', 'IGHA2' , 'IGHG2', 'IGHG1', 'IGHA1' ,'IGHG4', 'IGHG3' , 'IGHE']

    # df = df[df['clonotype'].isin(clonotype_list)]
    # clonotypes, df_filt = preprocess(df, roots,isotypes, cores=3, verbose=True )



    tr = Tribal(n_isotypes=len(isotypes), verbose=True, restarts=2, niter=2)
    shm_score, csr_likelihood, best_scores, transmat = tr.fit(clonotypes=clonotypes, mode="refinement", cores=3)

