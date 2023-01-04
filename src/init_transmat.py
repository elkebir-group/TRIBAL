
import argparse
import numpy as np
def gen_trans_mat(alpha,n ):
    trans = np.zeros((n,n))
    for i in range(n):
        j_vals = np.arange(i+1, n)
        for j in range(n):
            if i > j:
                trans[i,j] = 0
            elif i ==j:
                trans[i,j] = alpha
            else:
                trans[i,j] = (1-alpha)/j_vals.shape[0]
    

    trans[n-1, n-1] = 1
    # print(trans.sum(axis=1))
    return trans
    

def add_noise(transmat, rng, mu=0.05, sigma=0.005, min_prob=0.01):

        white_noise =rng.normal(mu, sigma, transmat.shape)
        tmat = transmat + white_noise 
        for i in range(transmat.shape[0]):
            for j in range(transmat.shape[1]):
                if i > j:
                    tmat[i,j]=0
                elif tmat[i,j] < 0:
                    tmat[i,j] = min_prob
        
        tmat =tmat/tmat.sum(axis=1).reshape(-1,1)
    
        return tmat 



    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alpha", type=float, default=0.75)
    parser.add_argument("-n", "--nisotypes", type=int, default=7)
    parser.add_argument("-o", "--outfile")
    parser.add_argument("--add_noise", action="store_true")
    parser.add_argument("--mu", type=float, default=0.05, help="mean of gaussian white noise to add for distortion")
    parser.add_argument("--sigma", type=float, default=0.006, help="std of gaussian white noise to add for distortion")
    parser.add_argument("--min_prob", type=float, default=0.01, help="minimum probability of any transition")
    parser.add_argument("-s", "--seed", type=int, help="random number seed", default=1026)
    args = parser.parse_args()
    # args = parser.parse_args([
    #     "-a", "0.8",
    #     "-n", "7",
    #     "-o", "/scratch/projects/tribal/real_data/test/test_tmat.txt",
    #     "--add_noise",
    #     "--mu", "0.1",
    #     "--sigma", "0.05",
    #     "--min_prob", "0.05",
    #     "-s", "12"
    # ])


    trans = gen_trans_mat(args.alpha, args.nisotypes)


    if args.add_noise:
        rng = np.random.default_rng(args.seed)
        trans = add_noise(trans, rng, mu=args.mu, sigma=args.sigma, min_prob=args.min_prob)
        # for i in range(trans.shape[0]):
        #     print(trans[i,:])

    if args.outfile is not None:
        np.savetxt(args.outfile, trans)



