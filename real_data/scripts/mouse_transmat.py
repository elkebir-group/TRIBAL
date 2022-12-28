
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
    
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alpha", type=float, default=0.75)
    parser.add_argument("-n", "--nisotypes", type=int, default=7)
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()


    trans = gen_trans_mat(args.alpha, args.nisotypes)

    if args.outfile is not None:
        np.savetxt(args.outfile, trans)



