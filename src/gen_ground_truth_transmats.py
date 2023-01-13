import init_transmat as tm
import numpy as np
import argparse
from draw_state_diagram import DrawStateDiag    




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--mu", type=float, default=0.07, help="mean of gaussian white noise to add for distortion")
    parser.add_argument("--sigma", type=float, default=0.05, help="std of gaussian white noise to add for distortion")
    parser.add_argument("-s", "--seed" , type=int, default=1026)
    parser.add_argument("-n", "--nreps", type=int,  default=1, help="number of reps")
    parser.add_argument("-p", "--path", type=str, required=True  )
    parser.add_argument("--n_isotypes", type=int, default=7, help="the number of isotypes states to use if isotype encoding file is not providied")
    parser.add_argument("--min_jump", type=float, default=0.1)
    parser.add_argument("--max_jump", type=float, default=0.25)
    parser.add_argument("--noise", action="store_true")
    parser.add_argument("--type", choices=["seq", "direct"], default="direct")
    parser.add_argument("--epsilon", type=float, default=0.01)



    args= parser.parse_args()

    

    possible_jump_probs = np.arange(args.min_jump, args.max_jump, step=0.05)



    rng = np.random.default_rng(args.seed)
    jps = rng.choice(possible_jump_probs, size=args.nreps)
    for i in range(args.nreps):
        if args.type == "direct":
            print(f"direct {i} jps: {jps[i]}")
            tmat = tm.gen_trans_mat(1-jps[i], args.n_isotypes)
            if args.noise:
                tmat = tm.add_noise(tmat, rng, args.mu,args.sigma)
        else:
            print(f"seq {i} jps: {jps[i]}")
            tmat = tm.gen_seq_mat(1-jps[i], args.n_isotypes, args.epsilon)
        np.savetxt(f"{args.path}/transmat{i+1}.txt", tmat)
        DrawStateDiag(tmat).heatmap(f"{args.path}/transmat{i+1}.png")
    print("DONE!")


