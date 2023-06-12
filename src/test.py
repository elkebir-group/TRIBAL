
import numpy as np
from utils import pickle_load
tmat_fname = "/scratch/projects/tribal/experimental_data/day_14/tribal/0.1/transmat.txt"

reps = [i+1 for i in range(20)]

pt ="/scratch/projects/tribal/benchmark_pipeline/sim_data/tmat_inf/seq/transmats"


def check_triangle_inequality(distance_matrix):
    n = distance_matrix.shape[0]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if distance_matrix[i, k] > distance_matrix[i, j] + distance_matrix[j, k]:
                    return False
    return True

for i in reps:
    tmat_fname = f"{pt}/transmat{i}.txt"



    tmat = np.loadtxt(tmat_fname)

    tmat = -1*np.log(tmat)
    is_valid = check_triangle_inequality(tmat)
    print(f"{i}: {is_valid}")
    # if is_valid:
    #     print("The distance matrix satisfies the triangle inequality.")
    # else:
    #     print("The distance matrix does not satisfy the triangle inequality.")


res1 = pickle_load("/scratch/projects/tribal/test/refine_ilp.pickle")[0]

res2 = pickle_load("/scratch/projects/tribal/test/refine_ilp_new.pickle")[0]

print(res1)


print(res2)


# Example distance matrix

tmat_fname = "/scratch/projects/tribal/experimental_data/day_14/tribal/0.1/transmat.txt"
tmat = -1*np.log(np.loadtxt(tmat_fname))
# Check if the distance matrix satisfies the triangle inequality
is_valid = check_triangle_inequality(tmat)

if is_valid:
    print("The distance matrix satisfies the triangle inequality.")
else:
    print("The distance matrix does not satisfy the triangle inequality.")

def compute_iso_score(res):
    tree = res.tree
    iso = res.isotypes
    score = 0
    score2 =0
    nodes = tree.preorder_traversal()
    for n in nodes:
        if n == "3_0":
            score2 = 0
        t= iso[n]
        for c in tree.children(n):
            s = iso[c]
            score += tmat[t,s]
            score2 += tmat[t,s]
    print(f"Score: {score}  {score2}")
    return score, score2


s1, s1_2 = compute_iso_score(res1)
s2, s2_2 = compute_iso_score(res2)





