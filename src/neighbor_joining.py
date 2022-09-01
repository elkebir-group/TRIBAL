

import numpy as np
import networkx as nx


def compute_q_matrix(dist_mat):
    q_matrix = np.zeros_like(dist_mat)
    n = dist_mat.shape[0]
    for i in range(n):
        for j in range(n):
            if i ==j:
                q_matrix[i,j] = np.Inf
            else:
                q_matrix[i,j] = (n-2)*dist_mat[i,j] - np.sum(dist_mat[i,:].sum()) - np.sum(dist_mat[:,j])
    return q_matrix


def update_distance_matrix(dmat, a,b, ids, next_node):

    new_row = np.zeros(len(ids)-2)
    k = 0
    for i, node in enumerate(ids):
        if node not in [ids[a], ids[b]]:

 
            new_row[k] = 0.5*(dmat[a,i] + dmat[b,i] - dmat[a,b])
            k += 1
    
    dmat2 = np.delete(dmat, [a,b], axis=0)
    dmat2 = np.delete(dmat2, [a,b], axis=1)
    dmat2 = np.vstack([dmat2, new_row.reshape(1,-1)])
 
    new_col = np.zeros(len(ids)-1)
    new_col[0:len(ids)-2] = new_row
    new_col[len(ids)-2] = 0
    dmat2 = np.hstack([dmat2, new_col.reshape(-1,1)])


    return dmat2


def neighbor_joining(dist_mat, ids):

    n = dist_mat.shape[0]
    next_node = n
    tree = nx.Graph()
    center = 2*n -3 
    for i in range(n):
        tree.add_edge(i, center)


    for i in range(n-3):
        q_matrix = compute_q_matrix(dist_mat)
        a, b = (np.unravel_index(q_matrix.argmin(), q_matrix.shape))
    
        tree.remove_edge(ids[a],center )
        tree.remove_edge(ids[b], center)
        tree.add_edge(next_node,ids[a])
        tree.add_edge(next_node, ids[b])
        tree.add_edge(center, next_node)
        dist_mat = update_distance_matrix(dist_mat, a,b, ids, next_node)

        to_remove = [ids[a], ids[b]]
        for val in to_remove:
            ids.remove(val)

        ids.append(next_node)
        # print(ids)

        next_node += 1


    return tree

def hamming_distance(s1, s2):

    return (np.array(s1) != np.array(s2)).sum()


def create_distance_matrix(alignment):
    
    dist_mat = np.array([[hamming_distance(s1, s2) for k2, s2 in alignment.items()] for k1, s1 in alignment.items()])
    return dist_mat.astype(float)


# s1 = ["a", "b", "c"]
# s2 = ["c", "b", "c"]
# s3 = ["a", "c", "c"]
# s4 = ["a", "a", "a"]
# alignment = {1: s1, 2: s2, 3 : s3, 4: s4}
# ids = list(alignment.keys())
# print(hamming_distance(s1,s2))
# dmat = create_distance_matrix(alignment)
# tree = neighbor_joining(dmat, ids)
# print(list(tree.edges()))
# print(dmat)
# dist = np.array([[0,5,9,9,8],  [5,0,10,10,9],  [9,10,0,8,7], [9,10,8,0,3],  [8,9,7,3,0]], dtype=float)

# ids = [i for i in range(dist.shape[0])]

# tree = neighbor_joining(dist, ids)
# print(list(tree.edges))