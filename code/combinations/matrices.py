
import numpy as np
import rustworkx as rx
from itertools import product, combinations_with_replacement as cwr


def matrices(betas, D=2):
    # Set the default probability parameter.
    p = lambda beta: np.exp(-beta)
    M = []

    # Construct the graph, add self-loops on all edges save for the root; get the
    # transition matrix as a mask for relabeling later. 
    T = rx.generators.full_rary_tree(2, 2**(D+1)-1)
    Tv = T.node_indices()
    for v in Tv[1:]: T.add_edge(v, v, None)
    mask = np.array(rx.adjacency_matrix(T))

    # Construct transition matrices for the set of temperatures.
    dimension = np.array(range(len(mask)))

    # Indices of root -> branch vertices.
    rootDown = np.array([[0, 0], [1, 2]])
    branchToRoot = rootDown[[1, 0]]

    # Indices of loops at branches.
    branchLoops = np.array([dimension[1:2**D-1], dimension[1:2**D-1]])

    # Indices of downward moves from branch -> branch vertices.
    branchDown = np.array([
        np.repeat(dimension[1:2**D-1], 2), dimension[3:2**(D+1)]
    ])

    # Indices of upward moves from branch -> branch (or branch -> root) vertices.
    branchUp = branchDown[[1, 0]][:, :2**(D+1)-2**D-4]

    # Indices of upward moves from leaf -> branch vertices.
    leafUp = np.array([
        dimension[2**(D)-1:],
        np.repeat(dimension[2**(D)-2**(D-1)-1:2**(D)-1], 2)
    ])

    # Indices of loops at leaves.
    leafLoops = np.array([dimension[2**D-1:], dimension[2**D-1:]])

    # Powers by depth. If there are none, 
    branchPowers = np.repeat(np.array(
        sum([[1/k]*2**(k-1) for k in range(1, D)], [])
    ), 2)

    leafPowers = np.array([1/D]*2**(D))

    # Set probabilities.
    for beta in betas:
        Q = mask.copy()
        Q[*rootDown] = 1/2
        Q[*branchToRoot] = p(beta)/3
        Q[*branchDown] = 1/3
        Q[*branchLoops] = (1-p(beta)**branchPowers)/3
        Q[*branchUp] = (p(beta)**branchPowers[2:])/3

        Q[*leafUp] = p(beta)**leafPowers
        Q[*leafLoops] = 1-p(beta)**leafPowers

        M.append(Q)
        # np.savetxt(f"./matrices/{D}-{beta}.npz", mask)

    return np.array(M)


def distributions(T, D=2):
    D = np.array([np.zeros(2**(D+1)-1) for _ in range(len(T))])
    for d in D: d[0] = 1

    return D


def weights(D=2):
    T = rx.generators.full_rary_tree(2, 2**(D+1)-1)

    W = np.array([1] + [2**-int(np.log2(j+1)) for j in range(1, len(T))])
    W[[2**k-1 for k in range(1, int(np.log2(len(T)))+1)]] *= 1/2

    return W


def swaps(T, W, M):
    """
    Args:
        T (Iterable): Iterable of temperatures.
        W (Iterable): Iterable of weights.

        Returns:
            A 3D matrix of transition probabilites.
    """
    S = np.zeros((len(T), len(W), len(W)))
    indices = list(product(range(len(W)), range(len(W))))

    for t in range(len(T)-1):
        lower, higher = T[t], T[t+1]

        for i, j in indices:
            S[t][i,j] = min(1, np.exp(-higher*W[j])/np.exp(-lower*W[i]))

        S[t] = S[t]/S[t].sum(axis=1, keepdims=True)

    # Set the last swap matrix to be the transition matrix.
    S[-1] = M[-1]
    return S

