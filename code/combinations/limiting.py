
import numpy as np
import networkx as nx
import scipy as scp
import matplotlib.pyplot as plt

betas = np.linspace(1, 200, 1000)
# betas = range(1, 50)

plt.rcParams["text.usetex"] = True
plt.rc("text.latex", preamble=r"\usepackage{nicefrac}")

fig, ax = plt.subplots()
size = size = "x-large"
ax.set_yticks([1/4, 1/2, 3/4, 1])
ax.set_yticklabels([r"$\nicefrac 14$", r"$\nicefrac 12$", r"$\nicefrac 34$", r"$1$"], fontsize=size)
ax.set_ylim(1/4-0.05, 1+0.05)

# ax.set_xticks([1, 25, 50, 75])
# ax.set_xlim(1-0.05, 75+0.05)
D = 32

for L in range(2, D+1, 2):
    X = []
    Y = []
    for beta in betas:
        p = np.exp(-beta/2)

        def V(r): return 2**r
        def rb(r): return 1
        def br(r): return 2*((1-p)/3)
        def bs(r): return V(r)*((1-p**(1/r))/3)
        def bd(r): return V(r+1)*(1/3)
        def bu(r): return V(r)*((p**(1/r))/3)
        def bb(r): return bs(r) + bu(r) + bd(r)
        def ll(r): return V(r)*(1-p**(1/r))
        def lb(r): return V(r)*(p**(1/r))

        np.set_printoptions(suppress=True)
        G = nx.path_graph(L+2)
        A = np.array(nx.adjacency_matrix(G).toarray()).astype(float)
        for k in range(1, len(G)): A[k,k] = 1

        # top branches
        for r in range(1, 2):
            nz = list(A[r].nonzero()[0])
            for n, q in zip(nz, [br(r), bs(r), bd(r)]):
                A[r,n] = q

        # middle branches
        for r in range(1, len(G)-1):
            nz = list(A[r].nonzero()[0])
            for n, q in zip(nz, [bu(r), bs(r), bd(r)]):
                A[r,n] = q

        # leaves
        r = len(G)-1
        for n, q in zip(list(A[r].nonzero()[0]), [lb(r), ll(r)]): A[r,n]=q

        # print(A)
        B = A/A.sum(axis=1, keepdims=1)

        # print(B)

        P = np.linalg.matrix_power(B, 2**(L+2))
        g = np.array([1]+[0]*(L-1+2))@P
        prob = g[-1]
        X.append(beta)
        Y.append(prob)
    
    ax.plot(X, Y, lw=1.5, color="k", alpha=((D+1)-L)/D)

plt.savefig("figures/finds-leaf-at-depth.png", dpi=600, bbox_inches="tight")
    