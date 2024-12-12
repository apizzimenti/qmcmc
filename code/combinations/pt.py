
from matrices import matrices, distributions, weights, swaps
import numpy as np

# Make arrays more legible.
np.set_printoptions(suppress=True)

# Construct matrices, initial distributions, weights, and a swap probability
# matrix.
depth = 4
lowest, highest, step = 1, 256, 16
betas = np.array(range(lowest, highest+1, step))
# betas = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

T = len(betas)
M = matrices(betas, D=depth)
D = distributions(betas, depth)
W = weights(depth)

# For some fixed number of iterations N and an interval n, iterate.
N = 2**(depth+2)
n = depth

S = swaps(betas, W, M)

# Save data?
Z = np.zeros((N//n, T, len(W)))

for k in range(N//n):
    # Update each replica.
    for j in range(T):
        D[j] = D[j]@np.linalg.matrix_power(M[j], n)

    # Do "replica swaps."
    for j in range(T):
        D[j] = D[j]@S[j]
    
    Z[k] = D

print(Z[-1].round(2))
np.savez("iterations", **{ str(k): Z[k] for k in range(N//n) })
