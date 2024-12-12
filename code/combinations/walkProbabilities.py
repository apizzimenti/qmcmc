
import numpy as np
import scipy as scp

betas = list(range(1, 25, 2))

for beta in betas:
    # r = (1/3)*np.exp(-beta/2)
    # s = ((1-np.exp(-beta/2))/3)
    # t = np.exp(-beta/4)

    r = (1/3)*np.exp(-beta/2)
    s = ((1-np.exp(-beta/2))/3)
    t = np.exp(-beta/4)

    P = np.array([
        [0, 1/2, 0, 0, 1/2, 0, 0],
        [r, s, 1/3, 1/3, 0, 0, 0],
        [0, t, 1-t, 0, 0, 0, 0],
        [0, t, 0, 1-t, 0, 0, 0],
        [r, 0, 0, 0, s, 1/3, 1/3],
        [0, 0, 0, 0, t, 1-t, 0],
        [0, 0, 0, 0, t, 0, 1-t],
    ])

    hr = (1/3)*np.exp(-beta*3/4)
    hs = ((1-np.exp(-beta*3/4))/3)
    ht = np.exp(-beta*3/8)

    HP = np.array([
        [0, 1/2, 0, 0, 1/2, 0, 0],
        [hr, hs, 1/3, 1/3, 0, 0, 0],
        [0, ht, 1-ht, 0, 0, 0, 0],
        [0, t, 0, 1-t, 0, 0, 0],
        [r, 0, 0, 0, s, 1/3, 1/3],
        [0, 0, 0, 0, t, 1-t, 0],
        [0, 0, 0, 0, t, 0, 1-t],
    ])

    print(scp.linalg.eig(HP))