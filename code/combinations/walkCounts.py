
import numpy as np
import pandas as pd
import math
from walkLinearProgrammingSolutions import sufficient
import scipy as scp

split = np.array([
    [0,1,0],
    [2,1,1],
    [0,2,1]
])

initial = np.array([
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 4]
])

adj = np.array([
    [0,1,0,0,1,0,0],
    [1,1,1,1,0,0,0],
    [0,1,1,0,0,0,0],
    [0,1,0,1,0,0,0],
    [1,0,0,0,1,1,1],
    [0,0,0,0,1,1,0],
    [0,0,0,0,1,0,1]
])

N = 17
records = []
leafrecords = []
betas = list(range(1, 25, 2))

def compute(k):
    walks = np.linalg.matrix_power(adj, k)
    total = (walks[0]).sum()
    root = walks[0,0]
    branch = walks[0,1]+walks[0,4]
    leaf = walks[0,2]+walks[0,3]+walks[0,5]+walks[0,6]
    # slns = { n: sufficient(k, n) for n in ["root", "branch", "leaf"] }
    return root, branch, leaf, total

solutions = {
    k: compute(k) for k in range(2, N)
}

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

    # for k in solutions:
    #     S, root, branch, leaf, total = solutions[k]
    #     print("count:")
    #     print("\ttotal:\t", total)
    #     print("\troot:\t", root)
    #     print("\tbrnch:\t", branch)
    #     print("\tleaf:\t", leaf)

    #     prob = np.linalg.matrix_power(P,k)
    #     proot = prob[0,0]
    #     pbranch = prob[0,1]+prob[0,4]
    #     pleaf = prob[0,2]+prob[0,3]+prob[0,5]+prob[0,6]
    #     ptotal = prob[0].sum()

    #     print("weighted:")
    #     print("\ttotal:\t", ptotal)
    #     print("\troot:\t", proot)
    #     print("\tbrnch:\t", pbranch)
    #     print("\tleaf:\t", pleaf)

    #     print()

    for k in solutions:
        root, branch, leaf, total = solutions[k]
        record = {"T": beta}

        prob = np.linalg.matrix_power(P,k)
        hprob = np.linalg.matrix_power(HP, k)

        proot = prob[0,0]
        pbranch = prob[0,1]+prob[0,4]
        pleaf = prob[0,2]+prob[0,3]+prob[0,5]+prob[0,6]
        ptotal = prob[0].sum()

        hroot = hprob[0,0]
        hbranch = hprob[0,1]+hprob[0,4]
        hleaf = hprob[0,2]+hprob[0,3]+hprob[0,5]+hprob[0,6]
        htotal = hprob[0].sum()

        for name, count, probability, hinted in zip(["root", "branch", "leaf"], [root, branch, leaf], [proot, pbranch, pleaf], [hroot, hbranch, hleaf, htotal]):
            record.update({
                f"P_{name}": probability,
                f"P_{name}_hinted": hinted,
                f"P_count_{name}": count/total
            })

        records.append(record)

for beta in betas:
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

    for k in solutions:
        record = { "T": beta }
        walks = np.linalg.matrix_power(adj, k)
        prob = np.linalg.matrix_power(P,k)
        hprob = np.linalg.matrix_power(HP, k)

        hintcount = walks[0,2]
        anycount = walks[0,3]+walks[0,5]+walks[0,6]
        count = hintcount + anycount

        phintleaf = prob[0,2]
        panyleaf = prob[0,3]+prob[0,5]+prob[0,6]

        hhintleaf = hprob[0,2]
        hanyleaf = hprob[0,3]+hprob[0,5]+hprob[0,6]

        record.update({
            f"P_count_hint": hintcount/count,
            f"P_count_any": anycount/count,
            f"P_ignorant_hint": phintleaf/(phintleaf+panyleaf),
            f"P_ignorant_any": panyleaf/(phintleaf+panyleaf),
            f"P_hint_hint": hhintleaf/(hhintleaf+hanyleaf),
            f"P_hint_any": hanyleaf/(hhintleaf+hanyleaf)
        })

        leafrecords.append(record)
        

probabilities = pd.DataFrame().from_records(records)
probabilities["gutcheck"] = probabilities["P_root"]+probabilities["P_branch"]+probabilities["P_leaf"]
probabilities["gutcheck_hinted"] = probabilities["P_root_hinted"]+probabilities["P_branch_hinted"]+probabilities["P_leaf_hinted"]
probabilities.to_csv("treecounts.csv", index=False)

probabilities = pd.DataFrame().from_records(leafrecords)
probabilities.to_csv("leafcounts.csv", index=False)