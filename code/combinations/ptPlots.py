
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Load data.
Z = np.load("iterations.npz")
Z = { int(k): Z[k] for k in Z.files }
d = 4
pos = sum(2**k for k in range(d))

plt.rcParams["text.usetex"] = True
plt.rc("text.latex", preamble=r"\usepackage{nicefrac}")
fig, ax = plt.subplots(figsize=(6,2))

lowest, highest, step = 1, 256, 16
betas = np.array(range(lowest, highest+1, step))
colors = mpl.colormaps["coolwarm"]
norm = mpl.colors.Normalize(0, 1)

# Plot the probability we find the hinted leaf over time.
for k in range(len(Z[0])-1):
	Y = [Z[t][k,pos] for t in range(len(Z))]
	ax.plot(np.arange(len(Y)), Y, color=colors(betas[k]/betas.max()), alpha=1/2 if k < len(Z[0])-2 else 1, zorder=10000)


# Re-set x-ticks.
# ax.set_xlim(0, len(Z))
xticks = list(ax.get_xticks())
xticklabels = [rf"${2+int(2*t)}$" for t in xticks]
ax.set_xticks(xticks[1:-1])
ax.set_xticklabels(xticklabels[1:-1])
ax.set_yticks([0, 1/4, 1/2, 3/4, 1])
ax.set_yticklabels([r"$0$", r"$\nicefrac 14$", r"$\nicefrac 12$", r"$\nicefrac 34$", r"$1$"], fontsize="medium")
ax.set_ylim(0,1)
ax.spines[['right', 'top']].set_visible(False)

plt.savefig("figures/pt.jpeg", bbox_inches="tight", dpi=600, pad_inches=0)