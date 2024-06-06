
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

walks = pd.read_csv("treecounts.csv")
temps = np.array(walks["T"].unique())
colors = mpl.colormaps["coolwarm"]
norm = mpl.colors.Normalize(0, 1)
size = "x-large"

plt.rcParams["text.usetex"] = True
plt.rc("text.latex", preamble=r"\usepackage{nicefrac}")

for ext in ["", "_hinted"]:
    for prob in ["root", "branch", "leaf"]:
        fig, ax = plt.subplots()
        for temp in temps:
            c = colors(1-temp/max(temps))
            probabilities = walks[walks["T"] == temp][f"P_{prob+ext}"]
            X = list(range(2, len(probabilities)+2))
            ax.plot(X, probabilities, lw=1.5, color=c)

        base = walks[walks["T"] == temps[0]][f"P_count_{prob}"]
        X = list(range(2, len(probabilities)+2))
        ax.plot(X, base, lw=1.5, ls=":", color="k")
        ax.set_yticks([0, 1/4, 1/2, 3/4, 1])
        ax.set_yticklabels([r"$0$", r"$\nicefrac 14$", r"$\nicefrac 12$", r"$\nicefrac 34$", r"$1$"], fontsize=size)
        ax.set_ylim(0, 1)

        plt.savefig(f"figures/{prob+ext}.jpeg", bbox_inches="tight", dpi=600)
        plt.close()



walks = pd.read_csv("leafcounts.csv")
temps = np.array(walks["T"].unique())

for ext in ["_hint", "_any"]:
    for prob in ["ignorant", "hint"]:
        fig, ax = plt.subplots()

        for temp in temps:
            c = colors(1-temp/max(temps))
            probabilities = walks[walks["T"] == temp][f"P_{prob+ext}"]
            X = list(range(2, len(probabilities)+2))
            ax.plot(X, probabilities, lw=1.5, color=c)

        base = walks[walks["T"] == temps[0]][f"P_count{ext}"]
        X = list(range(2, len(probabilities)+2))
        ax.plot(X, base, lw=1.5, ls=":", color="k")
        ax.set_yticks([0, 1/4, 1/2, 3/4, 1])
        ax.set_yticklabels([r"$0$", r"$\nicefrac 14$", r"$\nicefrac 12$", r"$\nicefrac 34$", r"$1$"], fontsize=size)
        ax.set_ylim(0, 1)

        plt.savefig(f"figures/conditional-{prob+ext}.png", transparent=True, bbox_inches="tight", dpi=600)
        plt.close()


fig, bar = plt.subplots(figsize=(1/3,5))
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colors), cax=bar, orientation="vertical")
cbar.ax.invert_yaxis()
T = 1-(temps/temps.max())
cbar.ax.set_yticks(T)
cbar.set_ticklabels([rf"${t}$" for t in temps], fontsize=size)
plt.savefig("figures/colorbar-vertical.jpeg", bbox_inches="tight", dpi=600)
plt.close()

fig, bar = plt.subplots(figsize=(5, 1/3))
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colors), cax=bar, orientation="horizontal")
cbar.ax.invert_xaxis()
T = 1-(temps/temps.max())
cbar.ax.set_xticks(T)
cbar.set_ticklabels([rf"${t}$" for t in temps], fontsize=size)
plt.savefig("figures/colorbar-horizontal.jpeg", bbox_inches="tight", dpi=600)
