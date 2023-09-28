
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.patches import Circle, Polygon
from gerrytools.plotting import latex

# Read in data and initialize the plots.
trials = pd.read_csv("output/runs/estimate-3200.csv").to_records()

plt.rcParams.update({ "text.usetex": True })
fig, (dartboard, estimate) = plt.subplots(1, 2, figsize=(6, 3))

# For each timestamp, create a plot and write it to file; we'll set an interval
# at ~25 just so we don't have too many plots.
interval = 20

# Set a tolerance for our estimation.
epsilon = 1/100

# Set up the plots appropriately.
dartboard.set_xlim(-1.1, 1.1)
dartboard.set_ylim(-1.1, 1.1)
dartboard.set_aspect("equal")
dartboard.add_patch(Polygon([(-1,-1), (1,-1), (1,1), (-1,1)], facecolor="None", edgecolor="k", lw=2))
dartboard.add_patch(Circle((0,0), 1, facecolor="None", edgecolor="k", lw=2, ls=":"))
dartboard.axis("off")

estimate.set_ylim(3-epsilon, 3.5+epsilon)
estimate.set_yticks([3, np.pi, 3.5], [r"$3$", r"$\pi$", r"$3.5$"])
estimate.axhspan(np.pi-epsilon, np.pi+epsilon, color=latex["Amber"], alpha=1/8)
estimate.axhline(y=np.pi, ls=":", lw=1, color=latex["Amethyst"])


with tqdm(total=len(trials)) as bar:
    for stamp in range(len(trials)):
        trial = trials[stamp]
        last = trials[stamp-1] if stamp > 0 else trial

        # Plot points on the dartboard...
        dartboard.plot(
            trial["x"], trial["y"],
            marker="x" if trial["hit"] else "o",
            markeredgecolor="r" if trial["hit"] else "b",
            markerfacecolor="None",
            markersize=3,
            alpha=1/2
        )

        # ... and our running estimate on the plot.
        estimate.plot(
            [last["t"], trial["t"]], [last["estimate"], trial["estimate"]],
            marker=".", markersize=2, color=latex["Alizarin"], lw=1
        )

        # Write to file if the interval's right.
        if not (stamp % interval):
            plt.savefig(f"output/figures/dartboards/{str(stamp).zfill(6)}.jpg", dpi=200, bbox_inches="tight")

        bar.update()
