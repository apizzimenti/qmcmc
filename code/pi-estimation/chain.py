

import numpy as np
import pandas as pd


def hit(x, y):
    """
    Determines whether the given coordinates specify a point within the unit
    circle.
    """
    return int(x**2 + y**2 <= 1)

# Choose the number of steps; we're going to use 10k for now.
N = 3200

# Distributions on the appropriate interval.
X = lambda: np.random.uniform(low=-1, high=1)
Y = lambda: np.random.uniform(low=-1, high=1)

# Initializing our estimator, and keeping some records for visualizations later.
hits = 0
records = []

# Chain-ing!
for t in range(1, N):
    x, y = X(), Y()
    hits += hit(x,y)
    records.append(dict(t=t, x=x, y=y, hit=hit(x,y), hits=hits, estimate=4*(hits/t)))

# Write records to file. 
df = pd.DataFrame.from_records(records)
df.to_csv(f"output/runs/estimate-{str(N)}.csv", index=False)
