from itertools import combinations
import math

import matplotlib as mpl
mpl.use("svg")
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

def hamming_distance(a, b):
    return sum(x == y for x, y in zip(a, b))

# read data
codebook = pd.read_table(snakemake.input.codebook, index_col=0, dtype=np.dtype(str))["codeword"]

cell = int(snakemake.wildcards.cell)

counts = pd.read_table(snakemake.input.counts, index_col=[0, 1]).loc[cell]
counts = counts["corrected"] + counts["exact"]

exprs = pd.read_table(snakemake.input.exprs, index_col=[0, 1])["expr_map"].loc[cell]
counts = counts.reindex(exprs.index, fill_value=0)

data = pd.DataFrame({"counts": counts, "exprs": exprs})

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.subplot(111, aspect="equal")
sns.stripplot(x=data=data
plt.scatter(counts, exprs, c="k", edgecolors="face", rasterized=True, alpha=0.5)

# plot neighborhood
codeword = codebook[snakemake.wildcards.feature]
neighbors = [feat for feat, w in codebook.iteritems() if hamming_distance(codeword, w)]

plt.scatter(counts[neighbors], exprs[neighbors], edgecolors="face", c="r")
plt.scatter(counts[snakemake.wildcards.feature], exprs[snakemake.wildcards.feature], edgecolors="face", c="r", s=2)

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
