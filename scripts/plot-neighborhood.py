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
counts = counts.reindex(codebook.index)
codeword = codebook[snakemake.wildcards.feature]
dist = [hamming_distance(codeword, codebook[feat]) for feat in counts.index]
counts = pd.DataFrame({"feat": counts.index, "count": counts, "dist": dist})

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
sns.stripplot(x="feat", y="count", hue="dist", data=counts)

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
