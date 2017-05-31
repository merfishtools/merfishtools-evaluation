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
    return sum(x != y for x, y in zip(a, b))

# read data
codebook = pd.read_table(snakemake.input.codebook, index_col=0, dtype=np.dtype(str))["codeword"]

if snakemake.wildcards.pred == "raw":
    counts = pd.read_table(snakemake.input.counts, index_col=1)
    counts["total"] = counts["corrected"] + counts["exact"]
else:
    counts = pd.read_table(snakemake.input.exprs, index_col=1)
    counts["total"] = counts["expr_map"]
counts = counts.loc[codebook.index]
counts = counts.reset_index()
codeword = codebook[snakemake.wildcards.feature]
counts["dist"] = [hamming_distance(codeword, codebook[feat]) for feat in counts["feat"]]

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
sns.stripplot(x="feat", y="total", hue="dist", data=counts, rasterized=True, color="red", size=2)
plt.xticks([])
plt.ylim((0, 500))
plt.legend().set_title("hamming distance to {}".format(snakemake.wildcards.feature))
plt.ylabel("raw counts")
plt.xlabel("genes")

sns.despine(offset=5)
plt.savefig(snakemake.output[0], bbox_inches="tight")
