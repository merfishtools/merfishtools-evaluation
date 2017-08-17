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

feature = snakemake.wildcards.feature
if feature == "maxexpr":
    feature = counts["feat"][counts["total"].argmax()]

codeword = codebook[feature]
counts["dist"] = [hamming_distance(codeword, codebook[feat]) for feat in counts["feat"]]

counts.sort_values("dist", inplace=True)

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
sns.stripplot(x="feat", y="total", hue="dist", data=counts, rasterized=True, color="red", size=2)
plt.xticks([])
plt.ylim((0, 500))
plt.legend().set_title("hamming distance to {}".format(feature))
plt.ylabel("raw counts")
plt.xlabel("genes")

sns.despine(offset=5)
plt.savefig(snakemake.output.strip, bbox_inches="tight")


plt.figure()
dist = 2 if 2 in counts["dist"].values else 4
for _, d in counts[counts["dist"] != dist].groupby("cell"):
    sns.kdeplot(np.log10(d["total"]), color="k", label=None, legend=False, linewidth=1, alpha=0.5)
for _, d in counts[counts["dist"] == dist].groupby("cell"):
    sns.kdeplot(np.log10(d["total"]), color="r", label=None, legend=False, linewidth=1, alpha=0.5)

if snakemake.wildcards.pred == "raw":
    plt.xlabel("log10 raw counts")
else:
    plt.xlabel("log10 posterior expression")

plt.ylabel("density")
sns.despine()

plt.savefig(snakemake.output.kde, bbox_inches="tight")
