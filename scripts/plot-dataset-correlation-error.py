import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
import numpy as np

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

small = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.small], axis=1)
large = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.large], axis=1)

small_counts = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.small_counts], axis=1)
large_counts = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input.large_counts], axis=1)

common_genes = small.index.intersection(large.index)
means = lambda estimates: estimates.loc[common_genes].mean(axis="columns")

small = means(small)
large = means(large)
small_counts = means(small_counts)
large_counts = means(large_counts)

errors = pd.concat([pd.DataFrame({"mhd4": small, "error": np.log2((large + 1) / (small + 1)), "approach": "posterior estimates"}),
                    pd.DataFrame({"mhd4": small, "error": np.log2((large_counts + 1) / (small_counts + 1)), "approach": "raw counts"})])

min_value = 0.1#min(small.min(), large.min())
max_value = 100#max(small.max(), large.max())

fig = plt.figure(figsize=snakemake.config["plots"]["figsize"])

sns.swarmplot(y="error", x="approach", data=errors, palette=["red", "black"], size=4)
#plt.scatter("mhd4", "error", c="black", s=3, data=errors[errors["approach"] == "raw counts"])
#plt.scatter("mhd4", "error", c="red", s=3, data=errors[errors["approach"] == "posterior estimates"])
#plt.scatter(small_counts, large_counts - small_counts, s=3, c="black", label="raw counts")

#plt.ylabel("MHD2 - MHD4")
plt.xlabel("")
plt.ylabel("log2 fold change")
plt.legend()

sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
