import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

data = []
for exprs, counts in zip(snakemake.input.estimates, snakemake.input.counts):
    exprs = pd.read_table(exprs, index_col=[0, 1])
    counts = pd.read_table(counts, index_col=[0, 1])

    d = pd.concat([exprs, counts], axis="columns")
    d["diff"] = d["expr_map"] - d["expr_naive"]
    data.append(d)
data = pd.concat(data)


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
fx, fy = snakemake.config["plots"]["figsize"]
plt.figure(figsize=[fx*2, fy])
colors = sns.xkcd_palette(["light red"])

hist = data["diff"].value_counts(normalize=True)
hist = hist.reindex(np.arange(hist.index.min(), hist.index.max(), dtype=np.int64), fill_value=0)

sns.barplot(hist.index, hist, color=colors[0], log=True)
plt.xlabel("MAP - naive")
plt.ylabel("fraction")
sns.despine()

plt.savefig(snakemake.output.hist, bbox_inches="tight")

plt.figure(figsize=[fx, fy])

plt.scatter(x="exact", y="corrected", c="diff", data=data, cmap="viridis", edgecolors="face", alpha=0.7, rasterized=True)

xlim = plt.xlim()
ylim = plt.ylim()
plt.xlim([0, xlim[1]])
plt.ylim([0, ylim[1]])
xlim = plt.xlim()
ylim = plt.ylim()

# helper line
vmin = min(xlim[1], ylim[1])
plt.plot([0, vmin], [0, vmin], "k--", linewidth=1, alpha=0.7)

plt.xlabel("exact")
plt.ylabel("corrected")
sns.despine()
cb = plt.colorbar()
cb.set_label("MAP - naive")

plt.savefig(snakemake.output.scatter, bbox_inches="tight")
