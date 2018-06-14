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

common_genes = small.index.intersection(large.index)
means = lambda estimates: estimates.loc[common_genes].mean(axis="columns")

small = means(small)
large = means(large)
smape = (large - small).abs().sum() / (large + small).sum()
#relative_accuracy = np.log10(large / small)
small = np.log10(small)
large = np.log10(large)

min_value = np.log10(0.1)#min(small.min(), large.min())
max_value = np.log10(100)#max(small.max(), large.max())

fig = plt.figure(figsize=snakemake.config["plots"]["figsize"])
fig.add_subplot(111, aspect='equal')

sns.regplot(small, large, line_kws={"color": "red"}, scatter_kws={"color": "black", "s": 3})
plt.plot([min_value, max_value], [min_value, max_value], "k--")

label = "raw counts" if snakemake.wildcards.type == "counts" else "posterior counts"
plt.xlabel("MHD4 " + label)
plt.ylabel("MHD2 " + label)
plt.xlim((min_value, max_value))
plt.ylim((min_value, max_value))
plt.text(x=0.5, y=1, s="SMAPE={:.0%}".format(smape), horizontalalignment="center", verticalalignment="top", transform=plt.gca().transAxes)

logticks = lambda ticks: ["$10^{{{:.0f}}}$".format(y) for y in ticks]
ax = plt.gca()
ax.locator_params(nbins=4)
ax.set_xticklabels(logticks(ax.get_xticks()))
ax.set_yticklabels(logticks(ax.get_yticks()))
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
