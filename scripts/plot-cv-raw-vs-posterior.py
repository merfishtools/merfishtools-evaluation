import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=[snakemake.config["plots"]["figsize"][1]] * 2)

def get_means(raw_counts):
    return raw_counts.sum(axis=1).groupby(level=1).mean()


def get_cvs(batch_means):
    g = batch_means.groupby(level=1)
    return g.std() / g.mean()


def read_raw_counts(f):
    return pd.read_table(f, index_col=[0, 1])


batch_means = [get_means(read_raw_counts(f)) for f in snakemake.input.raw_counts]
batch_means = pd.concat(batch_means, keys=range(len(batch_means)))

raw_cvs = get_cvs(batch_means)
posterior_cvs = pd.read_table(snakemake.input.diffexp, index_col=0)["cv_ev"]

plt.scatter(posterior_cvs, raw_cvs.loc[posterior_cvs.index], s=2, c="k", alpha=0.6, edgecolors="face")
limits = (0, max(plt.xlim()[1], plt.ylim()[1]))

plt.plot(limits, limits, "r--")

plt.xlim(limits)
plt.ylim(limits)
plt.xlabel("posterior CV")
plt.ylabel("raw CV")
sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")

print(raw_cvs)
