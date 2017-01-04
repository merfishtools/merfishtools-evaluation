import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


PSEUDOCOUNTS = 1


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=[snakemake.config["plots"]["figsize"][1]] * 2)


def get_means(raw_counts):
    return raw_counts.groupby(level=1).mean() + PSEUDOCOUNTS


def get_cvs(batch_means):
    g = batch_means.groupby(level=1)
    return g.std() / g.mean()


def read_raw_counts(f):
    m = pd.read_table(f, index_col=[0, 1]).sum(axis=1).unstack(fill_value=0)
    return m.stack()


def normalize(raw_counts):
    quartiles = np.array([counts.quantile(0.75) for counts in raw_counts])
    mean_quartile = quartiles.mean()
    print([mean_quartile / q for counts, q in zip(raw_counts, quartiles)])
    return [counts * (mean_quartile / q) for counts, q in zip(raw_counts, quartiles)]

raw_counts = normalize([read_raw_counts(f) for f in snakemake.input.raw_counts])
batch_means = [get_means(counts) for counts in raw_counts]
batch_means = pd.concat(batch_means, keys=range(len(batch_means)))

raw_cvs = get_cvs(batch_means)

#print(raw_cvs[raw_cvs < 0.1])

posterior_cvs = pd.read_table(snakemake.input.diffexp, index_col=0)[snakemake.wildcards.estimate]

plt.scatter(posterior_cvs, raw_cvs.loc[posterior_cvs.index], s=2, c="k", alpha=0.6, edgecolors="face")
#sns.regplot(posterior_cvs, raw_cvs.loc[posterior_cvs.index], color="k", scatter_kws={"s": 2, "alpha": 0.6})
limits = (0, min(plt.xlim()[1], plt.ylim()[1]))
limits = (0, 1)

plt.plot(limits, limits, "r--")

plt.xlim(limits)
plt.ylim(limits)
plt.xlabel("posterior CV" + (" (conservative)" if snakemake.wildcards.estimate == "cv_ci_lower" else ""))
plt.ylabel("raw CV")
sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")

print(raw_cvs)
