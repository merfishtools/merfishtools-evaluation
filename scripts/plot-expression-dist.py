import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])

means = []
stds = []
for f in snakemake.input:
    exprs = pd.read_table(f, index_col=0)
    for cell in exprs:
        expr = exprs[cell]
        means.append(expr.mean())
        stds.append(expr.std())
        expr = np.log10(expr + 1)
        sns.kdeplot(expr, color="black", clip=[0, expr.max()], alpha=0.1, lw=1, label="", kernel="gau", bw="scott")

means = np.array(means)
stds = np.array(stds)
cvs = stds / means
sns.kdeplot(np.log10(means + 1), color="red", clip=[0, means.max()], label="means")
sns.kdeplot(np.log10(stds + 1), color="red", linestyle="--", clip=[0, stds.max()], label="standard deviations")


m = means.mean()
cv = cvs.mean()
print(m, cv)
#def nbinom_params(mean, cv):
#    variance = (cv * mean) ** 2
#    p = 1 - ((variance - mean) / variance)
#    n = (mean * p) / (1 - p)
#    return n, p

#d = np.random.negative_binomial(*nbinom_params(m, cv), 5000)
#print(np.mean(d), np.std(d), np.std(d) / np.mean(d))
#sns.kdeplot(np.log10(d + 1), color="red", linestyle=":", clip=[0, d.max()], label="simulated example")


plt.xlim([0, plt.xlim()[1]])
plt.xlabel("log10 counts")
plt.ylabel("density")
plt.legend(loc="upper right")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
