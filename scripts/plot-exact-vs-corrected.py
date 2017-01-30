import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


sns.set(style="ticks", palette="colorblind", context="paper")
plt.figure(figsize=snakemake.config["plots"]["figsize"])


for f in snakemake.input:
    counts = pd.read_table(f, index_col=[0, 1])
    plt.scatter(counts["exact"], counts["corrected"], s=2, c="k", alpha=0.6, edgecolors="face", rasterized=True)

plt.xlabel("exact counts")
plt.ylabel("corrected counts")
plt.xlim((0, 1000))
plt.ylim((0, 1000))
plt.locator_params(nbins=4)
#ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
