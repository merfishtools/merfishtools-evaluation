import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
idx = pd.IndexSlice

import merfishtools

sns.set(style="ticks",
        palette="colorblind",
        context=snakemake.wildcards.context)
plt.figure(figsize=snakemake.config["plots"]["figsize"])
cdf = merfishtools.read_cdf(snakemake.input.expr).loc[idx[:, snakemake.wildcards.gene], :]
est = merfishtools.read_exp_estimates(snakemake.input.expr_est).loc[idx[:, snakemake.wildcards.gene], :].iloc[0]
print(est)
raw_counts = pd.read_table(
    snakemake.input.raw_counts,
    index_col=1).loc[snakemake.wildcards.gene]

count_exact = raw_counts["exact"]
count_corrected = raw_counts["corrected"]
count_total = count_exact + count_corrected
merfishtools.plot_pmf(
    cdf,
    map_value=est["expr_map"],
    credible_interval=est[["expr_ci_lower", "expr_ci_upper"]],
    legend=False)

# plot raw counts
ylim = plt.ylim()
ylim = [0, ylim[1]]
plt.vlines([count_total],
           *ylim,
           colors="grey",
           linestyles="--",
           label="total count",
           zorder=6)
plt.vlines([count_exact],
           *ylim,
           colors="grey",
           linestyles=":",
           label="exact count",
           zorder=6)
plt.ylim(ylim)

plt.xlabel("expression")
plt.ylabel("PMF")
if snakemake.wildcards.legend == "legend":
    #plt.legend(loc="upper left", bbox_to_anchor=(0.5, 1))
    plt.legend(loc="best")
plt.xlim([0, plt.xlim()[1]])

sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
