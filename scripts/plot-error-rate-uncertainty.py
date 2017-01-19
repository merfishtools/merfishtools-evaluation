import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

codebook = pd.read_table(snakemake.input.codebook, index_col=0)

errors = []

for uncertainty in [0, 10, 20, 30]:
    u = "err+{}%".format(uncertainty) if uncertainty > 0 else "default"
    for mean, posterior_counts, known_counts in zip(snakemake.params.means, snakemake.input.get(u), snakemake.input.known_counts):
        known_counts = pd.read_table(known_counts, index_col=[0, 1])
        known_counts = known_counts.reindex(codebook.index, level=1)

        posterior_estimates = pd.read_table(posterior_counts, index_col=[0, 1])
        posterior_estimates = posterior_estimates.reindex(known_counts.index, fill_value=0)

        errors.append(pd.DataFrame({"error": posterior_estimates["expr_map"] - known_counts["count"]), "mean": mean, "uncertainty": uncertainty})

errors = pd.concat(errors)
sns.violinplot(x="uncertainty", y="error", data=pred_errors, bw=1, inner="quartile", palette=colors, linewidth=1)
plt.plot(plt.xlim(), [0, 0], "-k", linewidth=1, zorder=-5)
sns.despine()

plt.xlabel("error rate underestimation (%)")
plt.ylabel("predicted - truth")

plt.savefig(snakemake.output[0], bbox_inches="tight")
