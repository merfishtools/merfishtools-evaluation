import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import seaborn as sns
from seaborn.palettes import blend_palette
from seaborn.utils import set_hls_values

epsilon = 0.001

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

errors = []
counts = []
for i, (mean, posterior_counts, raw_counts, known_counts) in enumerate(zip(
    snakemake.params.means,
    snakemake.input.posterior_counts,
    snakemake.input.raw_counts,
    snakemake.input.known_counts)):

    posterior_counts = pd.read_table(posterior_counts, index_col=[0, 1])["expr_ev"]
    raw_counts = pd.read_table(raw_counts, index_col=[0, 1])
    raw_counts = raw_counts["exact"] + raw_counts["corrected"]
    known_counts = pd.read_table(known_counts, index_col=[0, 1])
    codebook = "mhd4" if snakemake.wildcards.dist == "4" else "mhd2"
    known_counts = known_counts[known_counts[codebook]]

    raw_counts = raw_counts.reindex(known_counts.index, fill_value=0)
    posterior_counts = posterior_counts.reindex(known_counts.index, fill_value=0)
    biased = posterior_counts[posterior_counts > 70].index
    unbiased = posterior_counts[(posterior_counts > 49) & (posterior_counts < 51)].index

    if mean == 30:
        print(posterior_counts.loc[biased])
        print(raw_counts.loc[biased])
        print(known_counts.loc[biased])
        print("unbiased")
        print(posterior_counts.loc[unbiased])
        print(raw_counts.loc[unbiased])
        print(known_counts.loc[unbiased])


    #plt.plot(known_counts["count"], posterior_counts, "r.", label="conditional expectation" if i == 0 else "", zorder=1, alpha=0.01, rasterized=True)
    #plt.plot(known_counts["count"], raw_counts, "k.", label="raw counts" if i == 0 else "", zorder=0, alpha=0.01, rasterized=True)

    counts.append(pd.DataFrame({"known": known_counts["count"], "raw": raw_counts, "posterior": posterior_counts}))

    errors.append(pd.DataFrame({"error": raw_counts - known_counts["count"], "mean": mean, "type": "raw"}))
    errors.append(pd.DataFrame({"error": posterior_counts - known_counts["count"], "mean": mean, "type": "posterior"}))

counts = pd.concat(counts)
print(counts.describe())

def plot_hexbin(type, path, color):
    color_rgb = mpl.colors.colorConverter.to_rgb(color)
    colors = [set_hls_values(color_rgb, l=l) for l in np.linspace(1, 0, 12)]
    cmap = blend_palette(colors, as_cmap=True)

    plt.figure(figsize=0.75 * snakemake.config["plots"]["figsize"])
    plt.subplot(111, aspect="equal")
    #plt.scatter(counts["known"], counts[type], s=1, c="k", alpha=0.3, rasterized=True, edgecolors="face", marker="o")
    plt.hexbin(counts["known"], counts[type], cmap=cmap, gridsize=25, clip_on=True)

    maxv = max(plt.xlim()[1], plt.ylim()[1])

    plt.plot([0, maxv], [0, maxv], "--k")
    plt.xlim((0, maxv))
    plt.ylim((0,maxv))
    plt.ylabel("predicted")
    plt.xlabel("truth")
    sns.despine()
    plt.savefig(path, bbox_inches="tight")

colors = sns.xkcd_palette(["grey", "light red"])

plot_hexbin("raw", snakemake.output.scatter_raw, colors[0])
plot_hexbin("posterior", snakemake.output.scatter_posterior, colors[1])

errors = pd.concat(errors)
s = (errors["type"] == "posterior") & (errors["mean"] == 5)

x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x * 1.5, y))
sns.violinplot(x="mean", y="error", hue="type", data=errors, bw=1, split=True, inner="quartile", linewidth=1, palette=colors)
plt.plot(plt.xlim(), [0, 0], "-k", linewidth=1, zorder=-5)

plt.xlabel("mean expression")
plt.ylabel("predicted - truth")
plt.legend(loc="lower left")
sns.despine()

plt.savefig(snakemake.output.violin, bbox_inches="tight")
