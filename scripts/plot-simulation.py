import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import seaborn as sns
from seaborn.palettes import blend_palette
from seaborn.utils import set_hls_values


def ci_error(lower, upper, truth):
    below = np.maximum(lower - truth, 0)
    above = np.maximum(truth - upper, 0)
    return below + above


epsilon = 0.001

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)

codebook = pd.read_table(snakemake.input.codebook, index_col=0)

errors = []
counts = []
ci_errors = []
for i, (mean, posterior_counts, raw_counts, known_counts) in enumerate(zip(
    snakemake.params.means,
    snakemake.input.posterior_counts,
    snakemake.input.raw_counts,
    snakemake.input.known_counts)):

    posterior_estimates = pd.read_table(posterior_counts, index_col=[0, 1])
    raw_counts = pd.read_table(raw_counts, index_col=[0, 1])

    
    raw_counts = raw_counts["exact"] + raw_counts["corrected"]
    known_counts = pd.read_table(known_counts, index_col=[0, 1])
    known_counts = known_counts.reindex(codebook.index, level=1)
    # remove PTPN14 because we have artificially increased it's simulated expression
    # this will bias our plots
    known_counts = known_counts[known_counts.index.get_level_values(1) != "PTPN14"]

    raw_counts = raw_counts.reindex(known_counts.index, fill_value=0)
    posterior_estimates = posterior_estimates.reindex(known_counts.index, fill_value=0)

    dropouts = known_counts[(known_counts["count"] > 0) & (raw_counts == 0)]
    print("dropouts", dropouts)

    counts.append(pd.DataFrame({"known": known_counts["count"], "raw": raw_counts, "posterior": posterior_estimates["expr_map"]}))
    errors.append(pd.DataFrame({"error": raw_counts - known_counts["count"], "mean": mean, "type": "raw", "known": known_counts["count"]}))
    errors.append(pd.DataFrame({"error": posterior_estimates["expr_map"] - known_counts["count"], "mean": mean, "type": "posterior", "known": known_counts["count"]}))
    errors.append(pd.DataFrame({
        "error": ci_error(posterior_estimates["expr_ci_lower"], posterior_estimates["expr_ci_upper"], known_counts["count"]),
        "mean": mean,
        "type": "ci",
        "known": known_counts["count"]}))

counts = pd.concat(counts)

def plot_hexbin(type, path, color):
    color_rgb = mpl.colors.colorConverter.to_rgb(color)
    colors = [set_hls_values(color_rgb, l=l) for l in np.linspace(1, 0, 12)]
    cmap = blend_palette(colors, as_cmap=True)

    plt.figure(figsize=0.75 * np.array(snakemake.config["plots"]["figsize"]))
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

x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x * 1.5, y))

pred_errors = errors[(errors["type"] == "raw") | (errors["type"] == "posterior")]
#bins = pd.cut(pred_errors["known"], 
#              [0, 6, 11, 16, 21, 26, 30, 100000], 
#              right=False, 
#              labels=["0-5", "6-10", "11-15", "16-20", "21-25", "26-30", "≥30"])
#pred_errors["bin"] = bins


sns.violinplot(x="mean", y="error", hue="type", data=pred_errors, bw=1, split=True, inner="quartile", palette=colors, linewidth=1)
plt.plot(plt.xlim(), [0, 0], "-k", linewidth=1, zorder=-5)

plt.xlabel("mean expression")
plt.ylabel("predicted - truth")
plt.legend(loc="lower left")
sns.despine()

plt.savefig(snakemake.output.violin, bbox_inches="tight")


# plot CI errors
ci_errors = errors[errors["type"] == "ci"]["error"].astype(np.int64)

hist = ci_errors.value_counts(normalize=True)
hist = hist[[0, 1, 2]]

plt.figure(figsize=(0.2 * hist.size, y))
sns.barplot(hist.index, hist, color=colors[1])
plt.ylabel("fraction")
plt.xlabel("CI error")
sns.despine()

plt.savefig(snakemake.output.ci_errors, bbox_inches="tight")


# print RMSE
rmse = lambda errors: np.sqrt((errors ** 2).mean())
print("posterior RMSE", rmse(errors[errors["type"] == "posterior"]["error"]))
print("raw RMSE", rmse(errors[errors["type"] == "raw"]["error"]))

# write errors
errors.to_csv(snakemake.output.errors, sep="\t")

# RMSE
# rmse = lambda errors: np.sqrt((errors ** 2).mean())
# quantile = lambda errors: errors.quantile(0.98)
# maxval = lambda errors: errors.max()
# sd = lambda errors: errors.std()
#
# ci_errors = ci_errors["error"].abs()
# raw_errors = errors.loc[errors["type"] == "raw", "error"].abs()
# posterior_errors = errors.loc[errors["type"] == "posterior", "error"].abs()
# print(posterior_errors.max())
# all_errors = [ci_errors, raw_errors, posterior_errors]
#
# d = pd.DataFrame({
#     "rmse": list(map(rmse, all_errors)),
#     "upper_quantile": list(map(quantile, all_errors)),
#     "max": list(map(maxval, all_errors)),
#     "sd": list(map(sd, all_errors))},
#     index=["ci", "raw", "posterior"])
# d.to_csv(snakemake.output.rmse, sep="\t")
