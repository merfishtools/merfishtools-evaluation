import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import merfishtools

sns.set(style="ticks", context=snakemake.wildcards.context)
sns.set_palette("Paired", n_colors=7)

fig_x, fig_y = snakemake.config["plots"]["figsize"]
figsize = (fig_x * 4, fig_y * 4)
plt.figure(figsize=figsize)

exprs = [pd.read_table(f, index_col=0).stack(0) for f in snakemake.input.exprs]
exprs = pd.concat(exprs, keys=snakemake.params.expmnts)

diffexp = merfishtools.read_diffexp_estimates(snakemake.input.diffexp)
diffexp.sort_values("cv_map", inplace=True)

cdfs = merfishtools.read_cdf(snakemake.input.diffexp_cdf)

significant = diffexp[diffexp["diff_fdr"] <= 0.05]
if not significant.empty:
    significant = significant[~significant.index.str.startswith("blank") & ~significant.index.str.startswith("notarget")]

    exprs = exprs.loc[:, significant.index]

    n = exprs.index.levels[1].size
    for i, (gene, gene_exprs) in enumerate(exprs.groupby(level=1)):
        vmax = gene_exprs.quantile(0.95)
        ax = plt.subplot(4, 4, i + 1)
        for _, exp_exprs in gene_exprs.groupby(level=0):
            ax = sns.kdeplot(exp_exprs, ax=ax, linewidth=1, alpha=0.5, clip=(0, vmax), clip_on=True)
        if i == n - 1:
            plt.xlabel("gene expression")
        plt.ylabel(gene)
        plt.xlim((0, vmax))
        sns.despine()
        plt.yticks([])

        # plot cdf
        cdf_ax = inset_axes(ax, width="30%", height=0.5, loc=1)
        cdf = cdfs.loc[gene]
        est = significant.loc[gene]
        
        merfishtools.plot_cdf(cdf, map_value=est["cv_map"], credible_interval=est[["cv_ci_lower", "cv_ci_upper"]], legend=False)
        
        plt.setp(cdf_ax.get_xticklabels(), rotation=45, ha="right")
        cdf_ax.tick_params(pad=1)
        plt.locator_params(nbins=4)
        sns.despine()

    plt.tight_layout()
plt.savefig(snakemake.output[0], bbox_inches="tight")
