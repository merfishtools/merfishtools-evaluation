import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
fx, fy = snakemake.config["plots"]["figsize"]
plt.figure(figsize=[fx*2, fy])
colors = sns.xkcd_palette(["light red"])

estimates = pd.concat([pd.read_table(f) for f in snakemake.input.estimates])
diff = estimates["expr_map"] - estimates["expr_naive"]
hist = diff.value_counts(normalize=True)
hist = hist.reindex(np.arange(hist.index.min(), hist.index.max(), dtype=np.int64), fill_value=0)

sns.barplot(hist.index, hist, color=colors[0], log=True)
plt.xlabel("MAP - naive")
plt.ylabel("fraction")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
