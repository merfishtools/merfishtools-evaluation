import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.figure(figsize=[snakemake.config["plots"]["figsize"][1]] * 2)

estimates = pd.read_table(snakemake.input.estimates, index_col=[0, 1])
plt.scatter(estimates["expr_map"], estimates["expr_naive"], s=2, c="k", alpha=0.6, edgecolors="face")
sns.despine()
plt.xlabel("MAP")
plt.ylabel("naive")
plt.savefig(snakemake.output[0], bbox_inches="tight")
