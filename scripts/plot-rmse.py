import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


rmses = pd.DataFrame()
for m, f in zip(snakemake.params.ms, snakemake.input):
    rmse = pd.Series.from_csv(f, sep="\t")
    rmse["m"] = m
    rmses = rmses.append(rmse, ignore_index=True)

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
plt.plot(rmses["m"], rmses["raw"], ".k-", label="raw")
plt.plot(rmses["m"], rmses["posterior"], ".r-", label="posterior")
plt.plot(rmses["m"], rmses["ci"], ".r--", label="ci")

plt.xlabel("1-bits")
plt.ylabel("RMSE")
plt.xticks(snakemake.params.ms)
plt.ylim([0, 6])
plt.legend(loc="best")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
