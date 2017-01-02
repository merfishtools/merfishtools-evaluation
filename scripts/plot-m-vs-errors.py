import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


errors = []
for m, f in zip(snakemake.params.ms, snakemake.input):
    d = pd.read_table(f)
    d["m"] = m
    errors.append(d)
errors = pd.concat(errors)
errors["error"] = errors["error"].abs()
print(errors)

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
colors = sns.xkcd_palette(["grey", "light red", "red"])
g = sns.FacetGrid(errors, col="type")
g.map(sns.boxplot, "m", "error").despine()

#plt.xlabel("1-bits")
#plt.ylabel("absolute error")
#plt.legend(loc="best")
#sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
