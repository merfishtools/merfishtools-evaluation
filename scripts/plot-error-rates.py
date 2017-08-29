import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x, y))


error_rates = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input],
                        keys=snakemake.params.experiments)
error_rates = pd.DataFrame(error_rates.stack()).reset_index()
error_rates.columns = ["pos", "error", "rate"]

sns.barplot(x="pos", y="rate", hue="error")

plt.savefig(snakemake.output[0], bbox_inches="tight")
