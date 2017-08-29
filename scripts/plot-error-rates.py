import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x * 1.5, y))


error_rates = pd.concat([pd.read_table(f, index_col=0) for f in snakemake.input],
                        keys=snakemake.params.experiments)

error_rates = pd.DataFrame(error_rates.stack()).reset_index()
error_rates.columns = ["experiment", "pos", "error-rate", "rate"]
error_rates["pos"] += 1

sns.barplot(x="pos", y="rate", hue="error-rate", data=error_rates, errwidth=1, linewidth=0, palette=["red", "grey"], saturation=0.7)
plt.ylim((0.0, plt.ylim()[1] + 0.1))
plt.ylabel("")
plt.xlabel("position")
sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")

print("mean error rates per position")
means = error_rates.groupby(["error-rate", "pos"])["rate"].mean()
print("p0", list(means.loc["p0"].values))
print("p1", list(means.loc["p1"].values))
print("mean overall")
print(error_rates.groupby("error-rate")["rate"].mean())
