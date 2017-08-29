import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x, y))


pd.concat([pd.read_table(f) for f in snakemake.input], keys=snakemake.input)
