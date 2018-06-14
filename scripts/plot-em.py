import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("svg")
import matplotlib.pyplot as plt
import seaborn as sns

import re

pattern = re.compile(r"DEBUG: TYPE=EM-iteration, CELL=(?P<cell>\d+), x=(?P<expr>\[(\d+, )+\d+\])")


exprs = []

with open(snakemake.input[0]) as f:
    for l in f:
        m = pattern.match(l)
        if re.match(r"DEBUG: TYPE=EM-iteration, CELL=(?P<cell>\d+)", l):
            print("match")
        if m and m.group("cell") == "1":
            exprs.append(eval(m.group("expr")))

exprs = pd.DataFrame(exprs)


sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
x, y = snakemake.config["plots"]["figsize"]
plt.figure(figsize=(x * 1.5, y))

for feat, e in exprs.iteritems():
    curr = e[1:].reset_index(drop=True)
    last = e[:-1].reset_index(drop=True)
    #e = ((curr + 1) / (last + 1))
    plt.semilogy(np.arange(e.size), e, "k-", linewidth=1, alpha=0.5)

sns.despine()

plt.xlabel("EM iteration")
plt.ylabel("expression")
#plt.ylim((0, 200))
plt.savefig(snakemake.output[0], bbox_inches="tight")
