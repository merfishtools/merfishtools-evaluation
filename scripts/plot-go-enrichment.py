import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import networkx as nx
from networkx.algorithms import bipartite
import pandas as pd
import seaborn as sns

terms = pd.read_table(snakemake.input.terms, index_col=0)
genes = pd.read_table(snakemake.input.genes)
genes.index = genes["goterm"]

terms = terms[terms["adjPvalue"] <= 0.05]

genes = genes.loc[terms.index]
genes = genes.set_index(["goterm", "gene"])
matrix = ~genes["cv"].unstack().isnull()
matrix = matrix.loc[terms.index]
matrix.sort_index(axis=1, inplace=True)
plt.figure(figsize=(2,2))
sns.heatmap(matrix, cmap=plt.cm.Greys, cbar=False)
plt.ylabel("")
plt.xlabel("")
plt.savefig(snakemake.output[0], bbox_inches="tight")
