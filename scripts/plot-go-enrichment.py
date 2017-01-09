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

terms = terms[terms["Pvalue"] <= 0.01]

genes = genes.loc[terms.index]
genes = genes.set_index(["goterm", "gene"])
#matrix = ~genes["cv"].unstack().isnull()
matrix = genes["cv"].unstack()
matrix = matrix.loc[terms.index]
matrix.sort_index(axis=1, inplace=True)
print(matrix)

annot = genes["fdr"].unstack()
annot = annot.loc[terms.index]
annot.sort_index(axis=1, inplace=True)
print(annot)
significant = annot <= 0.05
annot[significant] = "•"
annot = annot.fillna(" ")
print(annot)
annot[~significant] = " "

plt.figure(figsize=(3,2))
sns.heatmap(matrix, cmap=plt.cm.Reds, annot=annot, fmt="s", cbar_kws={"label": "CV"})
plt.ylabel("")
plt.xlabel("")
plt.savefig(snakemake.output[0], bbox_inches="tight")
