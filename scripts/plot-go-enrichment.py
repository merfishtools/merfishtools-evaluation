import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import networkx as nx
from networkx.algorithms import bipartite
import pandas as pd

terms = pd.read_table(snakemake.input.terms, index_col=0)
genes = pd.read_table(snakemake.input.genes)
genes.index = genes["goterm"]

terms = terms[terms["adjPvalue"] <= 0.05]
G = nx.DiGraph()
for goterm in terms.index:
    for _, (goterm, gene, cv) in genes.loc[goterm].iterrows():
        G.add_node(goterm, bipartite=0)
        G.add_node(gene, cv=cv, bipartite=1)
        G.add_edge(gene, goterm)

genes = genes.loc[terms.index]

pos = dict()
pos.update((n, (1, i)) for i, n in enumerate(genes["gene"].unique()))
pos.update((n, (2, i)) for i, n in enumerate(genes["goterm"].unique()))

plt.figure(figsize=(4,4))
ax = nx.draw_networkx(G, pos, arrows=False, cmap=plt.cm.viridis)
plt.colorbar(ax)
plt.axis('off')
plt.savefig(snakemake.output[0], bbox_inches="tight")
