from itertools import combinations
import math

import matplotlib as mpl
mpl.use("svg")
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def hamming_distance(a, b):
    return sum(x == y for x, y in zip(a, b))

# read data
codebook = pd.read_table(snakemake.input.codebook, index_col=0, dtype=np.dtype(str))["codeword"]

cell = int(snakemake.wildcards.cell)

counts = pd.read_table(snakemake.input.counts, index_col=[0, 1]).loc[cell]
counts = counts["corrected"] + counts["exact"]

exprs = pd.read_table(snakemake.input.exprs, index_col=[0, 1])["expr_map"].loc[cell]

# generate graph
G = nx.Graph()
G.add_nodes_from(codebook.index)
for (a, codeword_a), (b, codeword_b) in combinations(codebook.iteritems(), 2):
    if hamming_distance(codeword_a, codeword_b):
        G.add_edge(a, b)

# layout graph
pos = nx.spring_layout(G, k=4/math.sqrt(140))

# draw graph
nx.draw_networkx(G, pos=pos, node_color=counts[G.nodes()], cmap="Reds", with_labels=False, node_size=10)
plt.axis('off')
plt.savefig(snakemake.output.counts, bbox_inches="tight")

nx.draw_networkx(G, pos=pos, node_color=counts[G.nodes()], cmap="Reds", with_labels=False, node_size=10)
plt.axis('off')
plt.savefig(snakemake.output.exprs, bbox_inches="tight")
