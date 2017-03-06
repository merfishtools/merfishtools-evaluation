import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def hamming_distance(a, b):
    return sum(x == y for x, y in zip(a, b))

# read data
codebook = pd.read_table(snakemake.input.codebook, index_col=0, dtype=np.dtype(str))["codeword"]

counts = pd.read_table(snakemake.input.counts, index_col=[0, 1])
counts = counts["corrected"] + counts["exact"]

exprs = pd.read_table(snakemake.input.exprs, index_col=[0, 1])["expr_map"]

# generate graph
G = nx.Graph()
G.add_nodes_from(codebook.index)
for (a, codeword_a), (b, codeword_b) in combinations(codebook.iteritems(), 2):
    if hamming_distance(codeword_a, codeword_b):
        G.add_edge(a, b)

# layout graph
pos = nx.spring_layout()

# draw graph
nx.draw_networkx(G, pos=pos, node_color=counts[G.nodes()], cmap="Reds")
plt.savefig(snakemake.output.counts)

nx.draw_networkx(G, pos=pos, node_color=counts[G.nodes()], cmap="Reds")
plt.savefig(snakemake.output.exprs)
