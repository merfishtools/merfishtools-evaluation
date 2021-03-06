import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing, manifold, decomposition
from itertools import combinations
from matplotlib import gridspec

experiments = np.array(list(range(1, len(snakemake.input) + 1)))
# load data
exprs = [pd.read_table(f, index_col=0) for f in snakemake.input.exprs]
exprs = pd.concat(exprs, axis="columns", keys=experiments, names=["expmnt", "cell"])
exprs = np.log10(1 + exprs.transpose())
exprs.index = exprs.index.set_levels(exprs.index.levels[1].astype(np.int64), level=1)

cellprops = [pd.read_table(f, index_col=0) for f in snakemake.input.cellprops]
cellprops = pd.concat(cellprops, keys=experiments, names=["expmnt", "cell"])
cellprops.columns = ["area", "cell_x", "cell_y"]
# reduce cell position dimensionality to 1 dimension
#pca = decomposition.PCA(n_components=1)
tsne = manifold.TSNE(n_components=1, random_state=2390)
pos = tsne.fit_transform(cellprops[["cell_x", "cell_y"]])[:, 0]
cellprops["pos"] = pos

# calculate t-SNE embedding
tsne = manifold.TSNE(random_state=21498)
X = preprocessing.scale(exprs)
embedding = pd.DataFrame(tsne.fit_transform(X), columns=["x", "y"])
embedding.index = exprs.index

# plot
sns.set(style="ticks", palette="colorblind", context=snakemake.wildcards.context)
width, height = snakemake.config["plots"]["figsize"]
fig = plt.figure(figsize=[0.75 * snakemake.config["plots"]["figsize"][0]] * 2)
plt.subplot(111, aspect="equal")

for expmnt, codebook in zip(experiments, snakemake.params.codebooks):
    embedding.loc[expmnt, "codebook"] = codebook
embedding["codebook"] = embedding["codebook"].astype("category")
embedding = pd.concat([embedding, cellprops], axis="columns")
embedding.reset_index(inplace=True)

if snakemake.wildcards.highlight == "expmnt":
    colors = sns.color_palette("Paired", len(experiments))
    highlight = [colors[e] for e in embedding["expmnt"]]
    cmap = None
elif snakemake.wildcards.highlight == "codebook":
    codebooks = embedding["codebook"].cat.categories
    colors = dict(zip(codebooks, sns.color_palette("muted", len(codebooks))))
    highlight = [colors[c] for c in embedding["codebook"]]
    cmap = None
elif snakemake.wildcards.highlight == "cellsize":
    highlight = "area"
    cmap = "viridis"
elif snakemake.wildcards.highlight == "cellpos":
    highlight = "pos"
    cmap = "viridis"
ax = plt.scatter("x", "y", c=highlight, s=5, data=embedding, cmap=cmap, alpha=0.7, edgecolors="face")


plt.axis("off")
from matplotlib.ticker import NullLocator
plt.gca().xaxis.set_major_locator(NullLocator())
plt.gca().yaxis.set_major_locator(NullLocator())

if snakemake.wildcards.highlight == "cellsize":
    cb = plt.colorbar(ax)
    cb.set_label("cell size in nm²")
plt.savefig(snakemake.output[0])
