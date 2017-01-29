import pandas as pd

cdf = pd.read_table(snakemake.input.cdf, index_col=0)
scales = pd.read_table(snakemake.input.scales, index_col=0, squeeze=True, header=None)
cdf["expr"] *= scales[int(snakemake.wildcards.experiment)]
cdf.to_csv(snakemake.output[0], sep="\t")
