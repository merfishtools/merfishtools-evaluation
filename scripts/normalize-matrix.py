import pandas as pd

expr = pd.read_table(snakemake.input.expr, index_col=0)
scales = pd.read_table(snakemake.input.scales, index_col=0, squeeze=True, header=None)
expr *= scales[int(snakemake.wildcards.experiment)]
expr.to_csv(snakemake.output[0], sep="\t")
