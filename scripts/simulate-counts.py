from collections import defaultdict
import csv

import pandas as pd
import numpy as np


# Fix the seed for a particular mean, so that results are reproducible.
np.random.seed(int(snakemake.wildcards.mean))

codebook_mhd4 = pd.read_table(snakemake.input.mhd4, index_col=0, dtype=np.dtype(str))["codeword"]
codebook_mhd2 = pd.read_table(snakemake.input.mhd2, index_col=0, dtype=np.dtype(str))["codeword"]
genes = set(codebook_mhd2.index) | set(codebook_mhd4.index)

with open(snakemake.output[0], "w") as known_out:
    known_out = csv.writer(known_out, delimiter="\t")
    known_out.writerow(["cell", "feat", "mhd2", "mhd4", "count"])
    for cell in range(snakemake.params.cell_count):
        random_counts = np.random.poisson(int(snakemake.wildcards.mean), len(genes))
        for gene, count in zip(genes, random_counts):
            if gene == "PTPN14":  # PTPN14 occurs in all codebooks.
                count = np.random.poisson(1500)
            known_out.writerow([cell, gene, gene in codebook_mhd2.index, gene in codebook_mhd4.index, count])
