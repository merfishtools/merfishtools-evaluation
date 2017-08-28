from collections import Counter, defaultdict
from itertools import chain
import csv

import pandas as pd
import numpy as np
from bitarray import bitarray


noise_rate = 0.5

# Fix the seed for a particular mean, so that results are reproducible.
np.random.seed(int(snakemake.wildcards.mean))

p0 = snakemake.params.ds["err01"]
p1 = snakemake.params.ds["err10"]


def sim_errors(word):
    p = np.random.uniform(0, 1, len(word))
    err10 = bitarray(list(p <= p1)) & word
    err01 = bitarray(list(p <= p0)) & ~word
    readout = (word ^ err10) ^ err01
    errs = err10.count(True) + err01.count(True)
    assert errs == 0 or word != readout
    return readout, errs


def hamming1_env(word):
    for i in range(len(word)):
        w = word.copy()
        w[i] ^= 1
        yield w


codebook = pd.read_table(snakemake.input.codebook, index_col=0, dtype=np.dtype(str))
codebook["codeword"] = codebook["codeword"].apply(bitarray)
codebook["expressed"] = codebook["expressed"] == "1"
noise_word = bitarray("0" * len(codebook["codeword"].iloc[0]))

known_counts = pd.read_table(snakemake.input.known_counts, index_col=[0, 1])

def simulate(codebook, counts_path, stats_path, has_corrected=True):
    lookup_exact = {word.tobytes(): gene for gene, word, _ in codebook.itertuples()}
    if has_corrected:
        lookup_corrected = {w.tobytes(): gene for gene, word, _ in codebook.itertuples() for w in hamming1_env(word)}
    else:
        lookup_corrected = {}

    with open(counts_path, "w") as sim_out:
        sim_out = csv.writer(sim_out, delimiter="\t")
        sim_out.writerow(["cell", "feat", "dist", "cell_x", "cell_y", "x", "y"])

        stats = []
        for cell in range(snakemake.params.cell_count):
            readouts = []
            words = []
            genes = []
            errors = []

            total_counts = 0
            for gene, word, expressed in codebook.itertuples():
                if not expressed:
                    # Skip entries that are marked as not expressed in the codebook.
                    # These can act as misidentification probes.
                    continue
                count = known_counts.loc[cell, gene]["count"]
                total_counts += count
                for _ in range(count):
                    readout, errs = sim_errors(word)
                    errors.append(errs)
                    readouts.append(readout)
                    words.append(word)
                    genes.append(gene)

            # add noise
            noise_count = int(noise_rate * total_counts)
            print("noise count: ", noise_count)
            for _ in range(noise_count):
                readout, errs = sim_errors(noise_word)
                if has_corrected:
                    if errs < 3:
                        continue
                elif errs < 4:
                    continue
                errors.append(errs)
                readouts.append(readout)
                words.append(noise_word)
                genes.append("noise")

            exact_counts = Counter()
            corrected_counts = Counter()
            exact_miscalls = Counter()
            corrected_miscalls = Counter()
            exact_noise_miscalls = Counter()
            corrected_noise_miscalls = Counter()
            for readout, errs, word, orig_gene in zip(readouts, errors, words, genes):
                try:
                    gene = lookup_exact[readout.tobytes()]
                    exact_counts[gene] += 1
                    if errs > 0:
                        if orig_gene == "noise":
                            exact_noise_miscalls[gene] += 1
                        else:
                            exact_miscalls[gene] += 1
                except KeyError:
                    try:
                        gene = lookup_corrected[readout.tobytes()]
                        corrected_counts[gene] += 1
                        if errs > 1:
                            if orig_gene == "noise":
                                corrected_noise_miscalls[gene] += 1
                            else:
                                corrected_miscalls[gene] += 1
                    except KeyError:
                        # readout is lost
                        pass

            for gene in set(chain(exact_counts, corrected_counts)):
                for _ in range(exact_counts[gene]):
                    sim_out.writerow([cell, gene, 0, 0, 0, 0, 0])
                for _ in range(corrected_counts[gene]):
                    sim_out.writerow([cell, gene, 1, 0, 0, 0, 0])

            for gene in exact_counts:
                known = known_counts.loc[cell, gene]["count"]
                counts = exact_counts[gene] + corrected_counts[gene]
                _exact_miscalls = exact_miscalls[gene]
                _corrected_miscalls = corrected_miscalls[gene]
                noise_miscalls = exact_noise_miscalls[gene] + corrected_noise_miscalls[gene]
                miscalls = _exact_miscalls + _corrected_miscalls + noise_miscalls
                missed = known - counts + miscalls
                stats.append([cell, gene, known, missed, counts, miscalls, noise_miscalls])

        stats = pd.DataFrame(stats)
        stats.columns = ["cell", "gene", "truth", "missed", "counts", "miscalls", "noise-miscalls"]
        stats["total"] = stats["truth"] + stats["miscalls"]
        stats["calls"] = stats["counts"] - stats["miscalls"]
        stats["call-rate"] = stats["calls"] / stats["total"]
        stats["miscall-rate"] = stats["calls"] / stats["total"]
        stats["missed-rate"] = stats["missed"] / stats["total"]
        stats.to_csv(stats_path, sep="\t")

        print("rates vs counts")
        print("miscalls", (stats["miscalls"] / stats["counts"]).mean())
        print("calls", (stats["calls"] / stats["counts"]).mean())
        print("missed", (stats["missed"] / stats["counts"]).mean())
        print("rates vs truth")
        print("miscalls", (stats["miscalls"] / stats["truth"]).mean())
        print("calls", (stats["calls"] / stats["truth"]).mean())
        print("missed", (stats["missed"] / stats["truth"]).mean())
        print("rates vs total")
        print("miscalls", (stats["miscalls"] / stats["total"]).mean())
        print("calls", (stats["calls"] / stats["total"]).mean())
        print("missed", (stats["missed"] / stats["total"]).mean())



print("Simulating dataset")
simulate(codebook, snakemake.output.sim_counts, snakemake.output.stats, has_corrected=snakemake.params.ds["has_corrected"])
