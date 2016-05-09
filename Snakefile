from itertools import combinations
import pandas as pd
sys.path.insert(0, "../")
import merfishtools as mt
import svgutils.transform as sg


"""
Analysis of MERFISHtools, using the data published
at http://zhuang.harvard.edu/merfish.
"""


configfile: "config.yaml"


merfishtools = "../target/release/merfishtools"
#merfishtools = "merfishtools"


contexts = ["paper"]
datasets = ["140genesData"]
types = ["expressions", "normalized_expressions"]


def experiments(dataset):
    return range(1, config["datasets"][dataset]["experiments"] + 1)


def matrices(dataset, type="expressions", settings="default"):
    suffix = "." if type == "counts" else ".{}.".format(settings)
    return expand("{type}/{dataset}.{experiment}.all{suffix}matrix.txt",
                  dataset=dataset,
                  type=type,
                  suffix=suffix,
                  experiment=experiments(dataset))


rule all:
    input:
        "figures/fig_example.pdf",
        "figures/fig_simulation.pdf",
        expand("figures/fig_{dataset}.{type}.clustering.pdf", dataset=datasets, type=types),
        expand("figures/fig_{dataset}.multidiffexp.pdf", dataset=datasets),
        expand("results/{context}/{dataset}.{type}.default.qqplot.pdf", context="paper", dataset=datasets, type=types),
        "figures/model.pdf", "figures/sketch.pdf"


#### handling raw data ####


rule format:
    input:
        "data/{dataset}.csv.bz2"
    output:
        "data/{dataset}.{experiment}.{group}.txt"
    script:
        "scripts/format-dataset.py"


rule cell_props:
    input:
        "data/{dataset}.{experiment}.all.txt"
    output:
        "cell_properties/{dataset}.{experiment}.all.txt"
    script:
        "scripts/cell-properties.py"


rule raw_counts:
    input:
        "data/{dataset}.{experiment}.{group}.txt"
    output:
        "counts/{dataset}.{experiment}.{group}.txt"
    script:
        "scripts/raw-counts.py"


rule count_matrix:
    input:
        "counts/{dataset}.{experiment}.{group}.txt"
    output:
        "counts/{dataset}.{experiment}.{group}.matrix.txt"
    script:
        "scripts/count-matrix.py"


#### estimating gene expressions ####


def get_codebook(wildcards):
    codebooks = config["codebooks"][wildcards.dataset]
    if isinstance(codebooks, str):
        return codebooks
    return codebooks[int(wildcards.experiment)]


rule expressions:
    input:
        data="data/{dataset}.{experiment}.{group}.txt",
        codebook=get_codebook
    output:
        pmf="expressions/{dataset}.{experiment,[0-9]+}.{group}.{settings}.txt",
        est="expressions/{dataset}.{experiment,[0-9]+}.{group}.{settings}.est.txt",
    params:
        dist=lambda wildcards: config["datasets"][wildcards.dataset]["dist"],
        bits=lambda wildcards: config["datasets"][wildcards.dataset]["N"]
    benchmark:
        "bench/exp/{dataset}.{settings}.txt"
    threads: 8
    shell:
        "{merfishtools} exp --codebook {input.codebook} --hamming-dist {params.dist} -N {params.bits} "
        "--estimate {output.est} -t {threads} < {input.data} > {output.pmf}"


rule expression_matrix:
    input:
        "expressions/{dataset}.{experiment}.{group}.{settings}.est.txt"
    output:
        "expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt"
    script:
        "scripts/expression-matrix.py"


#### normalization ####


rule scale_factors:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    params:
        experiments=lambda wildcards: experiments(wildcards.dataset)
    script:
        "scripts/normalize.py"


rule normalize_expression_matrix:
    input:
        expr="expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt",
        scales="normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    output:
        "normalized_expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt"
    run:
        expr = pd.read_table(input.expr, index_col=0)
        scales = pd.read_table(input.scales, index_col=0, squeeze=True, header=None)
        print(scales)
        expr *= scales[int(wildcards.experiment)]
        expr.to_csv(output[0], sep="\t")


rule normalize_cdf:
    input:
        cdf="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        scales="normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    output:
        "normalized_expressions/{dataset}.{experiment}.{group}.{settings,(default)}.txt"
    run:
        cdf = pd.read_table(input.cdf, index_col=0)
        scales = pd.read_table(input.scales, index_col=0, squeeze=True, header=None)
        cdf["expr"] *= scales[int(wildcards.experiment)]
        cdf.to_csv(output[0], sep="\t")


#### differential expression analysis ####


def diffexp_input(wildcards):
    expr = "expressions" if wildcards.experiment1 == wildcards.experiment2 else "normalized_expressions"
    return ["{expr}/{dataset}.{experiment1}.{group1}.{settings}.txt".format(expr=expr, **wildcards),
            "{expr}/{dataset}.{experiment2}.{group2}.{settings}.txt".format(expr=expr, **wildcards)]


rule diffexp:
    input:
        diffexp_input
    output:
        pmf="diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.txt",
        est="diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.est.txt"
    benchmark:
        "bench/diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.txt"
    threads: 8
    shell:
        "{merfishtools} diffexp -t {threads} --pseudocounts 1 "
        "--max-null-log2fc 1.0 --pmf {output.pmf} {input} "
        "> {output.est}"


def multidiffexp_input(wildcards, matrix=False):
    suffix = ".matrix" if matrix else ""
    return expand("normalized_expressions/{dataset}.{expmnt}.all.{settings}{suffix}.txt",
                  expmnt=experiments(wildcards.dataset),
                  suffix=suffix,
                  **wildcards)


rule multidiffexp:
    input:
        multidiffexp_input
    output:
        pmf="multidiffexp/{dataset}.{settings}.txt",
        est="multidiffexp/{dataset}.{settings}.est.txt"
    benchmark:
        "bench/multidiffexp/{dataset}.{settings}.txt"
    threads: 24
    shell:
        "{merfishtools} --debug multidiffexp --pseudocounts 1 "
        "-t {threads} --max-null-cv 0.5 "
        "--pmf {output.pmf} {input} > {output.est}"


rule enrichment:
    input:
        "multidiffexp/{dataset}.{settings}.est.txt"
    output:
        table="results/paper/{dataset}.{settings}.go_enrichment.txt"
    script:
        "scripts/go-enrichment.R"


#### simulation ####


rule simulate:
    input:
        mhd4=config["codebooks"]["simulated-MHD4"],
        mhd2=config["codebooks"]["simulated-MHD2"]
    output:
        sim_counts_mhd4="data/simulated-MHD4.{mean}.all.txt",
        sim_counts_mhd2="data/simulated-MHD2.{mean}.all.txt",
        known_counts="data/simulated.{mean}.known.txt",
        stats_mhd4="data/simulated-MHD4.{mean}.stats.txt",
        stats_mhd2="data/simulated-MHD2.{mean}.stats.txt"
    params:
        cell_count=100
    script:
        "scripts/simulate-counts.py"


#### plots ####


rule plot_multidiffexp:
    input:
       diffexp="multidiffexp/{dataset}.{settings}.est.txt",
       exprs=partial(multidiffexp_input, matrix=True),
       diffexp_cdf="multidiffexp/{dataset}.{settings}.txt"
    output:
       "results/{context}/{dataset}.{settings}.diffexp.pdf"
    params:
       expmnts=lambda wildcards: experiments(wildcards.dataset)
    script:
       "scripts/plot-multidiffexp.py"


rule plot_enrichment:
    input:
        "results/paper/{dataset}.{settings}.go_enrichment.txt"
    output:
        "results/{context}/{dataset}.{settings}.go_enrichment.svg"
    script:
        "scripts/plot-go-enrichment.py"


rule plot_expression_pmf:
    input:
        expr="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        expr_est="expressions/{dataset}.{experiment}.{group}.{settings}.est.txt",
        raw_counts="counts/{dataset}.{experiment}.{group}.txt"
    output:
        "results/{context}/expression_pmf/{dataset}.{experiment}.{group}.{gene}.{settings}.expression_pmf.{legend,(legend|nolegend)}.svg"
    script:
        "scripts/plot-expression-pmf.py"


rule plot_foldchange_cdf:
    input:
        fc="diffexp/{dataset}.{experiment}.{group1}-vs-{experiment}.{group2}.{settings}.txt",
        fc_est="diffexp/{dataset}.{experiment}.{group1}-vs-{experiment}.{group2}.{settings}.est.txt"
    output:
        "results/{context}/foldchange_cdf/{dataset}.{experiment}.{group1}-vs-{group2}.{gene}.{settings}.foldchange_cdf.{legend,(legend|nolegend)}.svg"
    script:
        "scripts/plot-foldchange-cdf.py"


rule plot_qq:
    input:
        lambda wildcards: matrices(wildcards.dataset, type=wildcards.type, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{type}.{settings}.qqplot.svg"
    params:
        experiments=lambda wildcards: experiments(wildcards.dataset)
    script:
        "scripts/plot-qq.py"


rule plot_expression_dist:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.expression_dist.svg"
    script:
        "scripts/plot-expression-dist.py"


rule plot_overdispersion:
    input:
        lambda wildcards: matrices(wildcards.dataset, type="normalized_expressions", settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.overdispersion.svg"
    script:
        "scripts/plot-overdispersion.py"


rule plot_correlation:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.correlation.svg"
    script:
        "scripts/plot-correlation.py"


rule plot_tsne:
    input:
        exprs=lambda wildcards: matrices(wildcards.dataset,
                                         type=wildcards.type,
                                         settings=wildcards.settings),
        cellprops=lambda wildcards: expand("cell_properties/{dataset}.{experiment}.all.txt",
                                           dataset=wildcards.dataset,
                                           experiment=experiments(wildcards.dataset)),
    output:
        "results/{context}/{dataset}.{type}.{settings}.{highlight,(expmnt|codebook|cellsize|cellpos)}.tsne.svg"
    params:
        codebooks=lambda wildcards: [config["codebooks"][wildcards.dataset][expmnt] for expmnt in experiments(wildcards.dataset)]
    script:
        "scripts/plot-tsne.py"


means = list(range(5, 40, 5))


rule plot_simulation:
    input:
        posterior_counts=expand("expressions/simulated-MHD{{dist}}.{mean}.all.{{settings}}.est.txt", mean=means),
        raw_counts=expand("counts/simulated-MHD{{dist}}.{mean}.all.txt", mean=means),
        known_counts=expand("data/simulated.{mean}.known.txt", mean=means)
    output:
        violin="results/{context}/simulation-MHD{dist}/MHD{dist}.error.{settings}.svg",
        scatter="results/{context}/simulation-MHD{dist}/MHD{dist}.scatter.{settings}.svg"
    params:
        means=means
    script:
        "scripts/plot-simulation.py"


rule plot_simulation_pmf:
    input:
        pmfs=expand("expressions/simulated-MHD{{dist}}.{mean}.all.{{settings}}.txt", mean=means),
        known_counts=expand("data/simulated.{mean}.known.txt", mean=means)
    output:
        "results/{context}/simulation-MHD{dist}/MHD{dist}.pmf.{settings}.svg"
    params:
        means=means
    script:
        "scripts/plot-simulation-pmf.py"


rule plot_dataset_correlation:
    input:
        small=matrices("140genesData"),
        large=matrices("1001genesData"),
        small_counts=matrices("140genesData", type="counts"),
        large_counts=matrices("1001genesData", type="counts")
    output:
        "results/{context}/{settings}.dataset_correlation.svg"
    script:
        "scripts/plot-dataset-correlation.py"


comparisons = list(combinations(experiments("140genesData"), 2))


rule plot_go_term_enrichment:
    input:
        expand("diffexp/140genesData.{experiment[0]}.all-vs-{experiment[1]}.all.default.est.txt",
               experiment=comparisons)
    output:
        "results/{context}/comparisons/140genesData.comparisons.svg",
    params:
        comparisons=comparisons
    script:
        "scripts/experiment-diffexp.py"


#### figures ####


load_svg = lambda path: return sg.fromfile(path).getroot()
label_plot = lambda x, y, label: sg.TextElement(x, y, label, size=12, weight="bold")


rule figure_example:
    input:
        a="results/paper/expression_pmf/140genesData.1.cell34.FLNC.default.expression_pmf.nolegend.svg",
        b="results/paper/expression_pmf/140genesData.1.cell0.FLNC.default.expression_pmf.nolegend.svg",
        c="results/paper/foldchange_cdf/140genesData.1.cell0-vs-cell34.FLNC.default.foldchange_cdf.nolegend.svg",
        d="results/paper/expression_pmf/140genesData.1.cell34.PRKCA.default.expression_pmf.nolegend.svg",
        e="results/paper/expression_pmf/140genesData.1.cell0.PRKCA.default.expression_pmf.nolegend.svg",
        f="results/paper/foldchange_cdf/140genesData.1.cell0-vs-cell34.PRKCA.default.foldchange_cdf.legend.svg"
    output:
        "figures/fig_example_{gene}.svg"
    run:
        fig = sg.SVGFigure("6.4in", "1.8in")
        a = load_svg(input.a)
        b = load_svg(input.b)
        c = load_svg(input.c)
        d = load_svg(input.d)
        e = load_svg(input.e)
        f = load_svg(input.f)
        b.moveto(190, 0)
        c.moveto(380, 0)
        d.moveto(0, 160)
        e.moveto(190, 160)
        f.moveto(380, 160)

        la = label_plot(0,10, "a")
        lb = label_plot(190,10, "b")
        lc = label_plot(380,10, "c")
        ld = label_plot(0,170, "c")
        le = label_plot(190,170, "c")
        lf = label_plot(380,170, "c")

        fig.append([a, b, c, la, lb, lc])
        fig.save(output[0])


rule figure_simulation:
    input:
        b="results/paper/simulation-MHD4/MHD4.scatter.default.svg",
        d="results/paper/simulation-MHD2/MHD2.scatter.default.svg",
        a="results/paper/simulation-MHD4/MHD4.error.default.svg",
        c="results/paper/simulation-MHD2/MHD2.error.default.svg"
    output:
        "figures/fig_simulation.svg"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("4.6in", "3.6in")
        a = load_svg(input.a)
        b = load_svg(input.b)
        c = load_svg(input.c)
        d = load_svg(input.d)
        b.moveto(250, 0)
        c.moveto(0, 160)
        d.moveto(250, 160)

        la = label_plot(0,10, "a")
        lb = label_plot(250,10, "b")
        lc = label_plot(0,170, "c")
        ld = label_plot(250,170, "d")

        fig.append([a, b, c, d, la, lb, lc, ld])
        fig.save(output[0])


rule figure_clustering:
    input:
        a="results/paper/{dataset}.{type}.default.cellsize.tsne.svg",
        b="results/paper/{dataset}.{type}.default.cellpos.tsne.svg",
        c="results/paper/{dataset}.{type}.default.codebook.tsne.svg",
        d="results/paper/{dataset}.{type}.default.expmnt.tsne.svg"
    output:
        "figures/fig_{dataset}.{type}.clustering.svg"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("23.8cm", "5.3cm")
        a = load_svg(input.a)
        b = load_svg(input.b)
        c = load_svg(input.c)
        d = load_svg(input.d)
        b.moveto(258, 0)
        c.moveto(453, 0)
        d.moveto(650, 0)

        la = label_plot(0,10, "a")
        lb = label_plot(258,10, "b")
        lc = label_plot(453,10, "c")
        ld = label_plot(650,10, "d")

        fig.append([a, b, c, d, la, lb, lc, ld])
        fig.save(output[0])


rule figure_multidiffexp:
    input:
        "results/paper/{dataset}.default.diffexp.pdf",
    output:
        "figures/fig_{dataset}.multidiffexp.pdf"
    shell:
        # nothing to be done here
        "cp {input} {output}"


#### utils ####


rule convert_svg:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.{fmt,(pdf|png)}"
    shell:
        "rsvg-convert -f {wildcards.fmt} {input} > {output}"
