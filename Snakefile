"""
Analysis of MERFISHtools, using the data published
at http://zhuang.harvard.edu/merfish.
"""


configfile: "config.yaml"

merfishtools = "merfishtools"

contexts = ["paper"]
datasets = ["140genesData", "1001genesData"]
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
        "figures/fig_simulation_supp.pdf",
        "figures/fig_simulation_ci_error.pdf",
        "figures/fig_exact_vs_corrected.pdf",
        expand("figures/fig_{dataset}.{type}.clustering.pdf", dataset=datasets, type=types),
        "results/paper/140genesData.default.go_enrichment.pdf",
        expand(["figures/fig_{dataset}.multidiffexp.pdf",
                "figures/fig_cv_raw_vs_posterior.pdf",
                "results/{dataset}.default.go_enrichment.terms.txt"], dataset=datasets),
        expand("results/{context}/{dataset}.{type}.default.qqplot.pdf", context="paper", dataset=datasets, type=types),
        expand("results/{context}/simulation-MHD2-{m}/MHD2-{m}.error.default.svg", m=[4,6,8], context="paper"),
        expand("results/{context}/simulation-MHD4-{m}/MHD4-{m}.error.default.svg", m=[4,6,8], context="paper"),
        #expand("results/{context}/simulation-MHD{dist}.rmse.default.svg", dist=2, context="paper"),
        "figures/fig_model.pdf"


#### handling raw data ####


rule format:
    input:
        "data/{dataset}.csv.bz2"
    output:
        "data/{dataset}.{experiment}.{group}.txt"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/format-dataset.py"


rule generate_mhd2_codebook:
    input:
        template="codebook/simulated-MHD2.txt"
    output:
        "codebook/{dataset}.txt"
    wildcard_constraints:
        dataset="simulated-MHD2-[0-9]"
    params:
        ds=lambda wildcards: config["datasets"][wildcards.dataset]
    conda:
        "envs/analysis.yml"
    shell:
        "cut -f1 {input.template} | tail -n+2 | merfishtools gen-mhd2 "
        "-N {params.ds[N]} -m {params.ds[m]} "
        "> {output}"


rule generate_mhd4_codebook:
    input:
        template="codebook/simulated-MHD2.txt"
    output:
        "codebook/{dataset}.txt"
    wildcard_constraints:
        dataset="simulated-MHD4-[0-9]"
    params:
        ds=lambda wildcards: config["datasets"][wildcards.dataset]
    conda:
        "envs/analysis.yml"
    shell:
        "cut -f1 {input.template} | tail -n+2 | merfishtools gen-mhd4 "
        "-m {params.ds[m]} "
        "> {output}"


rule cell_props:
    input:
        "data/{dataset}.{experiment}.all.txt"
    output:
        "cell_properties/{dataset}.{experiment}.all.txt"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/cell-properties.py"


rule raw_counts:
    input:
        "data/{dataset}.{experiment}.{group}.txt"
    output:
        "counts/{dataset}.{experiment}.{group}.txt"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/raw-counts.py"


rule count_matrix:
    input:
        "counts/{dataset}.{experiment}.{group}.txt"
    output:
        "counts/{dataset}.{experiment}.{group}.matrix.txt"
    conda:
        "envs/analysis.yml"
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
        ds=lambda wildcards: config["datasets"][wildcards.dataset]
    benchmark:
        "bench/exp/{dataset}.{settings}.txt"
    conda:
        "envs/analysis.yml"
    threads: 8
    shell:
        "{merfishtools} exp {input.codebook} --p0 {params.ds[err01]} "
        "--p1 {params.ds[err10]} -N {params.ds[N]} -m {params.ds[m]} "
        "--dist {params.ds[dist]} "
        "--estimate {output.est} -t {threads} "
        "< {input.data} > {output.pmf}"


rule expression_matrix:
    input:
        "expressions/{dataset}.{experiment}.{group}.{settings}.est.txt"
    output:
        "expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt"
    conda:
        "envs/analysis.yml"
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
    conda:
        "envs/analysis.yml"
    script:
        "scripts/normalize.py"


rule normalize_expression_matrix:
    input:
        expr="expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt",
        scales="normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    output:
        "normalized_expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt"
    conda:
        "envs/analysis.yml"
    run:
        import pandas as pd
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
    conda:
        "envs/analysis.yml"
    run:
        import pandas as pd
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
        cdf="diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.txt",
        est="diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.est.txt"
    benchmark:
        "bench/diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.txt"
    conda:
        "envs/analysis.yml"
    threads: 8
    shell:
        "{merfishtools} diffexp -t {threads} --pseudocounts 1 "
        "--max-null-log2fc 1.0 --cdf {output.cdf} {input} "
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
        cdf="multidiffexp/{dataset}.{settings}.txt",
        est="multidiffexp/{dataset}.{settings}.est.txt"
    benchmark:
        "bench/multidiffexp/{dataset}.{settings}.txt"
    conda:
        "envs/analysis.yml"
    threads: 24
    shell:
        "{merfishtools} multidiffexp --pseudocounts 1 "
        "-t {threads} --max-null-cv 0.5 "
        "--cdf {output.cdf} {input} > {output.est}"


rule enrichment:
    input:
        "multidiffexp/{dataset}.{settings}.est.txt"
    output:
        terms="results/{dataset}.{settings}.go_enrichment.terms.txt",
        genes="results/{dataset}.{settings}.go_enrichment.genes.txt"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/go-enrichment.R"


#### simulation ####


CELL_COUNT = 100


rule simulate_counts:
    input:
        mhd4=config["codebooks"]["simulated-MHD4"],
        mhd2=config["codebooks"]["simulated-MHD2"]
    output:
        "data/simulated.{mean}.known.txt"
    params:
        cell_count=CELL_COUNT
    conda:
        "envs/analysis.yml"
    script:
        "scripts/simulate-counts.py"


rule simulate:
    input:
        known_counts="data/simulated.{mean}.known.txt",
        codebook=lambda wildcards: config["codebooks"][wildcards.dataset]
    output:
        sim_counts="data/{dataset}.{mean}.all.txt",
        stats="data/{dataset}.{mean}.stats.txt"
    wildcard_constraints:
        dataset="simulated.+"
    params:
        cell_count=CELL_COUNT,
        ds=lambda wildcards: config["datasets"][wildcards.dataset]
    conda:
        "envs/analysis.yml"
    script:
        "scripts/simulate-dataset.py"


#### plots ####


rule plot_cv_raw_vs_posterior:
    input:
        raw_counts=lambda wildcards: expand("counts/{dataset}.{expmnt}.all.txt", dataset=wildcards.dataset, expmnt=experiments(wildcards.dataset)),
        diffexp="multidiffexp/{dataset}.{settings}.est.txt"
    output:
        "results/{context}/{dataset}.{settings}.cv_raw_vs_posterior.{estimate}.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-cv-raw-vs-posterior.py"


rule plot_enrichment:
    input:
        terms="results/{dataset}.{settings}.go_enrichment.terms.txt",
        genes="results/{dataset}.{settings}.go_enrichment.genes.txt"
    output:
        "results/{context}/{dataset}.{settings}.go_enrichment.pdf"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-go-enrichment.py"



rule plot_multidiffexp:
    input:
       diffexp="multidiffexp/{dataset}.{settings}.est.txt",
       exprs=partial(multidiffexp_input, matrix=True),
       diffexp_cdf="multidiffexp/{dataset}.{settings}.txt"
    output:
       "results/{context}/{dataset}.{settings}.diffexp.pdf"
    params:
       expmnts=lambda wildcards: experiments(wildcards.dataset)
    conda:
      "envs/analysis.yml"
    script:
       "scripts/plot-multidiffexp.py"



rule plot_expression_pmf:
    input:
        expr="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        expr_est="expressions/{dataset}.{experiment}.{group}.{settings}.est.txt",
        raw_counts="counts/{dataset}.{experiment}.{group}.txt"
    output:
        "results/{context}/expression_pmf/{dataset}.{experiment}.{group}.{gene}.{settings}.expression_pmf.{legend,(legend|nolegend)}.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-expression-pmf.py"


rule plot_foldchange_cdf:
    input:
        fc="diffexp/{dataset}.{experiment}.{group1}-vs-{experiment}.{group2}.{settings}.txt",
        fc_est="diffexp/{dataset}.{experiment}.{group1}-vs-{experiment}.{group2}.{settings}.est.txt"
    output:
        "results/{context}/foldchange_cdf/{dataset}.{experiment}.{group1}-vs-{group2}.{gene}.{settings}.foldchange_cdf.{legend,(legend|nolegend)}.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-foldchange-cdf.py"


rule plot_qq:
    input:
        lambda wildcards: matrices(wildcards.dataset, type=wildcards.type, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{type}.{settings}.qqplot.svg"
    params:
        experiments=lambda wildcards: experiments(wildcards.dataset)
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-qq.py"


rule plot_expression_dist:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.expression_dist.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-expression-dist.py"


rule plot_overdispersion:
    input:
        lambda wildcards: matrices(wildcards.dataset, type="normalized_expressions", settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.overdispersion.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-overdispersion.py"


rule plot_correlation:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.correlation.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-correlation.py"


rule plot_tsne:
    input:
        exprs=lambda wildcards: matrices(wildcards.dataset,
                                         type=wildcards.type,
                                         settings=wildcards.settings),
        cellprops=lambda wildcards: expand("cell_properties/{dataset}.{experiment}.all.txt",
                                           dataset=wildcards.dataset,
                                           experiment=experiments(wildcards.dataset))
    output:
        "results/{context}/{dataset}.{type}.{settings}.{highlight,(expmnt|codebook|cellsize|cellpos)}.tsne.svg"
    params:
        codebooks=lambda wildcards: [config["codebooks"][wildcards.dataset][expmnt] for expmnt in experiments(wildcards.dataset)]
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-tsne.py"


means = list(range(5, 40, 5))


rule plot_simulation:
    input:
        posterior_counts=expand("expressions/simulated-MHD{{dist}}.{mean}.all.{{settings}}.est.txt", mean=means),
        raw_counts=expand("counts/simulated-MHD{{dist}}.{mean}.all.txt", mean=means),
        known_counts=expand("data/simulated.{mean}.known.txt", mean=means),
        codebook="codebook/simulated-MHD{dist}.txt"
    output:
        violin="results/{context}/simulation-MHD{dist}/MHD{dist}.error.{settings}.svg",
        scatter_raw="results/{context}/simulation-MHD{dist}/MHD{dist}.scatter-raw.{settings}.svg",
        scatter_posterior="results/{context}/simulation-MHD{dist}/MHD{dist}.scatter-posterior.{settings}.svg",
        ci_errors="results/{context}/simulation-MHD{dist}/MHD{dist}.ci-errors.{settings}.svg",
        errors="results/{context}/simulation-MHD{dist}/MHD{dist}.errors.{settings}.txt"
    params:
        means=means
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-simulation.py"


rule plot_m_vs_errors:
    input:
        expand("results/{{context}}/simulation-MHD{{dist}}-{m}/MHD{{dist}}-{m}.errors.{{settings}}.txt", m=[4, 6, 8])
    output:
        "results/{context}/simulation-MHD{dist}.m-vs-errors.{settings}.svg"
    params:
        ms=[4, 6, 8]
    script:
        "scripts/plot-m-vs-errors.py"


rule plot_simulation_pmf:
    input:
        pmfs=expand("expressions/simulated-MHD{{dist}}.{mean}.all.{{settings}}.txt", mean=means),
        known_counts=expand("data/simulated.{mean}.known.txt", mean=means)
    output:
        "results/{context}/simulation-MHD{dist}/MHD{dist}.pmf.{settings}.svg"
    params:
        means=means
    conda:
        "envs/analysis.yml"
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
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-dataset-correlation.py"


#### figures ####


def load_svg(path):
    import svgutils.transform as sg
    return sg.fromfile(path).getroot()


def label_plot(x, y, label):
    import svgutils.transform as sg
    return sg.TextElement(x, y, label, size=12, weight="bold")


rule figure_example:
    input:
        a="results/paper/expression_pmf/140genesData.1.cell34.FLNC.default.expression_pmf.nolegend.svg",
        b="results/paper/expression_pmf/140genesData.1.cell0.FLNC.default.expression_pmf.nolegend.svg",
        c="results/paper/foldchange_cdf/140genesData.1.cell0-vs-cell34.FLNC.default.foldchange_cdf.nolegend.svg",
        d="results/paper/expression_pmf/140genesData.1.cell34.PRKCA.default.expression_pmf.legend.svg",
        e="results/paper/expression_pmf/140genesData.1.cell0.PRKCA.default.expression_pmf.nolegend.svg",
        f="results/paper/foldchange_cdf/140genesData.1.cell0-vs-cell34.PRKCA.default.foldchange_cdf.nolegend.svg"
    output:
        "figures/fig_example.svg"
    conda:
        "envs/analysis.yml"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("6.4in", "3.6in")
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

        la = label_plot(5,10, "a")
        lb = label_plot(195,10, "b")
        lc = label_plot(385,10, "c")
        ld = label_plot(5,170, "d")
        le = label_plot(195,170, "e")
        lf = label_plot(385,170, "f")

        fig.append([a, b, c, d, e, f, la, lb, lc, ld, le, lf])
        fig.save(output[0])


rule figure_simulation:
    input:
        d="results/paper/simulation-MHD4/MHD4.scatter-raw.default.svg",
        f="results/paper/simulation-MHD2-8/MHD2-8.scatter-raw.default.svg",
        c="results/paper/simulation-MHD4/MHD4.scatter-posterior.default.svg",
        e="results/paper/simulation-MHD2-8/MHD2-8.scatter-posterior.default.svg",
        a="results/paper/simulation-MHD4/MHD4.error.default.svg",
        b="results/paper/simulation-MHD2-8/MHD2-8.error.default.svg",
    output:
        "figures/fig_simulation.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-simulation.py"


rule figure_simulation_supp:
    input:
        d="results/paper/simulation-MHD4-4/MHD4-4.scatter-raw.default.svg",
        f="results/paper/simulation-MHD2/MHD2.scatter-raw.default.svg",
        c="results/paper/simulation-MHD4-4/MHD4-4.scatter-posterior.default.svg",
        e="results/paper/simulation-MHD2/MHD2.scatter-posterior.default.svg",
        a="results/paper/simulation-MHD4-4/MHD4-4.error.default.svg",
        b="results/paper/simulation-MHD2/MHD2.error.default.svg",
    output:
        "figures/fig_simulation_supp.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-simulation.py"


rule figure_simulation_ci_error:
    input:
        a="results/paper/simulation-MHD4/MHD4.ci-errors.default.svg",
        b="results/paper/simulation-MHD2/MHD2.ci-errors.default.svg",
        c="results/paper/simulation-MHD4-4/MHD4-4.ci-errors.default.svg",
        d="results/paper/simulation-MHD2-8/MHD2-8.ci-errors.default.svg"
    output:
        "figures/fig_simulation_ci_error.svg"
    script:
        "scripts/fig-ci-error.py"


rule figure_clustering:
    input:
        a="results/paper/{dataset}.{type}.default.cellsize.tsne.svg",
        b="results/paper/{dataset}.{type}.default.cellpos.tsne.svg",
        c="results/paper/{dataset}.{type}.default.codebook.tsne.svg",
        d="results/paper/{dataset}.{type}.default.expmnt.tsne.svg"
    output:
        "figures/fig_{dataset}.{type}.clustering.svg"
    conda:
        "envs/analysis.yml"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("4.8in", "1.1in")
        a = load_svg(input.a)
        b = load_svg(input.b)
        c = load_svg(input.c)
        d = load_svg(input.d)
        b.moveto(160, 0)
        c.moveto(250, 0)
        d.moveto(340, 0)

        la = label_plot(5,10, "a")
        lb = label_plot(165,10, "b")
        lc = label_plot(255,10, "c")
        ld = label_plot(345,10, "d")

        fig.append([a, b, c, d, la, lb, lc, ld])
        fig.save(output[0])


ruleorder: figure_multidiffexp > convert_svg


rule figure_multidiffexp:
    input:
        "results/paper/{dataset}.default.diffexp.pdf",
    output:
        "figures/fig_{dataset}.multidiffexp.pdf"
    conda:
        "envs/analysis.yml"
    shell:
        # nothing to be done here
        "cp {input} {output}"



def get_cv_raw_vs_posterior_input(dataset):
    return expand("results/paper/{dataset}.default.cv_raw_vs_posterior.{estimate}.svg", dataset=dataset, estimate=["cv_ev", "cv_ci_lower"])


rule figure_cv_raw_vs_posterior:
    input:
        mhd4=get_cv_raw_vs_posterior_input(datasets[0]),
        mhd2=get_cv_raw_vs_posterior_input(datasets[1])
    output:
        "figures/fig_cv_raw_vs_posterior.svg"
    conda:
        "envs/analysis.yml"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("7.5in", "1.9in")
        a = load_svg(input.mhd4[0])
        b = load_svg(input.mhd4[1])
        c = load_svg(input.mhd2[0])
        d = load_svg(input.mhd2[1])
        b.moveto(170, 0)
        c.moveto(340, 0)
        d.moveto(510, 0)

        la = label_plot(5, 10, "a")
        lb = label_plot(175, 10, "b")
        lc = label_plot(345, 10, "c")
        ld = label_plot(515, 10, "d")

        fig.append([a, b, c, d, la, lb, lc, ld])
        fig.save(output[0])


rule figure_exact_vs_corrected:
    input:
        expand("counts/140genesData.{experiment}.all.txt", experiment=experiments("140genesData"))
    output:
        "figures/fig_exact_vs_corrected.svg"
    script:
        "scripts/plot-exact-vs-corrected.py"

rule figure_model:
    input:
        "figures/sketch-small.svg",
        "figures/events.svg",
        "figures/urn-model.svg",
        "figures/model.svg"
    output:
        "figures/fig_model.svg"
    conda:
        "envs/analysis.yml"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("7.7in", "3.6in")
        a = load_svg(input[0])
        b = load_svg(input[1])
        c = load_svg(input[2])
        d = load_svg(input[3])
        a.moveto(10, 10, scale=0.6)
        b.moveto(20, 130, scale=0.36)
        c.moveto(40, 190, scale=0.6)
        d.moveto(230, 10, scale=0.8)


        la = label_plot(5, 10, "a")
        lb = label_plot(5, 130, "b")
        lc = label_plot(5, 190, "c")
        ld = label_plot(215, 10, "d")

        fig.append([a, b, c, d, la, lb, lc, ld])
        fig.save(output[0])


#### utils ####


rule convert_svg:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.{fmt,(pdf|png)}"
    conda:
        "envs/analysis.yml"
    shell:
        "cairosvg -f {wildcards.fmt} {input} -o {output}"
