"""
Analysis of MERFISHtools, using the data published
at http://zhuang.harvard.edu/merfish.
"""

__author__ = "Johannes Köster"
__license__ = "MIT"


from snakemake.utils import min_version
min_version("3.10.1")
shell.executable("bash")

from copy import deepcopy

configfile: "config.yaml"

contexts = ["paper"]
datasets = ["140genesData", "1001genesData"]
types = ["expressions", "normalized_expressions"]
means = list(range(5, 40, 5))


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
        "figures/fig_dataset_correlation.pdf",
        expand("figures/fig_{dataset}.{type}.clustering.pdf", dataset=datasets, type=types),
        "results/paper/140genesData.default.go_enrichment.pdf",
        expand(["figures/fig_{dataset}.multidiffexp.pdf",
                "figures/fig_cv_raw_vs_posterior.pdf",
                "results/{dataset}.default.go_enrichment.terms.txt"], dataset=datasets),
        expand("results/{context}/{dataset}.{type}.default.qqplot.pdf", context="paper", dataset=datasets, type=types),
        expand("results/{context}/simulation-MHD2-{m}/MHD2-{m}.error.default.svg", m=[4,6,8], context="paper"),
        expand("results/{context}/simulation-MHD4-{m}/MHD4-{m}.error.default.svg", m=[4,6,8], context="paper"),
        #expand("results/{context}/simulation-MHD{dist}.rmse.default.svg", dist=2, context="paper"),
        "figures/fig_model.pdf",
        "figures/fig_error_rate_uncertainty.pdf",
        # "results/paper/140genesData.default.naive-vs-map-scatter.pdf",
        expand("results/{context}/neighborhoods/{dataset}.{experiment}.{feature}.default.{pred}.neighborhood.svg",
               context="paper",
               dataset="140genesData",
               experiment=experiments("140genesData"),
               feature=["THBS1"],
               pred=["raw", "posterior"]),
        expand("results/{context}/neighborhoods/simulated-MHD{dist}.{mean}.{feature}.default.{pred}.neighborhood.svg",
               mean=means,
               context="paper",
               dist=["2", "4", "2-8", "4-4"],
               feature="maxexpr",
               pred=["raw", "posterior"]),
        expand("results/{context}/{dataset}.default.expression_dist.svg",
               context="paper",
               dataset=datasets + ["simulated-MHD4", "simulated-MHD2"]),
        "figures/fig_error_rates.svg",
        expand("results/{context}/em-algorithm/{dataset}.{experiment}.all.default.em-algorithm.svg",
               context="paper",
               dataset="140genesData",
               experiment=experiments("140genesData"))
               

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
        "envs/merfishtools.yml"
    shell:
        "cut -f1 {input.template} | tail -n+2 | merfishtools gen-mhd2 "
        "-N {params.ds[N]} -m {params.ds[m]} "
        "--not-expressed '(notarget|blank).+' "
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
        "envs/merfishtools.yml"
    shell:
        "cut -f1 {input.template} | tail -n+2 | merfishtools gen-mhd4 "
        "-m {params.ds[m]} --not-expressed '(notarget|blank).+' "
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


def get_expressions_params(wildcards):
    ds = deepcopy(config["datasets"][wildcards.dataset])
    settings = config["settings"][wildcards.settings]
    err = config["error-rates"][get_codebook(wildcards)]
    ds["err01"] = [p * settings["err01-factor"] for p in err["err01"]]
    ds["err10"] = [p * settings["err10-factor"] for p in err["err10"]]
    return ds


rule estimate_real_error_rates:
    input:
        readouts="data/{dataset}.readouts.txt",
        codebook=get_codebook
    output:
        "error-rates/{dataset}.{experiment}.error-rates.txt"
    wildcard_constraints:
        dataset="140genesData"
    conda:
        "envs/analysis.yml"
    shell:
        "grep -P '^{wildcards.experiment}\\t' {input.readouts} | cut -f2,3,4 | "
        "merfishtools est-error-rates {input.codebook} > {output}"


rule estimate_simulated_error_rates:
    input:
        readouts="data/{dataset}.{experiment}.readouts.txt",
        codebook=get_codebook
    output:
        "error-rates/{dataset}.{experiment}.error-rates.txt"
    wildcard_constraints:
        dataset="simulated-.+"
    conda:
        "envs/merfishtools.yml"
    shell:
        "merfishtools est-error-rates {input.codebook} "
        "< {input.readouts} > {output}"


rule expressions:
    input:
        data="data/{dataset}.{experiment}.{group}.txt",
        codebook=get_codebook
    output:
        pmf="expressions/{dataset}.{experiment,[0-9]+}.{group}.{settings}.txt",
        est="expressions/{dataset}.{experiment,[0-9]+}.{group}.{settings}.est.txt",
        stats="expressions/{dataset}.{experiment,[0-9]+}.{group}.{settings}.stats.txt"
    params:
        ds=get_expressions_params
    benchmark:
        "bench/exp/{dataset}.{experiment}.{group}.{settings}.txt"
    conda:
        "envs/merfishtools.yml"
    log:
        "logs/exp/{dataset}.{experiment}.{group}.{settings}.exp.log"
    threads: 8
    shell:
        "merfishtools -v exp {input.codebook} --p0 {params.ds[err01]} "
        "--p1 {params.ds[err10]} --seed 9823749 "
        "--estimate {output.est} -t {threads} "
        "--stats {output.stats} "
        "< {input.data} > {output.pmf} 2> {log}"


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
    script:
        "scripts/normalize-matrix.py"


rule normalize_cdf:
    input:
        cdf="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        scales="normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    output:
        "normalized_expressions/{dataset}.{experiment}.{group}.{settings,(default)}.txt"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/normalize-cdf.py"


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
        "envs/merfishtools.yml"
    threads: 8
    shell:
        "merfishtools diffexp -t {threads} --pseudocounts 1 "
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
        "envs/merfishtools.yml"
    threads: 24
    shell:
        "merfishtools multidiffexp --pseudocounts 1 "
        "-t {threads} --max-null-cv 0.5 "
        "--cdf {output.cdf} {input} > {output.est}"


rule enrichment:
    input:
        "multidiffexp/{dataset}.{settings}.est.txt"
    output:
        terms="results/{dataset}.{settings}.go_enrichment.terms.txt",
        genes="results/{dataset}.{settings}.go_enrichment.genes.txt"
    conda:
        "envs/enrichment.yml"
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
        readouts="data/{dataset}.{mean}.readouts.txt",
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

def get_experiments(wildcards):
    if wildcards.dataset.startswith("simulated"):
        return means
    elif wildcards.dataset == "140genesData":
        return [e for e in experiments(wildcards.dataset) 
                if config["codebooks"][wildcards.dataset][e] == 
                   "codebook/140genesData.{}.txt".format(wildcards.codebook)]
    else:
        return experiments(wildcards.dataset)
    


rule plot_error_rates:
    input:
        lambda wildcards: expand("error-rates/{dataset}.{experiment}.error-rates.txt", experiment=get_experiments(wildcards), dataset=wildcards.dataset)
    output:
        "results/{context}/{dataset}.{codebook}.error-rates.svg"
    params:
        experiments=get_experiments
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-error-rates.py"


rule plot_error_rate_uncertainty:
    input:
        codebook="codebook/simulated-MHD{dist}.txt",
        known_counts=expand("data/simulated.{mean}.known.txt", mean=means),
        raw_counts=expand("counts/simulated-MHD{{dist}}.{mean}.all.txt", mean=means),
        **{
            uncertainty: expand("expressions/simulated-MHD{{dist}}.{mean}.all.{settings}.est.txt", mean=means, settings=uncertainty)
            for uncertainty in ["default", "err-5%", "err-10%", "err-20%", "err-30%"]
        }
    output:
        "results/{context}/simulation-MHD{dist}/MHD{dist}.error-rate-uncertainty.svg"
    params:
        means=means
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-error-rate-uncertainty.py"


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


rule plot_naive_vs_map:
    input:
        estimates=lambda wildcards: expand("expressions/{dataset}.{expmnt}.all.{settings}.est.txt", expmnt=experiments(wildcards.dataset), **wildcards),
        counts=lambda wildcards: expand("counts/{dataset}.{expmnt}.all.txt", expmnt=experiments(wildcards.dataset), **wildcards)
    output:
        hist="results/{context}/{dataset}.{settings}.naive-vs-map-hist.svg",
        scatter="results/{context}/{dataset}.{settings}.naive-vs-map-scatter.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-naive-vs-map.py"


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
        lambda wildcards: matrices(wildcards.dataset, type="counts", settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.expression_dist.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-expression-dist.py"


rule plot_overdispersion:
    input:
        lambda wildcards: matrices(wildcards.dataset, type="counts", settings=wildcards.settings)
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
    conda:
        "envs/analysis.yml"
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
        small=lambda w: matrices("140genesData", type=w.type),
        large=lambda w: matrices("1001genesData", type=w.type)
    output:
        "results/{context}/{settings}.dataset_correlation.{type,counts|expressions}.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-dataset-correlation.py"


rule plot_dataset_correlation_error:
    input:
        small=lambda w: matrices("140genesData", type="expressions"),
        large=lambda w: matrices("1001genesData", type="expressions"),
        small_counts=lambda w: matrices("140genesData", type="counts"),
        large_counts=lambda w: matrices("1001genesData", type="counts")
    output:
        "results/{context}/{settings}.dataset_correlation_error.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-dataset-correlation-error.py"


rule plot_neighborhood:
    input:
        codebook=get_codebook,
        counts="counts/{dataset}.{experiment}.all.txt",
        exprs="expressions/{dataset}.{experiment}.all.{settings}.est.txt"
    output:
        strip="results/{context}/neighborhoods/{dataset}.{experiment}.{feature}.{settings}.{pred}.neighborhood.svg",
        kde="results/{context}/neighborhoods/{dataset}.{experiment}.{feature}.{settings}.{pred}.neighborhood.kde.svg",
    conda:
        "envs/analysis.yml"
    wildcard_constraints:
        pred="(raw|posterior)"
    script:
        "scripts/plot-neighborhood.py"


rule plot_em:
    input:
        "logs/exp/{dataset}.{experiment}.{group}.{settings}.exp.log"
    output:
        "results/{context}/em-algorithm/{dataset}.{experiment}.{group}.{settings}.em-algorithm.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-em.py"


#### figures ####


rule figure_thbs1_bias:
    input:
        a="results/paper/neighborhoods/140genesData.2.THBS1.default.raw.neighborhood.svg",
        b="results/paper/neighborhoods/140genesData.2.THBS1.default.posterior.neighborhood.svg",
    output:
        "figures/fig_thbs1_bias.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-thbs1-bias.py"



rule figure_error_rates:
    input:
        a="results/paper/140genesData.1.error-rates.svg",
        b="results/paper/140genesData.2.error-rates.svg"
    output:
        "figures/fig_error_rates.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-error-rates.py"


rule figure_error_rate_uncertainty:
    input:
        a="results/paper/simulation-MHD4/MHD4.error-rate-uncertainty.svg",
        b="results/paper/simulation-MHD2-8/MHD2-8.error-rate-uncertainty.svg"
    output:
        "figures/fig_error_rate_uncertainty.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-error-rate-uncertainty.py"


rule figure_naive_vs_map:
    input:
        a="results/paper/140genesData.default.naive-vs-map-scatter.svg",
        b="results/paper/140genesData.default.naive-vs-map-hist.svg"
    output:
        "figures/figure_naive_vs_map.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-naive-vs-map.py"



rule figure_example:
    input:
        "results/paper/expression_pmf/140genesData.1.cell16.NUP98.default.expression_pmf.nolegend.svg",
        "results/paper/expression_pmf/140genesData.1.cell17.NUP98.default.expression_pmf.nolegend.svg",
        "results/paper/expression_pmf/140genesData.1.cell0.PRKCA.default.expression_pmf.nolegend.svg",
    output:
        "figures/fig_example.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-example.py"


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
    conda:
        "envs/analysis.yml"
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
    script:
        "scripts/fig-clustering.py"


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
    return expand("results/paper/{dataset}.default.cv_raw_vs_posterior.{estimate}.svg", dataset=dataset, estimate=["cv_map", "cv_ci_lower"])


rule figure_cv_raw_vs_posterior:
    input:
        mhd4=get_cv_raw_vs_posterior_input(datasets[0]),
        mhd2=get_cv_raw_vs_posterior_input(datasets[1])
    output:
        "figures/fig_cv_raw_vs_posterior.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-cv-raw-vs-posterior.py"


rule figure_exact_vs_corrected:
    input:
        expand("counts/140genesData.{experiment}.all.txt", experiment=experiments("140genesData"))
    output:
        "figures/fig_exact_vs_corrected.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/plot-exact-vs-corrected.py"

rule figure_model:
    input:
        "figures/sketch-small.svg",
        "figures/events.svg",
        "figures/urn-model.svg",
        "resources/model.svg"
    output:
        "figures/fig_model.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-model.py"


rule figure_dataset_correlation:
    input:
        expand("results/paper/default.dataset_correlation.{type}.svg", type=["counts", "expressions"])
    output:
        "figures/fig_dataset_correlation.svg"
    conda:
        "envs/analysis.yml"
    script:
        "scripts/fig-dataset-correlation.py"


#### utils ####


rule convert_svg:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.{fmt,(pdf|png)}"
    conda:
        "envs/cairosvg.yml"
    shell:
        "cairosvg -f {wildcards.fmt} {input} -o {output}"
