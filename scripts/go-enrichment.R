library("GOstats")
library("biomaRt")
library("GO.db")
library("org.Hs.eg.db")
library("mutoss")

est <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)
rownames(est) <- est$feat

# translate gene names to entrez ids
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ids <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = est$feat, mart = ensembl)
# take the first entrez id in case of ambiguity
ids <- ids[!duplicated(ids$hgnc_symbol), ]
rownames(ids) <- ids$hgnc_symbol
# lookup entrez ids
est$entrez <- ids[est$feat, "entrezgene"]

# define foreground genes
significant <- est$diff_fdr <= 0.05
foreground <- est[significant, ]
if(nrow(foreground) > 0) {
    # define test parameters
    params <- new("GOHyperGParams",
                  geneIds = foreground$entrez,
                  universeGeneIds = est$entrez,
                  ontology = "BP",
                  annotation = "org.Hs.eg",
                  pvalueCutoff = 0.05,
                  conditional = TRUE,
                  testDirection = "over")
    results <- hyperGTest(params)
    print(results)
    goterms <- summary(results)

    # correct for multiple testing. We use Benjamini-Yekuteli here, because the performed tests are strongly dependent
    #by <- BY(goterms$Pvalue, 0.05)
    #goterms$adjPvalue <- by[["adjPValues"]]

    write.table(goterms, file = snakemake@output[["terms"]], row.names = FALSE, quote = FALSE, sep = "\t")

    # associate genes with go terms
    ids <- ids[!is.na(ids$entrezgene), ]
    rownames(ids) <- ids$entrezgene
    genes <- stack(geneIdUniverse(results))
    colnames(genes) <- c("entrezid", "goterm")
    genes <- genes[, c("goterm", "entrezid")]
    genes$gene <- ids[genes$entrezid, "hgnc_symbol"]
    genes$entrezid <- NULL
    genes$cv <- est[genes$gene, "cv_map"]
    genes$fdr <- est[genes$gene, "diff_fdr"]

    write.table(genes, file = snakemake@output[["genes"]], row.names = FALSE, quote = FALSE, sep = "\t")
} else {
    file.create(snakemake@output[["terms"]])
    file.create(snakemake@output[["genes"]])
}
