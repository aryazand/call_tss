library(rtracklayer)

min_score <- ifelse(is.null(snakemake@params[["min_score"]]), 1000, snakemake@params[["min_score"]])

adjusted_fold_enrichment <- import.bw(snakemake@input[[1]])
print(adjusted_fold_enrichment)
TSS <- adjusted_fold_enrichment[adjusted_fold_enrichment$score > min_score]

if (length(TSS) == 0) {
    warning("No TSS found with score above the minimum threshold.")
}

if(snakemake@wildcards[["direction"]] == "plus") {
    strand(TSS) <- "+"
} else {
    strand(TSS) <- "-"
}

names(TSS) <- rep(snakemake@wildcards["sample"], length(TSS))

export.bed(TSS, snakemake@output[[1]])