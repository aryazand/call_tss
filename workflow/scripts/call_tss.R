library(rtracklayer)

min_score <- ifelse(is.null(snakemake@params[["min_score"]]), 1000, snakemake@params[["min_score"]])
adjusted_fold_enrichment <- import.bw(snakemake@input[[1]], as = "Rle")

# Get max value within a window of size Kernal
max_per_window <- S4Vectors::runq(adjusted_fold_enrichment , k=snakemake@params[["kernal_size"]], i = snakemake@params[["kernal_size"]], endrule = "constant")

# Call TSS
TSS <- (adjusted_fold_enrichment > min_score) & (adjusted_fold_enrichment == max_per_window)
TSS <- as(TSS, "GRanges")
TSS <- TSS[TSS$score == TRUE]

if (length(TSS) == 0) {
    warning("No TSS found with score above the minimum threshold.")
}

if(snakemake@wildcards[["direction"]] == "plus") {
    strand(TSS) <- "+"
} else {
    strand(TSS) <- "-"
}

TSS$score <- NULL
names(TSS) <- rep(snakemake@wildcards["sample"], length(TSS))

export.bed(TSS, snakemake@output[[1]])