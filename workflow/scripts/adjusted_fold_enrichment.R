library(rtracklayer)

fg_bw <- import.bw(snakemake@input[["fg"]],as = "RleList")
fold_enrichment <- import.bw(snakemake@input[["fold_enrichment"]],as = "RleList")

# Calculate fold enrichment
adjusted_fold_enrichment <- fold_enrichment*fg_bw*1e6/sum(fg_bw)

# Export fold enrichment
export.bw(adjusted_fold_enrichment, snakemake@output[[1]])