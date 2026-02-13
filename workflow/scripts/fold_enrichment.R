library(rtracklayer)

fg_bw <- import.bw(snakemake@input[["fg"]],as = "RleList")
bg_bw <- import.bw(snakemake@input[["bg"]],as = "RleList")
pseudo_count <- ifelse(is.null(snakemake@params[["pseudo_count"]]), 1, snakemake@params[["pseudo_count"]])

print(pseudo_count)

# Calculate fold enrichment
fold_enrichment <- fg_bw / (bg_bw + pseudo_count)

# Export fold enrichment
export.bw(fold_enrichment, snakemake@output[[1]])