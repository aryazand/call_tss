library(rtracklayer)

# Produce Kernal 
kernal_size <- snakemake@params[["kernal_size"]]
kernal <- rep(1/(kernal_size-1), kernal_size)
kernal[ceiling(kernal_size/2)] <- 0

# Calculate background
fiveprime_coverage <- import.bw(snakemake@input[[1]],as = "RleList")
fiveprime_background <- runwtsum(fiveprime_coverage, k = kernal_size, wt = kernal, endrule = "constant")

# Export background
export.bw(fiveprime_background, snakemake@output[[1]])