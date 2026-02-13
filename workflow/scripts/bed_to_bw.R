library(rtracklayer)

# Load BED file
fiveprime_bed <- import.bed(snakemake@input$bed)

# Create SeqInfo Object
fai <- read.table(snakemake@input$fai, header = FALSE, stringsAsFactors = FALSE)
fai <- fai[fai[,1] %in% snakemake@params$region,]
seqinfo(fiveprime_bed) <- Seqinfo(seqnames = fai[,1], seqlengths = fai[,2])

# Create BigWig Object
w = ifelse(length(fiveprime_bed) > 0, fiveprime_bed$score, 0)
fiveprime_coverage <- coverage(fiveprime_bed, weight = w)
export.bw(fiveprime_coverage, snakemake@output[[1]])