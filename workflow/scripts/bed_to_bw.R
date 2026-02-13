library(rtracklayer)

##############################
# LOAD GRANGES AS BED FILE
##############################

# Import bed file 
fiveprime_bed <- import.bed(snakemake@input$bed)

# Obtain seq lengths so when converting from GRanges to Rle object, the Rle object is same length as chromosome
fai <- read.table(snakemake@input$fai, header = FALSE, stringsAsFactors = FALSE)
fai <- fai[fai[,1] %in% snakemake@params$region,]
seqinfo(fiveprime_bed) <- Seqinfo(seqnames = fai[,1], seqlengths = fai[,2])

##############################
# CONVERT GRANGES TO RLE
##############################

# Object scores from bed file to be used as weights in coverage function
w <- fiveprime_bed$score
if(is.null(w)) {
    w <- 0
}

fiveprime_coverage <- coverage(fiveprime_bed, weight = w)

##############################
# EXPORT AS BIGIWG 
##############################
export.bw(fiveprime_coverage, snakemake@output[[1]])