rule deeptools_5prime_coverage:
    input:
        bam=get_cram,
        bai=get_crai,
    output:
        "results/deeptools/coverage/{sample}_{direction}_fiveprime.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus",
    threads: 4
    params:
        extra=deeptools_fiveprime_extra,
    log:
        "results/deeptools/5prime_coverage/{sample}_{direction}_fiveprime.log",
    message:
        "generate normalized 5prime coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"
