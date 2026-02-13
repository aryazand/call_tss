rule deeptools_coverage:
    input:
        bam=get_cram,
        bai=get_crai,
    output:
        "results/deeptools/coverage/{sample}_{direction}.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus",
    threads: 4
    params:
        extra=lambda wildcards: config["mapping_stats"]["deeptools_coverage"]["extra"]
        + (
            " --filterRNAstrand forward"
            if wildcards.direction in ["for", "forward", "plus"]
            else " --filterRNAstrand reverse"
        ),
    log:
        "results/deeptools/coverage/{sample}_{direction}.log",
    message:
        "generate normalized coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"


rule deeptools_5prime_coverage:
    input:
        bam=get_cram,
        bai=get_crai,
    output:
        "results/deeptools/5prime_coverage/{sample}_{direction}.bw",
    wildcard_constraints:
        direction="for|forward|plus|rev|reverse|minus",
    threads: 4
    params:
        extra=lambda wildcards: config["mapping_stats"]["deeptools_coverage"]["extra"]
        + (
            " --filterRNAstrand forward"
            if wildcards.direction in ["for", "forward", "plus"]
            else " --filterRNAstrand reverse"
        ),
    log:
        "results/deeptools/5prime_coverage/{sample}_{direction}.log",
    message:
        "generate normalized 5prime coverage using deeptools"
    wrapper:
        "v7.0.0/bio/deeptools/bamcoverage"
