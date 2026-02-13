rule filter_alignment:
    input:
        aln=get_cram,
        idx=get_crai,
    output:
        temp("results/frags/{sample}.filtered.bam"),
    message:
        """--- Filtering alignments for proper pairs and mapping quality >= 30."""
    log:
        "results/frags/{sample}.filter.log",
    params:
        extra=config["filter_alignment"]["extra"],  # optional params string
        region=config["filter_alignment"]["region"],  # optional region string
    threads: 2
    wrapper:
        "v8.1.1/bio/samtools/view"


rule samtools_collate:
    input:
        rules.filter_alignment.output,
    output:
        temp("results/frags/{sample}.collated.bam"),
    log:
        "results/frags/{sample}.collate.log",
    message:
        """--- Collating alignments."""
    params:
        extra="-u",
    threads: 2
    wrapper:
        "v8.1.1/bio/samtools/collate"


rule alignment_to_bedpe:
    input:
        rules.samtools_collate.output,
    output:
        temp("results/frags/{sample}.bedpe"),
    message:
        """--- Converting alignments to BEDPE format."""
    log:
        "results/frags/{sample}.bedpe.log",
    params:
        extra="-bedpe",
    wrapper:
        "v8.1.1/bio/bedtools/bamtobed"


rule bedpe_to_bed:
    input:
        rules.alignment_to_bedpe.output,
    output:
        "results/frags/{sample}.bed",
    message:
        """--- Converting BEDPE to BED format."""
    log:
        "results/frags/{sample}.bed.log",
    shell:
        """        
        awk 'BEGIN{{FS=OFS="\\t"}} 
            $0 !~ /^#/ && NF>=10 && $1==$4 {{
                s = ($2 < $5 ? $2 : $5);
                e = ($3 > $6 ? $3 : $6);
                print $1, s, e, $7, $8, $10
            }}' {input} > {output}
        """


rule collapse_frags:
    input:
        rules.bedpe_to_bed.output,
    output:
        "results/frags/{sample}.collapsed.bed",
    message:
        """--- Collapsing fragments to unique positions."""
    log:
        "results/frags/{sample}.collapse.log",
    shell:
        """
        sort -k1,1V -k2n -k3n {input} | awk '{{print $1,$2,$3,$6}}' | uniq -c | awk 'BEGIN {{OFS="\\t"}} {{print $2,$3,$4,".",$1,$5}}' > {output}
        """


rule filter_frags:
    input:
        rules.collapse_frags.output,
    output:
        "results/frags/{sample}.filtered.bed",
    message:
        """--- Filtering collapsed fragments for count >= 2."""
    params:
        min_count=config["filter_frags"]["min_count"],  # optional count threshold
        min_width=config["filter_frags"]["min_width"],
        max_width=config["filter_frags"]["max_width"],
    log:
        "results/frags/{sample}.filter_collapsed.log",
    shell:
        """
        awk '$5 >= {params.min_count} && ($3-$2) >= {params.min_width} && ($3-$2) <= {params.max_width} {{print}}' {input} > {output}
        """


rule stranded_fiveprime_bed:
    input:
        rules.filter_frags.output,
    output:
        "results/fiveprime/{sample}_{direction}.fiveprime.bed",
    message:
        """--- Extracting 5' ends of fragments."""
    log:
        "results/fiveprime/{sample}_{direction}.fiveprime.bed.log",
    container:
        "docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
    shell:
        """
        # If {{direction}} is "plus", extract 5' ends on the plus strand; if "minus", extract 5' ends on the minus strand
        # Add scores columns together at each position
        if [[ "{wildcards.direction}" == "plus" ]]; then
            awk 'BEGIN{{OFS="\\t"}} $6 == "+" {{print $1, $2, $2+1, $4, $5, $6}}' {input} \
                | sort -k1,1 -k2,2n \
                | bedtools merge -s -c 4,5,6 -o count,sum,distinct > results/fiveprime/{wildcards.sample}_plus.fiveprime.bed
        else
            awk 'BEGIN{{OFS="\\t"}} $6 == "-" {{print $1, $3-1, $3, $4, $5, $6}}' {input} \
                | sort -k1,1 -k2,2n \
                | bedtools merge -s -c 4,5,6 -o count,sum,distinct > results/fiveprime/{wildcards.sample}_minus.fiveprime.bed
        fi
        """


rule filtered_fiveprime_bw:
    input:
        bed=rules.stranded_fiveprime_bed.output,
        fai=config["reference_fai"]
    output:
        "results/fiveprime/{sample}_{direction}.fiveprime.bw",
    params:
        region=config["filter_alignment"]["region"]
    message:
        """--- Convert 5' ends BED to BigWig format."""
    log:
        "results/fiveprime/{sample}_{direction}.fiveprime.bw.log",
    container:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.66.0--r44h15a9599_1"
    script:
        "../scripts/bed_to_bw.R"


rule filtered_fiveprime_background:
    input:
        rules.filtered_fiveprime_bw.output,
    output:
        "results/fiveprime/{sample}_{direction}.background.bw",
    params:
        kernal_size=9,  # optional kernel for smoothing
    message:
        """--- Extracting background regions from filtered fragments."""
    log:
        "results/fiveprime/{sample}_{direction}.background.log",
    container:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.66.0--r44h15a9599_1"
    script:
        "../scripts/bw_smoothing.R"


rule fold_enrichment:
    input:
        fg=rules.filtered_fiveprime_bw.output,
        bg=rules.filtered_fiveprime_background.output,
    output:
        "results/fiveprime/{sample}_{direction}.fold_enrichment.bw",
    params:
        pseudocount=1,  # optional pseudocount to avoid division by zero
    message:
        """--- Calculating fold enrichment of 5' ends over background."""
    log:
        "results/fiveprime/{sample}_{direction}.fold_enrichment.log",
    container:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.66.0--r44h15a9599_1"
    script:
        "../scripts/fold_enrichment.R"


rule adjusted_fold_enrichment:
    input:
        fg=rules.filtered_fiveprime_bw.output,
        fold_enrichment=rules.fold_enrichment.output,
    output:
        "results/fiveprime/{sample}_{direction}.adjusted_fold_enrichment.bw",
    message:
        """--- Calculating adjusted fold enrichment of 5' ends over background."""
    log:
        "results/fiveprime/{sample}_{direction}.adjusted_fold_enrichment.log",
    container:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.66.0--r44h15a9599_1"
    script:
        "../scripts/adjusted_fold_enrichment.R"


rule call_tss:
    input:
        rules.adjusted_fold_enrichment.output,
    output:
        "results/tss/{sample}_{direction}.tss.bed",
    message:
        """--- Calling TSS from adjusted fold enrichment signal."""
    log:
        "results/tss/{sample}_{direction}.tss.log",
    params:
        min_score=config["call_tss"]["min_score"],  # optional threshold for TSS calling
    container:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.66.0--r44h15a9599_1"
    script:
        "../scripts/call_tss.R"


rule consolidate_tss:
    input:
        expand(
            "results/tss/{sample}_{direction}.tss.bed",
            sample=samples.index,
            direction=["plus", "minus"],
        ),
    output:
        "results/tss/TSS.bed",
    message:
        """--- Consolidating TSS calls across all samples."""
    params:
        ## Add optional parameters
        extra="-s -c 4,5,6 -o distinct,max,distinct"
    log:
        "results/tss/TSS.log",
    wrapper:
        "v8.1.1/bio/bedtools/merge"