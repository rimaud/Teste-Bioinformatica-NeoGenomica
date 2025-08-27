rule variant_calling:
    input:
        bam = f"{OUTDIR}/bwa/marked/{{sample}}.bam",
        bai = f"{OUTDIR}/bwa/marked/{{sample}}.bam.bai",
        ref = REF_FASTA,
        targets = TARGET_BED
    output:
        vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.vcf"
    log:
        f"{LOG_DIR}/variant_calling/freebayes_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/freebayes/variants {LOG_DIR}/variant_calling
        freebayes \
            --fasta-reference {input.ref} \
            --targets {input.targets} \
            --bam {input.bam} \
            --min-mapping-quality 30 \
            --min-base-quality 20 \
            --min-alternate-count 5 \
            --min-alternate-fraction 0.2 \
            --ploidy 2 \
            --use-best-n-alleles 4 \
            > {output.vcf} 2> {log}
        """
rule filter_vcf:
    input:
        vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.vcf"
    output:
        filtered_vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.vcf"
    log:
        f"{LOG_DIR}/variant_calling/filter_vcf_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {LOG_DIR}/variant_calling
        bcftools filter \
            -e 'QUAL<=20 || INFO/DP<=10 || AF[0]<=0.2' \
            {input.vcf} \
            -o {output.filtered_vcf} 2> {log}
        """

rule normalize_vcf:
    input:
        vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.vcf",
        ref = REF_FASTA
    output:
        norm_vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.norm.vcf"
    log:
        f"{LOG_DIR}/variant_calling/normalize_vcf_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {LOG_DIR}/variant_calling
        bcftools norm \
            -f {input.ref} \
            -m -both \
            {input.vcf} \
            -o {output.norm_vcf} 2> {log}
        """
rule separate_snps:
    input:
        vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.norm.vcf"
    output:
        vcf_gz = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.norm.snps.vcf.gz"
    log:
        f"{LOG_DIR}/variant_calling/separate_snps_{{sample}}.log"
    conda: "../envs/environment.yaml"
    shell:
        """
        bcftools view \
            --types snps \
            -Oz \
            -o {output.vcf_gz} \
            {input.vcf} > {log} 2>&1
        """

rule separate_indels:
    input:
        vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.norm.vcf"
    output:
        vcf_gz = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.norm.indels.vcf.gz"
    log:
        f"{LOG_DIR}/variant_calling/separate_indels_{{sample}}.log"
    conda: "../envs/environment.yaml"
    shell:
        """
        bcftools view \
            --types indels,mnps \
            -Oz \
            -o {output.vcf_gz} \
            {input.vcf} > {log} 2>&1
        """