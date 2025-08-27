rule snpeff_annotate:
    input:
        vcf = f"{OUTDIR}/freebayes/variants/{{sample}}.filtered.norm.{{variant_type}}.vcf.gz"
    output:
        vcf = f"{OUTDIR}/annotation/{{sample}}.{{variant_type}}.ann.vcf.gz",
        html = f"{OUTDIR}/annotation/{{sample}}.{{variant_type}}.snpeff.html"
    params:
        jar = f"{SNPEFF_PATH}/snpEff.jar",
        genome = SNPEFF_GENOME,
        data_dir = SNPEFF_DATA
    log:
        f"{LOG_DIR}/snpeff/snpeff_{{sample}}_{{variant_type}}.log"
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/annotation {LOG_DIR}/snpeff

        java -Xmx8g -jar {params.jar} \
            -v \
            -stats {output.html} \
            -dataDir {params.data_dir} \
            {params.genome} \
            {input.vcf} \
            | bgzip -c > {output.vcf} 2> {log}
        """