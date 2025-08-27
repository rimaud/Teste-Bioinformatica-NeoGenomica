rule fastqc:
    input:
        r1 = "data/samples/{sample}_R1_001.fastq.gz",
        r2 = "data/samples/{sample}_R2_001.fastq.gz"
    output:
        r1_html = f"{OUTDIR}/qc/fastqc/{{sample}}_R1_001_fastqc.html",
        r1_zip = f"{OUTDIR}/qc/fastqc/{{sample}}_R1_001_fastqc.zip",
        r2_html = f"{OUTDIR}/qc/fastqc/{{sample}}_R2_001_fastqc.html",
        r2_zip = f"{OUTDIR}/qc/fastqc/{{sample}}_R2_001_fastqc.zip"
    params:
        extra= "--kmers 7",
    log:
        f"{LOG_DIR}/qc/fastqc_raw_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/qc/fastqc {LOG_DIR}
        fastqc -o {OUTDIR}/qc/fastqc {input.r1} {input.r2} &> {log}
        """

rule fastp:
    input:
        r1 = "data/samples/{sample}_R1_001.fastq.gz",
        r2 = "data/samples/{sample}_R2_001.fastq.gz"
    output:
        trimmed_r1  = f"{OUTDIR}/qc/trimmed/{{sample}}_R1_001.fastq.gz",
        trimmed_r2  = f"{OUTDIR}/qc/trimmed/{{sample}}_R2_001.fastq.gz",
        html = f"{OUTDIR}/qc/fastp/{{sample}}_fastp.html",
        json  = f"{OUTDIR}/qc/fastp/{{sample}}_fastp.json"
    log:
        f"{LOG_DIR}/qc/trimmed_{{sample}}.log"
    threads: THREADS
    params:
        extra="--detect_adapter_for_pe --trim_front1 9 --trim_front2 9"
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {LOG_DIR}/qc/trimmed_{{sample}}.log
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.trimmed_r1} \
            -O {output.trimmed_r2} \
            -h {output.html} \
            -j {output.json} \
            -w {threads} \
            --trim_front1 9 \
            --trim_front2 9 \
            &> {log}
        """


rule fastq_post_fastp:
    input:
        r1 = f"{OUTDIR}/qc/trimmed/{{sample}}_R1_001.fastq.gz",
        r2 = f"{OUTDIR}/qc/trimmed/{{sample}}_R2_001.fastq.gz"
    output:
        r1_html = f"{OUTDIR}/qc/fastqc_post_trim/{{sample}}_R1_001_fastqc.html",
        r1_zip = f"{OUTDIR}/qc/fastqc_post_trim/{{sample}}_R1_001_fastqc.zip",
        r2_html = f"{OUTDIR}/qc/fastqc_post_trim/{{sample}}_R2_001_fastqc.html",
        r2_zip = f"{OUTDIR}/qc/fastqc_post_trim/{{sample}}_R2_001_fastqc.zip"
    params:
        extra= "--kmers 7"
    log:
        f"{LOG_DIR}/qc/fastqc_post_trim_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/fastqc {LOG_DIR}
        fastqc -o {OUTDIR}/fastqc {input.r1} {input.r2} &> {log}
        """