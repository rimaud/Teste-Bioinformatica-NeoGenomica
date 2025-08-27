rule download_fastq_r1:
    output:
        "data/samples/{sample}_R1_001.fastq.gz"
    log:
        f"{LOG_DIR}/download_fastq_r1_{{sample}}.log"
    params:
        url="https://bioresources.s3.amazonaws.com/{sample}_R1_001.fastq.gz"
    shell:
        "mkdir -p data/samples {LOG_DIR} && "
        "wget -c -O {output} {params.url} &> {log}"

rule download_fastq_r2:
    output:
        "data/samples/{sample}_R2_001.fastq.gz"
    log:
        f"{LOG_DIR}/download_fastq_r2_{{sample}}.log"
    params:
        url="https://bioresources.s3.amazonaws.com/{sample}_R2_001.fastq.gz"
    shell:
        "mkdir -p data/samples {LOG_DIR} && "
        "wget -c -O {output} {params.url} &> {log}"


rule download_reference:
    output:
        fasta=f"{REF_DIR}/hg19.fasta",
        fai=f"{REF_DIR}/hg19.fasta.fai",
        amb=f"{REF_DIR}/hg19.fasta.amb",
        ann=f"{REF_DIR}/hg19.fasta.ann",
        bwt=f"{REF_DIR}/hg19.fasta.bwt",
        pac=f"{REF_DIR}/hg19.fasta.pac",
        sa=f"{REF_DIR}/hg19.fasta.sa"
    log:
        f"{LOG_DIR}/download_reference.log"
    shell:
        """
        mkdir -p {REF_DIR} {LOG_DIR} &&
        wget -c -O {output.fasta} https://bioresources.s3.amazonaws.com/hg19.fasta &> {log} &&
        wget -c -O {output.fai} https://bioresources.s3.amazonaws.com/hg19.fasta.fai &>> {log} &&
        wget -c -O {output.amb} https://bioresources.s3.amazonaws.com/hg19.fasta.amb &>> {log} &&
        wget -c -O {output.ann} https://bioresources.s3.amazonaws.com/hg19.fasta.ann &>> {log} &&
        wget -c -O {output.bwt} https://bioresources.s3.amazonaws.com/hg19.fasta.bwt &>> {log} &&
        wget -c -O {output.pac} https://bioresources.s3.amazonaws.com/hg19.fasta.pac &>> {log} &&
        wget -c -O {output.sa} https://bioresources.s3.amazonaws.com/hg19.fasta.sa &>> {log}
        """
