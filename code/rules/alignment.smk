rule bwa_mem:
    input:
        r1 = f"{OUTDIR}/qc/trimmed/{{sample}}_R1_001.fastq.gz",
        r2 = f"{OUTDIR}/qc/trimmed/{{sample}}_R2_001.fastq.gz",
        ref = REF_FASTA
    output:
        bam = f"{OUTDIR}/bwa/alignment/{{sample}}.bam"
    log:
        f"{LOG_DIR}/alignment/alignment_bwa_mem_{{sample}}.log"
    threads: THREADS
    wildcard_constraints:
        sample="[^.]+"
    shell:
        """
        RG="@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"

        mkdir -p {OUTDIR}/bwa/alignment {LOG_DIR}

        bwa mem -M -t {threads} -R "$RG" {input.ref} {input.r1} {input.r2} | \
        samtools view -h -q 20 -bS - > {output.bam} 2> {log}
        """

rule sort_bam:
    input:
        bam = f"{OUTDIR}/bwa/alignment/{{sample}}.bam"
    output:
        sorted_bam = f"{OUTDIR}/bwa/alignment/{{sample}}.sorted.bam"
    log:
        f"{LOG_DIR}/alignment/sort_bam_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam} 2> {log}
        """

rule mark_duplicates:
    input:
        bam = f"{OUTDIR}/bwa/alignment/{{sample}}.sorted.bam"
    output:
        bam = f"{OUTDIR}/bwa/marked/{{sample}}.bam",
        metrics = f"{OUTDIR}/bwa/marked/{{sample}}.metrics.txt"
    log:
        f"{LOG_DIR}/mark_duplicates/alignment_mark_duplicates_{{sample}}.log"
    threads: THREADS
    conda: "../envs/environment.yaml"
    shell:
        """
        mkdir -p {OUTDIR}/bwa/marked {LOG_DIR}/mark_duplicates
        picard MarkDuplicates \
            CREATE_INDEX=true \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            VALIDATION_STRINGENCY=STRICT \
            REMOVE_DUPLICATES=true \
            TMP_DIR=output/tmp \
            > {log} 2>&1
        """


rule index_bam:
    input:
        bam = f"{OUTDIR}/bwa/marked/{{sample}}.bam"
    output:
        bai = f"{OUTDIR}/bwa/marked/{{sample}}.bam.bai"
    log:
        f"{LOG_DIR}/index_bam/alignment_index_bam_{{sample}}.log"
    threads: 20
    conda: "../envs/environment.yaml"
    shell:
        """
        samtools index {input.bam} 2> {log}
        """