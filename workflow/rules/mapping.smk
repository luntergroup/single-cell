rule bwa_index:
    input:
        fa="data/references/{reference}.fa",
        fai="data/references/{reference}.fa.fai"
    output:
        "data/references/{reference}.fa.amb",
        "data/references/{reference}.fa.ann",
        "data/references/{reference}.fa.bwt",
        "data/references/{reference}.fa.pac",
        "data/references/{reference}.fa.sa"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.fa}"

localrules: bwa_index

def _get_fastqs(wildcards):
    if "sample_annotations" in config:
        if wildcards.sample in config["sample_annotations"]:
            if "layout" in config["sample_annotations"][wildcards.sample]:
                if config["sample_annotations"][wildcards.sample]["layout"] == "SINGLE":
                    return "data/reads/raw/" + wildcards.project + "/" + wildcards.sample + ".fastq.gz"
    return ["data/reads/raw/" + wildcards.project + "/" + wildcards.sample + ".R1.fastq.gz",
            "data/reads/raw/" + wildcards.project + "/" + wildcards.sample + ".R2.fastq.gz"]

rule bwa_map:
    input:
        rules.bwa_index.output,
        fa = "data/references/{reference}.fa",
        fai = "data/references/{reference}.fa.fai",
        fqs = _get_fastqs
    output:
        "data/reads/mapped/{project}/{sample}.{reference}.bwa.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{project}\tPU:Illumina",
        sort_threads=4,
        sort_memory_per_thread="4G"
    log:
        "logs/bwa/{project}/{sample}.{reference}.log"
    threads: 20
    conda:
        "../envs/bwa.yaml"
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.fa} {input.fqs} | \
          samtools view -bh | \
          samtools sort -@ {params.sort_threads} -m {params.sort_memory_per_thread} -o {output}) \
         2> {log}"

rule samtools_index_bam:
    input:
        "data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam"
    output:
        "data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

localrules: samtools_index_bam
