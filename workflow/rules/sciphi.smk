import os.path as path

rule samtools_mpileup_sciphi:
    input:
        reference = "data/references/{reference}.fa",
        bams = lambda wildcards: expand("data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam",
                        project=wildcards.project, sample=config["groups"][wildcards.group],
                        reference=wildcards.reference, mapper=wildcards.mapper),
        bais = lambda wildcards: expand("data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam.bai",
                        project=wildcards.project, sample=config["groups"][wildcards.group],
                        reference=wildcards.reference, mapper=wildcards.mapper),
        bed = "data/references/{reference}.chromosomes.bed"
    output:
        "results/calls/{project}/SCIPhi/{group}.{reference}.{mapper}.SCIPhi.mpileup"
    log:
        "logs/{project}/SCIPhi/{group}.{reference}.{mapper}.mpileup.log"
    threads: 2
    shell:
        "(samtools mpileup "
        " -f {input.reference} "
        " -B "
        " -d 1000 "
        " -q 40 "
        " -Q 30 "
        " -l {input.bed} "
        " {input.bams} "
        " > {output} ) >{log}"

rule create_sciphi_bam_file_list:
    input:
        lambda wildcards: expand("data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam",
                                project=wildcards.project, sample=config["groups"][wildcards.group],
                                reference=wildcards.reference, mapper=wildcards.mapper)
    output:
        txt = "results/calls/{project}/SCIPhi/{group}.{reference}.{mapper}.SCIPhi.filenames.txt"
    run:
        with open(output.txt, 'w') as out_file:
            for bam in input:
                out_file.write(bam + "\tCT\n")

localrules: create_sciphi_bam_file_list

rule sciphi:
    input:
        filenames = "results/calls/{project}/SCIPhi/{group}.{reference}.{mapper}.SCIPhi.filenames.txt",
        mpileup = "results/calls/{project}/SCIPhi/{group}.{reference}.{mapper}.SCIPhi.mpileup"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.SCIPhi.vcf"
    params:
        out = lambda wildcards, output: path.splitext(str(output))[0]
    log:
        "logs/{project}/SCIPhi/{group}.{reference}.{mapper}.log"
    threads: 2
    conda:
        "../envs/sciphi.yaml"
    shell:
        "sciphi "
        " -o {params.out} "
        " --in {input.filenames} "
        " --lz 1 "
        " {input.mpileup} 2>&1 | tee {log}"
