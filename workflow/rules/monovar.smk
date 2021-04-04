rule create_monovar_bam_file_list:
    input:
        bams = lambda wildcards: expand("data/reads/mapped/{{project}}/{sample}.{{reference}}.{{mapper}}.bam",
                         sample=config["groups"][wildcards.group]),
        bais = lambda wildcards: expand("data/reads/mapped/{{project}}/{sample}.{{reference}}.{{mapper}}.bam.bai",
                         sample=config["groups"][wildcards.group]),
    output:
        txt = "results/calls/{project}/MonoVar/{group}.{reference}.{mapper}.MonoVar.filenames.txt"
    run:
        with open(output.txt, 'w') as out_file:
            for bam in input.bams:
                out_file.write(bam + '\n')

localrules: create_monovar_bam_file_list

rule monovar_call:
    input:
        reference = "data/references/{reference}.fa",
        filenames = rules.create_monovar_bam_file_list.output,
        bed = "data/references/{reference}.chromosomes.bed"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.vcf"
    params:
        other = config["caller_options"]["MonoVar"]
    log:
        "logs/{project}/MonoVar/{group}.{reference}.{mapper}.log"
    benchmark:
        "results/benchmarks/{project}/MonoVar/{group}.{reference}.{mapper}.tsv"
    threads: 20
    conda:
        "../envs/monovar.yaml"
    shell:
        "(samtools mpileup \
            -l {input.bed} \
            -Q 0 \
            -q 40 \
            -d 10000 \
            -f {input.reference} \
            -b {input.filenames} | \
        monovar \
            -m {threads} \
            -f {input.reference} \
            -b {input.filenames} \
            -o {output.vcf} \
        )2> {log}"
