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

ruleorder: fix_monovar_vcf > bgzip_and_tabix

rule fix_monovar_vcf:
    input:
        vcf = rules.monovar_call.output
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        cat \
            <(bcftools view -h {input} | awk 'gsub(".{wildcards.reference}.{wildcards.mapper}.bam","")'1)\
            <(bcftools view -H {input})\
            | bgzip > {output.vcf}
        tabix {output.vcf}
        """

localrules: fix_monovar_header

import pysam as ps

def is_monovar_somatic(rec, samples):
    gt = rec.samples[samples[0]]["GT"]
    return any(rec.samples[sample]["GT"] != gt for sample in samples[1:])

rule get_monovar_somatic:
    input:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.vcf.gz.tbi"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.somatics.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.MonoVar.somatics.vcf.gz.tbi"
    run:
        in_vcf = ps.VariantFile(input.vcf)
        out_vcf = ps.VariantFile(output.vcf, 'wz', header=in_vcf.header)
        samples = list(in_vcf.header.samples)
        for rec in in_vcf:
            if is_monovar_somatic(rec, samples):
                out_vcf.write(rec)
        out_vcf.close()
        ps.tabix_index(output.vcf, preset="vcf")

localrules: get_monovar_somatic
