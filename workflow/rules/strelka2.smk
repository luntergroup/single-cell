from pathlib import Path

rule strelka2_call:
    input:
        reference = "data/references/{reference}.fa",
        bams = lambda wildcards: expand("data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam",
                        project=wildcards.project, sample=config["groups"][wildcards.group],
                        reference=wildcards.reference, mapper=wildcards.mapper),
        bais = lambda wildcards: expand("data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam.bai",
                        project=wildcards.project, sample=config["groups"][wildcards.group],
                        reference=wildcards.reference, mapper=wildcards.mapper),
        bed = "data/references/{reference}.chromosomes.bed.gz",
        bed_idx = "data/references/{reference}.chromosomes.bed.gz.tbi"
    output:
        directory("results/calls/{project}/{group}.{reference}.{mapper}.Strelka2.runDir")
    params:
        bams = lambda wildcards: expand("--bam data/reads/mapped/{project}/{sample}.{reference}.{mapper}.bam",
                        project=wildcards.project, sample=config["groups"][wildcards.group],
                        reference=wildcards.reference, mapper=wildcards.mapper)
    log:
        "logs/{project}/strelka2/{group}.{reference}.{mapper}.log"
    threads: 20
    conda:
        "../envs/strelka2.yaml"
    shell:
        "(configureStrelkaGermlineWorkflow.py \
         --referenceFasta {input.reference} \
         {params.bams} \
         --callRegions {input.bed} \
         --runDir {output} \
         && {output}/runWorkflow.py \
         -m local \
         -j {threads} \
         )2> {log}"

rule copy_strelka2_variants:
    input:
        rules.strelka2_call.output
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Strelka2.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Strelka2.vcf.gz.tbi"
    shell:
        """
        cp {input}/results/variants/variants.vcf.gz {output.vcf}
        cp {input}/results/variants/variants.vcf.gz.tbi {output.index}
        """

localrules: copy_strelka2_variants

def _get_control_normal_bulk_bam(wildcards):
    if "controls" in config["groups"] and "normals" in config["groups"] and "bulks" in config["groups"]:
        for bulk in config["groups"]["bulks"]:
            if bulk in config["groups"]["controls"] and bulk in config["groups"]["normals"]:
                return "data/reads/mapped/" + wildcards.project + "/" + bulk + "." + wildcards.reference + "." + wildcards.mapper + ".bam"
    return None

def _get_control_tumour_bulk_bam(wildcards):
    if "controls" in config["groups"] and "normals" in config["groups"] and "bulks" in config["groups"]:
        for bulk in config["groups"]["bulks"]:
            if bulk in config["groups"]["controls"] and bulk not in config["groups"]["normals"]:
                return "data/reads/mapped/" + wildcards.project + "/" + bulk + "." + wildcards.reference + "." + wildcards.mapper + ".bam"
    return None

rule strelka2_somatic_call:
    input:
        reference = "data/references/{reference}.fa",
        normal_bam = _get_control_normal_bulk_bam,
        tumour_bam =_get_control_tumour_bulk_bam,
        bed = "data/references/{reference}.chromosomes.bed.gz",
        bed_idx = "data/references/{reference}.chromosomes.bed.gz.tbi"
    output:
        directory("results/calls/{project}/controls.{reference}.{mapper}.Strelka2.somatics.runDir")
    log:
        "logs/{project}/strelka2/controls.{reference}.{mapper}.somatics.log"
    threads: 20
    conda:
        "../envs/strelka2.yaml"
    shell:
        "(configureStrelkaSomaticWorkflow.py \
         --referenceFasta {input.reference} \
         --normalBam {input.normal_bam} \
         --tumorBam {input.tumour_bam} \
         --callRegions {input.bed} \
         --runDir {output} \
         && {output}/runWorkflow.py \
         -m local \
         -j {threads} \
         )2> {log}"

rule merge_strelka2_somatic_variants:
    input:
        rules.strelka2_somatic_call.output
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Strelka2.somatics.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Strelka2.somatics.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools concat -a -Oz -o {output.vcf} {input}/results/variants/somatic.snvs.vcf.gz {input}/results/variants/somatic.indels.vcf.gz
        tabix {output.vcf}
        """

localrules: copy_strelka2_somatic_variants
