def _get_sample_dropout_option(wildcards):
    concentrations = []
    for sample in config["groups"][wildcards.group]:
        if sample not in config["groups"]["cells"]:
            concentrations.append(sample + "=50")
    if len(concentrations) > 0:
        return "--sample-dropout " + " ".join(concentrations)
    else:
        return ""

def _get_normal_samples(wildcards):
    result = []
    if "normals" in config["groups"]:
        for sample in config["groups"][wildcards.group]:
            if sample in config["groups"]["normals"]:
                result.append(sample)
    return result

def _get_normal_samples_option(wildcards):
    normals = _get_normal_samples(wildcards)
    if len(normals) > 0:
        return "--normal-samples " + " ".join(normals)
    else:
        return ""

def _get_calling_model(wildcards):
    if any(sample in config["groups"]["cells"] for sample in config["groups"][wildcards.group]):
        return "cell"
    elif len(_get_normal_samples(wildcards)) > 0:
        return "cancer"
    elif len(config["groups"][wildcards.group]) == 1:
        return "individual"
    else:
        return "population"

def _get_sequence_error_model(wildcards):
    if any(sample in config["groups"]["cells"] for sample in config["groups"][wildcards.group]):
        return "10X"
    else:
        return "PCRF"

rule octopus_call:
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
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Octopus.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Octopus.vcf.gz.tbi"
    params:
        calling_model = _get_calling_model,
        sequence_model = _get_sequence_error_model,
        normals = _get_normal_samples_option,
        sample_dropouts = _get_sample_dropout_option,
        forest = "/opt/octopus/resources/forests/cell.forest",
        other = config["caller_options"]["Octopus"]
    log:
        "logs/{project}/octopus/{group}.{reference}.{mapper}.log"
    benchmark:
        "results/benchmarks/{project}/octopus/{group}.{reference}.{mapper}.tsv"
    threads: 20
    container:
        "docker://dancooke/octopus:sc_paper"
    shell:
        "(octopus \
         -C {params.calling_model} \
         -R {input.reference} \
         -I {input.bams} \
         -t {input.bed} \
         -o {output.vcf} \
         --threads {threads} \
         --mask-inverted-soft-clipping \
         --mask-3prime-shifted-soft-clipped-heads \
         --sequence-error-model {params.sequence_model} \
         --forest {params.forest} \
         {params.normals} \
         {params.sample_dropouts} \
         {params.other}) \
         2> {log}"

rule get_octopus_somatic:
    input:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Octopus.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Octopus.vcf.gz.tbi"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Octopus.somatics.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Octopus.somatics.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools view -i 'SOMATIC=1' -Oz -o {output.vcf} {input.vcf} && tabix {output.vcf}
        """

localrules: get_octopus_somatic
