rule prosolo:
    input:
        reference = "data/references/{reference}.fa",
        cell_bam = "data/reads/mapped/{project}/{cell}.{reference}.{mapper}.bam",
        cell_bai = "data/reads/mapped/{project}/{cell}.{reference}.{mapper}.bam.bai",
        bulk_bam = "data/reads/mapped/{project}/{bulk}.{reference}.{mapper}.bam",
        bulk_bai = "data/reads/mapped/{project}/{bulk}.{reference}.{mapper}.bam.bai",
        bed = "data/references/{reference}.chromosomes.bed"
    output:
        "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.Prosolo.bcf"
    params:
        other = config["caller_options"]["Prosolo"]
    log:
         "logs/{project}/Prosolo/{cell}+{bulk}.{reference}.{mapper}.log"
    benchmark:
        "results/benchmarks/{project}/Prosolo/{cell}+{bulk}.{reference}.{mapper}.tsv"
    threads: 20
    conda:
        "../envs/prosolo.yaml"
    shell:
        "(workflow/scripts/prosolo.py \
            --reference {input.reference} \
            --cell {input.cell_bam} \
            --bulk {input.bulk_bam} \
            --regions {input.bed} \
            --output {output} \
            --threads {threads} \
            )2> {log}"

def decode_phred(phred):
    return 10**(-phred[0]/10)

from math import log10

def encode_phred(p):
    return -10 * log10(1 - p)

import pysam as ps

rule add_prosolo_gt:
    input:
        "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.Prosolo.bcf"
    output:
        vcf = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.Prosolo.vcf.gz",
        index = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.Prosolo.vcf.gz.tbi"
    run:
        in_vcf = ps.VariantFile(input[0])
        in_vcf.header.formats.add("GT","1","String","Genotype")
        in_vcf.header.formats.add("GQ","1","Float","Conditional genotype quality (phred-scaled)")
        out_vcf = ps.VariantFile(output.vcf, "wz", header=in_vcf.header)
        for rec in in_vcf:
            try:
                prob_hom_ref = decode_phred(rec.info["PROB_HOM_REF"]) + decode_phred(rec.info["PROB_ERR_ALT"])
                prob_hom_alt = decode_phred(rec.info["PROB_HOM_ALT"]) + decode_phred(rec.info["PROB_ERR_REF"])
                prob_het = decode_phred(rec.info["PROB_HET"]) + decode_phred(rec.info["PROB_ADO_TO_ALT"]) + decode_phred(rec.info["PROB_ADO_TO_REF"])
                rec.qual = encode_phred(prob_hom_alt + prob_het)
                if prob_het >= prob_hom_ref and prob_het >= prob_hom_alt:
                    gt = (0,1)
                    gq = encode_phred(prob_het)
                elif prob_hom_alt >= prob_hom_ref:
                    gt = (1,1)
                    gq = encode_phred(prob_hom_alt)
                else:
                    gt = (0,0)
                    gq = encode_phred(prob_hom_ref)
                for sample, fields in rec.samples.items():
                    fields["GT"] = gt
                    fields["GQ"] = gq
                out_vcf.write(rec)
            except:
                for sample, fields in rec.samples.items():
                    fields["GT"] = (None, None)
                    fields["GQ"] = None
                out_vcf.write(rec)
        in_vcf.close()
        out_vcf.close()
        ps.tabix_index(output.vcf, preset="vcf")

def _get_prosolo_vcfs(wildcards):
    cells, bulks = [], []
    for sample in config["groups"][wildcards.group]:
        if sample in config["groups"]["cells"]:
            cells.append(sample)
        else:
            bulks.append(sample)
    bulk_group = None
    for group, samples in config["groups"].items():
        if set(samples) == set(bulks):
            bulk_group = group
            break
    if bulk_group is None:
        normal_cells = []
        if "normals" in config["groups"]:
            for cell in cells:
                if cell in config["groups"]["normals"]:
                    normal_cells.append(cell)
            for group, samples in config["groups"].items():
                if set(samples) == set(normal_cells):
                    bulk_group = group
                    break
        if bulk_group is None:
            raise Exception("No valid bulk group")
    res = []
    for cell in cells:
        res.append("results/calls/" + wildcards.project + "/" + cell + "+" + bulk_group + "." + wildcards.reference + "." + wildcards.mapper + ".Prosolo.vcf.gz")
    return res

rule merge_prosolo_samples:
    input:
        bcfs = _get_prosolo_vcfs,
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Prosolo.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Prosolo.vcf.gz.tbi"
    params:
        min_qual = 10
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        bcftools merge --force-samples {input.bcfs} | \
            bcftools view -aUc1 | \
            bcftools norm -d all | \
            bcftools filter -e 'QUAL<{params.min_qual}' -mx -s q{params.min_qual} -Oz -o {output.vcf} -
        tabix {output.vcf}
        '''

localrules: merge_prosolo_samples

import pysam as ps

def is_prosolo_somatic(rec, samples):
    gt = rec.samples[samples[0]]["GT"]
    return any(rec.samples[sample]["GT"] != gt for sample in samples[1:])

rule get_prosolo_somatic:
    input:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Prosolo.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Prosolo.vcf.gz.tbi"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.Prosolo.somatics.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.Prosolo.somatics.vcf.gz.tbi"
    run:
        in_vcf = ps.VariantFile(input.vcf)
        out_vcf = ps.VariantFile(output.vcf, 'wz', header=in_vcf.header)
        samples = list(in_vcf.header.samples)
        for rec in in_vcf:
            if is_prosolo_somatic(rec, samples):
                out_vcf.write(rec)
        out_vcf.close()
        ps.tabix_index(output.vcf, preset="vcf")

localrules: get_prosolo_somatic
