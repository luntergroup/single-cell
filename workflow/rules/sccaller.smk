rule generate_sccaller_bulk_het_snps:
    input:
        "results/calls/{project}/{bulk}.{reference}.{mapper}.Strelka2.vcf.gz"
    output:
        "results/calls/{project}/SCcaller_tmp/{bulk}.{reference}.{mapper}.het_snps.vcf"
    conda:
        "../envs/sccaller.yaml"
    shell:
        "bcftools view {input} -v snps -f PASS -g het > {output}"

rule sccaller:
    input:
        reference = "data/references/{reference}.fa",
        cell_bam = "data/reads/mapped/{project}/{cell}.{reference}.{mapper}.bam",
        cell_bai = "data/reads/mapped/{project}/{cell}.{reference}.{mapper}.bam.bai",
        bulk_bam = "data/reads/mapped/{project}/{bulk}.{reference}.{mapper}.bam",
        bulk_bai = "data/reads/mapped/{project}/{bulk}.{reference}.{mapper}.bam.bai",
        het_snps = "results/calls/{project}/SCcaller_tmp/{bulk}.{reference}.{mapper}.het_snps.vcf"
    output:
        "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.vcf"
    params:
        other = config["caller_options"]["SCcaller"]
    log:
        "logs/{project}/SCcaller/{cell}+{bulk}.{reference}.{mapper}.log"
    benchmark:
        "results/benchmarks/{project}/SCcaller/{cell}+{bulk}.{reference}.{mapper}.tsv"
    threads: 20
    conda:
        "../envs/sccaller.yaml"
    shell:
        "(sccaller \
            --fasta {input.reference} \
            --bam {input.cell_bam} \
            --bulk {input.bulk_bam} \
            --output {output} \
            --snp_type hsnp \
            --snp_in {input.het_snps} \
            --cpu_num {threads}) \
            2> {log}"

import csv
import gzip

rule fix_sccaller:
    input:
        vcf = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.{caller}.vcf"
    output:
        vcf = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.{caller}.fixed.vcf.gz"
    run:
        with open(input.vcf, 'rt') as in_vcf:
            with gzip.open(output.vcf, 'wt') as out_vcf:
                vcfreader = csv.reader(in_vcf, delimiter='\t')
                for rec in vcfreader:
                    if rec[0].startswith("##"):
                        if rec[0].startswith("##FILTER=<ID=No.variants<4"):
                            rec[0] = rec[0].replace("ID=No.variants<4", "ID=No_variants4")
                        elif rec[0].startswith("##FORMAT=<ID=PL,Number=G"):
                            rec[0] = rec[0].replace("ID=PL,Number=G", "ID=CPL,Number=.") # PL invalid: https://github.com/biosinodx/SCcaller/issues/10
                    elif rec[0].startswith("#"):
                        rec[-1] = wildcards.cell # SCCaller reports incorrect sample name
                    else:
                        if "No.variants<4" in rec[6]: continue
                        ref, alt = rec[3], rec[4]
                        if any(b not in ['A', 'C', 'G', 'T', 'N'] for b in ref): continue
                        alts = alt.split(',')
                        if '-' in alt:
                            longest_del = ""
                            for alt_idx, alt in enumerate(alts):
                                if alt[0] == '-':
                                    allele = ""
                                    for b in reversed(alt):
                                        if b.isdigit():
                                            break
                                        allele = b + allele
                                    if len(allele) > len(longest_del):
                                        longest_del = allele
                                    alts[alt_idx] = ref + allele
                            ref += longest_del
                            for alt_idx, alt in enumerate(alts):
                                if ref.startswith(alt):
                                    alts[alt_idx] = ref[0] + ref[len(alt):]
                        if '+' in alt:
                            for alt_idx, alt in enumerate(alts):
                                if alt[0] == '+':
                                    allele = ""
                                    for b in reversed(alt):
                                        if b.isdigit():
                                            break
                                        allele = b + allele
                                    alts[alt_idx] = ref + allele
                        rec[3], rec[4] = ref, ','.join(alts)
                        if "No.variants<4" in rec[7]:
                            rec[7] = rec[7].replace("No.variants<4", "No_variants4")
                        format_fields = rec[-2].split(':')
                        if "PL" in format_fields:
                            format_fields[format_fields.index("PL")] = "CPL"
                            rec[-2] = ':'.join(format_fields)
                    out_vcf.write('\t'.join(rec) + '\n')

localrules: fix_sccaller

ruleorder: bcftools_sort_sccaller > bgzip_and_tabix

rule bcftools_sort_sccaller:
    input:
        "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.fixed.vcf.gz"
    output:
        vcf = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.vcf.gz",
        idx = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools sort -Oz -o {output.vcf} {input}
        tabix {output.vcf}
        """
    
localrules: bcftools_sort

def _get_sccaller_vcfs(wildcards):
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
        res.append("results/calls/" + wildcards.project + "/" + cell + "+" + bulk_group + "." + wildcards.reference + "." + wildcards.mapper + ".SCcaller.vcf.gz")
    return res

rule merge_sccaller_samples:
    input:
        bcfs = _get_sccaller_vcfs
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.SCcaller.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.SCcaller.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge -Oz -o {output.vcf} {input} && tabix {output.vcf}"

localrules: merge_sccaller_samples

rule get_sccaller_somatics:
    input:
        "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.vcf.gz"
    output:
        vcf = "results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.somatics.vcf.gz",
        index ="results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.SCcaller.somatics.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        bcftools view -g het -i 'FORMAT/SO=="True"&&SUM(FORMAT/AD)>20' -Oz -o {output.vcf} {input}
        tabix {output.vcf}
        """

localrules: get_sccaller_somatics

def _get_sccaller_somatic_vcfs(wildcards):
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
        res.append("results/calls/" + wildcards.project + "/" + cell + "+" + bulk_group + "." + wildcards.reference + "." + wildcards.mapper + ".SCcaller.somatics.vcf.gz")
    return res

rule merge_sccaller_somatics:
    input:
        bcfs = _get_sccaller_somatic_vcfs
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.SCcaller.somatics.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.SCcaller.somatics.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools merge -Oz -o {output.vcf} {input} && tabix {output.vcf}"
    
localrules: merge_sccaller_somatics
