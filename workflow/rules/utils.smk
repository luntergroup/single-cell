rule samtools_index_fasta:
	input:
		"{fasta}"
	output:
		 "{fasta}.fai"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools faidx {input}"

localrules: samtools_index_fasta

rule rtg_format:
	input:
		"{prefix}.fa"
	output:
		directory("{prefix}.sdf")
	conda:
		"../envs/rtg.yaml"
	shell:
		"rtg format {input} -o {output}"

localrules: rtg_format

rule samtools_merge:
	input:
		lambda wildcards: expand(
			"data/reads/mapped/{{project}}/{sample}.{{reference}}.{{mapper}}.bam",
			sample=config["groups"][wildcards.group])
	output:
		"data/reads/mapped/{project}/{group}.{reference}.{mapper}.bam"
	wildcard_constraints:
			group='|'.join([re.escape(s) for s,_ in config["groups"].items()])
	run:
		if len(input) > 1:
			shell("samtools merge {output} {input}")
		else:
			shell("ln -s {input} {output} && touch -h {output}")

rule bgzip_bed:
	input:
		"data/references/{name}.bed"
	output:
		"data/references/{name}.bed.gz"
	shell:
		"bgzip -c {input} > {output}"

rule index_bgzip_bed:
	input:
		"data/references/{name}.bed.gz"
	output:
		"data/references/{name}.bed.gz.tbi"
	shell:
		"tabix {input}"

localrules: bgzip_bed, index_bgzip_bed

rule bgzip_and_tabix:
    input:
        "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.vcf"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.vcf.gz.tbi"
    shell:
        "bgzip {input} && tabix {output.vcf}"

localrules: bgzip_and_tabix
