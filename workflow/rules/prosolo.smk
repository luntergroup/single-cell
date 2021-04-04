rule prosolo:
	input:
		reference = "data/references/{reference}.fa",
		cell_bam = "data/reads/mapped/{project}/{cell}.{reference}.{mapper}.bam",
		cell_bai = "data/reads/mapped/{project}/{cell}.{reference}.{mapper}.bam.bai",
		bulk_bam = "data/reads/mapped/{project}/{bulk}.{reference}.{mapper}.bam",
		bulk_bai = "data/reads/mapped/{project}/{bulk}.{reference}.{mapper}.bam.bai",
		bed = "data/references/{reference}.chromosomes.bed"
	output:
		"results/calls/{project}/{cell}+{bulk}.{reference}.{mapper}.prosolo.bcf"
	params:
		other = config["caller_options"]["prosolo"]
	log:
		 "logs/{project}/prosolo/{cell}+{bulk}.{reference}.{mapper}.log"
	benchmark:
		"results/benchmarks/{project}/prosolo/{cell}+{bulk}.{reference}.{mapper}.tsv"
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

def _get_prosolo_bcfs(wildcards):
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
		raise Exception("No valid bulk group")
	res = []
	for cell in cells:
		res.append("results/calls/" + wildcards.project + "/" + cell + "+" + bulk_group + "." + wildcards.reference + "." + wildcards.mapper + ".prosolo.bcf")
	return res

rule merge_prosolo_samples:
	input:
		bcfs = _get_prosolo_bcfs
	output:
		vcf = "results/calls/{project}/{group}.{reference}.{mapper}.prosolo.vcf.gz",
		index = "results/calls/{project}/{group}.{reference}.{mapper}.prosolo.vcf.gz.tbi"
	shell:
		""
