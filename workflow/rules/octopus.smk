rule octopus_call:
	input:
		reference = "data/references/{reference}.fa",
		bams = lambda wildcards: expand("data/reads/mapped/{sample}.{{library}}.{{reference}}.{{mapper}}.bam", sample=config["batches"][wildcards.batch]),
		bais = lambda wildcards: expand("data/reads/mapped/{sample}.{{library}}.{{reference}}.{{mapper}}.bam.bai", sample=config["batches"][wildcards.batch]),
		bed = "data/references/{reference}.chromosomes.bed"
	output:
		vcf = "results/calls/{batch}.{library}.{reference}.{mapper}.Octopus.vcf.gz",
		vcf_index = "results/calls/{batch}.{library}.{reference}.{mapper}.Octopus.vcf.gz.tbi"
	params:
		other = config["caller_options"]["Octopus"]
	log:
		"logs/octopus/{batch}.{library}.{reference}.{mapper}.log"
	benchmark:
		"results/benchmarks/octopus/{batch}.{library}.{reference}.{mapper}.tsv"
	threads: 20
	shell:
		"(octopus \
		 -R {input.reference} \
		 -I {input.bams} \
		 -t {input.bed} \
		 -o {output.vcf} \
		 --threads {threads} \
		 {params.other}) \
		 2> {log}"
