def _get_prosolo_

rule prosolo_cell_and_bulk:
	input:
		reference = "data/references/{reference}.fa",
		cell_bam = "data/reads/mapped/{cell}.{library}.{reference}.{mapper}.bam",
		cell_bai = "data/reads/mapped/{cell}.{library}.{reference}.{mapper}.bam.bai",
		cell_bam = "data/reads/mapped/{bulk}.{library}.{reference}.{mapper}.bam",
		cell_bai = "data/reads/mapped/{bulk}.{library}.{reference}.{mapper}.bam.bai",
		candidates = "results/calls/{bulk}+{cell}.prosolo_candidates.{library}.{reference}.{mapper}.bcf",
		bed = "data/references/{reference}.chromosomes.bed"
	output:
		"results/calls/{cell}+{bulk}.{library}.{reference}.{mapper}.prosolo.vcf.gz"
	params:
		other = config["caller_options"]["prosolo"]
	log:
		"logs/prosolo/{batch}.{library}.{reference}.{mapper}.log"
	benchmark:
		"results/benchmarks/prosolo/{batch}.{library}.{reference}.{mapper}.tsv"
	threads: 20
	shell:
		"(prosolo \
		single-cell-bulk \
		--omit-indels \
        --candidates {input.candidates} \
        --output {output.bcf} \
        {input.cell_bam} \
        {input.bulk_bam} \
        {input.ref} \
		 2> {log}"
