rule bwa_index:
	input:
		fa="data/references/{reference}.fa",
		fai="data/references/{reference}.fa.fai"
	output:
		"data/references/{reference}.fa.amb",
		"data/references/{reference}.fa.ann",
		"data/references/{reference}.fa.bwt",
		"data/references/{reference}.fa.pac",
		"data/references/{reference}.fa.sa"
	conda:
		"../envs/bwa.yaml"
	shell:
		"bwa index {input.fa}"

localrules: bwa_index

rule bwa_map:
	input:
		rules.bwa_index.output,
		fa="data/references/{reference}.fa",
		fai="data/references/{reference}.fa.fai",
		fq1="data/reads/raw/{sample}.{library}.R1.fastq.gz",
		fq2="data/reads/raw/{sample}.{library}.R2.fastq.gz"
	output:
		"data/reads/mapped/{sample}.{library}.{reference}.bwa.bam"
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{library}\tPU:Illumina",
		sort_threads=4,
		sort_memory_per_thread="4G"
	log:
		"logs/bwa/{sample}.{library}.{reference}.log"
	threads: 20
	conda:
		"../envs/bwa.yaml"
	shell:
		"(bwa mem -t {threads} -R '{params.rg}' {input.fa} {input.fq1} {input.fq2} | \
		  samtools view -bh | \
		  samtools sort -@ {params.sort_threads} -m {params.sort_memory_per_thread} -o {output}) \
		 2> {log}"

rule samtools_index_bam:
	input:
		"data/reads/mapped/{sample}.{library}.{reference}.{mapper}.bam"
	output:
		"data/reads/mapped/{sample}.{library}.{reference}.{mapper}.bam.bai"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools index {input}"

localrules: samtools_index_bam
