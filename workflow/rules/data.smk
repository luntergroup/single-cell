rule download_hs38DH:
	output:
		"data/references/hs38DH.fa",
		"data/references/hs38DH.fa.alt"
	conda:
		"../envs/bwakit.yaml"
	shell:
		"""
		run-gen-ref hs38DH
		mv hs38DH.fa* $(dirname {output[0]})
		"""

localrules: download_hs38DH

rule generate_hs38DH_chromosomes_bed:
	input:
		"data/references/hs38DH.fa.fai"
	output:
		"data/references/hs38DH.chromosomes.bed"
	shell:
		"head -23 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

localrules: download_hs38DH

rule download_sra_single:
	output:
		"data/reads/raw/{project}/{sample}.fastq.gz",
	conda:
		"../envs/sra-tools.yaml"
	params:
		accession = lambda wildcards: config["sra"][(wildcards.sample)]
	wildcard_constraints:
			sample='|'.join([re.escape(s) for s,_ in config["sra"].items()])
	shell:
		"""
		prefetch --max-size 200G {params.accession}
		fastq-dump {params.accession}
		rm -rf {params.accession}
		gzip -c {params.accession}.fastq > {output}
		rm {params.accession}.fastq
		"""

rule download_sra_paired:
	output:
		fq1 = "data/reads/raw/{project}/{sample}.R1.fastq.gz",
		fq2 = "data/reads/raw/{project}/{sample}.R2.fastq.gz"
	conda:
		"../envs/sra-tools.yaml"
	params:
		accession = lambda wildcards: config["sra"][(wildcards.sample)]
	wildcard_constraints:
			sample='|'.join([re.escape(s) for s,_ in config["sra"].items()])
	shell:
		"""
		prefetch --max-size 200G {params.accession}
		fastq-dump --split-3 {params.accession}
		rm -rf {params.accession}
		gzip -c {params.accession}_1.fastq > {output.fq1}
		gzip -c {params.accession}_2.fastq > {output.fq2}
		rm {params.accession}_1.fastq {params.accession}_2.fastq
		"""

localrules: download_sra_single, download_sra_paired
