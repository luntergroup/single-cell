rule generate_controls:
    input:
        vcf = "results/calls/{project}/controls.{reference}.bwa." + config["truth_caller"] + ".vcf.gz",
        index = "results/calls/{project}/controls.{reference}.bwa." + config["truth_caller"] + ".vcf.gz.tbi",
        bed = "data/references/{reference}.chromosomes.bed"
    output:
        vcf = "data/truth/{project}/controls.{reference}.vcf.gz",
        index = "data/truth/{project}/controls.{reference}.vcf.gz.tbi",
        bed = "data/truth/{project}/controls.{reference}.bed",
        inconsistent_bed = temp("data/truth/{project}/controls.{reference}.inconsistent.bed"),
        inconsistent_sorted_bed = temp("data/truth/{project}/controls.{reference}.inconsistent.sorted.bed")
    conda:
        "../envs/bedops.yaml"
    shell:
        """
		ln {input.vcf} {output.vcf}
		ln {input.index} {output.index}
        bcftools view -H -f PASS {output.vcf} | \
            awk -v FS='\t' -v OFS='\t' \
                '{{ \
                    for(i=11; i<=NF; i++) \
                        if(substr($10,0,3) != substr($i,0,3)) {{ \
                            print $1,$2-1,$2+length($4)+1; \
                            break; \
                        }} \
                }}' \
                > {output.inconsistent_bed}
        sort-bed {output.inconsistent_bed} > {output.inconsistent_sorted_bed}
        bedops -d {input.bed} {output.inconsistent_sorted_bed} > {output.bed}
		"""

localrules: generate_controls

rule vcfeval:
	input:
		reference = "data/references/{reference}.sdf",
		baseline_vcf = "data/truth/{project}/controls.{reference}.vcf.gz",
		baseline_vcf_index = "data/truth/{project}/controls.{reference}.vcf.gz.tbi",
        evaluation_bed = "data/truth/{project}/controls.{reference}.bed",
		calls_vcf = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.vcf.gz",
		calls_vcf_index = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.vcf.gz.tbi"
	output:
		directory("results/eval/{project}/{group}.{reference}.{mapper}.{caller}.{sample}.vs.{control}.{filter}.{match}.vcfeval")
	log:
		"logs/{project}/eval/{caller}/{group}.{reference}.{mapper}.{sample}.vs.{control}.{filter}.{match}.vcfeval.log"
	params:
		score_field = lambda wildcards: config["caller_score_fields"][wildcards.caller],
		all_records = lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
		squash_ploidy = lambda wildcards: "--squash-ploidy" if wildcards.match in ["AL", "CB"] else "",
		decompose = lambda wildcards: "--decompose" if wildcards.match=="AL" else "",
		output_mode = lambda wildcards: "combine" if wildcards.match=="CB" else "split",
	threads: 20
	resources:
		mem_gb=40
	conda:
		"../envs/rtg.yaml"
	shell:
		"(rtg \
			RTG_MEM={resources.mem_gb}g \
			vcfeval \
			-t {input.reference} \
			-b {input.baseline_vcf} \
            --evaluation-regions {input.evaluation_bed} \
			-c {input.calls_vcf} \
			-o {output} \
			--sample {wildcards.control},{wildcards.sample} \
			-f {params.score_field} \
			--ref-overlap \
			--output-mode {params.output_mode} \
			--threads {threads} \
			{params.all_records} \
			{params.squash_ploidy} \
			{params.decompose} \
			)2> {log}"

rule generate_somatic_controls:
    input:
        vcf = "results/calls/{project}/controls.{reference}.bwa." + config["truth_caller"] + ".somatics.vcf.gz",
        index = "results/calls/{project}/controls.{reference}.bwa." + config["truth_caller"] + ".somatics.vcf.gz.tbi"
    output:
        vcf = "data/truth/{project}/controls.{reference}.somatics.vcf.gz",
        index = "data/truth/{project}/controls.{reference}.somatics.vcf.gz.tbi"
    shell:
        """
		ln {input.vcf} {output.vcf}
		ln {input.index} {output.index}
		"""

localrules: generate_somatic_controls

rule extract_sample_somatics:
    input:
        "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.somatics.vcf.gz"
    output:
        vcf = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.{sample}.somatics.vcf.gz",
        index = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.{sample}.somatics.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools view -s {wildcards.sample} -Uac1 -Oz -o {output.vcf} {input} && tabix {output.vcf}"

localrules: extract_sample_somatics

rule vcfeval_somatics:
	input:
		reference = "data/references/{reference}.sdf",
		baseline_vcf = "data/truth/{project}/controls.{reference}.somatics.vcf.gz",
		baseline_vcf_index = "data/truth/{project}/controls.{reference}.somatics.vcf.gz.tbi",
		calls_vcf = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.{sample}.somatics.vcf.gz",
		calls_vcf_index = "results/calls/{project}/{group}.{reference}.{mapper}.{caller}.{sample}.somatics.vcf.gz.tbi"
	output:
		directory("results/eval/{project}/{group}.{reference}.{mapper}.{caller}.{sample}.{filter}.somatics.vcfeval")
	log:
		"logs/{project}/eval/{caller}/{group}.{reference}.{mapper}.{sample}.{filter}.somatics.vcfeval.log"
	params:
		score_field = lambda wildcards: config["caller_score_fields"][wildcards.caller],
		all_records = lambda wildcards: "--all-records" if wildcards.filter=="raw" else "",
	threads: 20
	resources:
		mem_gb=40
	conda:
		"../envs/rtg.yaml"
	shell:
		"(rtg \
			RTG_MEM={resources.mem_gb}g \
			vcfeval \
			-t {input.reference} \
			-b {input.baseline_vcf} \
			-c {input.calls_vcf} \
			-o {output} \
			--sample ALT \
            --squash-ploidy \
			-f {params.score_field} \
			--ref-overlap \
			--threads {threads} \
			{params.all_records} \
			)2> {log}"
