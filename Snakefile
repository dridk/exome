import os
from os.path import basename

configfile : "config.yaml"


rule all:
	input:	config["OUTPUT"]



# # ===== CREATE REFERENCE DICT 
# rule create_reference_dict:
# 	input:
# 		config["BWA_GENOM_REF"]
# 	output:
# 		"/DATAS/test.txt"



# 	shell:
# 		"picard CreateSequenceDictionary R={input} O={output}" 


# ================ RULE clean 
rule clean:
	input:
		forward =   config["RAW_FOLDER"] + "/{sample}_1.fastq.gz", 
		reverse =   config["RAW_FOLDER"] + "/{sample}_2.fastq.gz"
	output:
		forward =   "{sample}_1.clean.fastq.gz", 
		reverse =   "{sample}_2.clean.fastq.gz",
		nullf   =   temp("{sample}_1.null.fastq.gz"),
		nullr   =   temp("{sample}_2.null.fastq.gz")

	threads:
		128
		
	shell:
		"trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} {output.forward} {output.nullf} {output.reverse} {output.nullr} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:125"

# =============== RULE ALIGNEMENT 
# READ GROUP SHOULD BE REPLACED BY LANE 
rule alignement : 
	input:
		forward =   "{sample}_1.clean.fastq.gz", 
		reverse =   "{sample}_2.clean.fastq.gz"
	output:
		temp("{sample}.sam")
	threads: 100 
	shell:
		"bwa mem -M -R '@RG\tID:foo\tSM:bar' -t {threads} {config[BWA_GENOM_REF]} {input.forward} {input.reverse} > {output}"
		


rule sam_to_bam:
	input:
		"{filename}.sam"
	output:
		temp("{filename}.bam")
	threads: 100 
	shell:
		"sambamba view --sam-input -t {threads} -f bam -o {output} {input}"

rule sort_bam:
	input:
		"{filename}.bam"
	output:
		temp("{filename}.sorted.bam")
	threads: 100 
	shell:
 		"sambamba sort -t {threads} {input} -o {output}" 

# Picard alternative 
rule mark_duplicate:
	input:
		"{sample}.sorted.bam"
	output:
		"{sample}.sorted.dup.bam"
	threads: 100 
	shell:
		"sambamba markdup -t {threads} {input} {output}"

# BQSR 

rule bqsr_analyse_covariance:
	input:
		"{sample}.sorted.bam"
	output:
		"{sample}.recal_data.table"
	shell:
		"gatk -T BaseRecalibrator -R {config[BWA_GENOM_REF]} -I {input} -L 20  -knownSites {config[DB_SNP]} -o {output}"

# # BQSR recalibration 
# rule recalibration:
# 	input:
# 		"{sample}.sorted.dup.bam"
# 	output:
# 		""


# mpileup 
rule mpilup:
	input:
		"{sample}.sorted.dup.bam"
	output:
		temp("{sample}.mpilup")
	shell:
		"samtools mpileup -f {config[BWA_GENOM_REF]} {input} > {output}"

# Run varscan snp 
rule varscan_snp:
	input:
		"{sample}.mpilup"
	output:
		"{sample}.varscan.snp.vcf"
	shell:
		"varscan mpileup2snp {input} --output-vcf > {output}"

# Run varscan indel 
rule varscan_indel:
	input:
		"{sample}.mpilup"
	output:
		"{sample}.varscan.indel.vcf"
	shell:
		"varscan mpileup2indel {input} --output-vcf > {output}"

 # Run platypus varcalling 
rule platypus_calling:
	input:
		"{sample}.sorted.dup.bam"
	output:
		"{sample}.platypus.vcf"
	threads: 
		128
	shell:
		"platypus callVariants --bamFiles={input} --output {output} --refFile={config[BWA_GENOM_REF]} --nCPU={threads}"	


rule bgzip_tabix:
	input:
		"{filename}.vcf"
	output:
		"{filename}.vcf.gz"
	shell:
		"bgzip {input}; tabix {output}"
