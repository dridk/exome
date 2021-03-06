import os
import glob

#configfile : "config.yml"



# # ===== CREATE REFERENCE DICT 
# rule create_reference_dict:
# 	input:
# 		config["BWA_GENOM_REF"]
# 	output:
# 		"/DATAS/test.txt"



# 	shell:
# 		"picard CreateSequenceDictionary R={input} O={output}" 



# rule create_samplename:
# 	input:
# 		config["ILLUMINA_SAMPLESHEET"]
# 	output:
# 		"samples.txt"
# 	shell:
# 		"cat {input}|sed '1,/Data/d'|cut -f 1 -d','|tail -n +2 > {output}"


# rule all:
# 	input:
# 		dynamic("fastq/{sample}_{ss}_{lane}_{sens}_{end}.fastq.gz")



# rule bcl2fastq: 
# #bcl2fastq -R 170522_NB501647_0005_AHMHGVBGX2/ -o 170522_NB501647_0005_AHMHGVBGX2/output --sample-sheet 170522_NB501647_0005_AHMHGVBGX2/SampleSheet20170522MedExome.csv --create-fastq-for-index-reads 
# #CGH1678_S5_L004_R2_001.fastq.gz 
# 	input:
# 		folder = config["ILLUMINA_FOLDER"],
# 		config = config["ILLUMINA_SAMPLESHEET"]
# 	output:
# 		dynamic("fastq/{sample}_{ss}_{lane}_{sens}_{end}.fastq.gz")
# 	shell:
# 		"bcl2fastq -R {input.folder} -o fastq/ --sample-sheet {input.config} --create-fastq-for-index-reads"



rule all:
	input:
		expand("{sample}.sample.vcf.gz", sample = config["SAMPLES"])



# =============== MERGE LANES 
# rule mergelane:
# 	input: 
# 		lambda wildcards : sorted(glob.glob(config["RAW_FOLDER"]+"/"+wildcards.sample+"*_*_"+wildcards.sens+"*.fastq.gz"))
# 	output:
# 		config["RAW_FOLDER"]+"/{sample}_{sens}.all.fastq.gz"
# 	shell:
# 		"zcat {input}| gzip > {output}"


# ================ RULE clean 
rule clean:
	input:
		forward =   config["RAW_FOLDER"]+"/{sample}_R1.all.fastq.gz", 
		reverse =   config["RAW_FOLDER"]+"/{sample}_R2.all.fastq.gz"
	output:
		forward =   temp("{sample}_R1.clean.fastq.gz"), 
		reverse =   temp("{sample}_R2.clean.fastq.gz"),
		nullf   =   temp("{sample}_R1.null.fastq.gz"),
		nullr   =   temp("{sample}_R2.null.fastq.gz")

	threads:
		128
		
	shell:
		"trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} {output.forward} {output.nullf} {output.reverse} {output.nullr} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:125"

# =============== RULE ALIGNEMENT 
# READ GROUP SHOULD BE REPLACED BY LANE 
# CGH1678_S5_L001_I1_001.fastq.gz 
rule alignement : 
	input:
		forward = lambda wildcards : glob.glob(config["RAW_FOLDER"]+"/"+wildcards.sample+"_*_L00"+wildcards.lane+"_R1_001.fastq.gz"),
		reverse = lambda wildcards : glob.glob(config["RAW_FOLDER"]+"/"+wildcards.sample+"_*_L00"+wildcards.lane+"_R2_001.fastq.gz")
	output:
		temp("{sample}.L00{lane}.sam")
	threads: 4 
	shell:
		"bwa mem -M -R '@RG\tID:{wildcards.sample}.{wildcards.lane}\tLB:medexome\tSM:{wildcards.sample}\tPL:ILLUMINA' -t {threads} {config[BWA_GENOM_REF]} {input.forward} {input.reverse} > {output}"
		


rule mergebam:
	input:
		lambda wildcards: expand("{sample}.L00{num}.bam", sample = wildcards.sample, num = range(1,5))
	output:
		temp("{sample}.unsorted.bam")
	threads: 100
	shell:
		"samtools merge -@ {threads} {output} {input}"


rule sam_to_bam:
	input:
		"{filename}.sam"
	output:
		temp("{filename}.bam")
	threads: 100 
	shell:
		"sambamba view -h --sam-input -t {threads} -f bam -o {output} {input}"


rule sort_bam:
	input:
		"{sample}.unsorted.bam"
	output:
		temp("{sample}.sorted.bam")
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

# Run varscan snp 
rule varscan_snp:
	input:
		"{sample}.sorted.dup.bam"
	output:
		temp("{sample}.varscan.snp.vcf")
	shell:
		"samtools mpileup -f {config[BWA_GENOM_REF]} {input}|varscan pileup2snp --min-coverage 2 --output-vcf > {output}"

# Run samtools var calling 
rule samtools_calling:
	input:
		"{sample}.sorted.dup.bam"
	output:
		temp("{sample}.samtools.vcf")
	shell:
		"samtools mpileup -ugf {config[BWA_GENOM_REF]} {input}|bcftools call -vmO v -o {output}"

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

# Run freebayes varcalling:
rule freebayes_calling:
	input:
		"{sample}.sorted.dup.bam"
	output:
		temp("{sample}.freebayes.vcf")
	shell:
		"freebayes -f {config[BWA_GENOM_REF]} {input} > {output}"



# Run freebayes varcalling:
rule freebayes_global_calling:
	input:
		expand("{sample}.sorted.dup.bam", sample = config["SAMPLES"])
	output:
		temp("all.freebayes.vcf")
	shell:
		"freebayes -f {config[BWA_GENOM_REF]} {input} > {output}"


rule bgzip_tabix:
	input:
		"{filename}.vcf"
	output:
		"{filename}.vcf.gz"
	shell:
		"bgzip {input}; tabix {output}"


rule normalization:
	input:
		"{name}.vcf"
	output:
		"{name}.norm.vcf"
	shell:
		"vt normalize {input} -r {config[BWA_GENOM_REF]} -o {output}"



# Annotation ===============
# SnpEff 
rule snpEff : 
	input:
		"{filename}.vcf"
	output:
		temp("{filename}.snpEff.vcf")
	shell:
		"snpEff -Xmx4g -c {config[SNPEFF_CONFIG]} -v hg19 {input} > {output}"


rule anno_dbSNP : 
	input:
		"{filename}.snpEff.vcf"
	output:
		temp("{filename}.dbsnp.vcf")
	shell:
		"SnpSift -Xmx4g annotate -dbsnp -db {config[DB_SNP]} {input} > {output}"

rule anno_clinvar: 
	input:
		"{filename}.dbsnp.vcf"
	output:
		temp("{filename}.clinvar.vcf")
	shell:
		"SnpSift -Xmx4g annotate -clinvar -db {config[DB_CLINVAR]} {input} > {output}"

rule anno_dbNSFP: 
	input:
		"{filename}.clinvar.vcf"
	output:
		"{filename}.ann.vcf"
	shell:
		"SnpSift -Xmx4g dbnsfp -db {config[DB_NSFP]} -v {input} > {output}"



# Annotation ===============


rule vcf_filter: 
	input:
		"{filename}.vcf.gz"
	output:
		temp("{filename}.filter.vcf")
	shell:
		"bcftools filter -i 'QUAL >{config[VCF_QUAL]}' -R {config[TARGET_BED_FILE]} -Ov {input} > {output}"


# Split global VCF by family according pedigree file 
rule split_vcf_by_family:
	input:
		vcf = "all.freebayes.norm.filter.ann.vcf.gz",
		ped = config["PEDIGREE_FILE"]
	output:
		"{family}.family.vcf"
	shell:
		"SAMPLES=$(cat {input.ped}|grep -v '#'|awk '$1==\"{wildcards.family}\"{{print $2}}'|paste -s -d',');"
		"echo $SAMPLES found in {wildcards.family};"
		"bcftools view -c1 -Ov -s ${{SAMPLES}} -o {output} {input.vcf}"

# Split global VCF by Sample name according pedigree file 	
rule split_vcf_by_sample:
	input:
		vcf = "all.freebayes.norm.filter.ann.vcf.gz"
	output:
		"{sample}.sample.vcf"
	shell:
		"bcftools view -c1 -Ov -s {wildcards.sample} -o {output} {input};"

rule vcf_to_csv:
	input:
		"{filename}.vcf.gz"
	output:
		"{filename}.tsv"
	shell: 
		"SnpSift extractFields TRIO1.family.vcf.gz CHROM POS REF ALT ID "
		"dbNSFP_ExAC_AF dbNSFP_1000Gp1_AF "
		"dbNSFP_MutationTaster_pred "
		"dbNSFP_SIFT_pred "
		"dbNSFP_MetaSVM_pred "
		"MUT "
		"CLNSRC "
		"CLNDSDB "
		"CLNDSDBID "
		" GEN[*].GT "
		" ANN[0].EFFECT "
		"CLNHGVS > {output}" 
