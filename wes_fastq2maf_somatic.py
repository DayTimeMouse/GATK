#!/usr/bin/env python
# coding: utf-8

shell.executable("/bin/bash")
#configfile: "./fastq/config.yaml"

(SAMPLES,) = glob_wildcards("fastq/{sample}_1.fastq.gz")
# READS=["1","2"]

# reference='/home/ljz/year1/GATK/reference/gatk_hg38'
hg38_vcf='/home/ljz/year1/GATK/GATK_demo/GATK_bundle'
fasta='/home/ljz/year1/amlwes/script/ref/Homo_sapiens_assembly38.fasta'
bed ='/home/ljz/year1/amlwes/script/bed/hg38.exon.bed'
vep_data='/home/ljz/year1/GATK/GATK_demo/vep/grch38'
vep_path='/home/ljz/miniconda3/envs/vep/bin'

rule vep_maf:
    input:
        expand("7_vep/somatic_merged_vep.maf")

# 质控：Quality Control
rule fastp:
    input:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz"
    output:
        "2_fastp/{sample}_clean.1.fq.gz",
        "2_fastp/{sample}_clean.2.fq.gz",
        "2_fastp/{sample}_fastp.html",
        "2_fastp/{sample}.json"
    log:
        "log/{sample}_fastp.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/rnaseq
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -j {output[3]} \
        -w 12 -z 4 -q 20 -u 30 -n 10 -h {output[2]} 2> {log}
        """

#比对：alignment
rule bwa:
    input:
        "2_fastp/{sample}_clean.1.fq.gz",
        "2_fastp/{sample}_clean.2.fq.gz"
    output:
        "3_bwa/{sample}.sorted.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    log:
        "log/{sample}_bwa_align.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/rnaseq
        bwa mem -t 5 -M -R '{params.rg}' {fasta} {input[0]} {input[1]} | samtools sort -@ 8 -m 20G -O bam -o {output} - 2> {log}
        """

# 标记重复序列（MarkDuplicates）
rule MarkDuplicates:
    input:
       "3_bwa/{sample}.sorted.bam"
    output:
       "4_markdup/{sample}.sorted.markdup.bam",
       "4_markdup/{sample}.markdup_metrics.txt",
       "4_markdup/{sample}.sorted.markdup.bam.bai"
    log:
        "log/{sample}_markdup.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/gatk4
        gatk MarkDuplicates -I {input} \
        -M {output[1]}  \
        -O {output[0]} 2> {log}
        samtools index {output[0]}
        """

# 碱基质量分数重校准（BaseRecalibrator）
rule BaseRecalibrator:
    input:
        "4_markdup/{sample}.sorted.markdup.bam",
        "4_markdup/{sample}.sorted.markdup.bam.bai"
    output:
        "5_bqsr/{sample}.recal_data.table"
    log:
        "log/{sample}_bqsr.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/gatk4
        gatk --java-options '-Xmx30g' BaseRecalibrator \
        -I {input[0]} -R {fasta} \
        --known-sites {hg38_vcf}/1000G_phase1.snps.high_confidence.hg38.vcf.gz\
        --known-sites {hg38_vcf}/dbsnp_146.hg38.vcf.gz \
        --known-sites {hg38_vcf}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -L {bed} \
        -O {output} 2> {log}
        """

# 应用碱基质量分数重校准（ApplyBQSR）
rule ApplyBQSR:
    input:
        "4_markdup/{sample}.sorted.markdup.bam",
        "5_bqsr/{sample}.recal_data.table"
    output:
        "5_bqsr/{sample}.sorted.markdup.BQSR.bam"
    log:
        "log/{sample}_ApplyBQSR.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/gatk4
        gatk --java-options '-Xmx30g' ApplyBQSR \
        --bqsr-recal-file {input[1]} \
        -R {fasta} \
        -L {bed} \
        -I {input[0]} -O {output} 2> {log}
        """

rule Mutect2:
    input:
        "5_bqsr/{sample}.sorted.markdup.BQSR.bam"
    output:
        "6_vcf/{sample}.somatic.vcf.gz"
    log:
        "log/{sample}_Mutect2.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/gatk4
        gatk --java-options "-Xmx40g" Mutect2 \
        -R {fasta} \
        --native-pair-hmm-threads 10 \
        --germline-resource {hg38_vcf}/af-only-gnomad.hg38.vcf.gz \
        --panel-of-normals {hg38_vcf}/somatic-hg38_1000g_pon.hg38.vcf.gz \
        -I {input} \
        -O {output} \
        -L {bed} 2> {log}
        """

rule FilterMutectCalls:
    input:
        "6_vcf/{sample}.somatic.vcf.gz"
    output:
        "6_vcf/{sample}.somatic.filter.vcf.gz"
    log:
        "log/{sample}_FilterMutectCalls.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/gatk4
        gatk --java-options "-Xmx40g" FilterMutectCalls \
        -R {fasta} \
        -V {input} \
        -O {output} \
        -L {bed} 2> {log}
        """

rule gunzip_vcf:
    input:
        "6_vcf/{sample}.somatic.filter.vcf.gz"
    output:
        "6_vcf/{sample}.somatic.filter.vcf"
    shell:
        """
        gunzip {input}
        """

rule pass_vcf:
    input:
        "6_vcf/{sample}.somatic.filter.vcf"
    output:
        "6_vcf/{sample}.somatic.PASS.filter.vcf"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/gatk4
        awk '{{if($7 == "PASS" || /^#/)print}}' {input} > {output}
        """
        
rule vep:
    input:
        "6_vcf/{sample}.somatic.PASS.filter.vcf"
    output:
        "7_vep/{sample}.somatic.vep.maf"
    log:
        "log/{sample}_vep.log"
    shell:
        """
        source activate /home/ljz/miniconda3/envs/vep
        vcf2maf.pl --input-vcf {input} \
        --output-maf {output} \
        --ref-fasta {fasta} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --tumor-id {wildcards.sample} \
        --ncbi-build GRCh38 2> {log}
        """

rule merge_maf:
    input:
        expand("7_vep/{sample}.somatic.vep.maf", sample=SAMPLES)
    output:
        "7_vep/somatic_merged_vep.maf",
        "7_vep/tmp",
        "7_vep/header"
    log:
        "log/merge_maf.log"
    shell:
        """
        cat {input} | grep -v '^#' | grep -v '^Hugo_Symbol' > {output[1]}
        grep 'Hugo_Symbol' $(ls -1 {input} | head -n 1) >  {output[2]}
        cat {output[2]} {output[1]} > {output[0]} 2> {log}
        """

