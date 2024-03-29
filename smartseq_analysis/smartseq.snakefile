# =========================================================
# smartseq.snakefile
# pipeline for SMART-seq data analysis 
# Joy Pai
# Satpathy Lab, Stanford University
# May 20, 2019
# =========================================================

import glob
import os

# SET UP
CWD = os.getcwd()

RSEM='/seq/RSEM-1.3.0/' 			# path to RSEM 
STAR='/seq/STAR-2.5.2b/bin/Linux_x86_64/'	# path to STAR
FASTQC='/seq/fastqc-11.2/FastQC/fastqc'		# path to STAR

HUMAN_REF_FASTA='/seq/genomes/hg38/hg38.fa'
HUMAN_REF_GTF='/seq/genomes/hg38/hg38.gtf'

SAMPLES = [ s.split('/')[-1].split('_R1')[0] for s in glob.glob(CWD+"/*/*_R1_001.fastq.gz") ]
#L003_SAMPLES = [ s.split('/')[-1].split('_R1')[0] for s in glob.glob(CWD+"/*/*L003_R1_001.fastq.gz") ]


# RULES
rule all:
    input:
        rsem_files=expand('rsem/{sample}.genes.results', sample=SAMPLES),
        #fastqc_files=expand(CWD+'/fastqc/{sample}_R1_001_fastqc.html', sample=SAMPLES)

rule fastqc:
    input:
        fastq_r1=CWD+'/data/{sample}_R1_001.fastq.gz',
        fastq_r2=CWD+'/data/{sample}_R2_001.fastq.gz'
    params:
        fastqc=FastQC,
        outdir=CWD+'/fastqc'
    output: 
        CWD+'/fastqc/{sample}_R1_001_fastqc.html'
    log: 'fastqc/{sample}_fastqc.log'
    shell:
        """
        {params.fastqc} -f fastq {input.fastq_r1} {input.fastq_r2} --output {params.outdir} &> {log}
        """

rule build_ref:
    input: 
        fasta=HUMAN_REF_FASTA
    params:
        gtf=HUMAN_REF_GTF
        rsem=RSEM,
        star=STAR,
        refname='ref/hg38',
        threads=40,
    output: 'ref/hg38.idx.fa'
    # log: 'logs/star_ref.log'
    resources:
        mem_mb=100
    shell:
        """
        {params.rsem}/rsem-prepare-reference --gtf {params.gtf} \
                        --star \
                        --star-path {params.star} \
                        -p {params.threads} \
                        {input.fasta} {params.refname}
        """

rule quant:
    input: 
        fastq_r1=CWD+'/data/{sample}_R1_001.fastq.gz',
        fastq_r2=CWD+'/data/{sample}_R2_001.fastq.gz',
    params:
        ref=CWD+'/ref/hg38',
        rsem=RSEM,
        star=STAR,
        prefix='rsem/{sample}',
        threads=10,
    output: 'rsem/{sample}.genes.results'
    log: 'logs/{sample}.log'
    threads: 10
    resources:
        mem_mb=100
    shell:
        """
        {params.rsem}/rsem-calculate-expression -p {params.threads} --paired-end \
					-star --star-path {params.star} \
					--star-gzipped-read-file --estimate-rspd \
					--append-names \
					--output-genome-bam \
					{input.fastq_r1} {input.fastq_r2} \
					{params.ref} {params.prefix} &> {log}

        """
