# ==================================================================
# cellranger pipeline for paired scRNA-seq + TCR-seq (10x) analysis
# Joy Pai
# Satpathy Lab, Stanford University
# May 20, 2019
# ==================================================================

import glob
import os

# SET UP
CWD = os.curdir()

EXP_samples = glob.glob(CWD+"data/*-WTA")
TCR_samples = glob.glob(CWD+"data/*-TCR")


# RULES
rule all:
    input:
        bam_files=expand('{sample}/outs/possorted_genome.bam.bam', sample=EXP_samples)

rule count:
    input: 
        fastq_dir='{sample}'
    params:
        sample_name='{sample}',
        ref='/storage/joypai/software/..fillin../refdata-cellranger-GRCh38-3.0.0'
    output:
        '{params.sample_name}/outs/possorted_genome_bam.bam'
    log: 'cellranger_count_{params.sample_name}'
    threads: 10
    shell:
        """
        cellranger count --id={params.sample_name} \
                   --transcriptome={params.ref} \
                   --fastqs={input.fastq_dir} \
                   --sample={params.sample_name} \
                   --localcores={threads} &> {log}
        """

rule aggr:
    input: 
    params: 
        aggr_csv="MS_libraries.csv",
        cur_dir=CWD
    output: CWD+"/runs/MS/outs/aggregation.csv"
    run:
        with open({params.aggr_csv}, "w) as outf:
            outf.write("library_id,molecule_h5\n")

            for i in input:
                outf.write(i+",{params.cur_dir}/"+i+"/outs/molecule_info.h5\n")
        
        # run cellranger aggr
        import subprocess

        aggr_command = "cellranger aggr --id=MS \
                            --csv=MS_libraries.csv \
                            --normalize=mapped"

        process = subprocess.Popen(aggr_command, stdout=subprocess.PIPE)
        output, error = process.communicate()

rule vdj:
    input:
        fastq_dir='{sample}'
    params:
        vdj_ref='/storage/joypai/software/..fillin.../refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0'
    log: 'cellranger_vdj_{params.sample_name}.log'
    threads: 10
    shell:
        """
        cellranger vdj --id=sample345 \
                 --reference={params.vdj_ref} \
                 --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \
                 --sample=mysample \
                 --localcores={threads} &> {log}

        """