# =========================================================
# cellranger.snakefile
# pipeline for paired scRNA-seq + TCR-seq (10x) analysis
# Joy Pai
# Satpathy Lab, Stanford University
# May 20, 2019
# =========================================================

import glob
import os

# SET UP
CWD = os.getcwd()
CELLRANGER_DIR='/storage/joypai/software/cellranger-3.0.2'

EXP_samples = [ s.split('/')[-1] for s in glob.glob(CWD+"/data/*-WTA") ]
TCR_samples = [ s.split('/')[-1] for s in glob.glob(CWD+"/data/*-TCR") ]

print("exp samples:", EXP_samples)
print("tcr samples:", TCR_samples)


# RULES
rule all:
    input:
        h5_files=expand('{sample}/outs/molecule_info.h5', sample=EXP_samples),
        clone_files=expand('{sample}/outs/clonotypes.csv', sample=TCR_samples),
        aggr=CWD+"/runs/MS/outs/aggregation.csv"

rule count:
    input: 
        fastq_dir=CWD+'/data/{sample}'
    params:
        sample_name='{sample}',
        ref=CELLRANGER_DIR+'/refdata-cellranger-GRCh38-3.0.0'
    output: '{sample}/outs/molecule_info.h5'
    log: 'cellranger_count_{sample}.log'
    threads: 10
    shell:
        """
        cellranger count --id={params.sample_name}
                   --transcriptome={params.ref}
                   --fastqs={input.fastq_dir}
                   --sample={params.sample_name}
                   --localcores={threads} &> {log}
        """

rule aggr:
    input: expand("{sample}/outs/molecule_info.h5", sample=EXP_samples)
    params:
        aggr_csv="MS_libraries.csv",
        run_id="MS"
    output: CWD+"/runs/MS/outs/aggregation.csv"
    log: 'cellranger_aggr.log'
    threads: 1
    run:
        # write aggregation csv file
        with open({params.aggr_csv}, "w") as outf:
            outf.write("library_id,molecule_h5\n")

            for i in input:
                sample = i.split('/')[0]
                outf.write(sample+","+i+"\n")
        
        # run cellranger aggr
        shell("cellranger aggr --id={params.run_id} \
                --csv={params.aggr_csv} \
                --normalize=mapped &> {log}")

rule vdj:
    input:
        fastq_dir=CWD+'/data/{sample}'
    params:
        sample_name='{sample}',
        vdj_ref=CELLRANGER_DIR+'/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0'
    output: '{sample}/outs/clonotypes.csv'
    log: 'cellranger_vdj_{sample}.log'
    threads: 10
    shell:
        """
        cellranger vdj --id={params.sample_name} \
                 --reference={params.vdj_ref} \
                 --fastqs={input.fastq_dir} \
                 --sample={params.sample_name} \
                 --localcores={threads} &> {log}
        """