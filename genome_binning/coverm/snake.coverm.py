#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 13:52:17 2021

A snakemake for processing of bins after binning of metagenomic reads
To be developed further, currently only for abundance caculations

Run this snakemake using:
# Running snakemake pipeline
snakemake --snakefile /home/projects/cge/people/maloj/src/binprocessing/abundance/snake.coverm.py --cluster "qsub -l walltime={params.walltime} -l nodes=1:ppn={params.ppn} -l mem={params.mem} -A cge -W group_list=cge" --jobs 120 --latency-wait 60 --rerun-incomplete -k 2> snake.run && echo "I finished my run" | mailx -r maloj@food.dtu.dk -s "Snake-Abundance" maloj@food.dtu.dk

# Dryrun
snakemake --snakefile /home/projects/cge/people/maloj/src/binprocessing/abundance/snake.coverm.py --cluster "qsub -l walltime={params.walltime} -l nodes=1:ppn={params.ppn} -l mem={params.mem} -A cge -W group_list=cge" --jobs 120 --latency-wait 60 -n > snake.dryrun

@author: maloj
"""

CoverM_load = "module load coverm/0.6.1 minimap2/2.17r941 samtools/1.9 bwa-mem2/2.0"
coverage_calculation="/home/projects/cge/people/maloj/src/binprocessing/abundance/expected_coverage.R"
extract_abundances="/home/projects/cge/people/maloj/src/binprocessing/abundance/extract_rpkm_ab.R"

# Be in the project path directory when running script
PROJECT_PATH = "/home/projects/cge/data/projects/5001/Binning_vamb"

# read ID's for all samples
IDS, = glob_wildcards("/home/projects/cge/data/projects/5001/Binning_vamb/abundance/combined_abundance/all_reads/{id}.R1.trim.fq.gz") 
#IDS, = glob_wildcards("/home/projects/cge/data/projects/5001/maloj_symlinks/catted_files/{id}.R1.trim.fq.gz") 
rule all:
   input:
       expand("coverm/coverm.{sample}.tsv", sample=IDS),
       expand("coverm/coverage/{sample}.coverage.tsv", sample=IDS),
       expand("coverm/coverage/TPM/{sample}.TPM.tsv", sample=IDS),
       expand("coverm/coverage/read_count/{sample}.read_count.tsv", sample=IDS)
       
 
#  Calculating abundance using CoverM
rule CoverM:
   input:
      forward="/home/projects/cge/data/projects/5001/Binning_vamb/abundance/combined_abundance/all_reads/{sample,\w+}.R1.trim.fq.gz", #"/home/projects/cge/data/projects/5001/maloj_symlinks/catted_files/{sample,\w+}.R1.trim.fq.gz", #
      reverse="/home/projects/cge/data/projects/5001/Binning_vamb/abundance/combined_abundance/all_reads/{sample,\w+}.R2.trim.fq.gz", #"/home/projects/cge/data/projects/5001/maloj_symlinks/catted_files/{sample,\w+}.R2.trim.fq.gz", #
      NC_fastas="/home/projects/cge/data/projects/5001/Binning_vamb/abundance/combined_abundance/all_genomes" #"/home/projects/cge/data/projects/5001/Binning_vamb/best_genomes" #
   output:
      abundance_file="coverm/coverm.{sample,\w+}.tsv"
   params:
      walltime="864000", nodes="1", ppn="20", mem="90gb", project="cge", mem_string="90"
   threads: int("20")
   log:
      "log/coverm/{sample,\w+}.coverm.log"
   shell:
      "{CoverM_load}; coverm genome -t {threads} -d {input.NC_fastas} -1 {input.forward} -2 {input.reverse} -m relative_abundance rpkm tpm count covered_fraction length --min-covered-fraction 0 > {output.abundance_file} 2> {log}"
    
#  Calculating observed/expected coverage ratio
rule coverage:
   input:
      forward="/home/projects/cge/data/projects/5001/Binning_vamb/abundance/combined_abundance/all_reads/{sample,\w+}.R1.trim.fq.gz", #"/home/projects/cge/data/projects/5001/maloj_symlinks/catted_files/{sample,\w+}.R1.trim.fq.gz", #
      abundance_file="coverm/coverm.{sample,\w+}.tsv"
   output:
      coverage_file="coverm/coverage/{sample,\w+}.coverage.tsv"
   params:
      walltime="864000", nodes="1", ppn="1", mem="4gb", project="cge", mem_string="18"
   log:
      "log/coverage/{sample,\w+}.coverage.log"
   shell:
      "Rscript --vanilla {coverage_calculation} {input.forward} {input.abundance_file} {output.coverage_file} 2> {log}"

#  extracting rpkm and relative abundances based on expected coverage ratio (>=0.5)
rule extract_ab:
   input:
      coverage_file="coverm/coverage/{sample,\w+}.coverage.tsv"
   output:
      TPM_file="coverm/coverage/TPM/{sample,\w+}.TPM.tsv",
      read_count_file="coverm/coverage/read_count/{sample,\w+}.read_count.tsv"
   params:
      walltime="864000", nodes="1", ppn="1", mem="4gb", project="cge", mem_string="18"
   log:
      "log/coverage/{sample,\w+}.extract_ab.log"
   shell:
      "Rscript --vanilla {extract_abundances} {input.coverage_file} {output.TPM_file} {output.read_count_file} 2> {log}"
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      