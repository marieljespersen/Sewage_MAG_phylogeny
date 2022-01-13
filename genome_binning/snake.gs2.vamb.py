### dir med contigs: /home/projects/cge/analysis/global_sewage_2/
### links to read files: /home/projects/cge/data/cge_private/FOOD_newSeq/upload/2016_1/Symbolic_links/Trimmed_data

# set programs
SPADES = "module load anaconda3/4.0.0; module load SPAdes/3.9.0; spades.py"
LOAD_SAMTOOLS = "module load samtools/1.9"
SAMTOOLS = "samtools"
LOAD_MINIMAP2 = "module load minimap2/2.6"
MINIMAP2 = "minimap2"
LOAD_JGI = "module unload R/3.6.1 gcc/8.2.0; module load perl/5.24.0 metabat/2.10.2"
JGI = "jgi_summarize_bam_contig_depths"
NPZ = "/home/projects/cge/people/maloj/src/snakemakes/vamb/ab_to_npz.py"
TNF = "/home/projects/cge/people/maloj/src/snakemakes/vamb/calc_tnf.py"
FILTER_CONTIG_LENGTH = "module load anaconda3/4.0.0; python /home/projects/cge/people/maloj/src/snakemakes/vamb/filter_contig_length.py"
RENAME_CONTIGS = '/home/projects/cge/people/maloj/src/snakemakes/vamb/rename_contigs.sh'

# Be in the project path directory when running script
PROJECT_PATH = "/home/projects/cge/data/projects/5001/"

# read ID's for all samples
IDS, = glob_wildcards("/home/projects/cge/data/projects/5001/maloj_symlinks/dot_links/{id}.R1.trim.fq.gz")

rule all:
   input:
       expand("assembly/{sample}/contigs.fasta", sample=IDS),
       expand("assembly/{sample}/{sample}.renamed.contigs.min2000.fa", sample=IDS),
       "vamb/cat.min2000.fa",
       "vamb/cat.min2000.mmi",
       expand("vamb/mapped/{sample}/{sample}.bam", sample=IDS),
       expand("vamb/mapped/abundance/{sample}.jgi", sample=IDS),
       #"vamb/mapped/abundance/jgi.column1to3",
       expand("vamb/mapped/abundance/{sample}.column4to5", sample=IDS),
       "vamb/mapped/abundance/jgi.abundance.dat",
       "vamb/mapped/abundance/jgi.abundance.npz",
       "vamb/tnfs.npz", #Gøres først på contog katalog med den rigtige min threshold
       "vamb/names.npz",
       "vamb/lengths.npz"
       #"vamb/bins/"


# Assembling trimmed reads with SPAdes
rule spades_denovo:
   input:
      forward="/home/projects/cge/data/projects/5001/maloj_symlinks/{sample,\w+}_R1.trim.fq.gz",
      reverse="/home/projects/cge/data/projects/5001/maloj_symlinks/{sample,\w+}_R2.trim.fq.gz"
   output:
      folder="assembly/{sample}/contigs.fasta"
   params:
      walltime="864000", nodes="1", ppn="12", mem="450gb", project="cge", mem_string="450"
   threads: int("12")
   log:
      "log/assembly/{sample,\w+}/spades.log"
   shell:
      "{SPADES} --meta -t {threads} -1 {input.forward} -2 {input.reverse} -o {output.folder} -m {params.mem_string} 2> {log}"
      
# filter contigs to minimum length 
rule contigs_filter:
   input:
      fasta="assembly/{sample,\w+}/contigs.fasta"
   output:
      "assembly/{sample,\w+}/{sample}.contigs.min2000.fa"
   params:
      walltime="86400", nodes="1", ppn="1", mem="1gb", project="cge", length=2000
   log:
      "log/assembly/{sample,\w+}/contigs.filter.log"
   shell:
      "{FILTER_CONTIG_LENGTH} --i {input.fasta} --min {params.length} --o {output} 2> {log}"

# rename contigs to include samplename in header
rule contigs_rename:
   input:
      "assembly/{sample,\w+}/{sample}.contigs.min2000.fa"
   output:
      "assembly/{sample,\w+}/{sample}.renamed.contigs.min2000.fa"
   params:
      walltime="86400", nodes="1", ppn="1", mem="115gb", project="cge"
   log:
      "log/assembly/{sample,\w+}/contigs.rename.log"
   shell:
      "{RENAME_CONTIGS} {input} {output} 2> {log}"

# cat contigs into one file
rule contigs_cat:
   input:
      expand("assembly/{sample}/{sample}.renamed.contigs.min2000.fa", sample=IDS)
   output:
      "vamb/cat.min2000.fa"
   params:
      walltime="432000", nodes="1", ppn="1", mem="115gb", project="cge"
   log:
      "log/assembly/contigs.cat.log"
   shell:
      """cat {input} > {output} 2> {log}"""
   
# index fasta containing contigs
rule contigs_index:
   input:
      "vamb/cat.min2000.fa"
   output:
      "vamb/cat.min2000.mmi"
   params:
      walltime="864000", nodes="1", ppn="1", mem="500gb", project="cge"
   log:
      "log/vamb/contigs.index.log"
   shell:
      "{MINIMAP2} -I 35G -d {output} {input} 2> {log}"

# map all samples to contig file
rule mapping:
   input:
      index="vamb/cat.min2000.mmi",
      assembly="vamb/cat.min2000.fa",
      forward="/home/projects/cge/data/projects/5001/maloj_symlinks/dot_links/{sample,\w+}.R1.trim.fq.gz",
      reverse="/home/projects/cge/data/projects/5001/maloj_symlinks/dot_links/{sample,\w+}.R2.trim.fq.gz"
   output:
      "vamb/mapped/{sample}/{sample}.bam"
   params:
      walltime="864000", nodes="1", ppn="20", mem="90gb", project="cge"
   threads: 20
   log:
      general="log/vamb/{sample,\w+}.mapping.log",
      time="log/vamb/{sample,\w+}.mapping.timelog" #For timing include: /usr/bin/time -v -o {log.time}    
   shell: 
      "{LOAD_SAMTOOLS}; {LOAD_MINIMAP2}; {MINIMAP2} -t {threads} -N 50 -ax sr {input.index} {input.forward} {input.reverse} | {SAMTOOLS} view -F 3584 -b --threads 20 -T {input.assembly} > {output} 2> {log.general}"


"""
# Sorting bam files prior to calculating abundance
rule sort_bam: 
   input:
      "vamb/mapped/{sample,\w+}.bam"
   output:
      "vamb/mapped/{sample,\w+}.sorted.bam"
   params: 
      walltime="86400", nodes="1", ppn="1", mem="100gb", project="cge" 
   log:
      "log/vamb/{sample,\w+}.sort_bam.log"
   shell: 
      "{LOAD_SAMTOOLS}; {SAMTOOLS} sort {input} > {output} 2> {log}"
"""      
# Calculating abundances prior to running vamb
rule abundance: 
   input:
      "vamb/mapped/{sample}/{sample,\w+}.bam"
   output:
      "vamb/mapped/abundance/{sample,\w+}.jgi"
   params: 
      walltime="86400", nodes="1", ppn="1", mem="10gb", project="cge" 
   log:
      general="log/vamb/{sample,\w+}.abundance.log",
      time="log/vamb/{sample,\w+}.abundance.timelog" #For timing include: /usr/bin/time -v -o {log.time}
   shell: 
      "{LOAD_JGI};  {JGI} --noIntraDepthVariance --outputDepth {output} {input} 2> {log.general}"

"""
# cutting out general columns from the first sample file
rule cut_column1to3: 
   input:
      "vamb/mapped/abundance/DTU_2019_1041_1_MG_CM1802_L98.jgi" #### Hardcode et prøvenavn ind! ######
   output:
      "vamb/mapped/abundance/jgi.column1to3"
   params:
      walltime="86400", nodes="1", ppn="1", mem="1gb", project="cge"
   shell: 
      "cut -f1-3 {input} > {output}"     
"""
# cutting out columns with abundance data from each sample file
rule cut_column4to5:
   input:
      "vamb/mapped/abundance/{sample,\w+}.jgi"
   output:
      "vamb/mapped/abundance/{sample,\w+}.column4to5" 
   params:
      walltime="86400", nodes="1", ppn="1", mem="10gb", project="cge"
   shell: 
      "cut -f1-3 --complement {input} > {output}"
"""
# Combining abundace data from all samples
rule paste_abundances:
   input:
      column_names="vamb/mapped/abundance/jgi.column1to3",
   output:
      "vamb/mapped/abundance/jgi.abundance.dat" 
   params:
      walltime="86400", nodes="1", ppn="1", mem="1gb", project="cge"
   shell: 
      "paste {input.column_names} vamb/mapped/abundance/*.column4to5 > {output}" ### looper ikke over samples, fordi det skal gøres samtidig. Bruger * i stedet.
"""
# Converting abundance files to npz file format
rule ab_npz:
   input:
      "vamb/mapped/abundance/jgi.abundance.dat"
   output:
      "vamb/mapped/abundance/jgi.abundance.npz" 
   params:
      walltime="86400", nodes="1", ppn="1", mem="750gb", project="cge"
   log:
      general="log/vamb/npz.log",
      time="log/vamb/ab_npz.timelog" #For timing include: /usr/bin/time -v -o {log.time}
   shell: 
      "{NPZ} --i {input} --o {output} 2> {log.general}"  

# Calculating TNF prior to running vamb           
rule tnf_names_length_npz:
   input:
      "vamb/cat.min2000.fa",
   output:
      tnf="vamb/tnfs.npz",
      names="vamb/names.npz",
      lengths="vamb/lengths.npz"
   params:
      walltime="86400", nodes="1", ppn="1", mem="100gb", project="cge"
   log:
      general="log/vamb/tnf.log"
   shell: 
      "{TNF} --i {input} --tnf {output.tnf} --names {output.names} --lengths {output.lengths} 2> {log}"  

rule vamb:
   input:
      abundance="vamb/mapped/abundance/jgi.abundance.npz",
      tnf="vamb/tnfs.npz",
      names="vamb/names.npz",
      lengths="vamb/lengths.npz",
   output:
      directory("vamb/bins/")
   params:
      walltime="864000", nodes="1", ppn="8", mem="115gb", project="cge"
   log:
      "log/vamb/vamb.log"
   shell: 
      "vamb -o _NODE_ -a 0.01 --rpkm {input.abundance} --tnfs {input.tnf} --names {input.names} --lengths {input.lengths} --outdir {output}" 
