
### snakefile to align rnaseq reads in parallel from cells treated with jq1/dmso/naive and infected/mock
#   with sars2

SAMPLES, = glob_wildcards("/path/to/fastq_files/{id}_R1_001.fastq.gz")


rule run_all:
  input: expand("aligned/bams/{sample}.bam", sample=SAMPLES)
    

rule align_star:
    """
    Align sequencing reads using STAR.
    """
    input: fq1 = 'fastq/{sample}_R1_001.fastq.gz',
           fq2 = 'fastq/{sample}_R2_001.fastq.gz',
           ref = 'hg38scov2_index_star'
    output: bam = 'aligned/bams/{sample}.bam'  
    params:
      prefix = '{sample}' 
    threads: 8
    run:
        res_dir = 'aligned'
        
        shell("""
              #mkdir aligned/{params.prefix}
              STAR --runThreadN {threads}                 \
                   --runMode alignReads                   \
                   --outFileNamePrefix {res_dir}/{params.prefix}         \
                   --genomeDir {input.ref}  \
                   --readFilesIn {input.fq1} {input.fq2}  \
                   --readFilesCommand zcat  \
                   --outFilterType BySJout  \
                   --outFilterMultimapNmax 20   \
                   --alignSJoverhangMin 8   \
                   --outSJfilterOverhangMin 12 12 12 12   \
                   --outSJfilterCountUniqueMin 1 1 1 1  \
                   --outSJfilterCountTotalMin 1 1 1 1   \
                   --outSJfilterDistToOtherSJmin 0 0 0 0  \
                   --outFilterMismatchNmax 999  \
                   --outFilterMismatchNoverReadLmax 0.04  \
                   --scoreGapNoncan -4  \
                   --scoreGapATAC -4  \
                   --chimOutType WithinBAM HardClip   \
                   --chimScoreJunctionNonGTAG 0   \
                   --alignSJstitchMismatchNmax -1 -1 -1 -1  \
                   --alignIntronMin 20  \
                   --alignIntronMax 1000000   \
                   --alignMatesGapMax 1000000   \
                   --outSAMunmapped Within  \
                   --genomeLoad NoSharedMemory  \
                   --outSAMattributes NH HI AS NM MD  \
                   --outSAMtype BAM Unsorted  
              mv {res_dir}/{params.prefix}Aligned.out.bam {output.bam}    
              """)
