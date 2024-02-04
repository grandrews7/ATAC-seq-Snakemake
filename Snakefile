import pandas as pd
import os

baseDir = "/data/zusers/andrewsg/ATACseq-Snakemake/"
fastqDir = "/data/zusers/sheddn/iCREs/Calderon_ATAC/fastqs/"
trimmedFastqDir = baseDir + "1-Trimmed-Fastq/"
QCdir = baseDir + "2-QC/"
BAMdir = baseDir + "3-BAM/"
peaksDir = baseDir + "4-Peaks/"
bigWigDir = baseDir + "5-bigWig/"
cCREScoresDir = baseDir + "6-cCRE-Scores/"

samples = [_.replace("_R1_001.fastq.gz", "") for _ in os.listdir(fastqDir) if "_R1_001.fastq.gz" in _]
samples = samples[:1]

rule all:
    input:
        expand(cCREScoresDir + "{sample}.npy", sample=samples)


rule trim:
    input:
        fastqDir + "{sample}_R1_001.fastq.gz",
        fastqDir + "{sample}_R2_001.fastq.gz"

    output:
        trimmedFastqDir + "{sample}.R1.P.fastq.gz",
        trimmedFastqDir + "{sample}.R2.P.fastq.gz"
    threads: 12
    log: stdout="/home/andrewsg/logs/{sample}.trim.out", stderr="/home/andrewsg/logs/{sample}.trim.err"
    resources:
        mem_mb=40000,
        c=12,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"

    shell:
        """
        MINWAIT=10
        MAXWAIT=120
        sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
        tmpDir=/tmp/{wildcards.sample}
        rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
        echo "Copying raw reads"
        cp {input[0]} .
        cp {input[1]} .
        singularity exec -B /data:/data /home/andrewsg/bin/atac.sif java -jar /software/trimmomatic-0.38.jar PE -threads {threads} -phred33 {wildcards.sample}_R1_001.fastq.gz {wildcards.sample}_R2_001.fastq.gz {wildcards.sample}.R1.P.fastq.gz {wildcards.sample}.R1.U.fastq.gz {wildcards.sample}.R2.P.fastq.gz  {wildcards.sample}.R2.U.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        cp {wildcards.sample}.R1.P.fastq.gz {output[0]}
        cp {wildcards.sample}.R2.P.fastq.gz {output[1]}
        cd 
        rm -rf $tmpDir
        """


rule map:
    input: 
        trimmedFastqDir + "{sample}.R1.P.fastq.gz",
        trimmedFastqDir + "{sample}.R2.P.fastq.gz"

    output:
        QCdir + "{sample}.flagstat",
        BAMdir + "{sample}.bam"
    threads: 12
    resources:
        mem_mb=40000,
        c=12,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"

    shell:
        """
        MINWAIT=10
        MAXWAIT=120
        sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
        tmpDir=/tmp/{wildcards.sample}
        rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
        echo "Copying trimmed reads"
        cp {input[0]} .
        cp {input[1]} .

        # Map with Bowtie2
        singularity exec /home/andrewsg/bin/atac.sif bowtie2 -p {threads} -x /home/andrewsg/genome/hg38/hg38.fa -1 {wildcards.sample}.R1.P.fastq.gz -2 {wildcards.sample}.R2.P.fastq.gz > {wildcards.sample}.sam

        # sort
        singularity exec /home/andrewsg/bin/bioinformatics.sif samtools view -b -@ {threads} {wildcards.sample}.sam | singularity exec /home/andrewsg/bin/bioinformatics.sif samtools sort -@ {threads} -o {wildcards.sample}.bam -

        #Calculating mapping statistics
        singularity exec /home/andrewsg/bin/bioinformatics.sif samtools flagstat -@ {threads} {wildcards.sample}.bam > {wildcards.sample}.flagstat

        # Filter alignment and convert to BAM
        singularity exec /home/andrewsg/bin/bioinformatics.sif samtools view -@ {threads} -q 30 -F4 -O SAM {wildcards.sample}.bam | egrep -v chrM | singularity exec /home/andrewsg/bin/bioinformatics.sif samtools view -b -@ {threads} -o {wildcards.sample}.filtered.bam -T /home/andrewsg/genome/hg38/hg38.fa -

        cp {wildcards.sample}.flagstat {output[0]}
        cp {wildcards.sample}.filtered.bam {output[1]}
        cd
        rm -rf $tmpDir
        """

rule callpeaks:
    input: BAMdir + "{sample}.bam"
    output: peaksDir + "{sample}_peaks.narrowPeak"
    threads: 12
    resources:
        mem_mb=40000,
        c=12,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    shell:
        """
        MINWAIT=10
        MAXWAIT=120
        sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
        tmpDir=/tmp/{wildcards.sample}
        rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
        echo "Copying BAM"
        cp {input[0]} .
        singularity exec /home/andrewsg/bin/macs2.sif macs2 callpeak -t {wildcards.sample}.bam  --nolambda --nomodel -g hs --keep-dup all --call-summits -n {wildcards.sample}
        cp {wildcards.sample}_peaks.narrowPeak {output[0]}
        """

rule bamCoverage:
    input: BAMdir + "{sample}.bam"
    output: bigWigDir + "{sample}.bigWig"
    threads: 24
    resources:
        mem_mb=90000,
        c=24,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"
    shell:
        """
        MINWAIT=10
        MAXWAIT=120
        sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
        tmpDir=/tmp/{wildcards.sample}
        rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
        echo "Copying BAM"
        cp {input[0]} .
        echo "Indexing BAM"
        singularity exec /home/andrewsg/bin/bioinformatics.sif samtools index -@ {threads} {wildcards.sample}.bam
        echo "Obtaining signal track"
        singularity exec /home/andrewsg/bin/bioinformatics.sif bamCoverage -b {wildcards.sample}.bam -o {wildcards.sample}.bigWig --extendReads --ignoreForNormalization chrX -bs 10 -p {threads} --normalizeUsing RPKM
        cp {wildcards.sample}.bigWig {output[0]}
        """

rule scorecCREs:
    input: bigWigDir + "{sample}.bigWig"
    output: cCREScoresDir + "{sample}.npy"
    threads: 12
    resources:
        mem_mb=40000,
        c=12,
        runtime=240,
        nodes=1,
        slurm_partition="4hours"

    shell:
        """
        MINWAIT=10
        MAXWAIT=120
        sleep $((MINWAIT+RANDOM % (MAXWAIT-MINWAIT)))
        tmpDir=/tmp/{wildcards.sample}
        rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
        echo "Copying bigWig"
        cp {input[0]} .
        cp /data/zusers/andrewsg/ATACseq-Snakemake/getSignal.py .
        cp /data/projects/encode/Registry/V4/GRCh38/GRCh38-cCREs.bed .
        echo "Scoring cCREs"
        singularity exec /home/andrewsg/bin/bioinformatics.sif python3 getSignal.py GRCh38-cCREs.bed {wildcards.sample}.bigWig {wildcards.sample}.npy
        cp {wildcards.sample}.npy {output[0]}
        cd
        rm -rf $tmpDir
        """