#!/bin/bash

#BATCH --job-name=RNAalign      ## Name of the job.
#SBATCH -A ginalee_lab         ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --array=1-21         ## number of tasks to launch, given hint below wc -l $file is helpful
#SBATCH --cpus-per-task=2    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load samtools/1.10
module load hisat2/2.2.1

# or pass the file name to the shell script, how would I do this?
file="Bam/prefixesRNAalign.txt"
# is the file indexed for bwa?
ref="ref/GCA_000001405.15_GRCh38_full_analysis_set"
# here is a hint if you had a tab delimited input file
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1` 

hisat2 -p 2 -x $ref -1 ${prefix}_1.fq.gz -2 ${prefix}_2.fq.gz | samtools view -bS > Bam/${prefix}.bam
samtools sort Bam/${prefix}.bam -o Bam/${prefix}.sorted.bam
samtools index Bam/${prefix}.sorted.bam
