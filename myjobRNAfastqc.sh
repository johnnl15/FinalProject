#!/bin/bash

#BATCH --job-name=RNAfastq      ## Name of the job.
#SBATCH -A ginalee_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --array=1-42         ## number of tasks to launch, given hint below wc -l $file is helpful
#SBATCH --cpus-per-task=1    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load fastqc/0.11.9

file="FASTQC/prefixfastqc.txt"
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1` 

fastqc ${prefix}
