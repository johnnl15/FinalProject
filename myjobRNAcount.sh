#!/bin/bash
#SBATCH --job-name=RNAseq      ## Name of the job.
#SBATCH -A ginalee_lab          ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --array=1         ## number of tasks to launch, given hint below wc -l $file is helpful
#SBATCH --cpus-per-task=2    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load subread/2.0.1

ref="ref/GCA_000001405.15_GRCh38_full_analysis_set"

gtf="ref/Homo_sapiens.GRCh38.84.gtf"
myfile=`cat Bam/prefixesRNAalign.txt| sed 's:^:Bam/:'|sed 's/$/.bam/'|tr "\n" " "`
featureCounts -p -T 4 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o Bam/human_counts_RNA.txt $myfile