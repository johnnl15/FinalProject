My goal for this final project was to analyze my in-house bulk RNA seq data which meant I started from the beginnining with fastq files and ended with a volcano plot. 
First creating prefixes for fastqc runs. 

```
ls *.fq.gz > FASTQC/prefixfastqc.txt

```
Then I ran the following sbatch [myjobRNAfastqc.sh](myjobRNAfastqc.sh) 

Then created prefixes for my RNAalignment

```
ls *1.fq.gz | sed 's/_1.fq.gz//' > Bam/prefixesRNAalign.txt
```

For the reference genome, I first tried to download the fna.gz file from the ncbi. 

```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
```
Then i ran hisat2 on it. 
```
module load hisat2/2.2.1
hisat2-build GRCh38_latest_genomic.fna GRCh38_latest_genomic.fna
```

However this gave the error of out of memory so the job was killed.

So I downloaded the human genome reference that was already run with hisat2 on the ncbi and opened it with tar. 

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.hisat2_index.tar.gz
gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.hisat2_index.tar.gz
tar -zxvf GCA_000001405.15_GRCh38_full_analysis_set.fna.hisat2_index.tar.gz
```
Then I downloaded the gene annotated file for human
```
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz  
gzip -d Homo_sapiens.GRCh38.84.gtf.gz
```
Then I ran alignment with [myjobRNAalign.sh](myjobRNAalign.sh)
Then I ran [myjobRNAcount.sh](myjobRNAcount.sh)
Then I ran [DeSeq.R](DeSeq.R) with inputs shortRNAseq.txt and human_counts_RNA.txt however, due to privacy, I will not be providing that txt file. Some quality control files I received were [dispersion](disp.png), [histogram](histo.png), [MAplot](MAplot.png), and [HeatMap](HeatMap.png). I had to delete some informatoin from bottom axis from heatmap for privacy reasons. 

Below is the evidence for successful run:

```
drwxrws--- 2 johnnl15 som  564 Mar 15 15:15 Bam
drwxrws--- 5 johnnl15 som   49 Mar 15 14:30 .
-rw-rw-r-- 1 johnnl15 som  681 Mar 15 14:28 myjobRNAcount.sh
drwxrwsr-x 2 johnnl15 som   11 Mar 15 14:13 ref
-rw-rw-r-- 1 johnnl15 som  937 Mar 15 01:13 myjobRNAalign.sh
drwxrws--- 2 johnnl15 som   84 Mar 14 16:21 FASTQC
-rw-rw-r-- 1 johnnl15 som  521 Mar 14 15:40 myjobRNAfastqc.sh
-rw-rw---- 1 johnnl15 som 2.4G Mar 14 06:00 LD1_1.fq.gz
```

Finally, I want to do a volcano plot which Dr. Long had said would be appropriate for me since I was running a RNAseq analysis pipeline on my in-house RNA sequencing files. In my DeSeq.R file I also included this script: 

```
png(filename='LKLDVolcano.png', width=800, height=750)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
graphics.off()
```

and my product was this [Volcano.png](Volcano.png). 
