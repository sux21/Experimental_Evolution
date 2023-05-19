# Su_Xingyuan_May_to_August_2023
Bioinformatics project on Rhizobium leguminosarum 

```

aaa

```
# Step 0 - Background Information of Illumina Sequencing
Paired-End Sequencing: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html.

# Step 1 - Genome Assembly using SPAdes <br>
version: SPAdes genome assembler v3.15.2

## Before Assembly
### 1. Running FastQC (Practice on 2018 strains)
FastQC v0.11.5

Running FastQC: ``nohup fastqc -o /home/xingyuan/2018_strains/fastQC_raw_reads GSF2234-295A_S28_R2_001.fastq &`` 

### 2. Running MultiQC (Practice on 2018 strains) 
multiqc, version 1.9

Running MultiQC: ``multiqc .`` in the directory with the FastQC reports. 

### 3. Running Trimmomatic (Practice on 2018 strains)
version: trimmomatic-0.39.jar

Running for all files using for loop in shell script: <br>
``#!/bin/bash`` <br>
``for R1 in *R1*`` <br>
``do`` <br>
``R2=${R1//R1_001.fastq/R2_001.fastq}`` <br>
``R1_P=${R1//001.fastq/P_001.fq.gz}`` <br>
``R1_UP=${R1//001.fastq/UP_001.fq.gz}`` <br>
``R2_P=${R2//001.fastq/P_001.fq.gz}`` <br>
``R2_UP=${R2//001.fastq/UP_001.fq.gz}`` <br>

``java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/2018_strains/raw_reads/$R1 /home/xingyuan/2018_strains/raw_reads/$R2 /home/xingyuan/2018_strains/trim_2nd_attempt/$R1_P /home/xingyuan/2018_strains/trim_2nd_attempt/$R1_UP /home/xingyuan/2018_strains/trim_2nd_attempt/$R2_P /home/xingyuan/2018_strains/trim_2nd_attempt/$R2_UP ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15`` <br>
``done``

### 4. Repeat 1 and 2 (Practice on 2018 strains)

Run FastQC: ``nohup fastqc -o /home/xingyuan/2018_strains/fastQC_trim_2nd_attempt *_P_* &``. 

Run MultiQC: ``multiqc -f .`` in the directory with the FastQC reports.

## During Assembly 
https://www.melbournebioinformatics.org.au/tutorials/tutorials/assembly/spades/

## After Assembly 


### Codes
``nohup spades.py --pe1-1 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R1_001.fastq.gz --pe1-2 9_7_9_ACTTGTTATC-TCTAGGCGCG_L002_R2_001.fastq.gz -o 9_7_9-spades &``

``Bandage image assembly_graph_with_scaffolds.gfa assembly_graph_with_scaffolds.jpg``

``quast.py contigs.fasta``

``scp xingyuan@info.mcmaster.ca:/home/xingyuan/2018_strains/fastQC_trimmed_reads/GSF2234-105A_S2_R2_P_001_fastqc.html /Users/xingyuansu/Desktop``

## Error codes 

