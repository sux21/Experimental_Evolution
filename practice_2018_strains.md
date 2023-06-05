# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
### 1. Run FastQC for raw data (Practice on 2018 strains)
FastQC v0.11.5

Running FastQC: 
```
nohup fastqc -o /home/xingyuan/2018_strains/fastQC_raw_reads GSF2234-295A_S28_R2_001.fastq &
```

Input files in /home/xingyuan/2018_strains/raw_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_raw_reads.

### 2. Run MultiQC for raw data (Practice on 2018 strains) 
multiqc, version 1.9

Running MultiQC in the directory with the FastQC reports: 
```
multiqc .
```

Input files in /home/xingyuan/2018_strains/fastQC_raw_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_raw_reads.

### 3. Run Trimmomatic (Practice on 2018 strains)
version: trimmomatic-0.39.jar

Running for all files using "for loop" in shell script: <br>
```
#!/bin/bash 
for R1 in *R1* 
do 
R2=${R1//R1_001.fastq/R2_001.fastq} 
R1_P=${R1//001.fastq/P_001.fq.gz} 
R1_UP=${R1//001.fastq/UP_001.fq.gz} 
R2_P=${R2//001.fastq/P_001.fq.gz} 
R2_UP=${R2//001.fastq/UP_001.fq.gz} 

java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/2018_strains/raw_reads/$R1 /home/xingyuan/2018_strains/raw_reads/$R2 /home/xingyuan/2018_strains/trim_2nd_attempt/$R1_P /home/xingyuan/2018_strains/trim_2nd_attempt/$R1_UP /home/xingyuan/2018_strains/trim_2nd_attempt/$R2_P /home/xingyuan/2018_strains/trim_2nd_attempt/$R2_UP ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 
done
```

Input files in /home/xingyuan/2018_strains/fastQC_raw_reads. <br>
Output files in /home/xingyuan/2018_strains/trimmed_reads.

### 4. Run FastQC again for trimmed data (Practice on 2018 strains)
Run FastQC: 
```
nohup fastqc -o /home/xingyuan/2018_strains/fastQC_trim_2nd_attempt *_P_* &
```

Input files in /home/xingyuan/2018_strains/trimmed_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_trimmed_reads.

### 5. Run MultiQC again for trimmed data (Practice on 2018 strains)
Run MultiQC in the directory with the FastQC reports: 
```
multiqc -f . 
```
Input files in /home/xingyuan/2018_strains/fastQC_trimmed_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_trimmed_reads.

### 6. Steps 3, 4, 5 were repeated one time to remove the first 15bp of the read.

## During Assembly 
### 1. Running SPAdes (Practice on 2018 strains)
version: SPAdes genome assembler v3.15.2

#### 214C, 218A, 272A, 295A
command: ``spades -1 R1.fq.gz -2 R2.fq.gz --isolate -o spades-sample#`` <br>
command: ``quast -m 0 -o directory file``

contigs  
| Not filtered |           |                   |                 |              |
|:--------------:|:-----------:|:------:| :---:|:------:|
| sample | N50       | number of contigs |  largest contig | total length |
| 214C   |  225893   |    283            | 913776          |  7356909     |
| 218A   |  315517   |    306            | 739666          |  7202393     |
| 272A   |  435813   |    447            | 883349          |  7196064     |
| 295A   |  435863   |    615            | 652001          |  7226378     |

scaffolds     
| sample | N50       | number of scaffolds |  largest scaffold | total length |
|--------|-----------|---------------------|-------------------|--------------|
| 214C   | 273955    |         277         |      1490525      |    7357250   |
| 218A   | 559785    |         295         |      1032152      |    7203126   |
| 272A   | 521954    |         442         |      1690701      |    7196338   |
| 295A   | 521252    |         608         |      1101732      |    7226762   |


## After Assembly 
### 1. Run BWA
Version: 0.7.17-r1188 

#### 214C, 218A, 272A, 295A
1. Index the contigs/scaffolds
```
bwa index /home/xingyuan/2018_strains/trim_2nd_attempt/spades-295A/scaffolds.fasta     
```
2. align reads back to contigs/scaffolds
```
bwa mem -t 2 /home/xingyuan/2018_strains/trim_2nd_attempt/spades-295A/scaffolds.fasta /home/xingyuan/2018_strains/trim_2nd_attempt/GSF2234-295A_S28_R1_P_001.fq.gz /home/xingyuan/2018_strains/trim_2nd_attempt/GSF2234-295A_S28_R2_P_001.fq.gz > 295A_scaffolds.sam 
```
3. convert SAM file to BAM file (see commands [here](http://www.htslib.org/doc/samtools-view.html))
```
samtools view -bS 295A_scaffolds.sam > 295A_scaffolds.bam 
```
4. sort the BAM file (see commands [here](http://www.htslib.org/doc/samtools-sort.html))
```
samtools sort -o 295A_scaffolds.sorted.bam 295A_scaffolds.bam 
```
5. index the BAM file (see commands [here](http://www.htslib.org/doc/samtools-index.html))
```
samtools index 295A_scaffolds.sorted.bam 
```
6. obtain summary statistics (see commands [here](http://www.htslib.org/doc/samtools-flagstat.html))
```
samtools flagstat 295A_scaffolds.sorted.bam > mapping_summary 
```
7. Run qualimap for more information (see commands [here](http://qualimap.conesalab.org/doc_html/analysis.html))
```
qualimap bamqc -outdir bamqc -bam 295A_scaffolds.sorted.bam 
```

### Run Prokka 
prokka 1.12-beta

**sample 214C**
```
/home/sam/miniconda3/envs/pangenome/bin/prokka --force --outdir prokka-214C scaffolds.fasta
```
