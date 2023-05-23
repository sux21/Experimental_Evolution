# Su_Xingyuan_May_to_August_2023
Bioinformatics project on Rhizobium leguminosarum 

# Step 0 - Background Information of Illumina Sequencing
Paired-End Sequencing: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html.

# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
### 1. Running FastQC for raw data (Practice on 2018 strains)
FastQC v0.11.5

Running FastQC: 
```
nohup fastqc -o /home/xingyuan/2018_strains/fastQC_raw_reads GSF2234-295A_S28_R2_001.fastq &
```

Input files in /home/xingyuan/2018_strains/raw_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_raw_reads.

### 2. Running MultiQC for raw data (Practice on 2018 strains) 
multiqc, version 1.9

Running MultiQC in the directory with the FastQC reports: 
```
multiqc .
```

Input files in /home/xingyuan/2018_strains/fastQC_raw_reads. <br>
Output files in /home/xingyuan/2018_strains/fastQC_raw_reads.

### 3. Running Trimmomatic (Practice on 2018 strains)
version: trimmomatic-0.39.jar

Running for all files using for loop in shell script: <br>
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
https://www.melbournebioinformatics.org.au/tutorials/tutorials/assembly/spades/

#### For 101A <br>
**Run SPAdes on trimmed reads**: 
```
nohup spades.py --pe1-1 GSF2234-101A_S1_R1_P_001.fq.gz --pe1-2 GSF2234-101A_S1_R2_P_001.fq.gz -o spades-101A &
```

```
Quast output:
Assembly                   scaffolds
# contigs (>= 0 bp)        146      
# contigs (>= 1000 bp)     68       
Total length (>= 0 bp)     7138120  
Total length (>= 1000 bp)  7114449  
# contigs                  81       
Largest contig             553567   
Total length               7123407  
GC (%)                     60.85    
N50                        190187   
N75                        104963   
L50                        12       
L75                        25       
# N's per 100 kbp          8.70 
```

**Run SPAdes on trimmed reads - using --isolate**:
```
spades.py -1 GSF2234-101A_S1_R1_P_001.fq.gz -2 GSF2234-101A_S1_R2_P_001.fq.gz --isolate -o spades-101A-test
```

```
Quast results:
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        210      
# contigs (>= 1000 bp)     49       
Total length (>= 0 bp)     7155610  
Total length (>= 1000 bp)  7128413  
# contigs                  52       
Largest contig             564345   
Total length               7130378  
GC (%)                     60.85    
N50                        226578   
N75                        127550   
L50                        10       
L75                        20       
# N's per 100 kbp          8.55  
```
**Test the following command and compare with the above to see if there are differences**
```
spades.py --pe1-1 GSF2234-101A_S1_R1_P_001.fq.gz --pe1-2 GSF2234-101A_S1_R2_P_001.fq.gz --isolate -o spades-101A-test2 
```
```
Quast output:
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        210      
# contigs (>= 1000 bp)     49       
Total length (>= 0 bp)     7155610  
Total length (>= 1000 bp)  7128413  
# contigs                  52       
Largest contig             564345   
Total length               7130378  
GC (%)                     60.85    
N50                        226578   
N75                        127550   
L50                        10       
L75                        20       
# N's per 100 kbp          8.55 
```
**Test metaSPAdes**
```
spades.py -1 GSF2234-101A_S1_R1_P_001.fq.gz -2 GSF2234-101A_S1_R2_P_001.fq.gz --meta --isolate -o metaspades-101A-test
```

## After Assembly 

## Other Infomation
**Codes**
``spades.py -1 GSF2234-101A_S1_R1_P_001.fq.gz -2 GSF2234-101A_S1_R2_P_001.fq.gz --careful -o spades-101A-test``

``Bandage image assembly_graph_with_scaffolds.gfa assembly_graph_with_scaffolds.jpg``

``quast.py scaffolds.fasta``

``scp xingyuan@info.mcmaster.ca:/home/xingyuan/2018_strains/fastQC_trimmed_reads/GSF2234-105A_S2_R2_P_001_fastqc.html /Users/xingyuansu/Desktop``


