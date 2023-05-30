# Su_Xingyuan_Summer_2023
Bioinformatics project on Rhizobium leguminosarum 

Work done on INFO server, node info114. 

# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
Create symbolic link for raw reads:
```
% ln -s /2/scratch/batstonelab/RltEE2020-PE_reads/*gz .
```
Check the number of files:
```
% ls *.gz | wc -l
726
```

### 1. Run FastQC for raw data
Versions: FastQC v0.11.5

**Sample: 1_1_2**

Run FastQC:
```
% nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads 1_1_2* &
```

### 2. Run Trimmomatic 
Version: 0.39


**Sample: 1_1_2**
```
% java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_001.fastq.gz /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_UP_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_UP_001.fastq.gz ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 LEADING:3 TRAILING:3 MINLEN:36
```

### 3. Run FastQC on trimmed reads 
**Sample: 1_1_2**
Run FastQC: 
```
% fastqc -o /home/xingyuan/rhizo_ee/fastQC_trimmomatic_reads *_P_* 
```

## During Assembly 
### 1. Running SPAdes (Practice on 2018 strains)
Version: 3.15.2

**Sample: 1_1_2**
```
% spades.py -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2
```
spades.py -1 GSF2234-125A_S3_R1_P_001.fq.gz -2 GSF2234-125A_S3_R2_P_001.fq.gz --isolate -o spades-125A


## After Assembly 
### 1. Run Quast on scaffolds
Run ``quast.py scaffolds.fasta``: <br>
Version 3.1, build 29.08.2015 16:09
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

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
### 2.
prokka 1.12-beta

 

