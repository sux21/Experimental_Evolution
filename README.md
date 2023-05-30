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

## After Assembly 
### 1. Run Quast on scaffolds
Version: 5.2.0, 3d87c606

**Sample: 1_1_2**
Run Quast:
```
% quast.py scaffolds.fasta
```

 

