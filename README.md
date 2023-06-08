# Experimental Evoluntion Strains
Bioinformatics project on *Rhizobium leguminosarum*

Work done on INFO server, node info114. 

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
Create symbolic link for raw reads:
```
ln -s /2/scratch/batstonelab/RltEE2020-PE_reads/*gz .
```
Check the number of files:
```
ls *.gz | wc -l
```

### 1. Run FastQC for raw data
Versions: FastQC v0.11.5

**Sample: 1_1_2**

Run FastQC:
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads 1_1_2* &
```
**All 363 samples (726 files)**
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads *gz &
```

### 2. Run MultiQC for raw data
```
multiqc . 
```
### 3. Run Trimmomatic 
Version: 0.39


**Sample: 1_1_2**
```
nohup java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_001.fastq.gz /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_UP_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_UP_001.fastq.gz ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 CROP:130 LEADING:3 TRAILING:3 MINLEN:36 &
```
**All 363 samples (726 files)**
```
nano run_trimmomatic.sh
```
```
#!/bin/bash 
for R1 in *R1* 
do 
R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz} 
R1_P=${R1//001.fastq.gz/P_001.fastq.gz} 
R1_UP=${R1//001.fastq.gz/UP_001.fastq.gz} 
R2_P=${R2//001.fastq.gz/P_001.fastq.gz} 
R2_UP=${R2//001.fastq.gz/UP_001.fastq.gz} 

java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/rhizo_ee/raw_reads/$R1 /home/xingyuan/rhizo_ee/raw_reads/$R2 /home/xingyuan/rhizo_ee/trimmomatic_reads/$R1_P /home/xingyuan/rhizo_ee/trimmomatic_reads/$R1_UP /home/xingyuan/rhizo_ee/trimmomatic_reads/$R2_P /home/xingyuan/rhizo_ee/trimmomatic_reads/$R2_UP ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 CROP:130 LEADING:3 TRAILING:3 MINLEN:36
done
```
```
nohup bash run_trimmomatic.sh &
```

### 4. Run FastQC for trimmed reads 
**Sample: 1_1_2**

Run FastQC: 
```
fastqc -o /home/xingyuan/rhizo_ee/fastQC_trimmomatic_reads *_P_* 
```

**All 363 samples (726 files)**
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_trimmomatic_reads *_P_* &
```

### 5. Run MultiQC for trimmed reads
Version: MultiQC v1.9

**All 363 samples (726 files)**
```
multiqc . 
```

## During Assembly 
### 1. Run SPAdes 
Version: 3.15.2

**Sample: 1_1_2**
```
nohup spades.py --careful -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2 &
```

**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
vi running-spades.sh
```
```
#!/bin/bash 
for R1 in 10_1_8_*R1_P_* 10_1_9_*R1_P_* 10_7_6_*R1_P_* 11_4_2_*R1_P_* 11_4_4_*R1_P_* 11_5_6_*R1_P_* 13_4_1_*R1_P_* 14_4_6_*R1_P_* 14_5_3_*R1_P_* 15_4_4_*R1_P_* 15_4_6_*R1_P_* 16_1_6_*R1_P_* 16_1_7_*R1_P_* 16_1_8_*R1_P_* 16_4_2_*R1_P_* 16_4_3_*R1_P_* 16_6_6_*R1_P_* 17_2_1_*R1_P_* 17_2_8_*R1_P_* 17_2_9_*R1_P_* 18_1_4_*R1_P_* 18_1_5_*R1_P_* 19_1_1_*R1_P_* 19_5_8_*R1_P_* 2_2_5_*R1_P_* 2_5_2_*R1_P_* 2_6_4_*R1_P_* 3_1_5_*R1_P_* 3_2_1_*R1_P_* 3_2_3_*R1_P_* 3_2_6_*R1_P_* 3_2_7_*R1_P_* 3_3_5_*R1_P_* 3_3_7_*R1_P_* 3_3_9_*R1_P_* 4_1_2_*R1_P_* 4_1_4_*R1_P_* 4_2_1_*R1_P_* 6_4_5_*R1_P_* 6_4_7_*R1_P_* 6_7_5_*R1_P_* 7_1_2_*R1_P_* 7_1_5_*R1_P_* 7_6_3_*R1_P_* 7_6_9_*R1_P_* 7_7_2_*R1_P_* 7_7_3_*R1_P_* 8_4_10_*R1_P_* 8_4_4_*R1_P_* 9_3_7_*R1_P_* 9_7_6_*R1_P_* 9_7_9_*R1_P_* 
do
R2=${R1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}

spades.py --careful -1 $R1 -2 $R2 -o /home/xingyuan/rhizo_ee/spades_assembly/${R1%_*_L002_*gz}
done
```
```
nohup bash running-spades.sh &
```
manipulating string (useful for writing loop): https://mywiki.wooledge.org/BashFAQ/100, https://tldp.org/LDP/abs/html/string-manipulation.html

### 2. Run plasmidSPAdes
**Sample: 1_1_2**
```
nohup spades.py --plasmid -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2-plasmids &
```

## After Assembly 
### 1. Run Quast on scaffolds
Version: 5.2.0, 3d87c606

**Sample: 1_1_2** <br>
```
quast.py -m 0 scaffolds.fasta
```
**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
quast -m 0 scaffolds.fasta 
```
### 2. Run Spine on scaffolds
Version: 0.3.2

**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets** <br>
Copy scaffolds.fasta from ``rhizo_ee/spades_assembly`` to ``rhizo_ee/2008_2020_strains_comparison``:
```
#!/bin/bash 
for i in 10_1_8 13_4_1 15_4_6 16_4_2 17_2_8 19_1_1 2_6_4 3_2_6 3_3_9 6_4_5 7_1_5 7_7_3 9_7_6 10_1_9 11_4_2 14_4_6 16_1_6 16_4_3 17_2_9 19_5_8 3_1_5 3_2_7 4_1_2 6_4_7 7_6_3 8_4_10 9_7_9 10_7_6 11_4_4 14_5_3 16_1_7 16_6_6 18_1_4 2_2_5 3_2_1 3_3_5 4_1_4 6_7_5 7_6_9 8_4_4 11_5_6 15_4_4 16_1_8 17_2_1 18_1_5 2_5_2 3_2_3 3_3_7 4_2_1 7_1_2 7_7_2 9_3_7
do

cp /home/xingyuan/rhizo_ee/spades_assembly/$i/scaffolds.fasta /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/"$i-scaffolds.fasta"
done
```

Install Spine following instructions on https://github.com/Alan-Collins/Spine-Nucmer-SNPs:
```
conda config --append channels conda-forge
conda config --append channels bioconda

# Create conda environment called "sns" for installation and install conda package
conda create -n sns -c alan-collins spine-nucmer-snps

# activate sns
source activate sns
```

Run Spine:
```
ls *.fasta | awk 'BEGIN { FS="\t"; OFS="\t" } { print "/home/xingyuan/rhizo_ee/2008_2020_strains_comparison/"$1, $1, "fasta" }' > ./SPINE/config.txt

spine.pl -f /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/config.txt 
```







