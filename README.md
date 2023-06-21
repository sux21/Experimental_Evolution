# 2020 Experimental Evoluntion 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on info server. Compute canada server will be used if the info server cannot run a program. Results produced by compute canada server will be transferred to info server. 

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
### 1. Run FastQC for raw data
Versions: FastQC v0.11.5 <br>
Work done on info114

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
Version: MultiQC v1.9 <br>
Work done on info114

```
multiqc . 
```
### 3. Run Trimmomatic 
Version: 0.39 <br>
Work done on info114

**Sample: 1_1_2**
```
nohup java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_001.fastq.gz /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_UP_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_UP_001.fastq.gz ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 CROP:130 LEADING:3 TRAILING:3 MINLEN:36 &
```
**All 363 samples (726 files)**
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

### 4. Run FastQC for trimmed reads 
Versions: FastQC v0.11.5 <br>
Work done on info114

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
Version: MultiQC v1.9 <br>
Work done on info114

**All 363 samples (726 files)**
```
multiqc . 
```

## During Assembly 
### 1. Run SPAdes 
Version: 3.15.2 <br>
Work done on info114

**Sample: 1_1_2**
```
nohup spades.py --careful -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2 &
```

**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
#!/bin/bash 
for R1 in 10_1_8_*R1_P_* 10_1_9_*R1_P_* 10_7_6_*R1_P_* 11_4_2_*R1_P_* 11_4_4_*R1_P_* 11_5_6_*R1_P_* 13_4_1_*R1_P_* 14_4_6_*R1_P_* 14_5_3_*R1_P_* 15_4_4_*R1_P_* 15_4_6_*R1_P_* 16_1_6_*R1_P_* 16_1_7_*R1_P_* 16_1_8_*R1_P_* 16_4_2_*R1_P_* 16_4_3_*R1_P_* 16_6_6_*R1_P_* 17_2_1_*R1_P_* 17_2_8_*R1_P_* 17_2_9_*R1_P_* 18_1_4_*R1_P_* 18_1_5_*R1_P_* 19_1_1_*R1_P_* 19_5_8_*R1_P_* 2_2_5_*R1_P_* 2_5_2_*R1_P_* 2_6_4_*R1_P_* 3_1_5_*R1_P_* 3_2_1_*R1_P_* 3_2_3_*R1_P_* 3_2_6_*R1_P_* 3_2_7_*R1_P_* 3_3_5_*R1_P_* 3_3_7_*R1_P_* 3_3_9_*R1_P_* 4_1_2_*R1_P_* 4_1_4_*R1_P_* 4_2_1_*R1_P_* 6_4_5_*R1_P_* 6_4_7_*R1_P_* 6_7_5_*R1_P_* 7_1_2_*R1_P_* 7_1_5_*R1_P_* 7_6_3_*R1_P_* 7_6_9_*R1_P_* 7_7_2_*R1_P_* 7_7_3_*R1_P_* 8_4_10_*R1_P_* 8_4_4_*R1_P_* 9_3_7_*R1_P_* 9_7_6_*R1_P_* 9_7_9_*R1_P_* 
do
R2=${R1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}

spades.py --careful -1 $R1 -2 $R2 -o /home/xingyuan/rhizo_ee/spades_assembly/${R1%_*_L002_*gz}
done
```

### 2. Run plasmidSPAdes
**Sample: 1_1_2**
```
nohup spades.py --plasmid -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2-plasmids &
```

## After Assembly 
### 1. Run Quast on scaffolds
Version: 5.2.0, 3d87c606 <br>
Work done on info114 

**Sample: 1_1_2** <br>
```
quast.py -m 0 scaffolds.fasta
```
**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
quast -m 0 scaffolds.fasta 
```
**Samples: 28 samples from 2008 in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
#!/bin/bash 
for i in Rht*
do

/home/xingyuan/tools/quast-5.2.0/quast.py $i -m 0 -o /home/xingyuan/rhizo_ee/quast_2008_long_reads/${i%.fasta}
done
```

# Step 2 - Data Analysis
## Step 1 - Find the most related 2008 strain for each 2020 strain 
### Method 1: Spine-Nucmer-SNPs 
Commands are taken from https://github.com/Alan-Collins/Spine-Nucmer-SNPs. 

**Samples: 52 samples from 2020 + 28 samples from 2008 in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets** <br>

#### (1) Copy contigs.fasta from ``rhizo_ee/spades_assembly`` to ``rhizo_ee/2008_2020_strains_comparison``:
```
#!/bin/bash 
for i in 10_1_8 13_4_1 15_4_6 16_4_2 17_2_8 19_1_1 2_6_4 3_2_6 3_3_9 6_4_5 7_1_5 7_7_3 9_7_6 10_1_9 11_4_2 14_4_6 16_1_6 16_4_3 17_2_9 19_5_8 3_1_5 3_2_7 4_1_2 6_4_7 7_6_3 8_4_10 9_7_9 10_7_6 11_4_4 14_5_3 16_1_7 16_6_6 18_1_4 2_2_5 3_2_1 3_3_5 4_1_4 6_7_5 7_6_9 8_4_4 11_5_6 15_4_4 16_1_8 17_2_1 18_1_5 2_5_2 3_2_3 3_3_7 4_2_1 7_1_2 7_7_2 9_3_7
do

cp /home/xingyuan/rhizo_ee/spades_assembly/$i/contigs.fasta /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/"$i-contigs.fasta"
done
```

#### (2) Run Spine
Version: 0.3.2 <br>
Work done on info114. 

```
ls *.fasta | awk 'BEGIN { FS="\t"; OFS="\t" } { print "/home/xingyuan/rhizo_ee/2008_2020_strains_comparison/"$1, $1, "fasta" }' > ./SPINE/config.txt
```
```
spine.pl -f /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/SPINE/config.txt 
```

#### PAUSE (3) Run Nucmer 
Version: 3.1 <br>
Work done on info114

```
ls ../SPINE/*.core.fasta | while read i; do acc=${i%.core*}; acc=${acc#../SPINE/output.}; nucmer --prefix=${acc}_core ../SPINE/output.backbone.fasta $i; delta-filter -r -q ${acc}_core.delta > ${acc}_core.filter; show-snps -Clr ${acc}_core.filter > ${acc}_core.snps; done
```

#### PAUSE (4) Run snps2fasta.py
Work done on info114

```
python3 snps2fasta.py -r ../SPINE/output.backbone.fasta -f variant_core.fasta -p '(.*)_core\.snps' ../NUCMER/*.snps
```

#### PAUSE (5) Run fasta2diffmat.py 
Work done on info114

```
python fasta2diffmat.py -f variant_core.fasta -d diff_dict.pkl -t 5 
```
#### PAUSE (6) Run get_snps_support_MP.py
**Run BWA to map reads back to contigs, and create SAM files**
Version: 0.7.17-r1188 <br>
Work done on info114

```
#!/bin/bash 
for i in 10_1_8 13_4_1 15_4_6 16_4_2 17_2_8 19_1_1 2_6_4 3_2_6 3_3_9 6_4_5 7_1_5 7_7_3 9_7_6 10_1_9 11_4_2 14_4_6 16_1_6 16_4_3 17_2_9 19_5_8 3_1_5 3_2_7 4_1_2 6_4_7 7_6_3 8_4_10 9_7_9 10_7_6 11_4_4 14_5_3 16_1_7 16_6_6 18_1_4 2_2_5 3_2_1 3_3_5 4_1_4 6_7_5 7_6_9 8_4_4 11_5_6 15_4_4 16_1_8 17_2_1 18_1_5 2_5_2 3_2_3 3_3_7 4_2_1 7_1_2 7_7_2 9_3_7; do

bwa index $i-contigs.fasta

done
```
```
#!/bin/bash 
for x1 in /home/xingyuan/rhizo_ee/trimmomatic_reads/{10_1_8_*R1_P_*,10_1_9_*R1_P_*,10_7_6_*R1_P_*,11_4_2_*R1_P_*,11_4_4_*R1_P_*,11_5_6_*R1_P_*,13_4_1_*R1_P_*,14_4_6_*R1_P_*,14_5_3_*R1_P_*,15_4_4_*R1_P_*,15_4_6_*R1_P_*,16_1_6_*R1_P_*,16_1_7_*R1_P_*,16_1_8_*R1_P_*,16_4_2_*R1_P_*,16_4_3_*R1_P_*,16_6_6_*R1_P_*,17_2_1_*R1_P_*,17_2_8_*R1_P_*,17_2_9_*R1_P_*,18_1_4_*R1_P_*,18_1_5_*R1_P_*,19_1_1_*R1_P_*,19_5_8_*R1_P_*,2_2_5_*R1_P_*,2_5_2_*R1_P_*,2_6_4_*R1_P_*,3_1_5_*R1_P_*,3_2_1_*R1_P_*,3_2_3_*R1_P_*,3_2_6_*R1_P_*,3_2_7_*R1_P_*,3_3_5_*R1_P_*,3_3_7_*R1_P_*,3_3_9_*R1_P_*,4_1_2_*R1_P_*,4_1_4_*R1_P_*,4_2_1_*R1_P_*,6_4_5_*R1_P_*,6_4_7_*R1_P_*,6_7_5_*R1_P_*,7_1_2_*R1_P_*,7_1_5_*R1_P_*,7_6_3_*R1_P_*,7_6_9_*R1_P_*,7_7_2_*R1_P_*,7_7_3_*R1_P_*,8_4_10_*R1_P_*,8_4_4_*R1_P_*,9_3_7_*R1_P_*,9_7_6_*R1_P_*,9_7_9_*R1_P_*}
do
x2=${x1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}
R1=${x1#/home/xingyuan/rhizo_ee/trimmomatic_reads/}
R2=${x2#/home/xingyuan/rhizo_ee/trimmomatic_reads/}


bwa mem -t 5 ${R1%_*_L002_*gz}-contigs.fasta $x1 $x2 > /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/SAMS/${R1%_*_L002_*gz}.aln.pe.sam

done
```

#### (6) FastANI on core genomes produced by Spine in step (1)
Instructions are taken from https://github.com/ParBLiSS/FastANI.

Version: 1.32 <br>
Work done on graham.computecanada.ca

```
#!/bin/bash
#SBATCH --time=00-01:00:00
#SBATCH --account=def-batstone

module load fastani/1.32
fastANI --ql query_list --rl reference_list --matrix -o fastani.out # query_list contains core genomes output by Spine from 52 2020 samples, reference_list contains core genomes output by Spine from 28 2008 samples
```
```
sbatch fastani.sh
Submitted batch job 7238556
```

#### (7) FastANI on contigs
Version: 1.32 <br>
Work done on info2020

```
/usr/local/bin/fastANI --ql query_list --rl reference_list --matrix -o fastani.contigs.out # query_list contains contigs from 52 2020 samples, reference_list contains all sequences from 28 2008 samples
```

#### (8) FastANI on core genomes produced by Spine in step (1) (using first sequence of the 28 original strains (assume it is chromosome) as reference) 
Version: 1.32 <br>
Work done on info2020

```
/usr/local/bin/fastANI --ql query_list --rl reference_list_first_seq --matrix -o fastani.2020_contigs-2008_first_seq.out # query_list contains contigs from 52 2020 samples, reference_list contains only the first sequences from each 28 2008 samples
```

### Method 2: Prokka-Roary-
(1) Run Prokka
Version: 1.12-beta <br>
Work on info114

```
prokka --locustag 10_1_8.scaffolds --cpus 5 --outdir /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/PROKKA/10_1_8-scaffolds --prefix 10_1_8.scaffolds /home/xingyuan/rhizo_ee/spades_assembly/10_1_8/scaffolds.fasta
```

 (2) Run Roary 

```
roary -p 5 *.gff 
```

## Step 2: Presence and Absence of Genes
### Method 1 - Glimmer-InterProScan-CDHIT
Method is taken from Genome-Wide Association Analyses in the Model Rhizobium Ensifer meliloti. https://journals.asm.org/doi/10.1128/mSphere.00386-18. 
#### (1) Glimmer
http://ccb.jhu.edu/software/glimmer/index.shtml 

[glim302notes.pdf](https://github.com/sux21/2020_Experimental_Evoluntion/files/11793551/glim302notes.pdf)

https://academic.oup.com/bioinformatics/article/23/6/673/419055

Change directory paths in ``g3-iterated.csh``
```
Before changes:
set awkpath = /fs/szgenefinding/Glimmer3/scripts
set glimmerpath = /fs/szgenefinding/Glimmer3/bin
set elphbin = /nfshomes/adelcher/bin/elph

After changes:
set awkpath = /home/xingyuan/tools/glimmer3.02/scripts
set glimmerpath = /home/xingyuan/tools/glimmer3.02/bin
set elphbin = /home/xingyuan/tools/ELPH/bin/Linux-i386/elph
```

```
#!/bin/bash
for i in 10_1_8 10_1_9 10_7_6 11_4_2 11_4_4 11_5_6 13_4_1 14_4_6 14_5_3 15_4_4 15_4_6 16_1_6 16_1_7 16_1_8 16_4_2 16_4_3 16_6_6 17_2_1 17_2_8 17_2_9 18_1_4 18_1_5 19_1_1 19_5_8 2_2_5 2_5_2 2_6_4 3_1_5 3_2_1 3_2_3 3_2_6 3_2_7 3_3_5 3_3_7 3_3_9 4_1_2 4_1_4 4_2_1 6_4_5 6_4_7 6_7_5 7_1_2 7_1_5 7_6_3 7_6_9 7_7_2 7_7_3 8_4_10 8_4_4 9_3_7 9_7_6 9_7_9; do

if [[ ! -d $i ]]; then 
   mkdir $i
fi

/home/xingyuan/tools/glimmer3.02/scripts/g3-iterated.csh /home/xingyuan/rhizo_ee/spades_assembly/$i/contigs.fasta $i-contigs; mv `ls $i-contigs*` $i/
done

for j in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/ASSEMBLIES/*Rht*_?.fasta; do

k=${j#/home/xingyuan/rhizo_ee/2008_2020_strains_comparison/ASSEMBLIES/}

if [[ ! -d ${k%.fasta} ]]; then 
   mkdir ${k%.fasta}
fi

/home/xingyuan/tools/glimmer3.02/scripts/g3-iterated.csh $j ${k%.fasta}; mv `ls ${k%.fasta}*` ${k%.fasta}/
done
```









