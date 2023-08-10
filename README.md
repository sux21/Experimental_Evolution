# 2020 Experimental Evolution 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on **info server** (contact Brian Golding at golding@mcmaster.ca for an info account). **Compute canada server** will be used if the info server cannot run the program (register for an compute canada account at https://ccdb.alliancecan.ca/security/login). Results produced by compute canada server will be transferred to info server. 

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

# Samples used in this project
- **Original strains**: starting populations at the start of this experimental evolution experiment. Number of strains is 56. I receive complete genomes for these strains.
- **Experimentally evolved strains**: derived populations at the end of this experimental evolution experiment. Number of strains is 363. I receive Illumina paired-end reads for these strains.

# Step 1 - Genome Assembly of experimentally evolved strains
## Before Assembly
### 1. Run FastQC for raw data
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Versions: FastQC v0.11.5 <br>
Work done on info114

**All 363 samples (726 files)**
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads *gz &
```

### 2. Run MultiQC for raw data
https://multiqc.info/

Version: MultiQC v1.9 <br>
Work done on info114

```
multiqc . 
```

### 3. Run Trimmomatic 
http://www.usadellab.org/cms/?page=trimmomatic

Version: 0.39 <br>
Work done on info114

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
https://github.com/ablab/spades

Version: 3.15.2 <br>
Work done on info114

**52 samples**
```
#!/bin/bash 
for R1 in 10_1_8_*R1_P_* 10_1_9_*R1_P_* 10_7_6_*R1_P_* 11_4_2_*R1_P_* 11_4_4_*R1_P_* 11_5_6_*R1_P_* 13_4_1_*R1_P_* 14_4_6_*R1_P_* 14_5_3_*R1_P_* 15_4_4_*R1_P_* 15_4_6_*R1_P_* 16_1_6_*R1_P_* 16_1_7_*R1_P_* 16_1_8_*R1_P_* 16_4_2_*R1_P_* 16_4_3_*R1_P_* 16_6_6_*R1_P_* 17_2_1_*R1_P_* 17_2_8_*R1_P_* 17_2_9_*R1_P_* 18_1_4_*R1_P_* 18_1_5_*R1_P_* 19_1_1_*R1_P_* 19_5_8_*R1_P_* 2_2_5_*R1_P_* 2_5_2_*R1_P_* 2_6_4_*R1_P_* 3_1_5_*R1_P_* 3_2_1_*R1_P_* 3_2_3_*R1_P_* 3_2_6_*R1_P_* 3_2_7_*R1_P_* 3_3_5_*R1_P_* 3_3_7_*R1_P_* 3_3_9_*R1_P_* 4_1_2_*R1_P_* 4_1_4_*R1_P_* 4_2_1_*R1_P_* 6_4_5_*R1_P_* 6_4_7_*R1_P_* 6_7_5_*R1_P_* 7_1_2_*R1_P_* 7_1_5_*R1_P_* 7_6_3_*R1_P_* 7_6_9_*R1_P_* 7_7_2_*R1_P_* 7_7_3_*R1_P_* 8_4_10_*R1_P_* 8_4_4_*R1_P_* 9_3_7_*R1_P_* 9_7_6_*R1_P_* 9_7_9_*R1_P_* 
do
R2=${R1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}

spades.py --careful -1 $R1 -2 $R2 -o /home/xingyuan/rhizo_ee/spades_assembly/${R1%_*_L002_*gz}
done
```

**Remaining 311 samples**
```
#!/bin/bash
for R1 in *R1*; do
R2=${R1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}

spades.py --careful -1 $R1 -2 $R2 -o /home/xingyuan/rhizo_ee/spades_assembly/${R1%_*_L002_*gz}
done
```

## After Assembly 
### 1.1 Run Quast for 363 experimentally evolved strains (Contigs)
https://github.com/ablab/quast

Version: 5.2.0, 3d87c606 <br>
Work done on info114 

```
#!/bin/bash
for i in *; do

if [[ $i =~ ".sh" ]] || [[ $i =~ "out" ]]; then
   continue
fi

/home/xingyuan/tools/quast-5.2.0/quast.py $i -m 0 -t 5 -o /home/xingyuan/rhizo_ee/spades_assembly/quast_2020_short_reads/${i%.fasta}

done
```

### 1.2 Run Quast for 363 experimentally evolved strains (Scaffolds)
Version: 5.2.0, 3d87c606 <br>
Work done on info114 

```
#!/bin/bash
for i in *; do

if [[ $i =~ ".sh" ]] || [[ $i =~ "out" ]]; then
   continue
fi

/home/xingyuan/tools/quast-5.2.0/quast.py $i -m 0 -t 5 -o /home/xingyuan/rhizo_ee/spades_assembly/quast_2020_short_reads/${i%.fasta}

done
```

### 1.3 Run Quast for 56 original strains 
Version: v5.2.0, 3d87c606 <br>
Work done on info113

```
#!/bin/bash 
for i in Rht*
do

/home/xingyuan/tools/quast-5.2.0/quast.py $i -m 0 -t 5 -o /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/quast_2008_long_reads/${i%.fasta}

done
```

# Step 2 - Data Analysis
## Analysis 1 - Find the most probable ancestor for each experimentally evolved strain 

**Samples: 363 experimentally evolved strains + 56 original strains**

### (Optional) 1. Run Spine: find core genomes 
https://github.com/Alan-Collins/Spine-Nucmer-SNPs

Version: 0.3.2 <br>
Work done on info114. 

```
#!/bin/bash
# Copy contigs.fasta from ``rhizo_ee/spades_assembly`` to ``rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY``:

for i in *; do

if [[ $i =~ ".sh" ]] || [[ $i =~ "fasta" ]] || [[ $i =~ "quast" ]]; then
   continue
fi

cp /home/xingyuan/rhizo_ee/spades_assembly/$i/contigs.fasta /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/"$i-contigs.fasta"
done
```

```
# Create a config.txt file for Spine
ls | awk 'BEGIN { FS="\t"; OFS="\t" } { print "/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/"$1, $1, "fasta" }' > ../SPINE/config.txt
```
```
# Run Spine
nohup spine.pl -f /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/SPINE/config.txt &
```

### 2. Run FastANI
https://github.com/ParBLiSS/FastANI

Version: 1.32  <br>
Work done on info2020

Query = 363 experimentally evolved strains <br>
Reference = 56 original strains 

**Query to reference (default)**
```
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -o fastani.contigs.out &
```
**Query to reference (--fragLen 4000)**
```
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 4000 -o fastani.contigs4000.out &
```
**Query to reference (--fragLen 3500)**
```
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 3500 -o fastani.contigs3500.out &
```
**Query to reference (--fragLen 2500)**
```
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 2500 -o fastani.contigs2500.out &
```
**Query to reference (--fragLen 2000)**
```
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 2000 -o fastani.contigs2000.out &
```
**Query (core genome) to reference (core genome)**
```
ls *contigs.fasta.core.fasta > core_query_list
ls *_?.fasta.core.fasta > core_reference_list

nohup /usr/local/bin/fastANI --ql core_query_list --rl core_reference_list -t 5 -k 5 --fragLen 10 -o fastani.que_core_to_ref_core.out 
```

**Reference to reference**
```
ls Rht* > reference_list

nohup /usr/local/bin/fastANI --ql reference_list --rl reference_list -o fastani.ref_to_ref.out &
```

**Query to query**
```
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl contigs_query_list -o fastani.que_to_que.out &
```

**Query+reference to query+reference (exclude as5_2_4)**
```
ls *.fasta > all_samples_no_as5_2_4 (remove as_5_2_4 from the list)

nohup /usr/local/bin/fastANI --ql all_samples_no_as5_2_4 --rl all_samples_no_as5_2_4 -t 5 -o fastani.all_to_all.out &
```

**Get mapping diagram for these comparisons: 20_6_10 vs 511_N, 20_6_10 vs 415_C; 20_2_6 vs 596_N, 20_2_6 vs 116_N; 16_1_6 vs 173_C, 16_1_6 vs 016_N; 16_3_4 vs 003_C, 16_3_4 vs 706_C**

Rscript Version: 3.6.3 (2020-02-29) <br>
Work done on info2020

**Open R and install genoplotr using ``install.packages("genoPlotR")`` before running the following lines**

```
/usr/local/bin/fastANI -q 20_6_10-contigs.fasta -r Rht_511_N.fasta --visualize -o 20_6_10-Rht_511_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_6_10-contigs.fasta Rht_511_N.fasta 20_6_10-Rht_511_N.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 20_6_10-contigs.fasta -r Rht_415_C.fasta --visualize -o 20_6_10-Rht_415_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_6_10-contigs.fasta Rht_415_C.fasta 20_6_10-Rht_415_C.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 20_2_6-contigs.fasta -r Rht_596_N.fasta --visualize -o 20_2_6-Rht_596_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_2_6-contigs.fasta Rht_596_N.fasta 20_2_6-Rht_596_N.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 20_2_6-contigs.fasta -r Rht_116_N.fasta --visualize -o 20_2_6-Rht_116_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_2_6-contigs.fasta Rht_116_N.fasta 20_2_6-Rht_116_N.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 16_1_6-contigs.fasta -r Rht_173_C.fasta --visualize -o 16_1_6-Rht_173_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_1_6-contigs.fasta Rht_173_C.fasta 16_1_6-Rht_173_C.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 16_1_6-contigs.fasta -r Rht_016_N.fasta --visualize -o 16_1_6-Rht_016_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_1_6-contigs.fasta Rht_016_N.fasta 16_1_6-Rht_016_N.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 16_3_4-contigs.fasta -r Rht_003_C.fasta --visualize -o 16_3_4-Rht_003_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_3_4-contigs.fasta Rht_003_C.fasta 16_3_4-Rht_003_C.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 16_3_4-contigs.fasta -r Rht_706_C.fasta --visualize -o 16_3_4-Rht_706_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_3_4-contigs.fasta Rht_706_C.fasta 16_3_4-Rht_706_C.fastani.out.visual 
```

## Analysis 2 - Phylogeny for the original strains 

**Samples: 56 original strains**

### 1. Run Spine-Nucmer-SNPs
https://github.com/Alan-Collins/Spine-Nucmer-SNPs

Version: 0.3.2 <br>
Work done on info113

**(1)**
```
ls Rht* | awk 'BEGIN { FS="\t"; OFS="\t" } { print "/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/"$1, $1, "fasta" }' > ../phylogeny_2008/config.txt
```
```
nohup spine.pl -f config.txt -t 5 &
```

**(2)**
```
ls ../spine/*.core.fasta | while read i; do acc=${i%.core*}; acc=${acc#../spine/output.}; nucmer --prefix=${acc}_core ../spine/output.backbone.fasta $i; delta-filter -r -q ${acc}_core.delta > ${acc}_core.filter; show-snps -Clr ${acc}_core.filter > ${acc}_core.snps; done
```

**(3)**
```
python3 ~/tools/Spine-0.3.2/snps2fasta.py -r ../spine/output.backbone.fasta -f variant_core.fasta -whole -m snp_matrix.csv -d '\t' -p '(.*)_core\.snps' ../nucmer/*.snps
```

### 2. IQ-Tree
http://www.iqtree.org/doc/Tutorial

Version: 2.2.0 <br>
Work done on info114

```
# For C_only population (28 samples)
nohup iqtree2 -T 5 -s C_only.fasta -bb 1000 -wbt --seqtype DNA &

# For N_only population (28 samples)
nohup iqtree2 -T 5 -s N_only.fasta -bb 1000 -wbt --seqtype DNA &

# For mixed population (28 samples)
iqtree2 -T 5 -s mix.fasta -bb 1000 -wbt --seqtype DNA
```

## Analysis 3: Gene Presence Absence
### 1. Prokka
https://github.com/tseemann/prokka

Version: 1.12-beta <br>
Work done on info114

```
#!/bin/bash

for i in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/*contigs.fasta; do
a=${i#/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/}
b=${a%-contigs.fasta}

/usr/local/prokka/bin/prokka --cpus 5 --outdir $b --prefix $b-contigs --locustag $b $i 

done

for j in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/*Rht*_?.fasta; do
c=${j#/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/}
d=${c%.fasta}

/usr/local/prokka/bin/prokka --cpus 5 --outdir $d --prefix $d-contigs --locustag $d $j

done
```
### 1. PGAP
#### Filter sequences shorter than 200 bp (pgap only takes sequences equal or longer than 200 bp)
https://github.com/shenwei356/seqkit

Seqkit Version: 2.3.0 <br>
Work done on graham cluster

```
#!/bin/bash
#SBATCH --time=00-00:05
#SBATCH --account=def-batstone

module load seqkit/2.3.1

for i in *contigs.fasta; do
sample=${i%-contigs.fasta}

seqkit seq -m 200 "$sample"-contigs.fasta > "$sample"-contigs.filter.fasta
done
```

#### Download pgap.py
https://github.com/ncbi/pgap/tree/1126_Test

Version: 2023-05-17.build6771 <br>
Work done on graham cluster 

```
module load apptainer 
wget -O pgap.py https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
chmod +x pgap.py
./pgap.py -D apptainer 
```

#### Run pgap (Prepare 5 scripts and run them separately, because one sample takes about 1h.)

Version: 2023-05-17.build6771 <br>
Work done on graham cluster 

**script 1 (94 samples)**
```
#!/bin/bash
#SBATCH --time=04-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/{1..5}_*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%-contigs.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 2 (70 samples)**
```
#!/bin/bash
#SBATCH --time=03-12:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/{6..9}_*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%-contigs.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 3 (109 samples)**
```
#!/bin/bash
#SBATCH --time=05-00:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/1{0..5}_*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%-contigs.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 4 (90 samples)**
```
#!/bin/bash
#SBATCH --time=04-12:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/1{6..9}_*filter.fasta /project/6078724/sux21/rhizo_ee/genomes/20_*filter.fasta /project/6078724/sux21/rhizo_ee/genomes/as5_2_4-contigs.filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%-contigs.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 5 (56 samples)**
```
#!/bin/bash
#SBATCH --time=03-00:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/Rht*fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
#### Verify species taxonomy 
https://github.com/ncbi/pgap/wiki/Taxonomy-Check

Version: 2023-05-17.build6771 <br>
Work done on cedar cluster

**Evolved strains**
```
#!/bin/bash
#SBATCH --time=01-10:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL

module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%-contigs.filter.fasta}

/project/6078724/sux21/tools/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap_2023-05-17.build6771.sif -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40 --taxcheck-only
done
```

**Original Strains**
```
#!/bin/bash
#SBATCH --time=00-13:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL

module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/Rht*fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/}
sample=${j%.fasta}

/project/6078724/sux21/tools/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap_2023-05-17.build6771.sif -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40 --taxcheck-only
done
```
**Samples with INCONCLUSIVE status: Rht_061_N, Rht_173_C, Rht_209_N, Rht_231_N, Rht_717_N, Rht_773_N**

### 2. Roary
https://sanger-pathogens.github.io/Roary/

Version: 1.007001 <br>
Work done on info115

#### Run roary
```
nohup /usr/local/bin/roary -p 5 /home/xingyuan/rhizo_ee/genes_presence_absence/prokka/*.gff &
```


### 2. GenAPI
https://github.com/MigleSur/GenAPI

Version: 1.0 <br>
Work done on info2020

```
nohup genapi -p 5 -m *.gff &
```

### 3. Don't use clustering softwares
#### Find genes lost in evolved strains 

Bedtools Version: 2.19.1 <br>
Work done info114

**Using prokka annotation**
```
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/call_snps/step1_bwa_mem/*bam; do
j=${i#/home/xingyuan/rhizo_ee/call_snps/step1_bwa_mem/*-}
ref=${j%.bam}
k=${i#/home/xingyuan/rhizo_ee/call_snps/step1_bwa_mem/}
sample=${k%.bam}

bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/prokka/"$ref".bed -b "$i" -header -v > "$sample".genes_lost
done
```

**Using prokka annotation (using mapped reads extracted by samtools)**

Samtools Version: 1.11 (using htslib 1.11) <br>
Bedtools Version: 2.19.1 <br>
Work done info114

```
#!/bin/bash
for i in *bam; do
sample=${i/.bam/.mapped.bam}

samtools view -F 260 -o "$sample".mapped.bam "$i"
done
```
```
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/call_snps/step1_bwa_mem/*.mapped.bam; do
j=${i#/home/xingyuan/rhizo_ee/call_snps/step1_bwa_mem/*-}
ref=${j%.mapped.bam.mapped.bam}
k=${i#/home/xingyuan/rhizo_ee/call_snps/step1_bwa_mem/}
sample=${k%.mapped.bam.mapped.bam}

bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/prokka/"$ref".bed -b "$i" -header -v > "$sample".genes_absence
done
```

**Using pgap annotation**
```
bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/pgap/Rht_056_N/annot_with_genomic_fasta.bed -b 14_2_2-Rht_056_N.bam -header -v > 14_2_2-Rht_056_N.genes_lost
```
#### Find genes gained in evolved strains 
Samtools Version: 1.11 (using htslib 1.11) <br>
Work done info114

**Extract unmapped reads**
```
#!/bin/bash
for i in *_?.bam; do
sample=${i/.bam/.unmapped.bam}

samtools view -f 4 -o "$sample" "$i"
done
```

## Analysis 4: Call SNPS between each evolved strain and its most probable ancestor
https://github.com/rtbatstone/how-rhizobia-evolve/blob/master/Variant%20discovery/Variant_calling.md

**Install the latest version of picard at https://github.com/broadinstitute/picard** 
```
wget https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar
```

### 1. BWA
https://bio-bwa.sourceforge.net/

**Useful tools: decoding SAM flags: https://broadinstitute.github.io/picard/explain-flags.html**

BWA Version: 0.7.17-r1188 <br>
Samtools Version: 1.11 (using htslib 1.11) <br>
Work done on info114

#### Prepare a csv file as the following
Since the csv file is saved by excel, it may say ``dos`` when opening using vi text editor. If so, unix can not read the file correctly. Opne the csv file using vi, and type ``:set ff=unix``. The ``dos`` should be gone. 

```
#Format: evolved strain,the most probable ancestor (samples are separated by comma)
20_6_10,Rht_511_N
10_3_2,Rht_511_N
19_6_10,Rht_511_N
20_2_6,Rht_596_N
.
.
.
(362 lines)
```

#### Index 56 genomes of original strains for bwa
```
#!/bin/bash
for i in Rht*fasta; do
bwa index $i
done
```

#### Run bwa (bwa mem -t 5 -M -R), convert .sam to .bam (samtools view -huS), and produce alignment statistics (samtools stats)
Samtools stats: http://www.htslib.org/doc/samtools-stats.html
```
#!/bin/bash
# Usage: ./ThisScript file.csv
while IFS=',' read -r a b; do
r1=/home/xingyuan/rhizo_ee/trimmomatic_reads/"$a"_*R1_P*
r2=/home/xingyuan/rhizo_ee/trimmomatic_reads/"$a"_*R2_P*
ref=/home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$b".fasta

bwa mem -t 5 -M -R "@RG\tID:"$a"\tSM:"$a $ref $r1 $r2 | samtools view -huS -o $a-"$b".bam - && samtools stats $a-"$b".bam > $a-"$b".stats

done < $1
```

### 2. Use picard to manipulate the SAM files
Picard Version: 3.0.0 <br>
Samtools Version: 1.13 <br>
Work done one info2020

#### (1) Reorder reads in the BAM file to match the contig ordering in reference file
https://gatk.broadinstitute.org/hc/en-us/articles/360037426651-ReorderSam-Picard-

##### Create sequence dictionary file (.dict) for reference 
https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-
```
#!/bin/bash
for i in Rht*fasta; do
ref=${i%.fasta}
java -jar /home/xingyuan/tools/picard.jar CreateSequenceDictionary -R "$ref".fasta -O "$ref".dict
done
```

##### Run ReorderSam
```
#!/bin/bash
for i in *bam; do
sample=${i%.bam}
ref=${sample#*-}

java -jar /home/xingyuan/tools/picard.jar ReorderSam -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".fasta -I "$sample".bam -O "$sample".reordered.bam -SD /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".dict
done
```

#### (2) Assign all the reads in a file to a single new read-group
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
```
#!/bin/bash
for i in *.reordered.bam; do
sample=${i%.reordered.bam}

java -jar /home/xingyuan/tools/picard.jar AddOrReplaceReadGroups -I "$sample".reordered.bam -O "$sample".new_rg.bam -ID "$sample" -LB rhizo_ee -PL Illumina -PU 1 -SM "$sample" 
done
```

#### (3) Sort the input BAM file by coordinate
https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-
```
#!/bin/bash
for i in *.new_rg.bam; do
sample=${i%.new_rg.bam}

java -jar /home/xingyuan/tools/picard.jar SortSam -I "$sample".new_rg.bam -O "$sample".coordinate_sorted.bam -SO coordinate
done
```

#### (4) Identify duplicate reads and index the BAM files
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard- <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-
```
#!/bin/bash
for i in *.coordinate_sorted.bam; do
sample=${i%.coordinate_sorted.bam}

java -jar /home/xingyuan/tools/picard.jar MarkDuplicates -I "$sample".coordinate_sorted.bam -O "$sample".marked_duplicates.bam -M "$sample".marked_dup_metrics.txt && java -jar /home/xingyuan/tools/picard.jar BuildBamIndex -I "$sample".marked_duplicates.bam
done
```

### 3. Run HaplotypeCaller
Samtools Version: 1.13 (using htslib 1.13) <br>
GATK Version: 4.4.0.0 <br>
Work done on info2020

**Download the latest release of gatk at https://github.com/broadinstitute/gatk/releases**
```
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip gatk-4.4.0.0.zip
```

#### Create index file (``.fai``) of reference for GATK
http://www.htslib.org/doc/samtools-faidx.html
```
#!/bin/bash
for i in Rht*fasta; do
samtools faidx $i
done
```

#### HaplotypeCaller
https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller
```
#!/bin/bash
for i in *.marked_duplicates.bam; do
sample=${i%.marked_duplicates.bam}
ref=${sample#*-}

/home/xingyuan/tools/gatk-4.4.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".fasta -I "$i" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O "$sample".g.vcf.gz -ERC GVCF
done
```

### 4. CombineGVCFs, GenotypeGVCFs
CombineGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832710975771-CombineGVCFs <br>
GenotypeGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs

GATK Version: 4.4.0.0 <br>
Work done on info2020

**Create a list, named ``MPA.list``,  for the most probable ancestors (25 most probable ancestors, 25 lines)**
```
Rht_016_N
Rht_056_N
Rht_074_C
Rht_097_N
Rht_108_C
Rht_113_C
Rht_156_N
Rht_173_C
Rht_325_C
Rht_415_C
Rht_438_C
Rht_449_C
Rht_460_C
Rht_462_C
Rht_493_C
Rht_511_N
Rht_527_N
Rht_559_C
Rht_596_N
Rht_706_C
Rht_758_C
Rht_773_N
Rht_837_C
Rht_861_C
Rht_901_C
```
**Run CombineGVCFs and GenotypeGVCFs**
```
#!/bin/bash
while read line; do

# Group evolved strains with the same most probable ancestors to a list. This should create 25 lists.
find *"$line"*.gz > MPA-"$line".list &&

# Run CombineGVCFs to combine the vcf files in each list to one vcf files. This should create 25 cohort.g.vcf.gz files.
/home/xingyuan/tools/gatk-4.4.0.0/gatk CombineGVCFs -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$line".fasta --variant MPA-"$line".list -O "$line".cohort.g.vcf.gz &&

# Run GenotypeGVCFs on each 25 cohort.g.vcf.gz files
/home/xingyuan/tools/gatk-4.4.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$line".fasta -V "$line".cohort.g.vcf.gz -ploidy 1 -O genotype_"$line".vcf.gz -stand-call-conf 30

done < MPA.list
```

### 5. Filter SNPs
Vcftools Version: 0.1.16 <br>
Work done on info2020

**Download the latest release of vcftools at https://github.com/vcftools/vcftools. Follow instructions at https://vcftools.github.io/examples.html for installation.**

#### Run vcftools (with --max-meanDP 230)
https://vcftools.github.io/man_latest.html 
```
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/home/xingyuan/tools/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 230 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt1
done
```
#### Run vcftools (with --max-meanDP 150)
```
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/home/xingyuan/tools/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 150 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt2
done
```
#### Run vcftools (with --max-meanDP 200)
```
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/home/xingyuan/tools/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 200 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt3
done
```
#### Run vcftools (with --max-meanDP 250)
```
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/home/xingyuan/tools/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 250 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt4
done
```
**After this steps, Rht_016_N, Rht_074_C, Rht_156_N, Rht_173_C, Rht_325_C, Rht_462_C, Rht_527_N, Rht_559_C, Rht_773_N, Rht_861_C have no SNPs for analysis in all four options of max-meanDP (10 with no data)**

### 6. VariantsToTable
https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable

GATK Version: 4.4.0.0 <br>
Work done on info2020

```
#!/bin/bash
for i in genotype_Rht_*vcf*; do
j=${i%.filt*vcf}
ref=${j#genotype_}

/home/xingyuan/tools/gatk-4.4.0.0/gatk VariantsToTable -V "$i" -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".fasta -F CHROM -F POS -F REF -F ALT -F QUAL -F AF -F ANN -F DP -GF GT -O "$i".table
done
```

### 7. Find genes at the positions of SNPs
https://github.com/bedops/bedops

#### Covert gff and vcf to bed format
Bedops Version: 2.4.41 <br>
Work done on info2020

**Prokka gff files**
```
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/genes_presence_absence/prokka/*gff; do
j=${i#/home/xingyuan/rhizo_ee/genes_presence_absence/prokka/}
sample=${j%-contigs.gff}

gff2bed < "$i" > "$sample".bed
done
```

**Pgap gff files**
```
gff2bed < annot_with_genomic_fasta.gff > annot_with_genomic_fasta.bed
```
**Vcf files**
```
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/call_snps/*vcf; do
j=${i#/home/xingyuan/rhizo_ee/call_snps/}
sample=${j%.recode.vcf}

vcf2bed < "$i" > "$sample".bed
done
```

#### Run bedtools 
Bedtools Version: 2.19.1 <br>
Work done on info114


```
bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/pgap/Rht_056_N/annot_with_genomic_fasta.bed -b genotype_Rht_056_N.filt1.bed -header -wa >  genotype_Rht_056_N.filt1.genes.2
```

**Find genes based on prokka annotation**
```
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/call_snps/*bed; do
j=${i#/home/xingyuan/rhizo_ee/call_snps/genotype_}
ref=${j%.filt?.bed}
k=${i#/home/xingyuan/rhizo_ee/call_snps/}
sample=${k%.bed}

bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/prokka/"$ref".bed -b "$i" -header -wa > "$sample".prokka_genes
done
```






