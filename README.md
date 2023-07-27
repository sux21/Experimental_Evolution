# 2020 Experimental Evolution 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on info server. Compute canada server will be used if the info server cannot run the program. Results produced by compute canada server will be transferred to info server. 

# Monday meeting
- Discuss graphing the data in R. Show the script which genapi uses to graphs the gene presence-absence matrix.
- For roary, should I disable split paralogs because genapi does not split paralogs?

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

# Samples used in this project
- Original strains: starting populations at the start of this experimental evolution experiment. Number of strains is 56.
- Experimental evolved strains: derived populations at the end of this experimental evolution experiment. Number of strains is 363.
# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
### 1. Run FastQC for raw data
Versions: FastQC v0.11.5 <br>
Work done on info114

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
## Analysis 1 - Find the most related original strain for each experimentally evolved strain 
Commands in step 1 are taken from https://github.com/Alan-Collins/Spine-Nucmer-SNPs. 

**Samples: 363 experimentally evolved strains + 56 original strains**

### (Optional) 1. Run Spine: find core genomes 
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

## Analysis 2 - Phylogeny for the original strains 
Commands in step 1 are taken from https://github.com/Alan-Collins/Spine-Nucmer-SNPs. 

**Samples: 56 original strains**

### 1. Run Spine-Nucmer-SNPs
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

### (Compare the results with genapi) 2. Roary
https://sanger-pathogens.github.io/Roary/

Version: 1.007001 <br>
Work done on info114

#### Run roary
```
/usr/local/bin/roary -p 5 *.gff
```
#### Visualize in R
```
# Import gene_presence_absence.Rtab file
pav <- read.delim('~/Desktop/gene_presence_absence.Rtab', row.names = 1, header = TRUE)

library(pheatmap)
pheatmap(pav,col= c('green', 'blue'), legend_breaks = c(T,F), legend_labels = c("Present", "Absent"), treeheight_row = 0, treeheight_col = 0)
```
### (Compare the results with roary) 2. GenAPI
https://github.com/MigleSur/GenAPI

Version: 1.0 <br>
Work done on info2020

```
nohup genapi --threads 5 --matrix *.gff &
```

## Analysis 4: Call SNPS between each evolved strain and its most probable ancestor based on ANI values in Analysis 1
https://github.com/rtbatstone/how-rhizobia-evolve/blob/master/Variant%20discovery/Variant_calling.md

Install the latest version of picard at https://github.com/broadinstitute/picard. Use this command as an example to download: ``wget https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar``.

### 1. BWA
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

### 2. Reorder BAM input file according to the order of contigs in a second reference file
https://gatk.broadinstitute.org/hc/en-us/articles/360037426651-ReorderSam-Picard-

Picard Version: 3.0.0 <br>
Work done one info114

```
#!/bin/bash
for i in *bam; do
sample=${i%.bam}
ref=${sample#*-}
java -jar /home/xingyuan/tools/picard.jar ReorderSam INPUT="$sample".bam OUTPUT="$sample".reordered.bam REFERENCE="$ref".fasta
```

### 3. Assign all the reads in a file to a single new read-group
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
```
#!/bin/bash
for i in *.reordered.bam; do
sample=${i%.reordered.bam}
java -jar /home/xingyuan/tools/picard.jar AddOrReplaceReadGroups I="$sample".reordered.bam O="$sample".rg.bam RGLB=rhizo_ee RGPL=Illumina RGPU=1 RGSM="$sample"
done
```

### 4. Sort the input BAM file by coordinate
https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-
```
#!/bin/bash
for i in *.rg.bam; do
sample=${i%.rg.bam}
java -jar /home/xingyuan/tools/picard.jar SortSam I="$sample".rg.bam O="$sample".coordinate_sorted.bam SORT_ORDER=coordinate
done
```

### 5. Identify duplicate reads
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
```
#!/bin/bash
for i in *.coordinate_sorted.bam; do
sample=${i%.coordinate_sorted.bam}
java -jar /home/xingyuan/tools/picard.jar MarkDuplicates I="$sample".coordinate_sorted.bam O="$sample".marked_duplicates.bam M=marked_dup_metrics.txt
done
```

### 6. Generate BAM index ".bai" file
https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-
```
#!/bin/bash
for i in *.marked_duplicates.bam; do
java -jar /home/xingyuan/tools/picard.jar BuildBamIndex I=$i
done
```

### 7. Create dictionary file (``.dict``) and index file (``.fai``) for gatk tools.
https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format#:~:text=We%20use%20the%20faidx%20command,coordinates%20in%20the%20FASTA%20file.

#### Create index file (``.fai``) for reference using samtools
http://www.htslib.org/doc/samtools-faidx.html
```
#!/bin/bash
for i in Rht*fasta; do
samtools faidx $i
done
```
#### Create dictionary file (``.dict``) for reference using picard.jar 
https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard- <br> 

```
#!/bin/bash
for i in Rht*fasta; do
ref=${i%.fasta}
java -jar /home/xingyuan/tools/picard.jar CreateSequenceDictionary R="$ref".fasta O="$ref".dict
done
```





