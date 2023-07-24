# 2020 Experimental Evolution 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on info server. Compute canada server will be used if the info server cannot run the program. Results produced by compute canada server will be transferred to info server. 

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

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
Version: <br>
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

Methods are from https://github.com/microgenomics/tutorials/blob/master/pangenome.md.

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

### 2 (Method 1). Roary
https://sanger-pathogens.github.io/Roary/

Version: 1.007001 <br>
Work done on info114

```
/usr/local/bin/roary -p 5 *.gff
```

### 2 (Method 2): GenAPI
https://github.com/MigleSur/GenAPI

Version: 1.0 <br>
Work done on info2020

```
nohup genapi --threads 5 --matrix *.gff &
```

## Analysis 3 (Method 3): Presence and Absence of Genes

### 1. PGAP
https://github.com/ncbi/pgap/tree/1126_Test 

```
# Create YAML file 

fasta:
   class: File
   location: K_pneu.fasta
submol:
   class: File
   location: K_pneu_meta.yaml

# Create metadata file

topology: 'linear'
organism:
   genus_species: 'Klebsiella pneumoniae'
contact_info:
   last_name: 'NAME'
   first_name: 'NAME'
   email: 'EMAIL'
   organization: '...'
   department: '...'
   street: '...'
   city: '...'
   postal_code: '...'
   state: '...'
   country: '...'
authors:
- author:
   first_name: '...'
   last_name: '...'

# Annotate

./pgap.py -r -o output_dir file.yaml

```
### 2. BWA and Samtools
https://github.com/lh3/bwa

https://broadinstitute.github.io/picard/explain-flags.html

https://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/

BWA Version: 0.7.17-r1188 <br>
Samtools Version: 1.11 <br>
makeblastdb Version: 2.13.0+ <br>
blastn Version: 2.13.0+ <br>
Work done on info114

```
Reference genomes: 56 original strains
Query sequences: Illumina reads

# Index reference genomes
for i in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/Rht*fasta; do bwa index $i; done

# Run BWA, Samtools
bwa mem /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/Rht_108_C.fasta /home/xingyuan/rhizo_ee/trimmomatic_reads/6_4_5_ATATGCATGT-CCAGGCACCA_L002_R1_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/6_4_5_ATATGCATGT-CCAGGCACCA_L002_R2_P_001.fastq.gz | samtools view -S -b > Rht_108_C-6_4_5.bam

# Get unmapped reads
samtools fasta -f 4 Rht_108_C-6_4_5.bam > Rht_108_C_unmapped_reads.fasta

Check number of reads:
# Count total reads: samtools view -c file.bam
# Count unmapped reads: samtools view -f 4 -c file.bam
# Count mapped reads: samtools view -c -F 260 file.bam

# Make blast database for original strains
for i in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/Rht*fasta; do j=${i#/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/}; k=${j%.fasta}; /usr/local/blast/2.13.0+/bin/makeblastdb -in $i -title "$k" -dbtype nucl; done

# Align unmapped reads from Rht_108_C to other original strains using blastn
for i in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/Rht*fasta; do j=${i#/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/}; k=${j%.fasta}; /usr/local/blast/2.13.0+/bin/blastn -qcov_hsp_perc 70 -query Rht_108_C_unmapped_reads.fasta -out Rht_108_C_unmapped_reads_blast_to_$k.out -db $i; done

/usr/local/blast/2.13.0+/bin/blastn -query Rht_108_C_unmapped_reads.fasta -out Rht_108_C_unmapped_reads_blast_to_$k.out -db /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/Rht_108_C.fasta
```

## Analysis 4: Characterize plasmids
Methods are taken from [Symbiosis genes show a unique pattern of introgression and selection within a *Rhizobium leguminosarum* species complex](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7276703/). 

### 1. tblastn
Version: 2.9.0+ <br>
Work done on info114

The following sequences are copied and pasted to 8 separate files as query sequences for tblastn: <br>
Table S4: https://figshare.com/articles/dataset/Gene_alignments_and_SNP_matrices_of_a_Rhizobium_complex/11568894/5?file=21156690
<img width="1081" alt="Screenshot 2023-07-14 at 7 39 54 PM" src="https://github.com/sux21/2020_Experimental_Evoluntion/assets/132287930/24dce3c3-6f02-46a0-9a4d-bad414242dd3">

```
# Make assemblies as blast database

makeblastdb -in 10_1_1-contigs.fasta -title "10_1_1-contigs" -dbtype nucl 
```

```
# Run tblastn

for i in RepA-Rh0{1..8}.fasta; do
a=${i#RepA-}
b=${a%.fasta}

/usr/local/bin/tblastn -query $i -db /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/10_1_1-contigs.fasta -qcov_hsp_perc 70 -out 10_1_1-$b.blast
done
```

