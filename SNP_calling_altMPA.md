# SNP calling using alternative reference strains

info2020

## Align 56 reference strains to 56 reference strains

fastANI version 1.32
```bash
nohup /usr/local/bin/fastANI --ql reference.txt --rl reference.txt --threads 5 --matrix -o ref_to_ref.txt &
```

## SNP calling
BWA version: 0.7.17-r1188 <br>
Samtools version: 1.13 (using htslib 1.13) <br>
Picard version: 3.0.0 <br>
The Genome Analysis Toolkit (GATK) version: 4.4.0.0 <br>
HTSJDK version: 3.0.5
 
Index 6 reference genomes for bwa (Create files ending with .amb, .ann, .bwt, .pac, .sa)
```bash
for i in Rht*fasta; do
/usr/bin/bwa index $i
done
```

Run bwa (bwa mem -t 5 -M -R), convert .sam to .bam (samtools view -huS), and produce alignment statistics (samtools stats)
```bash
#!/bin/bash
# Usage: ./ThisScript Derived-AltMPA.csv
#First 5 lines of the Derived-AltMPA.csv looks like the following:
#-----------------------------------------
#10_3_2,Rht_415_C
#10_5_1,Rht_415_C
#10_5_8,Rht_415_C
#12_7_1,Rht_415_C
#12_7_4,Rht_415_C
#.
#.
#.
#-----------------------------------------
while IFS=',' read -r a b; do # a=derived_strain,b=most_probable_ancestor
r1=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$a"_*R1_P*
r2=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$a"_*R2_P*
ref=/home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$b".fasta

/usr/bin/bwa mem -t 5 -M -R "@RG\tID:"$a"\tSM:"$a $ref $r1 $r2 | /usr/local/bin/samtools view -huS -o $a-"$b".bam - 

done < $1
```

Run picard to manipulate the SAM files

Create sequence dictionary file (.dict) for the 6 references
```bash
for i in Rht*fasta; do
ref=${i%.fasta}
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar CreateSequenceDictionary -R "$ref".fasta -O "$ref".dict
done
```

Run ReorderSam
```bash
for i in *bam; do
sample=${i%.bam}
ref=${sample#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar ReorderSam -R "$ref".fasta -I "$sample".bam -O "$sample".reordered.bam -SD "$ref".dict
done
```

Assign all the reads in a file to a single new read-group
```bash
for i in *.reordered.bam; do
sample=${i%.reordered.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar AddOrReplaceReadGroups -I "$sample".reordered.bam -O "$sample".new_rg.bam -ID "$sample" -LB rhizo_ee -PL Illumina -PU 1 -SM "$sample" 
done
```

Sort the input BAM file by coordinate
```bash
for i in *.new_rg.bam; do
sample=${i%.new_rg.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar SortSam -I "$sample".new_rg.bam -O "$sample".coordinate_sorted.bam -SO coordinate
done
```

Identify duplicate reads and index the BAM files
```bash
for i in *.coordinate_sorted.bam; do
sample=${i%.coordinate_sorted.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar MarkDuplicates -I "$sample".coordinate_sorted.bam -O "$sample".marked_duplicates.bam -M "$sample".marked_dup_metrics.txt && /scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar BuildBamIndex -I "$sample".marked_duplicates.bam
done
```

Create index file (.fai) for the 6 references for GATK
```bash
for i in Rht*fasta; do
/usr/local/bin/samtools faidx $i
done
```

Run HaplotypeCaller
```bash
for i in *.marked_duplicates.bam; do
sample=${i%.marked_duplicates.bam}
ref=${sample#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar HaplotypeCaller -R "$ref".fasta -I "$i" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O "$sample".g.vcf.gz -ERC GVCF
done
```

Run CombineGVCFs and GenotypeGVCFs
```bash
#!/bin/bash
#The AltMPA_list.txt looks like the following:
#-----------------------------------------
#Rht_003_C
#Rht_061_N
#Rht_116_N
#Rht_415_C
#Rht_511_N
#Rht_596_N
#-----------------------------------------
while read Rht; do

# Group derived strains with the same alternative MPAs to a list. This should create 6 lists.
find *"$Rht"*.vcf.gz > MPA-"$Rht".list &&

# Run CombineGVCFs to combine the vcf.gz files in each list to one vcf.gz file. This should create 26 cohort.g.vcf.gz files.
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar CombineGVCFs -R "$Rht".fasta --variant MPA-"$Rht".list -O "$Rht".cohort.g.vcf.gz &&

# Run GenotypeGVCFs on each 26 cohort.g.vcf files
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar GenotypeGVCFs -R "$Rht".fasta -V "$Rht".cohort.g.vcf.gz -ploidy 1 -O genotype_"$Rht".vcf.gz -stand-call-conf 30

done < AltMPA_list.txt
```

Filter SNPs
```bash

```

# GC content across genomes 

info2020

Create .genome file with the size of each scaffold
```bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

faidx $i -i chromsizes > "$sample".genome
done
```

Create .bed file with 200 Kb window
```bash
for i in *genome; do
sample=${i%.genome}

/2/scratch/batstonelab/bin/bedtools2/bin/bedtools makewindows -g $i -w 200000 > "$sample".bed
done
```

Calculate GC content
```bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

/2/scratch/batstonelab/bin/bedtools2/bin/bedtools nuc -fi $i -bed "$sample".bed > "$sample"_gc.bed
done
```

# Identify plamid sequences in the draft genome assembly
https://github.com/HubertTang/PLASMe?tab=readme-ov-file

Note the program will move your genome file to its directory. Create a symbolic link of all the genome files for the program specifically

```bash
ln -s /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta /home/xingyuan/rhizo_ee/testing/plasmid_identification
```

```bash
conda activate plasme
```

```bash
#!/bin/bash
for x in *.fasta; do
dir=${x%.fasta}

python /home/xingyuan/tools/PLASMe/PLASMe.py --database /home/xingyuan/tools/PLASMe/DB --coverage 0.9 --identity 0.9 --probability 0.5 --thread 5 "$x" "$x".plasme.fna;

mkdir "$dir";

mv temp/ "$x"* "$dir"
done
```






