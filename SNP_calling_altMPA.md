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
HTSJDK version: 3.0.5 <br>
VCFtools version: 0.1.16
 
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

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar ReorderSam -R /home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$ref".fasta -I "$sample".bam -O "$sample".reordered.bam -SD /home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$ref".dict
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
#!/bin/bash
for i in *.marked_duplicates.bam; do
sample=${i%.marked_duplicates.bam}
ref=${sample#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar HaplotypeCaller -R /home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$ref".fasta -I "$i" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O "$sample".g.vcf.gz -ERC GVCF
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

# Run CombineGVCFs to combine the vcf.gz files in each list to one vcf.gz file. This should create 6 cohort.g.vcf.gz files.
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar CombineGVCFs -R /home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$Rht".fasta --variant MPA-"$Rht".list -O "$Rht".cohort.g.vcf.gz &&

# Run GenotypeGVCFs on each 6 cohort.g.vcf files
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar GenotypeGVCFs -R /home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$Rht".fasta -V "$Rht".cohort.g.vcf.gz -ploidy 1 -O genotype_"$Rht".vcf.gz -stand-call-conf 30

done < AltMPA_list.txt
```

Filter SNPs
```bash
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 150 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt1

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 200 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt2

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 250 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt3
done
```

VariantsToTable
```bash
#!/bin/bash
for i in genotype_Rht_*.recode.vcf; do
j=${i%.filt*vcf}
ref=${j#genotype_}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar VariantsToTable -V "$i" -R /home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$ref".fasta -F CHROM -F POS -F REF -F ALT -F QUAL -F AF -F ANN -F DP -GF GT -O "$i".table
done
```




