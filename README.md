# Experimental Evolution 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on **info server** (contact Brian Golding at golding@mcmaster.ca for an info account). **Compute canada server** will be used if the info server cannot run the program (register for an compute canada account at https://ccdb.alliancecan.ca/security/login). Results produced by compute canada server will be transferred to info server. 

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

# Samples used in this project
- **Original strains**: starting populations at the start of this experimental evolution experiment. Number of strains is 56. I receive complete genomes for these strains.
- **Derived strains**: derived populations at the end of this experimental evolution experiment. Number of strains is 363. I receive Illumina paired-end reads for these strains.

# Step 1 - Trim reads for derived strains
## 1. Run FastQC to check quality of raw reads 
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Versions: FastQC v0.12.1 <br>
Work done on info114

```bash
nohup /2/scratch/batstonelab/bin/FastQC/fastqc --outdir /home/xingyuan/rhizo_ee/fastQC_raw_reads --threads 5 /home/xingyuan/rhizo_ee/raw_reads/*gz &
```

## 2. Run MultiQC to combine all FastQC reports to a single file
https://multiqc.info/

Version: 1.20.dev0 <br>
Work done on info2020

```bash
/2/scratch/batstonelab/bin/multiqc /home/xingyuan/rhizo_ee/fastQC_raw_reads --outdir multiqc_raw_reads --verbose
```

Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, *Bioinformatics*, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354

## 3. Run fastp to trim the reads for all 363 derived strains
https://github.com/OpenGene/fastp

Version: 0.23.4 <br>
Work done on info114

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/raw_reads/*R1*; do 
R1=${i#/home/xingyuan/rhizo_ee/raw_reads/}
R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz} 
R1_P=${R1//001.fastq.gz/P_001.fastq.gz} 
R1_UP=${R1//001.fastq.gz/UP_001.fastq.gz} 
R2_P=${R2//001.fastq.gz/P_001.fastq.gz} 
R2_UP=${R2//001.fastq.gz/UP_001.fastq.gz}
MERGE=${R1//R1_001.fastq.gz/merged_001.fastq.gz} 
sample=${R1%_*_L002_*gz}

/2/scratch/batstonelab/bin/fastp --in1 /home/xingyuan/rhizo_ee/raw_reads/"$R1" --out1 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R1_P" --unpaired1 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R1_UP" --in2 /home/xingyuan/rhizo_ee/raw_reads/"$R2" --out2 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R2_P" --unpaired2 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R2_UP" --merge --merged_out /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$MERGE" --failed_out /home/xingyuan/rhizo_ee/fastp_results/fastp_failed_reads/"$sample".failed.fastq.gz --dont_overwrite --qualified_quality_phred 15 --unqualified_percent_limit 40  --detect_adapter_for_pe --cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --correction --trim_poly_g --poly_g_min_len 10 --trim_poly_x --poly_x_min_len 10 --overrepresentation_analysis --html /home/xingyuan/rhizo_ee/fastp_results/fastp_logs/"$sample".html --json /home/xingyuan/rhizo_ee/fastp_results/fastp_logs/"$sample".json --report_title "$sample" --thread 5 --verbose

done
```

Shifu Chen. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

## 4. Run FastQC to check quality of trimmed reads 
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Versions: FastQC v0.12.1 <br>
Work done on info114

```bash
nohup /2/scratch/batstonelab/bin/FastQC/fastqc --outdir /home/xingyuan/rhizo_ee/fastQC_trimmed_reads --threads 5 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/*gz &
```

## 5. Run MultiQC to combine all FastQC reports to a single file
Version: 1.20.dev0 <br>
Work done on info2020

```bash
/2/scratch/batstonelab/bin/multiqc /home/xingyuan/rhizo_ee/fastQC_trimmed_reads --outdir multiqc_trimmed_reads --verbose
```

Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, *Bioinformatics*, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354

# Step 2 - Genome assembly
## 1. Run SPAdes to assemble the trimmed reads into genomes for derived strains
https://github.com/ablab/spades

Version: v3.15.5 <br>
Work done on info2020

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/*merged*gz; do
merged=${i#/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/}
R1_P=${merged//merged_001.fastq.gz/R1_P_001.fastq.gz}
R1_UP=${merged//merged_001.fastq.gz/R1_UP_001.fastq.gz}
R2_P=${merged//merged_001.fastq.gz/R2_P_001.fastq.gz}
R2_UP=${merged//merged_001.fastq.gz/R2_UP_001.fastq.gz}
sample=${merged%_*_L002_*gz}

/2/scratch/batstonelab/bin/SPAdes-3.15.5-Linux/bin/spades.py -1 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R1_P" -2 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R2_P" --merged /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$merged" -s /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R1_UP" -s /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$R2_UP" --isolate --threads 5 -o /home/xingyuan/rhizo_ee/spades_genomes/"$sample"
done
```

Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes de novo assembler. Current Protocols in Bioinformatics, 70, e102. doi: 10.1002/cpbi.102


## 2. Run Quast to check the quality of scaffolds 
https://github.com/ablab/quast

Version: 5.2.0, 3d87c606 <br>
Work done on info2020

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/*scaffolds.fasta; do
file=${i#/2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/}
sample=${file%-scaffolds.fasta}

/2/scratch/batstonelab/bin/quast-5.2.0/quast.py $i --min-contig 0 --threads 5 --output-dir /home/xingyuan/rhizo_ee/quast_genomes_quality_check/"$sample"

done
```

# Step 3 - Find the most probable ancestor for each derived strain 
## 1. Run FastANI to calculate pairwise whole-genome average nucleotide identity between derived strains and original strains
https://github.com/ParBLiSS/FastANI

Version: 1.32 <br>
Work done on info2020

Query = 363 derived strains <br>
Reference = 56 original strains 

```bash
nohup /usr/local/bin/fastANI --ql query.txt --rl reference.txt --threads 5 --matrix -o most_prob_ancestors.txt &
```

# Step 4 - Find gene presence absence variations 
## 1. Annotate genome using prokka
https://github.com/tseemann/prokka

Prokka Version: 1.14.6 <br>
Work done on info2020

```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

/home/xingyuan/tools/miniconda3/bin/prokka --outdir "$sample" --prefix "$sample" --addgenes --locustag "$sample" --kingdom Bacteria --gcode 11 --proteins GCF_011604505.1_Strain23B_protein.gbk --evalue 1e-9 --rfam --cpus 5 "$i"
done
```

## 2. Create pangenome using roary
https://github.com/sanger-pathogens/Roary/blob/master/README.md

**Based on ANI values, remove 4_4_10, Rht_773_N, as5_2_4 since they have low ANI (below 90)**

Roary Version: 3.12.0 <br>
Work done on info2020

```bash
nohup /home/xingyuan/tools/miniconda3/envs/roary/bin/roary -v -p 5 -f roary_results -y *gff &
```

```bash
#Create graphs using instructions from https://github.com/microgenomics/tutorials/blob/master/pangenome.md
python roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv
```

# Step 5: Call SNPS between each derived strain and its most probable ancestor
This step is based on https://github.com/rtdoyle/how-rhizobia-evolve/blob/master/Variant%20discovery/Variant_calling.md

## 1. Run BWA 
https://bio-bwa.sourceforge.net/

**Useful tools: decoding SAM flags: https://broadinstitute.github.io/picard/explain-flags.html**

BWA Version: 0.7.17-r1188 <br>
Samtools Version: 1.13 (using htslib 1.13) <br>
Work done on info2020

### Prepare a csv file as the following: Derived_Strains,Most_Probable_Ancestor

```bash
10_1_1,Rht_460_C
10_1_5,Rht_462_C
10_1_8,Rht_460_C
.
.
.
```

### Index 56 genomes of original strains for bwa (Create files ending with .amb, .ann, .bwt, .pac, .sa)
```bash
#!/bin/bash
for i in Rht*fasta; do
/usr/bin/bwa index $i
done
```

### Run bwa (bwa mem -t 5 -M -R), convert .sam to .bam (samtools view -huS), and produce alignment statistics (samtools stats)
Samtools stats: http://www.htslib.org/doc/samtools-stats.html

```bash
#!/bin/bash
# Usage: ./ThisScript file.csv
#csv file has the following format:
#-----------------------------------------
#Derived_Strains,Most_Probable_Ancestor
#10_1_1,Rht_460_C
#10_1_5,Rht_462_C
#10_1_8,Rht_460_C
#.
#.
#.
#-----------------------------------------
while IFS=',' read -r a b; do # a=derived_strain,b=most_probable_ancestor
r1=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$a"_*R1_P*
r2=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$a"_*R2_P*
ref=/home/xingyuan/rhizo_ee/SNPS/"$b".fasta

/usr/bin/bwa mem -t 5 -M -R "@RG\tID:"$a"\tSM:"$a $ref $r1 $r2 | /usr/local/bin/samtools view -huS -o $a-"$b".bam - 

done < $1
```

## 2. Run picard to manipulate the SAM files
Picard Version: 3.0.0 <br>
Work done one info2020

### (1) Reorder reads in the BAM file to match the contig ordering in reference file
https://gatk.broadinstitute.org/hc/en-us/articles/360037426651-ReorderSam-Picard-

### Create sequence dictionary file (.dict) for the 56 original strains (Required for ReorderSam)
https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-

```bash
#!/bin/bash
for i in Rht*fasta; do
ref=${i%.fasta}
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar CreateSequenceDictionary -R "$ref".fasta -O "$ref".dict
done
```

### Run ReorderSam
```bash
#!/bin/bash
for i in *bam; do
sample=${i%.bam}
ref=${sample#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar ReorderSam -R "$ref".fasta -I "$sample".bam -O "$sample".reordered.bam -SD "$ref".dict
done
```

### (2) Assign all the reads in a file to a single new read-group
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-

```bash
#!/bin/bash
for i in *.reordered.bam; do
sample=${i%.reordered.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar AddOrReplaceReadGroups -I "$sample".reordered.bam -O "$sample".new_rg.bam -ID "$sample" -LB rhizo_ee -PL Illumina -PU 1 -SM "$sample" 
done
```

### (3) Sort the input BAM file by coordinate
https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-

```bash
#!/bin/bash
for i in *.new_rg.bam; do
sample=${i%.new_rg.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar SortSam -I "$sample".new_rg.bam -O "$sample".coordinate_sorted.bam -SO coordinate
done
```

### (4) Identify duplicate reads and index the BAM files
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard- <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-

```bash
#!/bin/bash
for i in *.coordinate_sorted.bam; do
sample=${i%.coordinate_sorted.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar MarkDuplicates -I "$sample".coordinate_sorted.bam -O "$sample".marked_duplicates.bam -M "$sample".marked_dup_metrics.txt && /scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar BuildBamIndex -I "$sample".marked_duplicates.bam
done
```

**Create alignment statistics**
Samtools Version: 1.9 (using htslib 1.9) <br>

```br
#!/bin/bash
for i in *.marked_duplicates.bam; do
sample=${i%.marked_duplicates.bam}

/home/xingyuan/tools/miniconda3/envs/samtools/bin/samtools stats $i > "$sample".stats
done
```

**Plot the results**

```br
#!/bin/bash
for i in *stats; do
sample=${i%.stats}

/home/xingyuan/tools/miniconda3/envs/samtools/bin/plot-bamstats -p "$sample".output "$i"
done
```

## 3. Run HaplotypeCaller
Samtools Version: 1.13 (using htslib 1.13) <br>
GATK Version: 4.4.0.0 <br>
Work done on info2020

### Create index file (.fai) for the 56 original strains for GATK
http://www.htslib.org/doc/samtools-faidx.html

```bash
#!/bin/bash
for i in Rht*fasta; do
/usr/local/bin/samtools faidx $i
done
```

### HaplotypeCaller
https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller

```bash
#!/bin/bash
for i in *.marked_duplicates.bam; do
sample=${i%.marked_duplicates.bam}
ref=${sample#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar HaplotypeCaller -R "$ref".fasta -I "$i" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O "$sample".g.vcf.gz -ERC GVCF
done
```

## 4. CombineGVCFs, GenotypeGVCFs
CombineGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832710975771-CombineGVCFs <br>
GenotypeGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs

GATK Version: 4.4.0.0 <br>
Work done on info2020

### Create a list, named ``MPA_list.txt``,  for the most probable ancestors (26 most probable ancestors)
```bash
Rht_016_N
Rht_056_N
Rht_074_C
Rht_097_N
Rht_108_C
Rht_113_C
Rht_116_N
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
### Run CombineGVCFs and GenotypeGVCFs
```bash
#!/bin/bash
while read Rht; do

# Group derived strains with the same most probable ancestors to a list. This should create 26 lists.
find *"$Rht"*.vcf.gz > MPA-"$Rht".list &&

# Run CombineGVCFs to combine the vcf.gz files in each list to one vcf.gz file. This should create 26 cohort.g.vcf.gz files.
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar CombineGVCFs -R "$Rht".fasta --variant MPA-"$Rht".list -O "$Rht".cohort.g.vcf.gz &&

# Run GenotypeGVCFs on each 26 cohort.g.vcf files
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar GenotypeGVCFs -R "$Rht".fasta -V "$Rht".cohort.g.vcf.gz -ploidy 1 -O genotype_"$Rht".vcf.gz -stand-call-conf 30

done < MPA_list.txt
```

## 5. Filter SNPs
https://vcftools.github.io/man_latest.html 

Vcftools Version: 0.1.16 <br>
Work done on info2020

```bash
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 150 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt1

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 200 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt2

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 250 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt3
done
```

### Check which samples have no data left for analysis
```
grep "No data left for analysis" *log
#Files with no data: Rht_016_N, Rht_074_C, Rht_113_C, Rht_156_N, Rht_173_C, Rht_325_C, Rht_462_C, Rht_527_N, Rht_559_C, Rht_773_N, Rht_861_C
```

## 6. VariantsToTable
https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable

GATK Version: 4.4.0.0 <br>
Work done on info2020

```bash
#!/bin/bash
for i in genotype_Rht_*.recode.vcf; do
j=${i%.filt*vcf}
ref=${j#genotype_}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar VariantsToTable -V "$i" -R "$ref".fasta -F CHROM -F POS -F REF -F ALT -F QUAL -F AF -F ANN -F DP -GF GT -O "$i".table
done
```

## 7. Find genes at the positions of SNPs
https://github.com/bedops/bedops

### Covert gff and vcf to bed format
Bedops Version: 2.4.41 <br>
Work done on info2020

*Note: Add the bin containing the gff2bed and vcf2bed to .bashrc file. Without this step, this error may appear "convert2bed: command not found".*

**Convert gff of original strains to bed format**
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_prokka/Rht*gff; do

/home/xingyuan/tools/bin/gff2bed < "$i" > "$i".bed
done
```

**Convert vcf to bed**
```bash
for i in /home/xingyuan/rhizo_ee/SNPS/final_vcf/*recode.vcf; do
j=${i#/home/xingyuan/rhizo_ee/SNPS/final_vcf/}
sample=${j%.recode.vcf}

/home/xingyuan/tools/bin/vcf2bed < "$i" > "$sample".bed
done
```

### Find genes at the positions of SNPs or indels  
Bedtools Version: v2.31.1 <br>
Work done on info2020

```bash
for i in /home/xingyuan/rhizo_ee/SNPS/final_vcf/*recode.vcf; do
j=${i#/home/xingyuan/rhizo_ee/SNPS/final_vcf/genotype_}
ref=${j%.filt?.recode.vcf}
out=${j%.recode.vcf}

/2/scratch/batstonelab/bin/bedtools2/bin/bedtools closest -a "$i" -b /home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_prokka/"$ref".gff.bed > "$out".genes
done
```

## 8. Do the following steps in R

**Load these libraries**
```r
library(utils)
library(stringr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(viridis)
```

**Put all the .recode.vcf.table from step 6 in a folder**

```r
#Set working directory to location of files
setwd("/Users/xingyuansu/Desktop/rhizo_ee/SNPS")
```

**Load data**
```r
#Load metadata
treatments <- read.csv("/Users/xingyuansu/Desktop/rhizo_ee/metadata_363EE2008isolates.csv") %>%
  select(sample_ID_seq, Treatment, soil_slurry_sample) %>%
  mutate(sample_ID_seq = gsub("[[:space:]]", "", sample_ID_seq)) #remove white space

#Load SNPs data
#3 different filters with different --max-meanDP (maximum mean depth values). Filt1: --max-meanDP 150, Filt2: --max-meanDP 200, Filt3: --max-meanDP 250

#Use --max-meanDP 150 for the following analysis

genotypes_filt1_list <- vector(mode = "list", length = length(list.files(path = ".", pattern = "filt1")))
names(genotypes_filt1_list) <- list.files(path = ".", pattern = "filt1")

#Read vcf tables and store them in the list  
for (i in names(genotypes_filt1_list)) {
  genotypes_filt1_list[[i]] <- read.table(i, header = TRUE, tryLogical = FALSE)
} 

#Find files with no SNPs, this means no SNPs between the original strain (most probable ancestor) and its derived isolates
empty_files <- c()
  
for (i in 1:length(names(genotypes_filt1_list))) {
  if (nrow(genotypes_filt1_list[[i]]) == 0){
    empty_files[i] <- names(genotypes_filt1_list)[i]
  }
}

#Check files with no SNPs
empty_files[!is.na(empty_files)]

#Remove the files that have no SNPs
genotypes_filt1_list2 <- genotypes_filt1_list[!names(genotypes_filt1_list) %in% empty_files[!is.na(empty_files)]]
```

**Format the data**
```r
#Reformat the vcf tables: the row name is name of the derived isolate, each column is a SNP
ReformatFunction <- function(df) { #write a function for formatting one data frame
  df2 <- df %>% 
    rowwise %>%
    summarise(across(.cols = -c(1:8),
                     .fns = ~paste0(CHROM, "/", REF, "-", POS, "-", .x)))
  as.data.frame(t(df2))
}

#Use a loop to format all data frames in the list
genotypes_filt1_list3 <- vector(mode = "list", length = length(genotypes_filt1_list2))

names(genotypes_filt1_list3) <- names(genotypes_filt1_list2)

for (i in 1:length(genotypes_filt1_list2)) {
  genotypes_filt1_list3[[i]] <- ReformatFunction(genotypes_filt1_list2[[i]])
}

genotypes_filt1_df <- genotypes_filt1_list3 %>% #merge all data frames in the list to one
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2), .) %>%
  rownames_to_column("derived_isolates")

#Count number of SNPs shared and not shared, remove non-SNPs, and further format the data
genotypes_filt1_df2 <- genotypes_filt1_df %>%
  pivot_longer(cols = c(2:ncol(genotypes_filt1_df)), names_to = "Old_colname", values_to = "SNP") %>%
  filter(!is.na(SNP)) %>%
  select(-Old_colname) %>%
  mutate(isolate = str_extract(derived_isolates, "[:digit:]+_[:digit:]+_[:alnum:]+"),
         reference = str_extract(derived_isolates, "Rht_[:digit:]+_."),
         SNP = paste0(reference, "/", SNP)) %>%
  mutate(V1 = str_extract(SNP, "(?<=/)[ATCG]+"),   #remove non-SNPs, e.g. cluster_001_consensus/C-4720654-C. Check this for regular expression: https://stringr.tidyverse.org/articles/regular-expressions.html
         V2 = str_extract(SNP, "[ATCG]+$")) %>%
  filter(V1 != V2) %>%
  select(-c(V1,V2)) %>%
  group_by(SNP) %>% 
  mutate(isolate_count = n(), .after = SNP) %>% #count number of times a SNP appears. If a SNP appears more than one time, it is a shared SNP
  mutate(Shared = case_when(isolate_count == 1 ~ "Not Shared",
                            isolate_count > 1 ~ "Shared")) 
```

**Create final data set with only shared SNPs**
```r
#Keep shared SNPs which are found in more than 1 derived isolates
genotypes_filt1_df3 <- genotypes_filt1_df2 %>%
  filter(isolate_count > 1) %>%
  left_join(treatments, by = c("isolate" = "sample_ID_seq")) %>%
  select(-c(derived_isolates, isolate, reference, Shared))
```

**Plotting**
```r
#plot treatment
genotypes_filt1_df3 %>%
  count(SNP, Treatment) %>%
  ggplot(aes(SNP, n, fill=Treatment)) +
  geom_col() +
  geom_text(aes(label=n), position = position_stack(vjust=0.5), colour="white", size=5) +
  scale_fill_viridis(discrete = T) +
  labs(title="Filt1: --max-meanDP 150", y="Number of derived isolates") +
  scale_y_continuous(breaks = seq(0,25,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.text = element_text(colour="black")) 

#plot pot
genotypes_filt1_df3 %>%
  count(SNP, soil_slurry_sample) %>%
  ggplot(aes(SNP, n, fill=soil_slurry_sample)) +
  geom_col() +
  geom_text(aes(label=n), position = position_stack(vjust=0.5), colour="white", size=5) +
  scale_fill_viridis(discrete = T) +
  labs(title="Filt1: --max-meanDP 150", y="Number of derived isolates") +
  scale_y_continuous(breaks = seq(0,25,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.text = element_text(colour="black")) 
```

**Put all .genes files from step 7 in a folder**

**Add a column of the closest gene next to the SNP to the final data set**
```r
#### Find closest genes next to the SNPs ####

setwd("/Users/xingyuansu/Desktop/rhizo_ee/Closest_genes_near_SNPs")

genes_filt1_list <- vector(mode = "list", length = length(list.files(path = ".", pattern = "filt1")))
names(genes_filt1_list) <- list.files(path = ".", pattern = "filt1")

#Read bed files and store them in the list  
for (i in names(genes_filt1_list)) {
  if (file.info(i)$size == 0) next 
  genes_filt1_list[[i]] <- read.table(i, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="", tryLogical = FALSE)  %>%
    select(1, 2, 4, 5, ncol(.)-8, ncol(.)-7, ncol(.)) %>%
    rename("CHROM"=1, "POS"=2, "REF"=3, "ALT"=4, "Start"=5, "End"=6, "Gene"=7) %>%
    mutate(SNP = paste0(CHROM, "/", REF, "-", POS, "-", ALT), .before=CHROM) %>%
    filter(grepl("product=", Gene))
} 

#remove empty elements in the list
genes_filt1_list2 <- genes_filt1_list[lapply(genes_filt1_list, length) > 0]

#merge all data frames in the list to one
genes_filt1_df <- genes_filt1_list2 %>% 
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2), .) %>%
  mutate(reference=str_extract(Gene, "Rht_[:digit:]+_[NC]"),
         SNP=paste0(reference, "/", SNP)) %>%
  select(-c(CHROM, POS, REF, ALT, reference, Start, End))


#combine the SNPs and genes data frames
SNPs_results <- genotypes_filt1_df3 %>%
  left_join(genes_filt1_df, by=c("SNP" = "SNP"))

#Import this data set
write.csv(SNPs_results, "~/Desktop/rhizo_EE_SNPs.csv", row.names = FALSE)
```
