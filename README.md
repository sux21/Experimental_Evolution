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

# Step 4 - Find gene presence absence variations (Method 1)
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

**Based on ANI values, remove 4_4_10, Rht_773_N, as_5_2_4 since they have low ANI (below 90)**

Roary Version: 3.12.0

```bash
nohup /2/scratch/batstonelab/bin/Roary-3.13.0/bin/roary -p 5 -f roary_results -y *gff &
```

# Step 4 - Find gene presence absence variations (Method 2)
## 1. Filter sequences shorter than 200 bp (pgap only takes sequences equal or longer than 200 bp)
https://github.com/shenwei356/seqkit

Seqkit Version: v2.7.0 <br>
Work done on info2020

```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
file_in=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
file_out=${file_in//.fasta/.filtered.fasta}

/2/scratch/batstonelab/bin/seqkit seq --min-len 200 --threads 5 /home/xingyuan/rhizo_ee/derived+original_genomes/"$file_in" > /home/xingyuan/rhizo_ee/genes_presence_absence_variation/"$file_out"
done
```

**W Shen**, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. ***PLOS ONE***. doi:10.1371/journal.pone.0163962.

## 2. Annotate genome using pgap (use 6 scripts because pgap is slow, and run each script in a different directory)

Version: 2023-05-17.build6771 <br>
Work done on graham cluster 

**script 1**
```bash
#!/bin/bash
#SBATCH --time=04-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/{1..5}_*filtered.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%-scaffolds.filtered.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 2**
```bash
#!/bin/bash
#SBATCH --time=04-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/{6..9}_*filtered.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%-scaffolds.filtered.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 3**
```bash
#!/bin/bash
#SBATCH --time=04-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/1{0..5}_*filtered.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%-scaffolds.filtered.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 4**
```bash
#!/bin/bash
#SBATCH --time=04-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/1{6..9}_*filtered.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%-scaffolds.filtered.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 5**
```bash
#!/bin/bash
#SBATCH --time=04-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/Rht*fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%.filtered.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 6**
```bash
#!/bin/bash
#SBATCH --time=02-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/20_*filtered.fasta /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/as5_2_4*filtered.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%-scaffolds.filtered.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```

### Find samples that don't have annot.gbk: 
```bash
find . -type d \! -exec test -e '{}/annot_with_genomic_fasta.gff' \; -print
#Output: 4_4_10, 2_4_11, 19_4_7, 19_1_9, Rht_773_N
```

### Run PGAP for the failed samples (4_4_10, 2_4_11, 19_4_7, 19_1_9, Rht_773_N) using the option `--ignore-all-errors`)
```bash
#!/bin/bash
#SBATCH --time=01-30:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/4_4_10-scaffolds.filtered.fasta /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/2_4_11-scaffolds.filtered.fasta /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/19_4_7-scaffolds.filtered.fasta /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/19_1_9-scaffolds.filtered.fasta /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/Rht_773_N.filtered.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}

if [[ "$j" =~ "Rht_773_N" ]]; then
 sample=${j%.filtered.fasta}
else
 sample=${j%-scaffolds.filtered.fasta}
fi

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample"_draft -g "$i" -s 'Rhizobium leguminosarum' -c 40 --ignore-all-errors
done
```

## 3. Reformat PGAP's gff files for roary
### Rename each annot_with_genomic_fasta.gff with sample names
```bash
for i in */annot_with_genomic_fasta.gff; do
sample=${i%/annot_with_genomic_fasta.gff}

ln -s $i "$sample"_annot_with_genomic_fasta.gff
done
```

### Reformat PGAP's annot_with_genomic_fasta.gff for AGAT 
```bash
#!/bin/bash
#There are two things needed to be fixed:
#----------------First Thing: This is done by steps 1 to 3--------------------------
##sequence-region  1 1046175

needs to be changed into:

##sequence-region NODE_1_length_1046175_cov_26.396254  1 1046175
#-----------------------------------------------------------------------------------
#----------------Second Thing-------------------------------------------------------

## FASTA
------------------
#Run this script as ./ThisScript INPUT.gff
# Step 1: Extract all lines starting with "NODE" 
while IFS=$'\t' read -r a b; do
if [[ "$a" =~ ^"NODE" ]]; then
  echo $a >> "$1".int1
fi
done < $1

### Step 2: Remove redundant lines
cat -n "$1".int1 | sort -k2 -k1n  | uniq -f1 | sort -nk1,1 | cut -f2- > "$1".int2

### Step 3: Reformat "##sequence-region  1 970887" into "##sequence-region NODE_1_length_970887_cov_33.311050  1 970887"
i=0; while read -r line; do
  if [[ $line =~ ^"##sequence-region" ]]; then
    let i=$i+1
    new_info=`sed -n "$i,$i p" "$1".int2`
    echo "$line" | sed "s/##sequence-region/& "$new_info"/g"
  else
    echo "$line"
  fi
done < $1 > "$1".int3

### Step 4: Reformat fasta header 
sed -e 's/lcl|//' -e 's/Rhizobium leguminosarum chromosome, whole genome shotgun sequence//' "$1".int3 > "$1".fixed

### Step 5: Remove intermediate files
rm -f "$1".int*

#References:
#Step 2 command was taken from https://unix.stackexchange.com/questions/194780/remove-duplicate-lines-while-keeping-the-order-of-the-lines. Step 3 commands were taken from  https://stackoverflow.com/questions/72293847/using-sed-in-order-to-change-a-specific-character-in-a-specific-line, https://stackoverflow.com/questions/67396536/sed-insert-whitespace-at-the-nth-position-after-a-delimiter, https://stackoverflow.com/questions/50971418/using-sed-to-replace-single-line-in-while-read-loop
```

### Reformat gff file using AGAT

AGAT Version: v1.4.0 <br>
Work done on info19

```bash
for i in *fixed; do
/2/scratch/batstonelab/bin/AGAT/bin/agat_convert_sp_gxf2gxf.pl --gff $i --output ${i%_annot_with_genomic_fasta.gff.fixed}.pgap.gff
done
```

## 4. Run roary
https://sanger-pathogens.github.io/Roary/

Roary version: <br>
Work done on info2020

```
nohup /2/scratch/batstonelab/bin/Roary-3.13.0/bin/roary -p 5 -f roary_results -y *gff &
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

/usr/bin/bwa mem -t 5 -M -R "@RG\tID:"$a"\tSM:"$a $ref $r1 $r2 | /usr/local/bin/samtools view -huS -o $a-"$b".bam - && /usr/local/bin/samtools stats $a-"$b".bam > $a-"$b".stats

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

## 3. Run HaplotypeCaller
Samtools Version: 1.13 (using htslib 1.13) <br>
GATK Version: 4.4.0.0 <br>
Work done on info2020

### Create index file (``.fai``) of reference for GATK
http://www.htslib.org/doc/samtools-faidx.html
```bash
#!/bin/bash
for i in Rht*fasta; do
samtools faidx $i
done
```

### HaplotypeCaller
https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller
```bash
#!/bin/bash
for i in *.marked_duplicates.bam; do
sample=${i%.marked_duplicates.bam}
ref=${sample#*-}

/home/xingyuan/tools/gatk-4.4.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".fasta -I "$i" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O "$sample".g.vcf.gz -ERC GVCF
done
```

## 4. CombineGVCFs, GenotypeGVCFs
CombineGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832710975771-CombineGVCFs <br>
GenotypeGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs

GATK Version: 4.4.0.0 <br>
Work done on info2020

**Create a list, named ``MPA.list``,  for the most probable ancestors (25 most probable ancestors, 25 lines)**
```bash
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
### Run CombineGVCFs and GenotypeGVCFs
```bash
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

## 5. Filter SNPs
Vcftools Version: 0.1.16 <br>
Work done on info2020

**Download the latest release of vcftools at https://github.com/vcftools/vcftools. Follow instructions at https://vcftools.github.io/examples.html for installation.**

https://vcftools.github.io/man_latest.html 

```bash
#!/bin/bash
for i in genotype*gz; do
out=${i%.vcf.gz}

/home/xingyuan/tools/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 230 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt1
done
```

**After this steps, Rht_016_N, Rht_074_C, Rht_156_N, Rht_173_C, Rht_325_C, Rht_462_C, Rht_527_N, Rht_559_C, Rht_773_N, Rht_861_C have no SNPs (10 with no data)**

## 6. VariantsToTable
https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable

GATK Version: 4.4.0.0 <br>
Work done on info2020

```bash
#!/bin/bash
for i in genotype_Rht_*vcf*; do
j=${i%.filt*vcf}
ref=${j#genotype_}

/home/xingyuan/tools/gatk-4.4.0.0/gatk VariantsToTable -V "$i" -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".fasta -F CHROM -F POS -F REF -F ALT -F QUAL -F AF -F ANN -F DP -GF GT -O "$i".table
done
```

## 7. Find genes at the positions of SNPs
https://github.com/bedops/bedops

#### Covert gff and vcf to bed format
Bedops Version: 2.4.41 <br>
Work done on info2020

**Convert gff of original strains to bed**
```bash
#!/bin/bash
for i in Rht*gff; do

/home/xingyuan/tools/bin/gff2bed < "$i" > "$i".bed
done
```

**Convert vcf to bed**
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/call_snps/*vcf; do
j=${i#/home/xingyuan/rhizo_ee/call_snps/}
sample=${j%.recode.vcf}

vcf2bed < "$i" > "$sample".bed
done
```

#### Find genes at the positions of SNPs or indels  
Bedtools Version: 2.19.1 <br>
Work done on info114

```bash
#!/bin/bash
for i in genotype*.filt1.bed; do
j=${i#genotype_}
ref=${j%.filt1.bed}

/usr/local/bedtools/2.19.1/bin/bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/pgap/"$ref".gff.bed -b "$i" -header -wa > "$i".genes
done
```

## 8. Find alternate alleles that occur in more than 1 evolved strains 

**Example: the alternate allele G occurs in more than one evolved strains. Convergent evolution?**
```bash
CHROM   POS     REF     ALT     QUAL    AF      ANN     DP      10_1_1-Rht_460_C.GT     10_1_8-Rht_460_C.GT     10_7_9-Rht_460_C.GT     11_3_10-Rht_460_C.GT 11_3_8-Rht_460_C.GT     12_2_8-Rht_460_C.GT     13_4_5-Rht_460_C.GT     14_4_7-Rht_460_C.GT     14_5_7-Rht_460_C.GT     14_5_9-Rht_460_C.GT     15_2_10-Rht_460_C.GT    15_4_4-Rht_460_C.GT     15_5_10-Rht_460_C.GT    16_1_10-Rht_460_C.GT    16_1_5-Rht_460_C.GT     16_1_7-Rht_460_C.GT     16_4_10-Rht_460_C.GT    16_6_8-Rht_460_C.GT     17_2_8-Rht_460_C.GT     17_2_9-Rht_460_C.GT     17_4_3-Rht_460_C.GT     18_1_2-Rht_460_C.GT     19_1_2-Rht_460_C.GT     1_5_2-Rht_460_C.GT      1_5_4-Rht_460_C.GT      20_2_1-Rht_460_C.GT     2_2_3-Rht_460_C.GT      2_2_5-Rht_460_C.GT      2_5_3-Rht_460_C.GT      2_6_4-Rht_460_C.GT      2_7_3-Rht_460_C.GT      3_2_10-Rht_460_C.GT     3_2_3-Rht_460_C.GT      3_3_1-Rht_460_C.GT      4_1_3-Rht_460_C.GT      4_1_5-Rht_460_C.GT      4_2_7-Rht_460_C.GT      4_6_8-Rht_460_C.GT      5_3_1-Rht_460_C.GT      5_3_7-Rht_460_C.GT      6_4_2-Rht_460_C.GT      6_4_7-Rht_460_C.GT      6_4_9-Rht_460_C.GT      6_6_2-Rht_460_C.GT      6_6_9-Rht_460_C.GT      7_7_2-Rht_460_C.GT      8_4_4-Rht_460_C.GT      9_1_7-Rht_460_C.GT      9_2_1-Rht_460_C.GT      9_3_1-Rht_460_C.GT 9_5_10-Rht_460_C.GT     9_5_9-Rht_460_C.GT
cluster_001_consensus   1135251 A       G       47102.7 0.365   NA      1415    A       G       A       A       A       A       A       G       G       A A       G       G       A       G       G       A       G       A       A       G       G       A       A       A       G       A       G       G A       A       A       A       A       A       A       A       A       A       A       G       A       G       G       G       G       G       A A       A       A       A
```
**Run this script for each tab-delimited table from step 6**
```bash
#!/bin/bash
#Usage: ./ThisScript input.file
{
IFS=$'\t' read A B C D E F G H I
while IFS=$'\t' read a b c d e f g h i; do 

# Count number of samples having the alternate allele
number=$(echo "$i" | grep -o $(printf '\t')$d | grep -c .)  

# Skip line when only one sample has the alternate allele
if [[ $number -eq "1" ]]; then
   continue
fi

# Print the results
echo Position: $b
echo "alternate allele:" $d, "sample's alleles:" $i "("$I")"

done
} < $1 > "$1"-alternate_alleles
```

**BUG: There is a bug in this script. See the situation below: the alternate allele G only occurs once in the samples, but it is picked up by the script.**
```bash
alternate allele: G, sample's alleles: GC G GC GC (11_3_9-Rht_837_C.GT 1_1_3-Rht_837_C.GT 1_3_8-Rht_837_C.GT 7_4_8-Rht_837_C.GT)
```
