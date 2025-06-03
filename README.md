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

**Notes: In the html files, quality encoding is Sanger / Illumina 1.9, which phred=ord(b)-33**

## 2. Run MultiQC to combine all FastQC reports to a single file
https://github.com/MultiQC/MultiQC

Version: 1.28 <br>
Work done on info2020

```bash
/home/xingyuan/tools/miniconda3/bin/multiqc /home/xingyuan/rhizo_ee/fastQC_raw_reads --outdir multiqc_raw_reads --verbose
```

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

## 4. Run FastQC to check quality of trimmed reads 
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Versions: FastQC v0.12.1 <br>
Work done on info114

```bash
nohup /2/scratch/batstonelab/bin/FastQC/fastqc --outdir /home/xingyuan/rhizo_ee/fastQC_trimmed_reads --threads 5 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/*gz &
```

## 5. Run MultiQC to combine all FastQC reports to a single file
Version: 1.28 <br>
Work done on info2020

**for paired reads**
```bash
/home/xingyuan/tools/miniconda3/bin/multiqc /home/xingyuan/rhizo_ee/fastQC_trimmed_reads/*_P_* --outdir multiqc_trimmed_paired --verbose
```

**for single unpaired reads**
```bash
/home/xingyuan/tools/miniconda3/bin/multiqc /home/xingyuan/rhizo_ee/fastQC_trimmed_reads/*_UP_* --outdir multiqc_trimmed_unpaired --verbose
```

**for merged reads**
```bash
/home/xingyuan/tools/miniconda3/bin/multiqc /home/xingyuan/rhizo_ee/fastQC_trimmed_reads/*_merged_* --outdir multiqc_trimmed_merged --verbose
```

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

**combine QUAST results to one file**
```bash
cat */transposed_report.tsv | sed '2,${/^Assembly*/d;}' > quast_assembly_check.tsv
```

## 3. Run CheckM to check completeness and contamination of the assembly
CheckM Version: 1.2.3 <br>
Work done on info2020

**Following installing instruction from https://github.com/Ecogenomics/CheckM/wiki/Installation#system-requirements (If pysam fails to install, try https://anaconda.org/bioconda/pysam).**

```bash
#inform checkm where is the downloaded reference data
/home/xingyuan/tools/miniconda3/bin/checkm data setRoot /home/xingyuan/tools/checkm_data_2015_01_16
```

**Run CheckM with the recommended Lineage-specific Workflow using lineage-specific marker sets (https://github.com/Ecogenomics/CheckM/wiki/Workflows)** 
```
nohup /home/xingyuan/tools/miniconda3/bin/checkm lineage_wf -x fasta -t 1 /2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes 363EEgenomes_CheckM_Results &
```

## 4. Split 19_1_9 and 19_4_7 genomes into 2 because they likely contain two different genomes based on bimodal GC content distribution and CheckM 100% completeness and >= 100% contamination
SemiBin2 Version: 2.2.0 <br>
Bowtie2 Version: 2.5.4 <br>
Samtools Version: 1.13 <br> 
Work done on info2020

**Map reads to scaffolds and sort** <br>
Follow steps of "Mapping using bowtie2" from https://semibin.readthedocs.io/en/latest/generate/. 

```bash
for i in 19_1_9 19_4_7; do
mkdir -p -v index/"$i"

/home/xingyuan/tools/bowtie2-2.5.4-linux-x86_64/bowtie2-build -f /2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/"$i"*fasta index/"$i"/"$i"

done
```

```bash
for i in 19_1_9 19_4_7; do
/home/xingyuan/tools/bowtie2-2.5.4-linux-x86_64/bowtie2 -q --fr -x index/"$i"/"$i" -1 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$i"*R1_P_*fastq.gz -2 /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$i"*R2_P_*fastq.gz -S "$i".sam -p 4
done
```

```bash
for i in 19_1_9 19_4_7; do
/usr/local/bin/samtools view -h -b -S "$i".sam | samtools view -b -F 4 - | samtools sort - -o "$i".mapped.sorted.bam

/usr/local/bin/samtools index "$i".mapped.sorted.bam
done
```

**Run SemiBin2 with single-sample binning and self-supervised mode**
```bash
for i in 19_1_9 19_4_7; do

/home/xingyuan/tools/miniconda3/bin/SemiBin2 single_easy_bin --self-supervised -i /2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/"$i"*fasta -b "$i".mapped.sorted.bam -o "$i"_output --threads 5 --engine cpu --compression none --random-seed 12345

done
```


Results: both 19_1_9 and 19_4_7 scaffolds were separated to 4 bins
```
filename        nbps    nr_contigs      N50     L50
19_1_9_output/output_bins/SemiBin_0.fa  7134363 34      373995  6
19_1_9_output/output_bins/SemiBin_1.fa  1903653 31      122512  6
19_1_9_output/output_bins/SemiBin_2.fa  2948191 49      104221  10
19_1_9_output/output_bins/SemiBin_3.fa  360297  16      52740   3

######################################################################

filename        nbps    nr_contigs      N50     L50
19_4_7_output/output_bins/SemiBin_0.fa  4468896 19      463747  4
19_4_7_output/output_bins/SemiBin_1.fa  6916484 21      596258  6
19_4_7_output/output_bins/SemiBin_2.fa  474384  7       258841  1
19_4_7_output/output_bins/SemiBin_3.fa  358859  14      35161   3
```

**Run CheckM on binned scaffolds** <br>
CheckM Version: 1.2.3 <br>
Work done on info2020

Make a new directory for reconstructed bins 
```bash
mkdir reconstructed_bins
cd reconstructed_bins
```

Create symbolic links to the reconstructed bins
```bash
for i in /home/xingyuan/rhizo_ee/split_genomes/*_output/output_bins/*fa; do
j=${i#/home/xingyuan/rhizo_ee/split_genomes/}
filename=${j/output\/output_bins\/}

ln -s $i $filename
done
```

Run CheckM
```bash
nohup /home/xingyuan/tools/miniconda3/bin/checkm lineage_wf -x fa -t 1 /home/xingyuan/rhizo_ee/split_genomes/reconstructed_bins reconstructed_bins_CheckM_Results &
```

Results: 19_4_7 contained Rhizobiaceae and *Bacillus*. 19_1_9 contained Rhizobiaceae, Bacillaceae, and Bacteria. 19_4_7_SemiBin_1 and 19_1_9_SemiBin_0 would be used as true isolates.
```
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                   Marker lineage        # genomes   # markers   # marker sets    0     1    2   3   4   5+   Completeness   Contamination   Strain heterogeneity  
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  19_4_7_SemiBin_1   f__Rhizobiaceae (UID3591)       50         838           326         1    835   2   0   0   0       99.94            0.34               0.00          
  19_1_9_SemiBin_0   f__Rhizobiaceae (UID3591)       50         838           326         1    835   2   0   0   0       99.94            0.34               0.00          
  19_4_7_SemiBin_0      g__Bacillus (UID856)        101         682           227         3    671   7   1   0   0       98.97            2.96               0.00          
  19_1_9_SemiBin_2    f__Bacillaceae (UID829)       128         559           183        295   260   4   0   0   0       50.23            0.96               0.00          
  19_1_9_SemiBin_1      k__Bacteria (UID203)        5449        104            58         73    31   0   0   0   0       48.28            0.00               0.00          
  19_4_7_SemiBin_3          root (UID1)             5656         56            24         56    0    0   0   0   0        0.00            0.00               0.00          
  19_4_7_SemiBin_2          root (UID1)             5656         56            24         56    0    0   0   0   0        0.00            0.00               0.00          
  19_1_9_SemiBin_3          root (UID1)             5656         56            24         56    0    0   0   0   0        0.00            0.00               0.00          
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

## 2. Check the most probable ancestor of 19_4_7_SemiBin_1 and 19_1_9_SemiBin_0
Query = 19_4_7_SemiBin_1, 19_1_9_SemiBin_0 <br>
Reference = 56 original strains 

```bash
nohup /usr/local/bin/fastANI --ql 19_X_X_query.txt --rl reference.txt --threads 5 --matrix -o 19_X_X_most_prob_ancestors.txt &
```

# Step 4 - Find gene presence absence variations 

## 1. Annotate genome using Bakta
https://github.com/oschwengers/bakta?tab=readme-ov-file#database-download

Bakta Version: 1.11.0 <br>
Work done on info20

### Bakta requires python version >=3.9 and <3.12
```bash 
conda create -n py39 python=3.9 #create conda environment for python 3.9
conda activate py39 #activate the environment
python --version #check python version. I have Python 3.9.22
conda install -c conda-forge -c bioconda bakta=1.11 #install latest version of bakta
```

### Download full database
```bash
/home/xingyuan/tools/miniconda3/envs/py39/bin/bakta_db download --output /home/xingyuan/tools --type full #database version: 6.0
```

### Run Bakta for all samples except as5_2_4 since it may not be *Rhizobium leguminosarum* (Remember to run ``conda activate py39``)
```bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

/home/xingyuan/tools/miniconda3/envs/py39/bin/bakta --db /home/xingyuan/tools/db --prefix 10_1_1 --output 10_1_1 --genus Rhizobium --species leguminosarum  --keep-contig-headers --verbose --threads 5 /home/xingyuan/rhizo_ee/derived+original_genomes/10_1_1-scaffolds.fasta
done
```

## 2. Gene presence absence analysis using Panaroo
https://gthlab.au/panaroo/#/

**Based on ANI values, remove 4_4_10, Rht_773_N, as5_2_4 since they have low ANI (below 90), also remove 19_1_9 and 19_4_7 for now until the cleaned genomes are annotated - 414 total strains**

Panaroo Version: 1.5.2 <br>
Work done on info2020

```bash
#create a new directory for the results
mkdir panaroo_results

#run panaroo 
nohup /home/xingyuan/tools/miniconda3/bin/panaroo -i *gff -o panaroo_results --clean-mode strict --threshold 0.9 &

#filter potential pseudo genes, genes with unusual lengths, and fragmented genes
#/home/xingyuan/tools/miniconda3/bin/panaroo-filter-pa -i ./gene_presence_absence.csv -o ./ --type pseudo,length
```

## 3. Align DNA sequences of genes gained to each of the ancestral strains

Blastn Version: 2.16.0 <br>
Work done on info2020

Previous steps are done in R.

```bash
#specify format of output: qacc (Query accession), sacc (Subject accession), evalue (Expect value), bitscore (Bit score), length (Alignment length), pident (Percentage of identical matches), nident (Number of identical matches), gapopen (Number of gap openings), gaps (Total number of gaps), qcovs (Query Coverage Per Subject)

for i in /home/xingyuan/rhizo_ee/derived+original_genomes/Rht*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

/home/xingyuan/tools/ncbi-blast-2.16.0+/bin/blastn -query prokka_genes_gain.fasta -subject $i -outfmt "10 qacc sacc evalue bitscore length pident nident gapopen gaps qcovs" > "$sample"_blast.csv

done

#add variable names for the file
for i in *blast.csv; do
sed -i '1s/^/qacc,sacc,evalue,bitscore,length,pident,nident,gapopen,gaps,qcovs\n/' $i
done
```


## 4. Separate gene presence absence analysis for 4_4_10 and Rht_773_N

Panaroo Version: 1.5.2 <br>
Work done on info2020

```bash
#create a new directory for 4_4_10 and Rht_773_N
mkdir pav_for_4410and773
cd pav_for_4410and773

#rename files
for i in */*.gff; do
sample=${i%/*.gff};
ln -s $i "$sample".gff
done

#run panaroo 
nohup /home/xingyuan/tools/miniconda3/bin/panaroo -i *gff -o panaroo_results --clean-mode strict --threshold 0.9 &
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

### Find nearby genes of SNPs or indels  
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

## 8. Align raw long reads to corresponding MPA's genome assembly 

Minimap2 Version: 2.17-r941 <br>
Samtools Version: 1.13 <br>
Work done on info2020

**Align PacBio reads to corresponding MPA's genome (-x map-hifi for PacBio HiFi/CCS genomic reads (v2.19+))**
```bash
#!/bin/bash
for i in *fastq; do
sample=${i%.fastq}

/home/xingyuan/tools/minimap2-2.29_x64-linux/minimap2 -a -x map-hifi /home/xingyuan/rhizo_ee/derived+original_genomes/"$sample".fasta $i > "$sample"_aln.sam
done
```

```bash
chmod u+x run_minimap2.sh
nohup ./run_minimap2.sh &
```

**Convert SAM format to BAM format and create corresponding index file (.bai)**

Converting SAM to BAM
```bash
for i in *sam; do
sample=${i%_aln.sam}

/usr/local/bin/samtools view -S -b $i > "$sample".bam
done
```

Order alignments based upon their alignment coordinates on each chromosome
```bash
for i in *bam; do
sample=${i%.bam}

/usr/local/bin/samtools sort $i -o "$sample".sorted.bam
done
```

Index the BAM files
```bash
for i in *.sorted.bam; do
/usr/local/bin/samtools index $i
done
```

**Visualization in IGV**

-load reference genome (fasta formatted file): Genomes > Load Genome from File <br>
-load BAM file with its corresponding index file (.bam and .bai): File > Load from File






