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

## 2. Genome assembly QC
### Run Quast to check the quality of scaffolds 
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
cat */transposed_report.tsv | sed '2,${/^Assembly*/d;}' > quast_assembly_stats.tsv
```

### Run CheckM to check completeness and contamination of the assembly
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

For citations:

Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043–1055.

Matsen FA, Kodner RB, Armbrust EV. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: doi:10.1186/1471-2105-11-538.

Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28: 2223–2230.

http://hmmer.org/

### Run PGAP taxonomy check
Version: 2025-05-06.build7983 <br>
Work done on info20

Remove sequences less than 200 bp because PGAP does not accept these
```bash
#!/bin/bash
for i in /2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/*scaffolds.fasta; do
j=${i#/2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/}
output=${j//.fasta/.filtered.fasta}

/2/scratch/batstonelab/bin/seqkit seq --min-len 200 --threads 5 "$i" > ./"$output"
done
```

Run PGAP taxonomy check
```bash
#!/bin/bash
for i in *.filtered.fasta; do
sample=${i%.filtered.fasta}

/home/xingyuan/tools/pgap.py -D /usr/bin/apptainer --container-path /home/xingyuan/tools/pgap_2025-05-06.build7983.sif --report-usage-false -o "$sample" -g "$i" -s "Rhizobium leguminosarum" --cpu 8 --no-self-update --taxcheck-only
done
```

Extract ANI statuses to a tab separated file
```bash
#!/bin/bash
# this script will extract ANI statuses from ani-tax-report.txt for all 363 derived isolates. Output is a tab separated file.

#first line of the file consists of column names
head -n 6 10_1_1-scaffolds/ani-tax-report.txt | cut -d : -f 1 | tr "\n" "\t" > pgap_ani_status.txt

echo "" >> pgap_ani_status.txt # add newline character to the end of line

#remaining 363 lines consist of ANI status for each isolate
for i in *-scaffolds/ani-tax-report.txt; do
  head -n 6 "$i" | sed 's/: /\t/'| cut -f 2 | tr "\n" "\t" >> pgap_ani_status.txt

  echo "" >> pgap_ani_status.txt # add newline character
done
```

Extract ANI statistics to a tab separated file
```bash
#!/bin/bash

# this script will extract ANI tables from ani-tax-report.txt for all 363 derived isolates. Output is a tab separated file.

#first line is the column names

sed -n '18p' 10_1_1-scaffolds/ani-tax-report.txt | sed 's/./\t/8' | sed 's/./\t/22' | sed 's/./\t/31' | sed 's/./\t/40' | sed 's/./\t/50' | sed 's/./\t/54' | sed 's/ //g' | sed "s/^/ANI report for assembly\t/" > pgap_ani_table.txt

#remaining lines are ANI tables of all 363 derived isolates.

for i in *-scaffolds/ani-tax-report.txt; do
  file=`head -n 1 "$i" | sed 's/: /:/' | cut -d : -f 2`

  tail -n +20 $i | sed 's/./\t/8' | sed 's/./\t/22' | sed 's/./\t/31' | sed 's/./\t/40' | sed 's/./\t/50' | sed 's/./\t/54' | sed "s/^/$file\t/" >>  pgap_ani_table.txt
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

## 2. Contamination resolution and ancestral verification for 19_1_9 and 19_4_7 genomes which show bimodal GC content distribution and CheckM 100% completeness and >= 100% contamination
SemiBin2 Version: 2.2.0 <br>
Bowtie2 Version: 2.5.4 <br>
Samtools Version: 1.13 <br> 
Work done on info2020

### Use SemiBin2 to group DNA sequences from the same genome together
Map reads to scaffolds and sort <br>
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

Run SemiBin2 with single-sample binning and self-supervised mode.
```bash
for i in 19_1_9 19_4_7; do

/home/xingyuan/tools/miniconda3/bin/SemiBin2 single_easy_bin --self-supervised -i /2/scratch/batstonelab/N_adaptation_Rhizobium/2020_derived_strains_genomes/"$i"*fasta -b "$i".mapped.sorted.bam -o "$i"_output --threads 5 --engine cpu --compression none --random-seed 12345

done
```

Check binning results
```
(base) xingyuan@info20:~/rhizo_ee/split_genomes$ cat 19_1_9_output/recluster_bins_info.tsv
filename        nbps    nr_contigs      N50     L50
19_1_9_output/output_bins/SemiBin_0.fa  7134363 34      373995  6
19_1_9_output/output_bins/SemiBin_1.fa  1903653 31      122512  6
19_1_9_output/output_bins/SemiBin_2.fa  2948191 49      104221  10
19_1_9_output/output_bins/SemiBin_3.fa  360297  16      52740   3
(base) xingyuan@info20:~/rhizo_ee/split_genomes$ cat 19_4_7_output/recluster_bins_info.tsv
filename        nbps    nr_contigs      N50     L50
19_4_7_output/output_bins/SemiBin_0.fa  4468896 19      463747  4
19_4_7_output/output_bins/SemiBin_1.fa  6916484 21      596258  6
19_4_7_output/output_bins/SemiBin_2.fa  474384  7       258841  1
19_4_7_output/output_bins/SemiBin_3.fa  358859  14      35161   3
```

Both 19_1_9 and 19_4_7 scaffolds were separated to 4 bins


### Check which bin contains the true Rhizobium species (CheckM method)
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


### Check which bin contains the true Rhizobium species (PGAP Taxonomy Check method)
Version: 2025-05-06.build7983 <br>
Work done on info20

```bash
#!/bin/bash
for i in *.fa; do
sample=${i%.fa}

/home/xingyuan/tools/pgap.py -D /usr/bin/apptainer --container-path /home/xingyuan/tools/pgap_2025-05-06.build7983.sif --report-usage-false -o PGAP_taxo_check_"$sample" -g "$i" -s "Rhizobium leguminosarum" --cpu 8 --no-self-update --taxcheck-only
done
```

Extract ANI statuses to a tab separated file
```bash
#!/bin/bash
# this script will extract and combine ANI statuses from ani-tax-report.txt. Output is a tab separated file.

#first line of the file consists of column names
head -n 6 PGAP_taxo_check_19_1_9_SemiBin_0/ani-tax-report.txt | cut -d : -f 1 | tr "\n" "\t" > contamination_resolution_pgap_ani_status.txt

echo "" >> contamination_resolution_pgap_ani_status.txt # add newline character to the end of line

#remaining 363 lines consist of ANI status for each isolate
for i in PGAP_taxo_check*/ani-tax-report.txt; do
  head -n 6 "$i" | sed 's/: /\t/'| cut -f 2 | tr "\n" "\t" >> contamination_resolution_pgap_ani_status.txt

  echo "" >> contamination_resolution_pgap_ani_status.txt # add newline character
done
```

Extract ANI statistics to a tab separated file
```bash
!/bin/bash

# this script will extract and combine ANI tables from ani-tax-report.txt. Output is a tab separated file.

#first line is the column names

sed -n '18p' PGAP_taxo_check_19_1_9_SemiBin_0/ani-tax-report.txt | sed 's/./\t/8' | sed 's/./\t/22' | sed 's/./\t/31' | sed 's/./\t/40' | sed 's/./\t/50' | sed 's/./\t/54' | sed 's/ //g' | sed "s/^/ANI report for assembly\t/" > contamination_resolution_pgap_ani_table.txt

#remaining lines are ANI tables of all 363 derived isolates.

for i in PGAP_taxo_check*/ani-tax-report.txt; do
  file=`head -n 1 "$i" | sed 's/: /:/' | cut -d : -f 2`

  tail -n +20 $i | sed 's/./\t/8' | sed 's/./\t/22' | sed 's/./\t/31' | sed 's/./\t/40' | sed 's/./\t/50' | sed 's/./\t/54' | sed "s/^/$file\t/" >>  contamination_resolution_pgap_ani_table.txt
done
```

### Check the most probable ancestor of 19_4_7_SemiBin_1 and 19_1_9_SemiBin_0
Query = 19_4_7_SemiBin_1, 19_1_9_SemiBin_0 <br>
Reference = 56 original strains 

```bash
nohup /usr/local/bin/fastANI --ql 19_X_X_query.txt --rl reference.txt --threads 5 --matrix -o 19_X_X_most_prob_ancestors.txt &
```

# Step 4 - Gene gain and loss analysis
## 1A. Annotate genome using PGAP
### Use cleaned genomes for 19_1_9 and 19_4_7

```bash
ln -s /home/xingyuan/rhizo_ee/split_genomes/19_1_9_output/output_bins/SemiBin_0.fa 19_1_9_SemiBin_0.fasta
ln -s /home/xingyuan/rhizo_ee/split_genomes/19_4_7_output/output_bins/SemiBin_1.fa 19_4_7_SemiBin_1.fasta
```

### Filter sequences shorter than 200 bp (pgap only takes sequences equal or longer than 200 bp)
https://github.com/shenwei356/seqkit

Seqkit Version: v2.7.0 <br>
Work done on info20

```bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
file_in=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
file_out=${file_in//.fasta/.filtered.fasta}

/2/scratch/batstonelab/bin/seqkit seq --min-len 200 --threads 5 /home/xingyuan/rhizo_ee/derived+original_genomes/"$file_in" > /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/"$file_out"
done
```

**W Shen**, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. ***PLOS ONE***. doi:10.1371/journal.pone.0163962.

### Run pgap (This will take a long time, about 4 to 5 hours per sample)

Version: 2025-05-06.build7983 <br>
Work done on info20

```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/*.filtered.fasta; do
j=${i#/home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/}
sample=${j%.filtered.fasta}

/home/xingyuan/tools/pgap.py -D /usr/bin/apptainer --container-path /home/xingyuan/tools/pgap_2025-05-06.build7983.sif --report-usage-false -o "$sample" --prefix "$sample" -g "$i" -s "Rhizobium leguminosarum" --cpu 6 --no-self-update
done
```

Which samples failed to be annotated (don't have the *_with_genomic_fasta.gff file). Code taken from https://askubuntu.com/questions/196960/find-directories-that-dont-contain-a-file/196966#196966. 
```bash
find . -type d '!' -exec sh -c 'ls "{}"|egrep -i -q "_with_genomic_fasta.gff"' ';' -print
```

The folllowing samples failed to be annotated.
```bash
.
./Rht_156_N
./Rht_493_C
./Rht_328_N
./Rht_726_C
./Rht_773_N
./4_4_10-scaffolds
./2_4_11-scaffolds
./Rht_861_C
```

Rht_156_N, Rht_493_C, Rht_328_N, Rht_726_C, Rht_861_C all failed at the step of "Find_Best_Evidence_Alignments", which failed at the last sequence (cluster_005_consensus which is a plasmid). Remove the last sequence and annotate the first 4 sequences.

Split FASTA entry
```bash
#!/bin/bash
# A script to split each fasta entry
# run as './split_fasta.sh filename'
i=0
while read line; do
    if [[ $line =~ ">" ]]; then
        let i=$i+1
    fi
    if [[ $i -eq $i ]]; then
        echo $line >> "${1%.fasta}"."$i".fasta
    fi
done < $1
```

Run for all 5 samples
```bash
for i in Rht_156_N Rht_493_C Rht_328_N Rht_726_C Rht_861_C; do
  ./split_fasta.sh /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/"$i"*fasta
done
```

Concatenate the first 4 sequences.
```bash
for i in Rht_156_N Rht_493_C Rht_328_N Rht_726_C Rht_861_C; do
  cat /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/"$i".filtered.{1..4}.fasta > /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/"$i".filtered.1-4.fasta
done
```

Run PGAP again.
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/*.filtered.1-4.fasta; do
j=${i#/home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/}
sample=${j%.filtered.1-4.fasta}_partial

/home/xingyuan/tools/pgap.py -D /usr/bin/apptainer --container-path /home/xingyuan/tools/pgap_2025-05-06.build7983.sif --report-usage-false -o "$sample" --prefix "$sample" -g "$i" -s "Rhizobium leguminosarum" --cpu 6 --no-self-update
done
```

## 1B. Annotate genome using Bakta
https://github.com/oschwengers/bakta?tab=readme-ov-file#database-download

Bakta Version: 1.11.4 <br>
Work done on info20

### Bakta requires python version >=3.9 and <3.12
```bash 
conda create -n py39 python=3.9 #create conda environment for python 3.9
conda activate py39 #activate the environment
python --version #check python version. I have Python 3.9.22
conda install -c conda-forge -c bioconda bakta=1.11.4 #install latest version of bakta
```

### Download full database
```bash
/home/xingyuan/tools/miniconda3/envs/py39/bin/bakta_db download --output /home/xingyuan/tools --type full #database version: 6.0
```

### Run Bakta for all samples except as5_2_4 since it is not *Rhizobium leguminosarum*

```bash
conda activate py39
```

```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

#skip as5_2_4
if [[ "$sample" =~ "as5_2_4" ]]; then
continue
fi

/home/xingyuan/tools/miniconda3/envs/py39/bin/bakta --db /home/xingyuan/tools/db --prefix "$sample" --output "$sample" --genus Rhizobium --species leguminosarum  --keep-contig-headers --verbose --threads 8 "$i"
done
```

## 2B. Gene presence absence analysis using Panaroo
https://gthlab.au/panaroo/#/

**Based on ANI values, remove 4_4_10, Rht_773_N, as5_2_4 since they have low ANI (below 90), also remove 19_1_9 and 19_4_7 for now until the cleaned genomes are annotated - 414 total strains**

Panaroo Version: 1.5.2 <br>
Work done on info2020

```bash
#!/usr/bin/Rscript
#Run this script: /usr/bin/Rscript get_MPA_descendents_names.R

#load ANI data output by fastANI
ANI_dat <- read.table("E:/rhizobia_exp_evo_genetics/second_manuscript/MPA_data/most_prob_ancestors.txt")

#for each derived isolate, find the maximum ANI value
max_ANI <- aggregate(V3 ~ V1, data = ANI_dat, FUN = max)

#subset the ANI data to only include these maximum ANI values - these ancestral isolates are the most probable ancestors (MPAs)
MPA_dat <- subset(ANI_dat, paste0(V1, V3) %in% paste0(max_ANI$V1, max_ANI$V3))

#rename isolates
MPA_dat$V1 <- sub("-scaffolds.fasta", "", MPA_dat$V1, fixed = T)

MPA_dat$V2 <- sub(".fasta", "", MPA_dat$V2, fixed = T)

#get a list of MPA names
MPA_names <- data.frame(MPA_name=unique(MPA_dat$V2))

#output the list of MPA names
write.table(MPA_names, "E:/rhizobia_exp_evo_genetics/second_manuscript/MPA_names.csv", 
            row.names=F, quote=F, col.names=F)


#next, create 26 files containing descendents of each MPA

descendents <- vector(mode="list", length(MPA_names$MPA_name))

names(descendents) <- MPA_names$MPA_name

for (i in seq_along(MPA_names$MPA_name)) {
  descendents[[i]] <- subset(MPA_dat, V2 %in% MPA_names$MPA_name[i])[,1] # a vector of descendent names
  
  descendents[[i]] <- c(MPA_names$MPA_name[i], descendents[[i]]) #add MPA name to the vector
  
  #output as each name on one line
  
  file = paste0(MPA_names$MPA_name[i], "_descendents.txt")
  
  writeLines(c(descendents[[i]]), con = file)
}
```

**create a directory containing each MPA and its corresponding derived isolates (26 MPAs, 26 directories)**
```bash
for mpa in ; do
for i in `cat MPA_"$mpa".txt`; do
ls "$i".fasta
echo $mpa
done
done

```

```bash
#create a new directory for the results
mkdir panaroo_results

#run panaroo 
nohup /home/xingyuan/tools/miniconda3/bin/panaroo -i *gff -o panaroo_results --clean-mode strict --threshold 0.9 &

#filter potential pseudo genes, genes with unusual lengths, and fragmented genes
#/home/xingyuan/tools/miniconda3/bin/panaroo-filter-pa -i ./gene_presence_absence.csv -o ./ --type pseudo,length
```






# Step 5: Call SNPS between each derived strain and its most probable ancestor
This step is based on https://github.com/rtdoyle/how-rhizobia-evolve/blob/master/Variant%20discovery/Variant_calling.md

## 0. Filter the reads for the two derived isolates (19_4_7, 19_1_9). Recall the two derived isolates, 19_4_7 and 19_1_9, are contaminated. Align the trimmed reads to the split genomes and only use the aligned reads for subsequent steps.  
BWA version: 0.7.17-r1188 <br>
Samtools version: 1.13 (using htslib 1.13) <br>
bedtools version:2.31.1
Work done on info2020

Index the decontaminated genomes of 19_1_9 and 19_4_7 for bwa (Create files ending with .amb, .ann, .bwt, .pac, .sa).
```bash
for i in l{19_1_9,19_4_7}*SemiBin*fasta; do
/usr/bin/bwa index "$i" 
done
```

Use BWA to align reads to the split genomes.
```bash
for i in 19_1_9 19_4_7; do
r1=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$i"_*R1_P*
r2=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$i"_*R2_P*
ref=/home/xingyuan/rhizo_ee/derived+original_genomes/"$i"*SemiBin*fasta

/usr/bin/bwa mem -t 5 $ref $r1 $r2 | /usr/local/bin/samtools view -huS -o "$i".bam - 
done
```

Only keep reads that are mapped and paired. 
```bash
/usr/local/bin/samtools flags

#About: Convert between textual and numeric flag representation
#Usage: samtools flags FLAGS...
#
#Each FLAGS argument is either an INT (in decimal/hexadecimal/octal) representing
#a combination of the following numeric flag values, or a comma-separated string
#NAME,...,NAME representing a combination of the following flag names:
#
#   0x1     1  PAIRED         paired-end / multiple-segment sequencing technology
#   0x2     2  PROPER_PAIR    each segment properly aligned according to aligner
#   0x4     4  UNMAP          segment unmapped
#   0x8     8  MUNMAP         next segment in the template unmapped
#  0x10    16  REVERSE        SEQ is reverse complemented
#  0x20    32  MREVERSE       SEQ of next segment in template is rev.complemented
#  0x40    64  READ1          the first segment in the template
#  0x80   128  READ2          the last segment in the template
# 0x100   256  SECONDARY      secondary alignment
# 0x200   512  QCFAIL         not passing quality controls or other filters
# 0x400  1024  DUP            PCR or optical duplicate
# 0x800  2048  SUPPLEMENTARY  supplementary alignment
```

```bash
for i in *bam; do
/usr/local/bin/samtools view -u -f 1 -F 12 -o ${i/.bam/_mapped_and_paired.bam} $i #remove unmapped reads and singletons (one read maps and the other does not map)
done
```

Verify the reads are 100% mapped and 0% singletons
```bash
/usr/local/bin/samtools flagstat 19_1_9_mapped_and_paired.bam
#859371 + 0 in total (QC-passed reads + QC-failed reads)
#859206 + 0 primary
#0 + 0 secondary
#165 + 0 supplementary
#0 + 0 duplicates
#0 + 0 primary duplicates
#859371 + 0 mapped (100.00% : N/A)
#859206 + 0 primary mapped (100.00% : N/A)
#859206 + 0 paired in sequencing
#429603 + 0 read1
#429603 + 0 read2
#856918 + 0 properly paired (99.73% : N/A)
#859206 + 0 with itself and mate mapped
#0 + 0 singletons (0.00% : N/A)
#910 + 0 with mate mapped to a different chr
#874 + 0 with mate mapped to a different chr (mapQ>=5)

/usr/local/bin/samtools flagstat 19_4_7_mapped_and_paired.bam
#2386764 + 0 in total (QC-passed reads + QC-failed reads)
#2386370 + 0 primary
#0 + 0 secondary
#394 + 0 supplementary
#0 + 0 duplicates
#0 + 0 primary duplicates
#2386764 + 0 mapped (100.00% : N/A)
#2386370 + 0 primary mapped (100.00% : N/A)
#2386370 + 0 paired in sequencing
#1193185 + 0 read1
#1193185 + 0 read2
#2379888 + 0 properly paired (99.73% : N/A)
#2386370 + 0 with itself and mate mapped
#0 + 0 singletons (0.00% : N/A)
#1392 + 0 with mate mapped to a different chr
#1324 + 0 with mate mapped to a different chr (mapQ>=5)
```

Sort the reads by read name (-n option). Required by ``bedtools bamtofastq``. 
```bash
for i in *mapped_and_paired.bam; do
/usr/local/bin/samtools sort -n -o ${i/mapped_and_paired.bam/mapped_and_paired.qsort.bam} $i
done
```

Convert bam to fastq.
```bash
for i in *mapped_and_paired.qsort.bam; do
/home/xingyuan/tools/miniconda3/bin/bedtools bamtofastq -i $i -fq ${i/.bam/.R1.fastq} -fq2 ${i/.bam/.R2.fastq}
done
```

Verify the number of reads are correct.
```bash
grep -c "@A00419" 19_1_9_mapped_and_paired.qsort.R1.fastq
#429603
grep -c "@A00419" 19_1_9_mapped_and_paired.qsort.R2.fastq
#429603
grep -c "@A00419" 19_4_7_mapped_and_paired.qsort.R1.fastq
#1193185
grep -c "@A00419" 19_4_7_mapped_and_paired.qsort.R2.fastq
#1193185
```

Make a directory with symbolic link to the cleaned reads and other reads for BWA.
```bash
ln -s /home/xingyuan/rhizo_ee/fastp_results/fastp_reads/*_P_* /home/xingyuan/rhizo_ee/snp_indel/trimmed_paired_reads
rm -f 19_1_9* 19_4_7* #remove the old reads of the two contaminated samples
ln -s /home/xingyuan/rhizo_ee/fastp_results/clean_reads_for_19_X_X/*fastq /home/xingyuan/rhizo_ee/snp_indel/trimmed_paired_reads
ls 19_1_9* 19_4_7* #verify the reads are correct for these two samples
```

## 1. Align the trimmed reads of derived isolate to genomes of its most probable ancestor
https://bio-bwa.sourceforge.net/

BWA Version: 0.7.17-r1188 <br>
Samtools Version: 1.13 (using htslib 1.13) <br>
Work done on info2020

Prepare a CSV file with a pair of derived isolate and its most probable ancestor.
```r
#load FastANI results and metadata
ani_values <- read.table("./MPA_data/most_prob_ancestors.txt")

#select the ancestral isolate with the highest ANI (defined as the most probable ancestor in this study)
MPA <- ani_values %>%
  rename(isolate=V1, MPA=V2, ANI=V3, match=V4, total=V5) %>% #rename variable
  mutate(isolate = str_extract(isolate, "[:digit:]+_[:digit:]+_[:alnum:]+"),
         MPA = str_extract(MPA, "Rht_[:digit:]+_(N|C)"))  %>% #rename isolates
  group_by(isolate) %>%
  slice_max(tibble(ANI), n = 1) #select the top 1 ANI values for each derived isolate (MPA)

#prepare a csv file as the following: Derived_Strains,Most_Probable_Ancestor
write.table(MPA[,1:2], file="./MPA_data/derived_mpa.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

Index 56 genomes of ancestral strains for bwa (Create files ending with .amb, .ann, .bwt, .pac, .sa).
```bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/Rht*fasta; do
/usr/bin/bwa index "$i" 
done
```

**Run bwa (bwa mem -t 5 -M -R), convert SAM to BAM.**
```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
while IFS=',' read -r a b; do # a=derived_strain,b=most_probable_ancestor
r1=/home/xingyuan/rhizo_ee/snp_indel/trimmed_paired_reads/"$a"_*R1*
r2=/home/xingyuan/rhizo_ee/snp_indel/trimmed_paired_reads/"$a"_*R2*
ref=/home/xingyuan/rhizo_ee/derived+original_genomes/"$b".fasta

/usr/bin/bwa mem -t 5 -M -R "@RG\tID:"$a"-"$b"\tSM:"$a $ref $r1 $r2 | /usr/local/bin/samtools view -huS -o "$a"-"$b".bam - 

done < ../derived_mpa.csv
```

## 2. Run picard to manipulate the SAM files
Picard Version: 3.0.0 <br>
Work done one info2020

**Create sequence dictionary file (.dict) for the 56 original strains (Required for ReorderSam).** <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/Rht*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample_name=${j%.fasta}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar CreateSequenceDictionary -R $i -O /home/xingyuan/rhizo_ee/derived+original_genomes/"$sample_name".dict
done
```

**Reorder reads in the BAM file to match the contig ordering in reference file.** <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037426651-ReorderSam-Picard-

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/snp_indel/bwa_output/*bam; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/bwa_output/}
base_name=${j%.bam}
mpa_name=${base_name#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar ReorderSam -R /home/xingyuan/rhizo_ee/derived+original_genomes/"$mpa_name".fasta -I $i -O /home/xingyuan/rhizo_ee/snp_indel/reorderSAM_output/"$base_name".reordered.bam -SD /home/xingyuan/rhizo_ee/derived+original_genomes/"$mpa_name".dict
done
```

**Assign all the reads in a file to a single new read-group. Read group fields are required by GATK.** <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard- 


Obtain the flowcell ID and lane number from FASTQ for all reads. The format of the first line of FASTQ is as following (https://help.basespace.illumina.com/files-used-by-basespace/fastq-files): 
```
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
```

```bash
(base) [xingyuan@info2020 raw_reads]$ zcat *_R1_*fastq.gz | grep "^@" | cut -d ":" -f 3,4 | uniq
HTH5JDRX2:2
(base) [xingyuan@info2020 raw_reads]$ zcat *_R2_*fastq.gz | grep "^@" | cut -d ":" -f 3,4 | uniq
HTH5JDRX2:2
```
This indicates the flowcell ID is HTH5JDRX2 and the lane is 2 for all reads. Read more about read groups on https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups.

Add the following read group fields: 
* ``ID``: read group identifier, use combination of names of the derived isolate and its most probable ancestor (ID must be unique among all read groups as described in SAM Format Specification document)
* ``PU``: platform unit, use HTH5JDRX2.2 (flowcell ID and lane)
* ``SM``: sample, use the name of the derived isolate
* ``PL``: platform/technology used to produce the reads, use ILLUMINA
* ``LB``: library, use the name of the derived isolate

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/snp_indel/reorderSAM_output/*.reordered.bam; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/reorderSAM_output/}
base_name=${j%.reordered.bam}
derived_isolate_name=${base_name%-*}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar AddOrReplaceReadGroups -I $i -O /home/xingyuan/rhizo_ee/snp_indel/AddOrReplaceReadGroups_output/"$base_name".new_rg.bam -ID "$base_name" -LB "$derived_isolate_name" -PL ILLUMINA -PU HTH5JDRX2.2 -SM "$derived_isolate_name"
done
```

Verify read group is added correctly.
```bash
for i in *new_rg.bam; do 
samtools view -H $i | grep "^@RG"
done
```

**Sort the input BAM file by coordinate.** <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/snp_indel/AddOrReplaceReadGroups_output/*.new_rg.bam; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/AddOrReplaceReadGroups_output/}
base_name=${j%.new_rg.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar SortSam -I $i -O /home/xingyuan/rhizo_ee/snp_indel/sortSAM_output/"$base_name".coordinate_sorted.bam -SO coordinate
done
```

**Identify duplicate reads and index the BAM files. BuildBamIndex requires BAM files to be sorted by coordinate.** <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard- <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-

Duplicate reads are the reads from the same fragment of DNA. Optical duplicates occur when a single amplification cluster is detected incorrectly by the optical sensor to be multiple clusters. Duplicate types (DT): library/PCR-generated duplicates (LB), sequencing-platform artifact duplicates (SQ). When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates.

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/snp_indel/sortSAM_output/*.coordinate_sorted.bam; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/sortSAM_output/}
base_name=${j%.coordinate_sorted.bam}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar MarkDuplicates -I $i -O /home/xingyuan/rhizo_ee/snp_indel/markduplicate_and_index_output/"$base_name".marked_duplicates.bam -M /home/xingyuan/rhizo_ee/snp_indel/markduplicate_and_index_output/"$base_name".marked_dup_metrics.txt && /scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/picard.jar BuildBamIndex -I /home/xingyuan/rhizo_ee/snp_indel/markduplicate_and_index_output/"$base_name".marked_duplicates.bam
done
```

## 3A. Call SNPs and indels using GATK HaplotypeCaller
https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller

Samtools Version: 1.13 (using htslib 1.13) <br>
GATK Version: 4.4.0.0 <br>
Work done on info2020

**Create index file (.fai) for the 56 ancestral isolates for GATK.** <br>
http://www.htslib.org/doc/samtools-faidx.html

```bash
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/Rht*fasta; do
/usr/local/bin/samtools faidx $i
done
```

**Run HaplotypeCaller.**
```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/snp_indel/markduplicate_and_index_output/*.marked_duplicates.bam; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/markduplicate_and_index_output/}
base_name=${j%.marked_duplicates.bam}
mpa_name=${base_name#*-}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar HaplotypeCaller -R /home/xingyuan/rhizo_ee/derived+original_genomes/"$mpa_name".fasta -I "$i" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O /home/xingyuan/rhizo_ee/snp_indel/haplotypecaller_output/"$base_name".g.vcf.gz -ERC GVCF
done
```

**Merge SNPs and indels results and perform joint genotyping.** <br>
CombineGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832710975771-CombineGVCFs <br>
GenotypeGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs

Create a list with MPA name on each row.
```bash
cat derived_mpa.csv | cut -d "," -f 2 | sort | uniq > MPA_list
```

Run CombineGVCFs and GenotypeGVCFs (input file name for ``--variant`` must end in ".list").
```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
while read Rht; do

# Group derived strains with the same most probable ancestors to a list. This should create 26 lists.
find /home/xingyuan/rhizo_ee/snp_indel/haplotypecaller_output/*"$Rht"*.vcf.gz > /home/xingyuan/rhizo_ee/snp_indel/genotypegvcfs_output/MPA-"$Rht".list &&

# Run CombineGVCFs to combine the vcf.gz files in each list to one vcf.gz file. This should create 26 cohort.g.vcf.gz files.
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar CombineGVCFs -R /home/xingyuan/rhizo_ee/derived+original_genomes/"$Rht".fasta --variant MPA-"$Rht".list -O /home/xingyuan/rhizo_ee/snp_indel/genotypegvcfs_output/"$Rht".cohort.g.vcf.gz &&

# Run GenotypeGVCFs on each 26 cohort.g.vcf files
/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar GenotypeGVCFs -R /home/xingyuan/rhizo_ee/derived+original_genomes/"$Rht".fasta -V /home/xingyuan/rhizo_ee/snp_indel/genotypegvcfs_output/"$Rht".cohort.g.vcf.gz -ploidy 1 -O /home/xingyuan/rhizo_ee/snp_indel/genotypegvcfs_output/genotype_"$Rht".vcf.gz -stand-call-conf 30

done < ../MPA_list
```

**Filter SNPs.**
https://vcftools.github.io/man_latest.html 

Vcftools Version: 0.1.16 <br>
Work done on info2020

Only keep sites that meet all of the following thresholds: 
* bi-allelic sites
* sites with mean depth values (over all included individuals) greater than or equal to 20 and less than or equal to 230
* sites with Quality value above 30
* sites only allowed to have less than or equal to 10% missing data (set ``--max-missing`` to 0.9 since it is defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed)

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
for i in /home/xingyuan/rhizo_ee/snp_indel/genotypegvcfs_output/genotype*gz; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/genotypegvcfs_output/}
out=${j%.vcf.gz}

/2/scratch/batstonelab/bin/vcftools-0.1.16/bin/vcftools --gzvcf "$i" --min-meanDP 20 --max-meanDP 230 --minQ 30 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$out".filt

done
```

**Extract fields from a VCF file to a tab-delimited table.** <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable

See descriptions of each field on https://samtools.github.io/hts-specs/VCFv4.2.pdf (match the documentation version to your VCF file format version. VCF file format version here is VCFv4.2). 

```bash
for i in /home/xingyuan/rhizo_ee/snp_indel/vcftools_output/genotype_Rht_*.recode.vcf; do
j=${i#/home/xingyuan/rhizo_ee/snp_indel/vcftools_output/}
base_name=${j%.filt.recode.vcf}
mpa_name=${base_name#genotype_}

/scratch/batstonelab/bin/apps/jdk-21.0.2/bin/java -jar /scratch/batstonelab/bin/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar VariantsToTable -V "$i" -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F DP -GF GT -GF AD -GF DP -GF GQ -GF PL -O /home/xingyuan/rhizo_ee/snp_indel/variantstotable/"$base_name".table
done
```

**Download all 26 tab-delimited table to local computer to be analyzed in R.**

## 3B. Call SNPs and indels using snippy - To be continued
https://github.com/tseemann/snippy

snippy version: 4.6.0 <br>
Work done on info2020

```bash
#!/bin/bash
#Usage: nohup ./ThisScript &
while IFS=',' read -r a b; do # a=derived_strain,b=most_probable_ancestor
r1=/home/xingyuan/rhizo_ee/snp_indel/trimmed_paired_reads/"$a"_*R1*
r2=/home/xingyuan/rhizo_ee/snp_indel/trimmed_paired_reads/"$a"_*R2*
ref=/home/xingyuan/rhizo_ee/derived+original_genomes/"$b".fasta

/home/xingyuan/tools/snippy/bin/snippy --cpus 10 --ram 8 --outdir "$a"-"$b" --reference $ref --R1 $r1 --R2 $r2 --rgid "$a"-"$b" --mapqual 60 --basequal 30
done < ../derived_mpa.csv
```

Merge single-sample VCF file to multi-samples VCF files
```bash
/usr/local/bin/bcftools merge
```

## 5. Find genes at and near the positions of SNPs - To be continued after annotations are done
https://github.com/bedops/bedops

### Covert gff and vcf to bed format
Bedops Version: 2.4.41 <br>
Work done on info2020

*Note: Add the bin containing the gff2bed and vcf2bed to .bashrc file. Without this step, this error may appear "convert2bed: command not found".*

**Convert gff of original strains to bed format**
```bash
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






