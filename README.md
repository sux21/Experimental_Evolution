# 2020 Experimental Evolution 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on **info server** (contact Brian Golding at golding@mcmaster.ca for an info account). **Compute canada server** will be used if the info server cannot run the program (register for an compute canada account at https://ccdb.alliancecan.ca/security/login). Results produced by compute canada server will be transferred to info server. 

#
Try this to convert pgap's gff files to gtf files (https://agat.readthedocs.io/en/latest/gxf.html)
  
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
```bash
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads *gz &
```

### 2. Run MultiQC for raw data
https://multiqc.info/

Version: MultiQC v1.9 <br>
Work done on info114

```bash
multiqc . 
```

### 3. Run Trimmomatic 
http://www.usadellab.org/cms/?page=trimmomatic

Version: 0.39 <br>
Work done on info114

**All 363 samples (726 files)**
```bash
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
```bash
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_trimmomatic_reads *_P_* &
```

### 5. Run MultiQC for trimmed reads
Version: MultiQC v1.9 <br>
Work done on info114

**All 363 samples (726 files)**
```bash
multiqc . 
```

## During Assembly 
### 1. Run SPAdes 
https://github.com/ablab/spades

Version: 3.15.2 <br>
Work done on info114

**52 samples**
```bash
#!/bin/bash 
for R1 in 10_1_8_*R1_P_* 10_1_9_*R1_P_* 10_7_6_*R1_P_* 11_4_2_*R1_P_* 11_4_4_*R1_P_* 11_5_6_*R1_P_* 13_4_1_*R1_P_* 14_4_6_*R1_P_* 14_5_3_*R1_P_* 15_4_4_*R1_P_* 15_4_6_*R1_P_* 16_1_6_*R1_P_* 16_1_7_*R1_P_* 16_1_8_*R1_P_* 16_4_2_*R1_P_* 16_4_3_*R1_P_* 16_6_6_*R1_P_* 17_2_1_*R1_P_* 17_2_8_*R1_P_* 17_2_9_*R1_P_* 18_1_4_*R1_P_* 18_1_5_*R1_P_* 19_1_1_*R1_P_* 19_5_8_*R1_P_* 2_2_5_*R1_P_* 2_5_2_*R1_P_* 2_6_4_*R1_P_* 3_1_5_*R1_P_* 3_2_1_*R1_P_* 3_2_3_*R1_P_* 3_2_6_*R1_P_* 3_2_7_*R1_P_* 3_3_5_*R1_P_* 3_3_7_*R1_P_* 3_3_9_*R1_P_* 4_1_2_*R1_P_* 4_1_4_*R1_P_* 4_2_1_*R1_P_* 6_4_5_*R1_P_* 6_4_7_*R1_P_* 6_7_5_*R1_P_* 7_1_2_*R1_P_* 7_1_5_*R1_P_* 7_6_3_*R1_P_* 7_6_9_*R1_P_* 7_7_2_*R1_P_* 7_7_3_*R1_P_* 8_4_10_*R1_P_* 8_4_4_*R1_P_* 9_3_7_*R1_P_* 9_7_6_*R1_P_* 9_7_9_*R1_P_* 
do
R2=${R1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}

spades.py --careful -1 $R1 -2 $R2 -o /home/xingyuan/rhizo_ee/spades_assembly/${R1%_*_L002_*gz}
done
```

**Remaining 311 samples**
```bash
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

```bash
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

```bash
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

```bash
#!/bin/bash 
for i in Rht*
do

/home/xingyuan/tools/quast-5.2.0/quast.py $i -m 0 -t 5 -o /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/quast_2008_long_reads/${i%.fasta}

done
```

# Step 2 - Generating Data
## Analysis 1 - Find the most probable ancestor for each experimentally evolved strain 

**Samples: 363 experimentally evolved strains + 56 original strains**

### (Optional) 1. Run Spine: find core genomes 
https://github.com/Alan-Collins/Spine-Nucmer-SNPs

Version: 0.3.2 <br>
Work done on info114. 

```bash
#!/bin/bash
# Copy contigs.fasta from ``rhizo_ee/spades_assembly`` to ``rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY``:

for i in *; do

if [[ $i =~ ".sh" ]] || [[ $i =~ "fasta" ]] || [[ $i =~ "quast" ]]; then
   continue
fi

cp /home/xingyuan/rhizo_ee/spades_assembly/$i/contigs.fasta /home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/"$i-contigs.fasta"
done
```

```bash
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
```bash
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -o fastani.contigs.out &
```
**Query to reference (--fragLen 4000)**
```bash
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 4000 -o fastani.contigs4000.out &
```
**Query to reference (--fragLen 3500)**
```bash
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 3500 -o fastani.contigs3500.out &
```
**Query to reference (--fragLen 2500)**
```bash
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 2500 -o fastani.contigs2500.out &
```
**Query to reference (--fragLen 2000)**
```bash
ls Rht* > reference_list
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl reference_list -t 5 --fragLen 2000 -o fastani.contigs2000.out &
```
**Query (core genome) to reference (core genome)**
```bash
ls *contigs.fasta.core.fasta > core_query_list
ls *_?.fasta.core.fasta > core_reference_list

nohup /usr/local/bin/fastANI --ql core_query_list --rl core_reference_list -t 5 -k 5 --fragLen 10 -o fastani.que_core_to_ref_core.out 
```

**Reference to reference**
```bash
ls Rht* > reference_list

nohup /usr/local/bin/fastANI --ql reference_list --rl reference_list -o fastani.ref_to_ref.out &
```

**Query to query**
```bash
ls *contigs.fasta > contigs_query_list

nohup /usr/local/bin/fastANI --ql contigs_query_list --rl contigs_query_list -o fastani.que_to_que.out &
```

**Query+reference to query+reference (exclude as5_2_4)**
```bash
ls *.fasta > all_samples_no_as5_2_4 (remove as_5_2_4 from the list)

nohup /usr/local/bin/fastANI --ql all_samples_no_as5_2_4 --rl all_samples_no_as5_2_4 -t 5 -o fastani.all_to_all.out &
```

**Get mapping diagram for these comparisons: 20_6_10 vs 511_N, 20_6_10 vs 415_C; 20_2_6 vs 596_N, 20_2_6 vs 116_N; 16_1_6 vs 173_C, 16_1_6 vs 016_N; 16_3_4 vs 003_C, 16_3_4 vs 706_C**

Rscript Version: 3.6.3 (2020-02-29) <br>
Work done on info2020

**Open R and install genoplotr using ``install.packages("genoPlotR")`` before running the following lines**

```bash
/usr/local/bin/fastANI -q 20_6_10-contigs.fasta -r Rht_511_N.fasta --visualize -o 20_6_10-Rht_511_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_6_10-contigs.fasta Rht_511_N.fasta 20_6_10-Rht_511_N.fastani.out.visual 
```
```bash
/usr/local/bin/fastANI -q 20_6_10-contigs.fasta -r Rht_415_C.fasta --visualize -o 20_6_10-Rht_415_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_6_10-contigs.fasta Rht_415_C.fasta 20_6_10-Rht_415_C.fastani.out.visual 
```
```bash
/usr/local/bin/fastANI -q 20_2_6-contigs.fasta -r Rht_596_N.fasta --visualize -o 20_2_6-Rht_596_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_2_6-contigs.fasta Rht_596_N.fasta 20_2_6-Rht_596_N.fastani.out.visual 
```
```bash
/usr/local/bin/fastANI -q 20_2_6-contigs.fasta -r Rht_116_N.fasta --visualize -o 20_2_6-Rht_116_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 20_2_6-contigs.fasta Rht_116_N.fasta 20_2_6-Rht_116_N.fastani.out.visual 
```
```
/usr/local/bin/fastANI -q 16_1_6-contigs.fasta -r Rht_173_C.fasta --visualize -o 16_1_6-Rht_173_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_1_6-contigs.fasta Rht_173_C.fasta 16_1_6-Rht_173_C.fastani.out.visual 
```
```bash
/usr/local/bin/fastANI -q 16_1_6-contigs.fasta -r Rht_061_N.fasta --visualize -o 16_1_6-Rht_061_N.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_1_6-contigs.fasta Rht_061_N.fasta 16_1_6-Rht_061_N.fastani.out.visual 
```
```bash
/usr/local/bin/fastANI -q 16_3_4-contigs.fasta -r Rht_003_C.fasta --visualize -o 16_3_4-Rht_003_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_3_4-contigs.fasta Rht_003_C.fasta 16_3_4-Rht_003_C.fastani.out.visual 
```
```bash
/usr/local/bin/fastANI -q 16_3_4-contigs.fasta -r Rht_706_C.fasta --visualize -o 16_3_4-Rht_706_C.fastani.out && ~/tools/R/bin/Rscript visualize.R 16_3_4-contigs.fasta Rht_706_C.fasta 16_3_4-Rht_706_C.fastani.out.visual 
```

## Analysis 2 - Phylogeny for the original strains 

**Samples: 56 original strains**

### 1. Run Spine-Nucmer-SNPs
https://github.com/Alan-Collins/Spine-Nucmer-SNPs

Version: 0.3.2 <br>
Work done on info113

**(1)**
```bash
ls Rht* | awk 'BEGIN { FS="\t"; OFS="\t" } { print "/home/xingyuan/rhizo_ee/2008_2020_strains_comparison_All/ASSEMBLY/"$1, $1, "fasta" }' > ../phylogeny_2008/config.txt
```
```bash
nohup spine.pl -f config.txt -t 5 &
```

**(2)**
```bash
ls ../spine/*.core.fasta | while read i; do acc=${i%.core*}; acc=${acc#../spine/output.}; nucmer --prefix=${acc}_core ../spine/output.backbone.fasta $i; delta-filter -r -q ${acc}_core.delta > ${acc}_core.filter; show-snps -Clr ${acc}_core.filter > ${acc}_core.snps; done
```

**(3)**
```bash
python3 ~/tools/Spine-0.3.2/snps2fasta.py -r ../spine/output.backbone.fasta -f variant_core.fasta -whole -m snp_matrix.csv -d '\t' -p '(.*)_core\.snps' ../nucmer/*.snps
```

### 2. IQ-Tree
http://www.iqtree.org/doc/Tutorial

Version: 2.2.0 <br>
Work done on info114

```bash
# For C_only population (28 samples)
nohup iqtree2 -T 5 -s C_only.fasta -bb 1000 -wbt --seqtype DNA &

# For N_only population (28 samples)
nohup iqtree2 -T 5 -s N_only.fasta -bb 1000 -wbt --seqtype DNA &

# For mixed population (28 samples)
iqtree2 -T 5 -s mix.fasta -bb 1000 -wbt --seqtype DNA
```

## Analysis 3: Gene Presence Absence
### 1. Genome annotation - Prokka
#### Prokka 
https://github.com/tseemann/prokka

Version: 1.12-beta <br>
Work done on info114

```
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/spades_assembly/*/scaffolds.fasta; do
  j=${i#/home/xingyuan/rhizo_ee/spades_assembly/}
  sample=${j%/scaffolds.fasta}

  /usr/local/prokka/bin/prokka "$i" --outdir "$sample" --prefix "$sample" --locustag "$sample" --cpus 6
done

for i in /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/Rht*_?.fasta; do
  j=${i#/home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/}
  sample=${j%.fasta}

  /usr/local/prokka/bin/prokka "$i" --outdir "$sample" --prefix "$sample" --locustag "$sample" --cpus 6
done
```

#### Gene presence absence - Roary
https://sanger-pathogens.github.io/Roary/

Version: 1.007001 <br>
Work done on info114

**Exclude this sample: as5_2_4 (maybe Paenibacillus)**

```bash
nohup roary -p 6 -y */*gff &
```

### 2. Genome annotation - PGAP
#### Transfer scaffolds file from info to graham
```
# Re-name each scaffolds.fasta file with its sample name
for i in */scaffolds.fasta; do
sample_name=${i%/scaffolds.fasta}
cp $i "$sample_name"-scaffolds.fasta
done

# Put all renamed files into a new directory and compress the directory (Also create a checksum for the compressed directory to check file integrity before and after transfer)
mkdir rhizo_ee.scaffolds
mv *scaffolds.fasta rhizo_ee.scaffolds
tar -zcvf rhizo_ee.scaffolds.tar.gz /home/xingyuan/rhizo_ee/spades_assembly/rhizo_ee.scaffolds
md5sum rhizo_ee.scaffolds.tar.gz > md5sums.txt

# Transfer rhizo_ee.scaffolds.tar.gz and md5sums.txt from info computers to graham computers
(base) [xingyuan@infoserv spades_assembly]$ scp rhizo_ee.scaffolds.tar.gz sux21@graham.computecanada.ca:/home/sux21/2023_summer_coop/rhizo_ee/genomes
(base) [xingyuan@infoserv spades_assembly]$ scp md5sums.txt sux21@graham.computecanada.ca:/home/sux21/2023_summer_coop/rhizo_ee/genomes

# Verify the checksum
[sux21@gra-login1 genomes]$ md5sum -c md5sums.txt

# Extract rhizo_ee.scaffolds.tar.gz
[sux21@gra-login1 genomes]$ tar -xf rhizo_ee.scaffolds.tar.gz
```

#### Filter sequences shorter than 200 bp (pgap only takes sequences equal or longer than 200 bp)
https://github.com/shenwei356/seqkit

Seqkit Version: 2.3.0 <br>
Work done on graham cluster

```bash
#!/bin/bash
#SBATCH --time=00-00:10
#SBATCH --mem-per-cpu=1024M
#SBATCH --account=def-batstone

module load seqkit/2.5.1

for i in *scaffolds.fasta; do
sample=${i%-scaffolds.fasta}

seqkit seq -m 200 "$sample"-scaffolds.fasta > "$sample"-scaffolds.filter.fasta
done
```

#### Download pgap.py
https://github.com/ncbi/pgap/tree/1126_Test

Version: 2023-05-17.build6771 <br>
Work done on graham cluster 

```bash
module load apptainer 
wget -O pgap.py https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
chmod +x pgap.py
./pgap.py -D apptainer 
```

#### Verify species taxonomy 
https://github.com/ncbi/pgap/wiki/Taxonomy-Check

Version: 2023-05-17.build6771 <br>
Work done on cedar cluster

**Run one script at a time. When it finishes, then run the next one. This is because pgap creates yaml file when it is running. When two pgap jobs run at the same time, one of the jobs may use the yaml produced from the other job. Creating  an error like the following:**
```bash
Original command: /project/6078724/sux21/tools/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap_2023-05-17.build6771.sif -r -o 10_5_8 -g /project/6078724/sux21/rhizo_ee/genomes/10_5_8-contigs.filter.fasta -s Rhizobium leguminosarum -c 40 --taxcheck-only

Docker command: /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/apptainer/1.1.8/bin/apptainer exec --bind /home/sux21/.pgap/input-2023-05-17.build6771:/pgap/input:ro --bind /project/6078724/sux21/rhizo_ee/taxcheck:/pgap/user_input --bind /project/6078724/sux21/rhizo_ee/taxcheck/pgap_input_39ew9cbl.yaml:/pgap/user_input/pgap_input.yaml:ro --bind /tmp:/tmp:rw --bind /project/6078724/sux21/rhizo_ee/taxcheck/10_5_8:/pgap/output:rw --pwd /pgap /project/6078724/sux21/tools/pgap_2023-05-17.build6771.sif /bin/taskset -c 0-39 cwltool --timestamps --debug --disable-color --preserve-entire-environment --outdir /pgap/output pgap/taxcheck.cwl /pgap/user_input/pgap_input.yaml

--- Start YAML Input --- 
fasta:
    class: File
    location: /project/6078724/sux21/rhizo_ee/genomes/Rht_061_N.fasta
submol:
    class: File
    location: pgap_submol_bdp39o2k.yaml
supplemental_data: { class: Directory, location: /pgap/input }
report_usage: true
--- End YAML Input ---
```

**Script 1: Evolved strains**
```bash
#!/bin/bash
#SBATCH --time=00-20:00
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

**Script 2: Original Strains**
```bash
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
**Samples with INCONCLUSIVE status: 16_1_6, 4_4_10, as5_2_4, Rht_061_N, Rht_173_C, Rht_209_N, Rht_231_N, Rht_717_N, Rht_773_N**

#### Run pgap (Prepare 5 scripts and run each script in a different directory)

Version: 2023-05-17.build6771 <br>
Work done on graham cluster 

**script 1 (94 samples)**
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

for i in /project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/{1..5}_*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/}
sample=${j%-scaffolds.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 2 (70 samples)**
```bash
#!/bin/bash
#SBATCH --time=03-12:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/{6..9}_*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/}
sample=${j%-scaffolds.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 3 (109 samples)**
```bash
#!/bin/bash
#SBATCH --time=05-00:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/1{0..5}_*filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/}
sample=${j%-scaffolds.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 4 (90 samples)**
```bash
#!/bin/bash
#SBATCH --time=04-12:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/1{6..9}_*filter.fasta /project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/20_*filter.fasta /project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/as5_2_4-scaffolds.filter.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/rhizo_ee.scaffolds/}
sample=${j%-scaffolds.filter.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-05-17.build6771.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40
done
```
**script 5 (56 samples)**
```bash
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

 Find directories that don't have annot.gbk: ``find . -type d \! -exec test -e '{}/annot.gbk' \; -print`` <br>
**Samples that do not have gbk produced: 14_5_5, 11_4_4, 7_4_2, 2_3_4, 15_5_1, 6_3_2, 4_4_10, 2_5_9, 15_3_1, 1_1_3, 18_6_7, 2_6_3, 15_2_1, 2_4_11, 5_3_6, 3_1_3, 8_4_7, 4_1_4, 19_4_7, 19_1_9, 5_3_9, 3_3_5, 17_2_7, Rht_773_N (24 samples, 419-24=395 samples remain)**

#### Convert PGAP's gff files for roary

**Rename each annot_with_genomic_fasta.gff with sample names**
```bash
for i in */annot_with_genomic_fasta.gff; do
sample=${i%/annot_with_genomic_fasta.gff}

cp $i "$sample"_annot_with_genomic_fasta.gff
done
```

**Reformat PGAP's annot_with_genomic_fasta.gff for AGAT**
```bash
#!/bin/bash
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

**Standardise gff file using AGAT**

AGAT Version: 1.2.1 <br>
Work done on info19

```bash
for i in *fixed; do
/2/scratch/batstonelab/bin/AGAT/bin/agat_convert_sp_gxf2gxf.pl -g $i -o ${i%_annot_with_genomic_fasta.gff.fixed}.pgap.gff
done
```

#### Run roary
https://sanger-pathogens.github.io/Roary/

Roary version: 1.007001 <br>
Work done on info114

```
roary -p 6 *gff
```

### 2. Glimmer-Interproscan
#### 1. Glimmer (Do not assume sequence is circular for evolved strains, assume circular for original strains)
http://ccb.jhu.edu/software/glimmer/index.shtml

Version: glimmer3.02 <br>
Work done on info2020

**Change directory paths in ``g3-iterated.csh``**
```bash
Before changes:
set awkpath = /fs/szgenefinding/Glimmer3/scripts
set glimmerpath = /fs/szgenefinding/Glimmer3/bin
set elphbin = /nfshomes/adelcher/bin/elph

After changes:
set awkpath = /home/xingyuan/tools/glimmer3.02/scripts
set glimmerpath = /home/xingyuan/tools/glimmer3.02/bin
set elphbin = /home/xingyuan/tools/ELPH/bin/Linux-i386/elph
```
**Evolved strains (Don't assume circular sequences because they are not completed genomes):** 

**Copy ``g3-iterated.csh`` and create ``g3-iterated_evolved.csh``.**
```bash
cp g3-iterated.csh g3-iterated_evolved.csh
```

**In ``g3-iterated_evolved.csh``, Add linear options to glimmer and long-orfs**
```bash
Before changes:
# add/change glimmer options here
set glimmeropts = "-o50 -g110 -t30"
...
step1:
# Find long, non-overlapping orfs to use as a training set
echo "Step 1 of ${numsteps}:  Finding long orfs for training"
$glimmerpath/long-orfs -n -t 1.15 $genome $tag.longorfs

After changes:
# add/change glimmer options here
set glimmeropts = "--linear -o50 -g110 -t30"
...
step1:
# Find long, non-overlapping orfs to use as a training set
echo "Step 1 of ${numsteps}:  Finding long orfs for training"
$glimmerpath/long-orfs --linear -n -t 1.15 $genome $tag.longorfs
```
**Run glimmer for evolved strains**
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/*-contigs.fasta; do
j=${i#/home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/}
sample=${j%-contigs.fasta}

if [[ ! -d "$sample" ]]; then 
   mkdir "$sample"
fi

/home/xingyuan/tools/glimmer3.02/scripts/g3-iterated_evolved.csh "$i" "$sample"-contigs; mv `ls "$sample"-contigs*` "$sample"/
done
```
**Run glimmer for original strains (Assume circular sequences because they are completed genomes)**
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/Rht*_?.fasta; do
j=${i#/home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/}
sample=${j%.fasta}

if [[ ! -d "$sample" ]]; then 
   mkdir "$sample"
fi

/home/xingyuan/tools/glimmer3.02/scripts/g3-iterated.csh "$i" "$sample"; mv `ls "$sample"*` "$sample"/
done
```
#### Re-format ``.predict`` file
```bash
#!/bin/bash
# Format glimmer .predict file for the multi-extract program  
# Usage: ./ThisScript .predict

for (( c=1; c<=`grep -c ">" $1`; c++ )); do

i=0
let y=$c+1
while read line; do
  if [[ $line =~ ">" ]]; then
    let i=$i+1
  fi
  if [ $i -eq $c ] && [[ $line =~ ">" ]]; then
    tag=${line#>}
  fi
  if [ $i -eq $c ] && [[ $line =~ "orf" ]]; then
    id=${line:0:8}
    string=${line:8} 
    echo "$tag"_"$id" "$tag" "$string"
  fi
  if [ $i -eq $y ]; then
    break 
  fi
done < $1

done
```
```
#!/bin/bash
for i in */*-contigs.predict */Rht*_?.predict; do
j=${i#*/}
sample=${j%.predict}

./re-format.sh "$i" > "$sample".re_formatted.predict
done
```
#### Run multi-extract to extract genes from ``.predict`` file
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/genes_presence_absence/glimmer/*.re_formatted.predict; do
j=${i#/home/xingyuan/rhizo_ee/genes_presence_absence/glimmer/}
sample=${j%.re_formatted.predict}

/home/xingyuan/tools/glimmer3.02/bin/multi-extract /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$sample".fasta "$i" > "$sample".glimmer_genes.fasta 
done
```

#### Run interproscan
Interproscan version: 5.63-95.0 <br>
Work done on cedar cluster

**Install local lookup service using instructions at https://interproscan-docs.readthedocs.io/en/latest/LocalLookupService.html (The lookup service is about 1.4T).**



### 3. Don't use clustering softwares
#### Find genes lost in evolved strains 
Samtools Version: 1.11 (using htslib 1.11) <br>
Bedtools Version: 2.19.1 <br>
Work done info114

**Extract mapped reads**
```bash
#!/bin/bash
for i in *bam; do
sample=${i/.bam/.mapped.bam}

samtools view -F 260 -o "$sample" "$i"
done
```

**Find genomic regions where the most probable ancestor has but its evolved strains do not**
```bash
bedtools intersect -a /home/xingyuan/rhizo_ee/genes_presence_absence/pgap/Rht_056_N/annot_with_genomic_fasta.bed -b 14_2_2-Rht_056_N.bam -header -v > 14_2_2-Rht_056_N.genes_lost
```

#### Find genes gained in evolved strains 
Samtools Version: 1.11 (using htslib 1.11) <br>
Seqkit Version: 2.5.1 <br>
Blast Version: 2.13.0, build Sep 13 2022 22:19:14 <br>
Work done info114

**Extract unmapped reads**
```bash
#!/bin/bash
for i in *_?.bam; do
r1=${i/.bam/.unmapped.r1.fastq}
r2=${i/.bam/.unmapped.r2.fastq}

samtools fastq -f 4 -1 "$r1" -2 "$r2" "$i"
done
```

**Rename fasta header for the 56 original strains with sample name, and combine all genomes into one file**
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/Rht*fasta; do
j=${i#/home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/}
sample=${j/.fasta/_}
out=${j/.fasta/.rename.fasta}

/home/xingyuan/tools/seqkit replace -p ^ -r "$sample" "$i" > "$out"
done
```
```bash
cat *.rename.fasta > 56_original.combine.fasta
```

## Analysis 4: Call SNPS between each evolved strain and its most probable ancestor
https://github.com/rtbatstone/how-rhizobia-evolve/blob/master/Variant%20discovery/Variant_calling.md

**Install the latest version of picard at https://github.com/broadinstitute/picard** 
```bash
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

```bash
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
```bash
#!/bin/bash
for i in Rht*fasta; do
bwa index $i
done
```

#### Run bwa (bwa mem -t 5 -M -R), convert .sam to .bam (samtools view -huS), and produce alignment statistics (samtools stats)
Samtools stats: http://www.htslib.org/doc/samtools-stats.html
```bash
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
```bash
#!/bin/bash
for i in Rht*fasta; do
ref=${i%.fasta}
java -jar /home/xingyuan/tools/picard.jar CreateSequenceDictionary -R "$ref".fasta -O "$ref".dict
done
```

##### Run ReorderSam
```bash
#!/bin/bash
for i in *bam; do
sample=${i%.bam}
ref=${sample#*-}

java -jar /home/xingyuan/tools/picard.jar ReorderSam -R /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".fasta -I "$sample".bam -O "$sample".reordered.bam -SD /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$ref".dict
done
```

#### (2) Assign all the reads in a file to a single new read-group
https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
```bash
#!/bin/bash
for i in *.reordered.bam; do
sample=${i%.reordered.bam}

java -jar /home/xingyuan/tools/picard.jar AddOrReplaceReadGroups -I "$sample".reordered.bam -O "$sample".new_rg.bam -ID "$sample" -LB rhizo_ee -PL Illumina -PU 1 -SM "$sample" 
done
```

#### (3) Sort the input BAM file by coordinate
https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-
```bash
#!/bin/bash
for i in *.new_rg.bam; do
sample=${i%.new_rg.bam}

java -jar /home/xingyuan/tools/picard.jar SortSam -I "$sample".new_rg.bam -O "$sample".coordinate_sorted.bam -SO coordinate
done
```

#### (4) Identify duplicate reads and index the BAM files
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard- <br>
https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-
```bash
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
```bash
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip gatk-4.4.0.0.zip
```

#### Create index file (``.fai``) of reference for GATK
http://www.htslib.org/doc/samtools-faidx.html
```bash
#!/bin/bash
for i in Rht*fasta; do
samtools faidx $i
done
```

#### HaplotypeCaller
https://gatk.broadinstitute.org/hc/en-us/articles/13832687299739-HaplotypeCaller
```bash
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
**Run CombineGVCFs and GenotypeGVCFs**
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

### 5. Filter SNPs
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

### 6. VariantsToTable
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

### 7. Find genes at the positions of SNPs
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

### 8. Find alternate alleles that occur in more than 1 evolved strains 

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
