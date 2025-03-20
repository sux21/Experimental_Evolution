# Run Trimmomatic 
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

# Run Spine: find core genomes 
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
# Analysis - Phylogeny for the original strains 

**Samples: 56 original strains**

## Run Spine-Nucmer-SNPs
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

## IQ-Tree
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

# Genome annotation - Prokka
## Prokka 
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

## Gene presence absence - Roary
https://sanger-pathogens.github.io/Roary/

Version: 1.007001 <br>
Work done on info114

**Exclude this sample: as5_2_4 (maybe Paenibacillus)**

```bash
nohup roary -p 6 -y */*gff &
```

# Glimmer-Interproscan
## Glimmer (Do not assume sequence is circular for evolved strains, assume circular for original strains)
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
## Run multi-extract to extract genes from ``.predict`` file
```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/genes_presence_absence/glimmer/*.re_formatted.predict; do
j=${i#/home/xingyuan/rhizo_ee/genes_presence_absence/glimmer/}
sample=${j%.re_formatted.predict}

/home/xingyuan/tools/glimmer3.02/bin/multi-extract /home/xingyuan/rhizo_ee/find_most_probable_ancestors_all/ASSEMBLY/"$sample".fasta "$i" > "$sample".glimmer_genes.fasta 
done
```

## Run interproscan
Interproscan version: 5.63-95.0 <br>
Work done on cedar cluster

**Install local lookup service using instructions at https://interproscan-docs.readthedocs.io/en/latest/LocalLookupService.html (The lookup service is about 1.4T).**



# Don't use clustering softwares
## Find genes lost in evolved strains 
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

## Find genes gained in evolved strains 
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
#### Transfer scaffolds file from info to graham and cedar
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

# Transfer rhizo_ee.scaffolds.tar.gz and md5sums.txt from info computers to graham computers
(base) [xingyuan@infoserv spades_assembly]$ scp rhizo_ee.scaffolds.tar.gz sux21@cedar.computecanada.ca:/home/sux21/2023_summer_coop/rhizo_ee/genomes
(base) [xingyuan@infoserv spades_assembly]$ scp md5sums.txt sux21@cedar.computecanada.ca:/home/sux21/2023_summer_coop/rhizo_ee/genomes

# Verify the checksum
[sux21@gra-login1 genomes]$ md5sum -c md5sums.txt

# Extract rhizo_ee.scaffolds.tar.gz
[sux21@gra-login1 genomes]$ tar -xf rhizo_ee.scaffolds.tar.gz
```


#### Download pgap.py
https://github.com/ncbi/pgap/tree/1126_Test

Version:  <br>
Work done on 

```bash
module load apptainer 
wget -O pgap.py https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
chmod +x pgap.py
./pgap.py -D apptainer 
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
#Run this script as ./ThisScript INPUT.gff

#This script will make the following change:
#------------------------------------------------------------------
##sequence-region  1 1046175

#will change into

##sequence-region NODE_1_length_1046175_cov_26.396254  1 1046175
#------------------------------------------------------------------

# Step 1: Extract all lines starting with "NODE" 
while IFS=$'\t' read -r a b; do
if [[ "$a" =~ ^"Rht_" ]]; then #Change "NODE" to "cluster" for Rht* samples, and change to "Rht_" for Rht_262_N, Rht_643_N
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
done < $1 > "$1".fixed

### Step 4: Reformat fasta header (currently not needed)
#sed -e 's/lcl|//' -e 's/Rhizobium leguminosarum chromosome, whole genome shotgun sequence//' "$1".int3 > "$1".fixed

### Step 4: Remove intermediate files
rm -f "$1".int*

#References:
#Step 2 command was taken from https://unix.stackexchange.com/questions/194780/remove-duplicate-lines-while-keeping-the-order-of-the-lines. Step 3 commands were taken from  https://stackoverflow.com/questions/72293847/using-sed-in-order-to-change-a-specific-character-in-a-specific-line, https://stackoverflow.com/questions/67396536/sed-insert-whitespace-at-the-nth-position-after-a-delimiter, https://stackoverflow.com/questions/50971418/using-sed-to-replace-single-line-in-while-read-loop
```

### Reformat gff file using AGAT

AGAT Version: v1.2.1 <br>
Work done on info19

```bash
#riboswitch is not taken into account, add this feature for AGAT
#follow instructions from https://agat.readthedocs.io/en/latest/troubleshooting.html
#steps
agat levels --expose # a file called feature_levels.yaml will appear in your current directory
vi feature_levels.yaml # open the file with text editor for editing, Add "riboswitch: 1" to level 1
```

```bash
#!/bin/bash
for i in *fixed; do
/2/scratch/batstonelab/bin/AGAT-1.4.0/bin/agat_convert_sp_gxf2gxf.pl --gff $i --output ${i%_annot_with_genomic_fasta.gff.fixed}.pgap.gff
done
```

```bash
#Check the files are fixed by AGAT properly
ls -hlS *pgap.gff #files with smaller sizes will appear at the bottom of the output. Some sizes are 0 or smaller than the common size. These files may have no output or no DNA sequences at the bottom. Re-run AGAT for these files.
#To check whether the number of sequences in the gff match with the number of sequences in the *filtered.fasta
grep -c ^">" *pgap.gff > gff_seq_number.txt #Count and output the number of DNA sequences in the gff files fixed by AGAT. Use vi to edit each line to "Sample_name:number_of_sequences" to make the files comparable
grep -c ">" *filtered.fasta > fasta_seq_number.txt #Count and output the number of DNA sequences in the fasta files which are used as input files for PGAP. Use vi to edit each line to "Sample_name:number_of_sequences" to make the files comparable
comm -3 <(sort gff_seq_number.txt) <(sort fasta_seq_number.txt) #It is expected no output if the number of sequeneces match
```

## 4. Run roary
https://sanger-pathogens.github.io/Roary/

Roary version: 3.12.0 <br>
Work done on info2020

**Based on ANI values, remove 4_4_10, Rht_773_N, as5_2_4 since they have low ANI (below 90)**

```
nohup /home/xingyuan/tools/miniconda3/envs/roary/bin/roary -v -p 5 -f roary_results -y *pgap.gff &
```

```bash
#Create graphs using instructions from https://github.com/microgenomics/tutorials/blob/master/pangenome.md
python roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv
```



