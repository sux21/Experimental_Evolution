# Find gene presence absence variations (Annotation by PGAP and gene presence absence analysis by Panaroo)
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
