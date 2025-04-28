# Alternative method of finding gene presence absence variations (Annotation by PGAP and gene presence absence analysis by Panaroo)
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

## 2. Send results back to info cluster
```bash
#copy and rename each annot.gbk with sample names
for i in */annot.gbk; do
sample=${i%/annot.gbk}
cp $i "$sample"_annot.gbk
done

#compress results
tar -czvf pgap_annot.tar.gz *gbk

#verify data integrity with MD5
md5sum pgap_annot.tar.gz > pgap_annot.md5

#files transfer
scp pgap_annot.tar.gz xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap
scp pgap_annot.md5 xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap

#verify file integrity (do it on info)
md5sum -c pgap_annot.md5
```

## 3. Gene presence absence analysis using Panaroo
https://gthlab.au/panaroo/#/

**Based on ANI values, remove 4_4_10, Rht_773_N, as5_2_4 since they have low ANI (below 90) - 416 total strains**

Panaroo Version: 1.5.2 <br>
Work done on info2020

```bash
#create a new directory for the results
mkdir panaroo_results

#create a list with eack gbk file on one line
for i in *gbk; do  echo $i ; done > input_files.txt

#run panaroo 
nohup panaroo -i input_files.txt -o panaroo_results --clean-mode strict --remove-invalid-genes &
```

