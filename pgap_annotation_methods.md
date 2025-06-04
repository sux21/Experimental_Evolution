# Re-do gene presence absence analysis: results are in this new directory ``genes_pav``.

## 1. Use cleaned genomes for 19_1_9 and 19_4_7

```bash
ln -s /home/xingyuan/rhizo_ee/split_genomes/19_1_9_output/output_bins/SemiBin_0.fa 19_1_9_SemiBin_0.fasta
ln -s /home/xingyuan/rhizo_ee/split_genomes/19_4_7_output/output_bins/SemiBin_1.fa 19_4_7_SemiBin_1.fasta
```

## 2. Filter sequences shorter than 200 bp (pgap only takes sequences equal or longer than 200 bp)
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

## 3. Annotate genome using pgap (very slow, so include a if statement for restarting at where it stops) 

Version: 2025-05-06.build7983 <br>
Work done on info20

```
/home/xingyuan/tools/pgap.py -r -o mg37_results -g $HOME/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s "Mycoplasmoides genitalium"
```

```bash
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/*.filtered.fasta; do
j=${i#/home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/}
sample=${j%.filtered.fasta}

#delete directory not completed



/home/xingyuan/tools/pgap.py -D /home/xingyuan/tools/bin/apptainer --container-path /home/xingyuan/tools/pgap_2025-05-06.build7983.sif --report-usage-false -o "$sample" --prefix "$sample" -g "$i" -s "Rhizobium leguminosarum" --cpu 6 --no-self-update
done
```

## Below are old commands

**Find samples that don't have annot.gbk:**
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

### Run PGAP for 19_1_9_SemiBin_0 and 19_4_7_SemiBin_1.fasta
Version: 2023-10-03.build7061 <br>
Work done on minigraham cluster 

**transfer from info to graham**
```bash
scp xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/split_genomes/19_1_9_output/output_bins/SemiBin_0.fa 19_1_9_SemiBin_0.fasta
scp xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/split_genomes/19_4_7_output/output_bins/SemiBin_1.fa 19_4_7_SemiBin_1.fasta
```

**Run PGAP**
```bash
#!/bin/bash
#SBATCH --time=01-00:00
#SBATCH --account=def-batstone
#SBATCH --mem=32G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
 
module load apptainer

export APPTAINER_BIND=/project

for i in /project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/19*SemiBin*.fasta; do
j=${i#/project/6078724/sux21/rhizo_ee/genomes/genes_presence_absence_variation/}
sample=${j%.fasta}

/project/6078724/sux21/tools/pgap/pgap.py -D apptainer --container-path /project/6078724/sux21/tools/pgap/pgap_2023-10-03.build7061.sif --no-internet --no-self-update -r -o "$sample" -g "$i" -s 'Rhizobium leguminosarum' -c 40

done
```

## 2. Send results back to info cluster
```bash
#compress results
#!/bin/bash
#SBATCH --time=00-3:00:00
#SBATCH --account=def-batstone
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=sux21@mcmaster.ca
#SBATCH --mail-type=ALL
tar -czvf pgap-scaffolds.tar.gz pgap-scaffolds

#verify data integrity with MD5
md5sum pgap-scaffolds.tar.gz > pgap-scaffolds.md5

#files transfer
scp pgap-scaffolds.tar.gz xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap
scp pgap-scaffolds.md5 xingyuan@info.mcmaster.ca:/home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap

#verify file integrity (do it on info)
md5sum -c pgap-scaffolds.md5

#extract results
tar -xzvf pgap-scaffolds.tar.gz
```

## 3. Gene presence absence analysis using Panaroo
https://gthlab.au/panaroo/#/

Panaroo Version: 1.5.2 <br>
Work done on info2020

**Using R, prepare a file with each MPA on one line (``mpa.txt``). Prepare 26 files with MPA and its derived isolates each on one line (``MPA_"$mpa".txt``)**

```r
#open R on info2020
R

#load library
library(tidyverse)

#find most probable ancestor for each derived isolate
MPA <- read.table("/home/xingyuan/rhizo_ee/derived+original_genomes/most_prob_ancestors.txt") %>%
  group_by(V1) %>%
  slice_max(tibble(V3), n = 1) %>% #select the highest ANI value for each derived isolate
  select(V1, V2) %>%
  rename(isolate=V1, MPA=V2) %>% #rename variable
  mutate(isolate = str_extract(isolate, "[:digit:]+_[:digit:]+_[:alnum:]+"),
         MPA = str_extract(MPA, "Rht_[:digit:]+_(N|C)")) #rename strain name

#create a file with each MPA on one line
MPA_per_line <- data.frame(MPA=unique(MPA$MPA))
write.table(MPA_per_line, "mpa.txt", row.names=F, quote=F, col.names=F)

#create files with MPA and its derived isolates each on one line (26 files total)
ancestors <- unique(MPA$MPA)
output <- vector(mode="list", length(ancestors))
names(output) <- ancestors

for (i in seq_along(ancestors)) {
  output[[i]] <- filter(MPA, MPA==ancestors[i]) %>%
    select(1) 
  
  colnames(output[[i]]) <- names(output)[[i]]
  
  output[[i]] <- rbind(names(output[[i]]), as.data.frame(output[[i]]))
  
  write.table(output[[i]], paste0("MPA_", names(output[[i]]), ".txt"),
              row.names=F, quote=F, col.names=F)
}

#exit R
q()
```

**Make 26 working directories for Panaroo one for each MPA**

```bash
for mpa in `cat mpa.txt`; do
for i in `cat MPA_"$mpa".txt`; do

#skip files not existed
if [ ! -d "$i" ]; then
continue
fi

#make directory for each MPA
if [ ! -d "$mpa"-pav-analysis ]; then
mkdir "$mpa"-pav-analysis
fi

#make symbolic links for GFF files
ln -s /home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap/"$i"/annot_with_genomic_fasta.gff "$mpa"-pav-analysis/"$i".gff

#make symbolic links for FASTA files
ln -s /home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap/"$i"/*.filtered.fasta "$mpa"-pav-analysis/"$i".fasta

done
done
```

**Make input file for Panaroo in each of the 26 working directories**

```bash
#each line is a GFF file and a FASTA file of each isolate
for mpa in *-pav-analysis; do
for sample_gff in "$mpa"/*gff; do
#sample_gff=${i#"$mpa"/}
sample_fasta=${sample_gff/.gff/.fasta}

echo "$sample_gff" "$sample_fasta";
done > "$mpa"/input_files.txt

done
```

**Run Panaroo**

```bash
for i in *pav*; do
/home/xingyuan/tools/miniconda3/bin/panaroo -i "$i"/input_files.txt -o "$i"/panaroo_results --clean-mode strict --remove-invalid-genes --merge_paralogs
done
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

/home/xingyuan/tools/ncbi-blast-2.16.0+/bin/blastn -query pgap_genes_gain.fasta -subject $i -outfmt "10 qacc sacc evalue bitscore length pident nident gapopen gaps qcovs" > pgap_"$sample"_blast.csv

done

#add variable names for the file
for i in *blast.csv; do
sed -i '1s/^/qacc,sacc,evalue,bitscore,length,pident,nident,gapopen,gaps,qcovs\n/' $i
done
```

