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

## 3. Annotate genome using pgap

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

## -------------- Below are old commands ---------------

**Find samples that don't have annot.gbk:**
```bash
find . -type d \! -exec test -e '{}/annot_with_genomic_fasta.gff' \; -print
#Output: 4_4_10, 2_4_11, 19_4_7, 19_1_9, Rht_773_N
```

## -------------- Above are old commands ---------------

## 4. Gene presence absence analysis using Panaroo
https://gthlab.au/panaroo/#/

Panaroo Version: 1.5.2 <br>
Work done on info2020

**Using R, prepare a file with each MPA on one line (``mpa.txt``). Prepare 26 files with MPA and its derived isolates each on one line (``MPA_"$mpa".txt``)**

R version 4.3.1
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

**Run Panaroo for each of the 26 groups of MPA and its descendents**

```bash
for i in *pav*; do
/home/xingyuan/tools/miniconda3/bin/panaroo -i "$i"/input_files.txt -o "$i"/panaroo_results --clean-mode strict --remove-invalid-genes --merge_paralogs
done
```

## 5. Which ancestral isolates are the genes gained originated

Blastn Version: 2.16.0 <br>
Work done on info2020

Align DNA sequences of genes gained to all 56 ancestral genomes
```bash
for i in *_gene_gain.fasta; do
mpa_group=${i%_gene_gain.fasta}

for reference in /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/Rht*; do
j=${reference#/home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/}
ancestal_isolate=${j%.filtered.fasta}

/home/xingyuan/tools/ncbi-blast-2.16.0+/bin/blastn -query $i -subject $reference -outfmt "10 qacc qlen sacc slen qstart qend sstart send evalue bitscore length pident nident mismatch gapopen gaps sstrand qcovs qcovhsp qcovus" > MPA_"$mpa_group"_subject_"$ancestal_isolate"_blast.csv

done
done
```

Add variable names for the file
```bash
for i in *blast.csv; do
printf 'qacc,qlen,sacc,slen,qstart,qend,sstart,send,evalue,bitscore,length,pident,nident,mismatch,gapopen,gaps,sstrand,qcovs,qcovhsp,qcovus\n' | cat - $i > temp && mv -f temp $i
done
```
