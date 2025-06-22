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

## 3. Gene presence absence analysis using Panaroo
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

## 4. gene gain analysis (for each 26 pangenome)
Work done on info2020

```r
#Open R
R

#install.packages("tidyverse", lib = "/home/xingyuan/tools/R_library", repos = "http://cran.us.r-project.org") #lib = where to install the packages

#load packages
library(tidyverse, lib.loc = "/home/xingyuan/tools/R_library")

#remove all variables existed in the working directory
rm(list=ls())

setwd("/home/xingyuan/rhizo_ee/Genes_PAV/genome_annotation_pgap/")

#make a function to find genes gained for each most probable ancestor and derived isolate group
find_gene_gain <- function(mpa) {
  gene_pre_abs_tab <- read.table(paste0(mpa, "-pav-analysis/panaroo_results/gene_presence_absence.Rtab"), header = TRUE)

  gene_pre_abs_csv <- read.csv(paste0(mpa, "-pav-analysis/panaroo_results/gene_presence_absence.csv"), header = TRUE) 
  
  gene_pre_abs_csv2 <- gene_pre_abs_csv %>%
    pivot_longer(cols = c(4:ncol(gene_pre_abs_csv)), names_to = "derived", values_to = "derived_annotation_id") %>%
    filter(!grepl("Rht", derived, fixed = FALSE)) %>% #remove rows with ancestral isolate
    mutate(derived = str_extract(derived, "[:digit:]+_[:digit:]+_[:alnum:]+")) #rename isolate
  
  gene_data_csv <- read.csv(paste0(mpa, "-pav-analysis/panaroo_results/gene_data.csv"), header = TRUE)
  
  derived_subset <- gene_pre_abs_tab %>%
    pivot_longer(cols = c(2:ncol(gene_pre_abs_tab)), names_to = "derived", values_to = "derived_pres_abs") %>%
    filter(!grepl("Rht", derived, fixed = FALSE)) %>% #remove "Rht" samples and make a subset for the derived isolates
    mutate(derived = str_extract(derived, "[:digit:]+_[:digit:]+_[:alnum:]+")) #rename isolate
  
  ancestral_subset <- gene_pre_abs_tab %>%
    pivot_longer(cols = c(2:ncol(gene_pre_abs_tab)), names_to = "MPA", values_to = "MPA_pres_abs") %>%
    filter(grepl("Rht", MPA, fixed = FALSE)) #select "Rht" samples and make a subset for ancestral strains
  
  #join derived_subset and ancestral_subset
  pres_abs_data <- left_join(derived_subset, ancestral_subset, by = "Gene") %>%
    left_join(gene_pre_abs_csv[,c("Gene","Non.unique.Gene.name","Annotation")], by = "Gene")
  
  #obtain genes gained (1 in isolate, 0 in MPA)
  genes_gained <- filter(pres_abs_data, derived_pres_abs == 1 & MPA_pres_abs == 0) %>%
    left_join(gene_pre_abs_csv2[,c("Gene", "derived", "derived_annotation_id")], by = c("Gene", "derived")) %>% #add annotation id
    distinct(Gene, .keep_all = TRUE) %>% #only keep unique genes
    mutate(derived_annotation_id = str_remove(derived_annotation_id, "_pseudo")) %>% #remove "_pseudo" from anntation id as gene_data_csv does not has this
    left_join(gene_data_csv, by = c("derived"="gff_file", "derived_annotation_id"="annotation_id")) %>% #add DNA sequence
    mutate(derived_annotation_id2 = paste(derived, derived_annotation_id, sep = "_"), .after = derived_annotation_id) #include derived isolate name in derived_annotation_id
  
  #write genes gained to a FASTA file
  genes_gained$derived_annotation_id2 <- paste0(">", genes_gained$derived_annotation_id2)
  
  genes_gained_fasta <- c(rbind(genes_gained$derived_annotation_id2, genes_gained$dna_sequence))
  
  write(x = genes_gained_fasta, file = paste0(mpa, "_gene_gain.fasta"))
}

#loop through all 26 groups of derived isolates and their most probable ancestors

mpa_names <- list.files(path = ".", pattern = "-pav-analysis")
mpa_names <- unique(str_extract(mpa_names, ".*(?=-pav)"))
mpa_names

for (i in 1:length(mpa_names)) {
  find_gene_gain(mpa_names[i])
}
```


## 3. Which ancestral strain mostly had the genes gained?

Blastn Version: 2.16.0 <br>
Work done on info2020

Align DNA sequences of genes gained to all 56 ancestral genomes
```bash
specify format of output: qacc (Query accession), qlen (Query sequence length), sacc (Subject accession), sstart (Start of alignment in subject), send (End of alignment in subject), evalue (Expect value), bitscore (Bit score), length (Alignment length), pident (Percentage of identical matches), nident (Number of identical matches), gapopen (Number of gap openings), gaps (Total number of gaps), qcovs (Query Coverage Per Subject)

for genes_gain_file in *_gene_gain.fasta; do
gene_gain=${genes_gain_file%.fasta}

#skip files with no DNA sequences
content=$(cat -A $genes_gain_file | head -1) 
if [ $content = "$" ]; then
  continue
fi

for reference_file in /home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/Rht*; do
j=${reference_file#/home/xingyuan/rhizo_ee/genes_pav/pgap_method/input_sequences/}
reference=${j%.filtered.fasta}

/home/xingyuan/tools/ncbi-blast-2.16.0+/bin/blastn -query $genes_gain_file -subject $reference_file -outfmt "10 qacc qlen sacc sstart send evalue bitscore length pident nident gapopen gaps qcovs" > pgap_"$gene_gain"_"$reference"_blast.csv

/home/xingyuan/tools/ncbi-blast-2.16.0+/bin/blastn -query $genes_gain_file -subject $reference_file > pgap_"$gene_gain"_"$reference"_blast.out

done
done
```

Add variable names for the file
```bash
for i in *blast.csv; do
sed -i '1s/^/qacc,qlen,sacc,sstart,send,evalue,bitscore,length,pident,nident,gapopen,gaps,qcovs\n/' $i
done
```

Investigate the BLAST results in R
```R
#Open R
R

#install.packages("tidyverse", lib = "/home/xingyuan/tools/R_library", repos = "http://cran.us.r-project.org") #lib = where to install the packages

#load packages
library("tidyverse", lib.loc	= "/home/xingyuan/tools/R_library")

#load BLAST output
blast_out <- vector(mode = "list", length = length(list.files(path = ".", pattern = "_blast.csv")))

names(blast_out) <- list.files(path = ".", pattern = "_blast.csv")

for (i in names(blast_out)) {
  blast_out[[i]] <- read.csv(i, header = TRUE)
}

#add ancestral strain name as a column to each data frame
pgap_blast_results2 <- mapply(cbind, pgap_blast_results, "reference_genome" = names(pgap_blast_results), SIMPLIFY = F)


#merge all data frames in the list to one
blast_results_df <- pgap_blast_results2 %>% 
  Reduce(function(dtf1,dtf2) bind_rows(dtf1,dtf2), .) %>%
  mutate(reference_genome = str_extract(reference_genome, "Rht_[:digit:]+_(N|C)"),
         derived_isolate = str_extract(qacc, "[:digit:]+_[:digit:]+_[:alnum:]+"))

#load metadata of ancestral isolate community and MPA of the derived isolate
Klinger_Rhizobium_2008strains <- readxl::read_xlsx("/Users/xingyuansu/Desktop/rhizo_ee/Klinger_Rhizobium_2008strains.xlsx") %>%
  mutate(ancestral_isolate = str_extract(sample_name, "Rht_[:digit:]+_(N|C)")) %>%
  select(ancestral_isolate, community) %>%
  drop_na(community) %>%
  mutate(number = formatC(as.numeric(str_extract(ancestral_isolate, "[:digit:]+")), width = 3, flag = "0"),
         group = str_extract(ancestral_isolate, "(N|C)"),
         ancestral_isolate = paste0("Rht", "_", number, "_", group)) %>%
  select(-c(number, group))

MPA <- read.table("/Users/xingyuansu/Desktop/rhizo_ee/most_prob_ancestors.txt") %>%
  group_by(V1) %>%
  slice_max(tibble(V3), n = 1) %>% #select the highest ANI value for each derived isolate
  select(V1, V2) %>%
  rename(derived_isolate=V1, MPA=V2) %>% #rename variable
  mutate(derived_isolate = str_extract(derived_isolate, "[:digit:]+_[:digit:]+_[:alnum:]+"),
         MPA = str_extract(MPA, "Rht_[:digit:]+_(N|C)")) #rename strain name

blast_results_df2 <- left_join(blast_results_df, MPA, by="derived_isolate", relationship = "many-to-many") 

```








