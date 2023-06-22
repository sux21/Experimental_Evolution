```
### Take out the first scaffold from each fasta file, and put all scaffolds into a single file ordered from largest to smallest.

### Take out the largest scaffolds
#!/bin/bash
for x in 10_1_8 10_1_9 10_7_6 11_4_2 11_4_4 11_5_6 13_4_1 14_4_6 14_5_3 15_4_4 15_4_6 16_1_6 16_1_7 16_1_8 16_4_2 16_4_3 16_6_6
 17_2_1 17_2_8 17_2_9 18_1_4 18_1_5 19_1_1 19_5_8 2_2_5 2_5_2 2_6_4 3_1_5 3_2_1 3_2_3 3_2_6 3_2_7 3_3_5 3_3_7 3_3_9 4_1_2 4_1_4
 4_2_1 6_4_5 6_4_7 6_7_5 7_1_2 7_1_5 7_6_3 7_6_9 7_7_2 7_7_3 8_4_10 8_4_4 9_3_7 9_7_6 9_7_9; do

echo "Start: $x" 

i=0
while read line; do
  if [[ $line =~ ">" ]]; then  
    let i=$i+1
  fi
  if [ $i -eq 1 ]; then
    sampleid=${line//>/> $x-}
    echo $sampleid  
  fi

done < /home/xingyuan/rhizo_ee/spades_assembly/$x/scaffolds.fasta > $x-largest_scaffold.fasta

echo "Finish: $x"

done 


### put into order
cat `ls -1S *largest_scaffold.fasta` > ordered_largest_scaffolds.fasta 
```

```
### Find a gene in a genome

### Make the genome the database 
#!/bin/bash
for i in *.fasta; do
echo "makeblastdb -in $i -title "${i%.fasta}" -dbtype nucl"
done

### download the gene fasta file from NCBI and use it as a query 
#!/bin/bash
for i in /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/ASSEMBLIES/*.fasta; do
x=${i#/home/xingyuan/rhizo_ee/2008_2020_strains_comparison/ASSEMBLIES/}

if [[ $x =~ "contigs" ]]; then
    y=${x%-contigs.fasta}
fi

if [[ $x =~ "Rht" ]]; then
    y=${x%.fasta}
fi

blastn -query NZ_CP050089.1[121783..123303].fa -out $y-nifD.blast -db $i

done
```
```
### put the gene in a fasta file
blastn -query NZ_CP050089.1[121783..123303].fa -outfmt "6 sseqid length sseq" -out 10_1_8-nifD -db /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/ASSEMBLIES/10_1_8-contigs.fasta; grep "NODE" 10_1_8-nifD | awk -F ' BEGIN {OFS="\n"}{print ">gene:nifD, sample:10_1_8, which_contig_it_locates:"$1 ", " "gene_length:"$2" bp",$3}' > 10_1_8-nifD.fasta; rm -f 10_1_8-nifD
```
