# SNP calling using alternative reference strains

info2020

Align 56 reference strains to 56 reference strains

fastANI version 1.32
```bash
nohup /usr/local/bin/fastANI --ql reference.txt --rl reference.txt --threads 5 --matrix -o ref_to_ref.txt &
```

Prepare a csv file as the following: que_name, ref_name
```bash
10_3_2,Rht_415_C
10_5_1,Rht_415_C
12_7_4,Rht_415_C
16_3_4,Rht_003_C
19_6_10,Rht_415_C
19_6_2,Rht_415_C
19_6_8,Rht_415_C
20_6_10,Rht_415_C
5_6_1,Rht_415_C
8_1_6,Rht_415_C
8_1_9,Rht_415_C
8_2_5,Rht_511_N
8_2_9,Rht_511_N
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






