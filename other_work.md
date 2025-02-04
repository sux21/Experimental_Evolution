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






