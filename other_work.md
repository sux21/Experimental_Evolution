# GC content across genomes 

Create .genome file with the size of each scaffold

info2020
```
for i in /home/xingyuan/rhizo_ee/derived+original_genomes/*fasta; do
j=${i#/home/xingyuan/rhizo_ee/derived+original_genomes/}
sample=${j%.fasta}

faidx $i -i chromsizes > "$sample".genome
done
```

Create .bed file with 200 Kb window
```
bedtools makewindows -g .genome -w 200000 > .bed
```




