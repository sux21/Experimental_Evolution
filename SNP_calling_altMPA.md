# SNP calling using alternative reference strains

info2020

## Align 56 reference strains to 56 reference strains

fastANI version 1.32
```bash
nohup /usr/local/bin/fastANI --ql reference.txt --rl reference.txt --threads 5 --matrix -o ref_to_ref.txt &
```

## SNP calling
BWA version: 0.7.17-r1188 <br>
Samtools version: 1.13 (using htslib 1.13)
 
Index 6 reference genomes for bwa (Create files ending with .amb, .ann, .bwt, .pac, .sa)
```bash
for i in Rht*fasta; do
/usr/bin/bwa index $i
done
```

Run bwa (bwa mem -t 5 -M -R), convert .sam to .bam (samtools view -huS), and produce alignment statistics (samtools stats)
```bash
#!/bin/bash
# Usage: ./ThisScript Derived-AltMPA.csv
#First 5 lines of the Derived-AltMPA.csv looks like the following:
#-----------------------------------------
#10_3_2,Rht_415_C
#10_5_1,Rht_415_C
#10_5_8,Rht_415_C
#12_7_1,Rht_415_C
#12_7_4,Rht_415_C
#.
#.
#.
#-----------------------------------------
while IFS=',' read -r a b; do # a=derived_strain,b=most_probable_ancestor
r1=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$a"_*R1_P*
r2=/home/xingyuan/rhizo_ee/fastp_results/fastp_reads/"$a"_*R2_P*
ref=/home/xingyuan/rhizo_ee/snps_with_different_MPAs/references/"$b".fasta

/usr/bin/bwa mem -t 5 -M -R "@RG\tID:"$a"\tSM:"$a $ref $r1 $r2 | /usr/local/bin/samtools view -huS -o $a-"$b".bam - 

done < $1
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






