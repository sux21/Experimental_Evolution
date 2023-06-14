# 2020 Experimental Evoluntion 
Bioinformatics project on *Rhizobium leguminosarum*

Work done on INFO server, node info114 and my local computer, MacBook-Pro. Bash shell is used on both computers. 

# Key questions in this project
1. How did standing genetic variation change according to EE selective treatments (high-N, no plant; low-N, no-plant; high-N, plus plant; low-N, plus plant)
2. What genetic changes occurred throughout EE to each isolate (de novo mutation, small sequence variants (indels)
3. Can we detect HGT by examining presence/absence variation? 

# Step 1 - Genome Assembly using SPAdes <br>

## Before Assembly
### 1. Run FastQC for raw data
Versions: FastQC v0.11.5 <br>
Work done on info114

**Sample: 1_1_2**

Run FastQC:
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads 1_1_2* &
```
**All 363 samples (726 files)**
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_raw_reads *gz &
```

### 2. Run MultiQC for raw data
Version: MultiQC v1.9 <br>
Work done on info114

```
multiqc . 
```
### 3. Run Trimmomatic 
Version: 0.39 <br>
Work done on info114

**Sample: 1_1_2**
```
nohup java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_001.fastq.gz /home/xingyuan/rhizo_ee/raw_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_UP_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz /home/xingyuan/rhizo_ee/trimmomatic_reads/1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_UP_001.fastq.gz ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 CROP:130 LEADING:3 TRAILING:3 MINLEN:36 &
```
**All 363 samples (726 files)**
```
nano run_trimmomatic.sh
```
```
#!/bin/bash 
for R1 in *R1* 
do 
R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz} 
R1_P=${R1//001.fastq.gz/P_001.fastq.gz} 
R1_UP=${R1//001.fastq.gz/UP_001.fastq.gz} 
R2_P=${R2//001.fastq.gz/P_001.fastq.gz} 
R2_UP=${R2//001.fastq.gz/UP_001.fastq.gz} 

java -jar /usr/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE /home/xingyuan/rhizo_ee/raw_reads/$R1 /home/xingyuan/rhizo_ee/raw_reads/$R2 /home/xingyuan/rhizo_ee/trimmomatic_reads/$R1_P /home/xingyuan/rhizo_ee/trimmomatic_reads/$R1_UP /home/xingyuan/rhizo_ee/trimmomatic_reads/$R2_P /home/xingyuan/rhizo_ee/trimmomatic_reads/$R2_UP ILLUMINACLIP:/usr/local/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE HEADCROP:15 CROP:130 LEADING:3 TRAILING:3 MINLEN:36
done
```
```
nohup bash run_trimmomatic.sh &
```

### 4. Run FastQC for trimmed reads 
Versions: FastQC v0.11.5 <br>
Work done on info114

**Sample: 1_1_2**

Run FastQC: 
```
fastqc -o /home/xingyuan/rhizo_ee/fastQC_trimmomatic_reads *_P_* 
```

**All 363 samples (726 files)**
```
nohup fastqc -o /home/xingyuan/rhizo_ee/fastQC_trimmomatic_reads *_P_* &
```

### 5. Run MultiQC for trimmed reads
Version: MultiQC v1.9 <br>
Work done on info114

**All 363 samples (726 files)**
```
multiqc . 
```

## During Assembly 
### 1. Run SPAdes 
Version: 3.15.2 <br>
Work done on info114

**Sample: 1_1_2**
```
nohup spades.py --careful -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2 &
```

**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
vi running-spades.sh
```
```
#!/bin/bash 
for R1 in 10_1_8_*R1_P_* 10_1_9_*R1_P_* 10_7_6_*R1_P_* 11_4_2_*R1_P_* 11_4_4_*R1_P_* 11_5_6_*R1_P_* 13_4_1_*R1_P_* 14_4_6_*R1_P_* 14_5_3_*R1_P_* 15_4_4_*R1_P_* 15_4_6_*R1_P_* 16_1_6_*R1_P_* 16_1_7_*R1_P_* 16_1_8_*R1_P_* 16_4_2_*R1_P_* 16_4_3_*R1_P_* 16_6_6_*R1_P_* 17_2_1_*R1_P_* 17_2_8_*R1_P_* 17_2_9_*R1_P_* 18_1_4_*R1_P_* 18_1_5_*R1_P_* 19_1_1_*R1_P_* 19_5_8_*R1_P_* 2_2_5_*R1_P_* 2_5_2_*R1_P_* 2_6_4_*R1_P_* 3_1_5_*R1_P_* 3_2_1_*R1_P_* 3_2_3_*R1_P_* 3_2_6_*R1_P_* 3_2_7_*R1_P_* 3_3_5_*R1_P_* 3_3_7_*R1_P_* 3_3_9_*R1_P_* 4_1_2_*R1_P_* 4_1_4_*R1_P_* 4_2_1_*R1_P_* 6_4_5_*R1_P_* 6_4_7_*R1_P_* 6_7_5_*R1_P_* 7_1_2_*R1_P_* 7_1_5_*R1_P_* 7_6_3_*R1_P_* 7_6_9_*R1_P_* 7_7_2_*R1_P_* 7_7_3_*R1_P_* 8_4_10_*R1_P_* 8_4_4_*R1_P_* 9_3_7_*R1_P_* 9_7_6_*R1_P_* 9_7_9_*R1_P_* 
do
R2=${R1//R1_P_001.fastq.gz/R2_P_001.fastq.gz}

spades.py --careful -1 $R1 -2 $R2 -o /home/xingyuan/rhizo_ee/spades_assembly/${R1%_*_L002_*gz}
done
```
```
nohup bash running-spades.sh &
```
manipulating string (useful for writing loop): https://mywiki.wooledge.org/BashFAQ/100, https://tldp.org/LDP/abs/html/string-manipulation.html

### 2. Run plasmidSPAdes
**Sample: 1_1_2**
```
nohup spades.py --plasmid -1 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R1_P_001.fastq.gz -2 1_1_2_GAACTGAGCG-CGCTCCACGA_L002_R2_P_001.fastq.gz -o /home/xingyuan/rhizo_ee/spades_assembly/1_1_2-plasmids &
```

## After Assembly 
### 1. Run Quast on scaffolds
Version: 5.2.0, 3d87c606 <br>
Work done on info114 

**Sample: 1_1_2** <br>
```
quast.py -m 0 scaffolds.fasta
```
**Samples: 52 samples from 2020 strains in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
quast -m 0 scaffolds.fasta 
```
**Samples: 28 samples from 2008 in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets**
```
#!/bin/bash 
for i in Rht*
do

/home/xingyuan/tools/quast-5.2.0/quast.py $i -m 0 -o /home/xingyuan/rhizo_ee/quast_2008_long_reads/${i%.fasta}
done
```

# Step 2 - Data Analysis
## Question 1 - How closely related a subset of the experimentally evolved isolates are to the 2008 strain set?
### Method 1: Spine-Nucmer-SNPs 
Commands are taken from https://github.com/Alan-Collins/Spine-Nucmer-SNPs. 

**Samples: 52 samples from 2020 + 28 samples from 2008 in Rhizobium_leguminosarum_EE2021-Single_strain_experiment Google sheets** <br>

Copy contigs.fasta from ``rhizo_ee/spades_assembly`` to ``rhizo_ee/2008_2020_strains_comparison``:

```
#!/bin/bash 
for i in 10_1_8 13_4_1 15_4_6 16_4_2 17_2_8 19_1_1 2_6_4 3_2_6 3_3_9 6_4_5 7_1_5 7_7_3 9_7_6 10_1_9 11_4_2 14_4_6 16_1_6 16_4_3 17_2_9 19_5_8 3_1_5 3_2_7 4_1_2 6_4_7 7_6_3 8_4_10 9_7_9 10_7_6 11_4_4 14_5_3 16_1_7 16_6_6 18_1_4 2_2_5 3_2_1 3_3_5 4_1_4 6_7_5 7_6_9 8_4_4 11_5_6 15_4_4 16_1_8 17_2_1 18_1_5 2_5_2 3_2_3 3_3_7 4_2_1 7_1_2 7_7_2 9_3_7
do

cp /home/xingyuan/rhizo_ee/spades_assembly/$i/contigs.fasta /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/"$i-contigs.fasta"
done
```

#### (1) Run Spine
Version: 0.3.2 <br>
Work done on info114. 

```
ls *.fasta | awk 'BEGIN { FS="\t"; OFS="\t" } { print "/home/xingyuan/rhizo_ee/2008_2020_strains_comparison/"$1, $1, "fasta" }' > ./SPINE/config.txt

spine.pl -f /home/xingyuan/rhizo_ee/2008_2020_strains_comparison/SPINE/config.txt 
```

#### (2) Run Nucmer 
Version: 3.1 <br>
Work done on info114

```
ls ../SPINE/*.core.fasta | while read i; do acc=${i%.core*}; acc=${acc#../SPINE/output.}; nucmer --prefix=${acc}_core ../SPINE/output.backbone.fasta $i; delta-filter -r -q ${acc}_core.delta > ${acc}_core.filter; show-snps -Clr ${acc}_core.filter > ${acc}_core.snps; done
```

#### (3) Run snps2fasta.py
Work done on info114
```
python3 snps2fasta.py -r ../SPINE/output.backbone.fasta -f variant_core.fasta -p '(.*)_core\.snps' ../NUCMER/*.snps

-r reference backbone fasta file
-f output fasta file
-p regex pattern to capture genome ID from filename for use as fasta header and rowname in snp matrix
```

#### (4) Run fasta2diffmat.py
Work done on info114

```
python fasta2diffmat.py -f variant_core.fasta -d diff_dict.pkl -t 5 -g SNP_dist_hist.png -c under_2500_SNP_dist_hist.png -ct 2500

-f input aligned multifasta file produced by snps2fasta.py script
-d output filename with .pkl extension for pickled distance dict of format {(Genome1, Genome2) : #SNPs_between_them}
-t number of threads to use. Default: 1
-g filename for histogram of snp distances
-c filename for histogram of realtively highly related assemblies snp distances
-ct threshold for histogram of realtively highly related assemblies snp distances. Plots a histogram of only the SNP distances below this level

ERROR:
[xingyuan@info114 SNPS]$ python fasta2diffmat.py -f variant_core.fasta -d diff_dict.pkl -t 5 -g SNP_dist_hist.png -c under_2500_SNP_dist_hist.png -ct 2500
reading in fasta
Traceback (most recent call last):
  File "fasta2diffmat.py", line 141, in <module>
    plt.yscale('log', nonpositive='clip')
  File "/home/xingyuan/.local/lib/python3.5/site-packages/matplotlib/pyplot.py", line 3084, in yscale
    return gca().set_yscale(value, **kwargs)
  File "/home/xingyuan/.local/lib/python3.5/site-packages/matplotlib/axes/_base.py", line 3704, in set_yscale
    ax.yaxis._set_scale(value, **kwargs)
  File "/home/xingyuan/.local/lib/python3.5/site-packages/matplotlib/axis.py", line 767, in _set_scale
    self._scale = mscale.scale_factory(value, self, **kwargs)
  File "/home/xingyuan/.local/lib/python3.5/site-packages/matplotlib/scale.py", line 569, in scale_factory
    return _scale_mapping[scale](axis, **kwargs)
  File "/home/xingyuan/.local/lib/python3.5/site-packages/matplotlib/scale.py", line 249, in __init__
    "{!r}".format(kwargs))
ValueError: provided too many kwargs, can only pass {'basex', 'subsx', nonposx'} or {'basey', 'subsy', nonposy'}.  You passed {'nonpositive': 'clip'}
```
#### (5) Run get_snps_support_MP.py
```
paste <(ls ../NUCMER/*.snps) <(ls ../SPINE/*.core_coords.txt) <(ls ../ASSEMBLIES/*.fasta) <(ls ../SAMS/*.sam) <(ls ../SAMS/) | awk '{gsub("../SAMS/","",$5)}1 {gsub(".sam","",$5)}1' | sed 's/ /\t/g' > config.txt
```

## Presence-Absence Variation analysis
### 1. Gene prediction by Glimmer 
https://academic.oup.com/bioinformatics/article/23/6/673/419055

### 2. Protein function prediction by InterProScan (omit genes greater than or equal to 5000bp)
https://academic.oup.com/bioinformatics/article/30/9/1236/237988

### 3. Cluster genes by CD-HIT (minimum identity set to 90%, -aL and -AL set to 70%
https://academic.oup.com/bioinformatics/article/17/3/282/189639?login=false 









