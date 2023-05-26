# 2018 Strains
For genome infomation, see [Rhizobium leguminosarum bv. trifolii]
## 214C
Plasmids SPAdes command:
```
nohup spades.py -1 GSF2234-214C_S15_R1_P_001.fq.gz -2 GSF2234-214C_S15_R2_P_001.fq.gz --plasmid -o plasmids_spades-214C &
```
Quast command:
```
quast.py scaffolds.fasta 
```
result:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        15       
# contigs (>= 1000 bp)     7        
Total length (>= 0 bp)     783681   
Total length (>= 1000 bp)  781059   
# contigs                  9        
Largest contig             312200   
Total length               782300   
GC (%)                     60.19    
N50                        281438   
N75                        281438   
L50                        2        
L75                        2        
# N's per 100 kbp          25.57
```
