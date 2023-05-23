
# Run SPAdes for samples with long read data
## Samples ID provided by Health lab: 051, 125, 144D, 153C, 154B, 164A, 177, 214C, 218A, 225A, 272A, 295A, 298A, 301D, 336, 338A, 372, 377, 377A, 391, 524D.

### 125A
commands:
```
spades.py -1 GSF2234-125A_S3_R1_P_001.fq.gz -2 GSF2234-125A_S3_R2_P_001.fq.gz --isolate -o spades-125A
```

results:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        234      
# contigs (>= 1000 bp)     59       
Total length (>= 0 bp)     7857267  
Total length (>= 1000 bp)  7824375  
# contigs                  65       
Largest contig             1431316  
Total length               7828714  
GC (%)                     60.50    
N50                        281131   
N75                        184128   
L50                        8        
L75                        16       
# N's per 100 kbp          5.36 
```

### 153C
commands
```
spades.py -1 GSF2234-153C_S7_R1_P_001.fq.gz -2 GSF2234-153C_S7_R2_P_001.fq.gz --isolate -o spades-153C
```

results:
```


