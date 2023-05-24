
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

results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        225      
# contigs (>= 1000 bp)     35       
Total length (>= 0 bp)     7477155  
Total length (>= 1000 bp)  7440926  
# contigs                  46       
Largest contig             1176701  
Total length               7448904  
GC (%)                     60.74    
N50                        564110   
N75                        200370   
L50                        5        
L75                        11       
# N's per 100 kbp          8.32  
```
### 154B
commands
```
spades.py -1 GSF2234-154B_S8_R1_P_001.fq.gz -2 GSF2234-154B_S8_R2_P_001.fq.gz --isolate -o spades-154B
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        300      
# contigs (>= 1000 bp)     64       
Total length (>= 0 bp)     7589416  
Total length (>= 1000 bp)  7540765  
# contigs                  75       
Largest contig             879135   
Total length               7548301  
GC (%)                     60.60    
N50                        498124   
N75                        131695   
L50                        6        
L75                        13       
# N's per 100 kbp          6.76
```
### 164A
commands
```
spades.py -1 GSF2234-164A_S9_R1_P_001.fq.gz -2 GSF2234-164A_S9_R2_P_001.fq.gz --isolate -o spades-164A
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        266      
# contigs (>= 1000 bp)     69       
Total length (>= 0 bp)     7408350  
Total length (>= 1000 bp)  7371337  
# contigs                  81       
Largest contig             1047071  
Total length               7378442  
GC (%)                     60.71    
N50                        375436   
N75                        229011   
L50                        7        
L75                        13       
# N's per 100 kbp          9.89 
```
### 177C
commands
```
spades.py -1 GSF2234-177C_S11_R1_P_001.fq.gz -2 GSF2234-177C_S11_R2_P_001.fq.gz --isolate -o spades-177C
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        276      
# contigs (>= 1000 bp)     69       
Total length (>= 0 bp)     7342297  
Total length (>= 1000 bp)  7302780  
# contigs                  87       
Largest contig             965700   
Total length               7315114  
GC (%)                     60.65    
N50                        275845   
N75                        147112   
L50                        8        
L75                        16       
# N's per 100 kbp          6.01 
```
### 214C
command
```
nohup spades.py -1 GSF2234-214C_S15_R1_P_001.fq.gz -2 GSF2234-214C_S15_R2_P_001.fq.gz --isolate -o spades-214C &
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        277      
# contigs (>= 1000 bp)     47       
Total length (>= 0 bp)     7357250  
Total length (>= 1000 bp)  7315949  
# contigs                  53       
Largest contig             1490525  
Total length               7320607  
GC (%)                     60.66    
N50                        273955   
N75                        163034   
L50                        7        
L75                        16       
# N's per 100 kbp          5.74 
```
### 218A
command
```
nohup spades.py -1 GSF2234-218A_S16_R1_P_001.fq.gz -2 GSF2234-218A_S16_R2_P_001.fq.gz --isolate -o spades-218A &
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        295      
# contigs (>= 1000 bp)     35       
Total length (>= 0 bp)     7203126  
Total length (>= 1000 bp)  7154097  
# contigs                  41       
Largest contig             1032152  
Total length               7157999  
GC (%)                     60.80    
N50                        559785   
N75                        164838   
L50                        5        
L75                        10       
# N's per 100 kbp          10.20  
```
### 272A
command
```
nohup spades.py -1 GSF2234-272A_S22_R1_P_001.fq.gz -2 GSF2234-272A_S22_R2_P_001.fq.gz --isolate -o spades-272A &
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        442      
# contigs (>= 1000 bp)     28       
Total length (>= 0 bp)     7196338  
Total length (>= 1000 bp)  7106057  
# contigs                  33       
Largest contig             1690701  
Total length               7109725  
GC (%)                     60.85    
N50                        521954   
N75                        321218   
L50                        5        
L75                        9        
# N's per 100 kbp          4.50  
```
### 295A
command
```
nohup spades.py -1 GSF2234-295A_S28_R1_P_001.fq.gz -2 GSF2234-295A_S28_R2_P_001.fq.gz --isolate -o spades-295A &
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        608      
# contigs (>= 1000 bp)     29       
Total length (>= 0 bp)     7226762  
Total length (>= 1000 bp)  7102816  
# contigs                  34       
Largest contig             1101732  
Total length               7106221  
GC (%)                     60.86    
N50                        521252   
N75                        243938   
L50                        5        
L75                        10       
# N's per 100 kbp          6.05  
```
### 298A
command 
```
nohup spades.py -1 GSF2234-298A_S30_R1_P_001.fq.gz -2 GSF2234-298A_S30_R2_P_001.fq.gz --isolate -o spades-298A &
```
results
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   scaffolds
# contigs (>= 0 bp)        237      
# contigs (>= 1000 bp)     48       
Total length (>= 0 bp)     7280941  
Total length (>= 1000 bp)  7249696  
# contigs                  55       
Largest contig             1186861  
Total length               7254498  
GC (%)                     60.69    
N50                        273623   
N75                        175800   
L50                        8        
L75                        17       
# N's per 100 kbp          10.06   
```
