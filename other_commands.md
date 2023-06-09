# Lists of other commands 
**Task: for 52 scaffolds.fasta files, copy the largest scaffolds, rename it with its sample id, and store all of them in a single fasta file**

**Command:**
```
#!/bin/bash
for i in 
read firstline; echo $firstline
while read line;
do
  if [[ $line =~ ">" ]]; then
  sampleid=${line//>/> 10_1_8}
  exit;
  fi
echo $sampleid 
echo $line
done
```
