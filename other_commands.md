# Lists of other commands 

**use a script (filename.sh) to extract the first entry and rename the first line with sample ID**
```
#!/bin/bash
read firstline;
sampleid=${firstline//>/> 10_1_8-}
echo $sampleid 
while read line;
do
  if [[ $line =~ ">" ]]; then
  exit;
  fi
echo $line      
done 
```
