## Purpose:
This was developed specifically for the Ames Panel genotypic data file ```ZeaGBSv27_publicSamples_rawGenos_AGPv3_20170206.h5``` (downloaded from ```/iplant/home/shared/panzea/genotypes/GBS/v27/```)

Certain Ames Panel lines were genotyped multiple times as part of the original project, therefore they were combined to form a 'consensus' call for each genotype within the Ames Panel.

## Commandline code to enter

```
python consensus_call.py -i Example_Input_File.vcf -o Example_Output -f Example_Find_File.txt -c False
```

If that doesn't work call up python3 from command line


### If your input file doesn't work, you might need to try running this command on your input and find files 
```
sed -i 's/\r//g' Consensus_Script/Example_Input_File.vcf
```

This will remove the /r (carriage return) which windows adds and not sure if the script likes it or not. Better to remove it.


## Script specific for Ames Panel

```
python consensus_call_AmesDP.py
```

This specific script corrects some genotypes that were labeled differently:

1. PI543841
2. PI537101
3. PI543842
4. PI543845

These four lines have multiple copies (ending in -1...-2...-3) and were treated as the same sample genotyped multiple times
 
