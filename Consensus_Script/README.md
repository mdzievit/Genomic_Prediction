##Commandline code to enter

```
python consensus_call.py -i Example_Input_File.vcf -o Example_Output -f Example_Find_File.txt -c False
```

If that doesn't work call up python3 from command line


###If your input file doesn't work, you might need to try running this command on your input and find files 
```
sed -i 's/\r//g' Consensus_Script/Example_Input_File.vcf
```

This will remove the /r (carriage return) which windows adds and not sure if the script likes it or not. Better to remove it.


##Script specific for Ames Panel

```
python consensus_call_AmesDP.py
```

This specific script corrects some genotypes that were labeled differently:

1. PI543841
2. PI537101
3. PI543842
4. PI543845

These four lines have multiple copies (ending in -1...-2...-3) and were treated as the same sample genotyped multiple times
 