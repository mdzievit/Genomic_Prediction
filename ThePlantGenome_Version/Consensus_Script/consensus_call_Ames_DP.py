import argparse

parser = argparse.ArgumentParser(description='Use a line file to find duplicate entries to \
combine and make a consensus call')

parser.add_argument('--input_file', '-i', dest = 'input_file',
                    help='Ames_DP file that needs to create a consensus call for')
parser.add_argument('--output_file', '-o', dest = 'output_file',
                    help='Pefix for the output file of the consensus file, stored as .vcf')
parser.add_argument('--amesdp_lines', '-f', dest = 'amesdp_lines',
                    help='List of the Ames Panel lines that need to create a consensus call')
parser.add_argument('--printComments', '-c', dest = 'printComments',
                    help='True | False to print out comments')
args = parser.parse_args()

##Reads in the line file and then removes the last line if it is a blank
#with open("Ames_Panel-Data_Analysis/Genotype_Imputation/SNP_Uniquer/Unique_Yu_Panel_lines.txt",'r') as lines:
with open(args.amesdp_lines, 'r') as lines:
    ames_lines = [x.strip('\n') for x in lines.readlines()]

##Need to get the header row for the genotype names
with open(args.input_file,  'r') as lines:
#with open("Ames_Panel-Data_Analysis/Genotype_Imputation/SNP_Uniquer/small_head_chr10_amesdp.vcf",  'r') as lines:
    holder = lines.readlines()
    for lines in holder:
        if(lines[0:2] == '#C'):
            gen_header = lines.strip('/n').split('\t')
            

##Creates a dictionary that contains the column number for the specific GBS ID, use
##this to find GBS ids associated with all name copies
line_key = {}
for i in range(len(gen_header)):
    name = gen_header[i].lower()
    line_key[name] = i

##Creates a dictionary of unique pedigree names and their corresponding GBS ID,
##designed to use with the line_key
pedigree_key = {}
for i in range(len(gen_header)):
    spot = gen_header[i].find(':')
    name = gen_header[i][0:spot].lower()
    ##These 3 lines have increments added to the end, but the metadata file (panzea.org) says 
    ##that the are the exact same line, so in this situation, we are combining them
    if (name == 'pi543841-1' or name == 'pi543841-2' or name == 'pi543841-3'):
        name = 'pi543841'
    elif (name == 'pi537101-1' or name == 'pi537101-2' or name == 'pi537101-3'):
        name = 'pi537101'
    elif (name == 'pi543842-1' or name == 'pi543842-2' or name == 'pi543842-3'):
        name = 'pi543842'
    elif (name == 'pi543845-1' or name == 'pi543845-2' or name == 'pi543845-3'):
        name = 'pi543845'
    if (spot == -1):
        pass
    elif(bool(pedigree_key.get(name) == False)):
        pedigree_key[name] = gen_header[i].lower()
    else:
        pedigree_key.setdefault(name,[]).append(gen_header[i].lower())
        
##Creates a dictionary of the ames lines we want to find and the corresponding column number
##from the pedigree_key and line_key dictionary
found_key = {}
for i in range(len(ames_lines)):
    name = ames_lines[i].lower()
    if(bool(found_key.get(name) == False)):
        found_key[name] = line_key[pedigree_key.get(name)]
    else:
        found_key.setdefault(name,[]).append(pedigree_key.get(name))
        
fileName = args.output_file + '.vcf'
outputFile = open(fileName,'w')
#outputFile = open("Ames_Panel-Data_Analysis/Genotype_Imputation/SNP_Uniquer/chr10_consensus.vcf",'w')

##Want to ignore all lines starting with ## and only print the header row (#C).
##Write the first 10 columns of each line as they are the VCF file info
##Then want to calculate the consensus call
#inFile = open("Ames_Panel-Data_Analysis/Genotype_Imputation/SNP_Uniquer/small_head_chr10_amesdp.vcf",  'r')
inFile = open(args.input_file,  'r')
line = inFile.readline()
while line:
    current = line.strip().split('\t')
    if (current[0][0:2] == '##' and args.printComments == 'True'):
        outputFile.write('\t'.join(current) + '\n')
    elif (current[0][0:2] == '##' and args.printComments == 'False'):
        pass
    elif (current[0][0:2]== '#C'):
        out = current[0:9]
        for key in found_key:
            out.append(key.upper())
        outputFile.write('\t'.join(out) + '\n')
    else:
        out = current[0:9]
        for key in found_key:
            current_ids = pedigree_key.get(key)
            gen_calls = []
            for id in current_ids:
                gen_calls.append(current[line_key.get(id)])
            if len(gen_calls) == 1:
                out.append(gen_calls[0])
            elif len(gen_calls) == 2:
                if gen_calls[0] == gen_calls[1]:
                    out.append(gen_calls[0])
                elif gen_calls[0] == './.' and gen_calls[1] != './.':
                    out.append(gen_calls[1])
                elif gen_calls[1] == './.' and gen_calls[0] != './.':
                    out.append(gen_calls[0])
                else:
                    out.append('./.')
            else:
                a = 0
                b = 0
                n = 0
                het1 = 0
                het2 = 0
                for snp in gen_calls:
                    if snp == '1/1':
                        a += 1
                    elif snp == '0/0':
                        b += 1
                    elif snp == './.':
                        n += 1
                    elif snp == '1/0':
                        het1 += 1
                    else:
                        het2 += 1
                total = a + b + het1 + het2
                if total == 0:
                    out.append('./.')
                elif a/total > .5:
                    out.append('1/1')
                elif b/total > .5:
                    out.append('0/0')
                elif het1/total > .5:
                    out.append('1/0')
                elif het2/total > .5:
                    out.append('0/1')
                else:
                    out.append('./.')
        outputFile.write('\t'.join(out) + '\n')
    line = inFile.readline()

inFile.close()
outputFile.close()
