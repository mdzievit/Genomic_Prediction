### Workflow on how the genotypic data was processed in Dzievit et al.

1) Use tassel to convert input file ```ZeaGBSv27_publicSamples_rawGenos_AGPv3_20170206.h5``` (downloaded from 
```/iplant/home/shared/panzea/genotypes/GBS/v27/``` into a VCF file:

    tassel-5-standalone/./run_pipeline.pl -Xms512m -Xmx10g -h5 ZeaGBSv27_publicSamples_rawGenos_AGPv3_20170206.h5 -export OutFile -exportType VCF

2) Run the custom python script ```consensus_call_AmesDP``` to combine samples genotyped multiple times and to split into training and prediction sets:

    consensus_call_AmesDP -i OutFile.vcf -o Training_Out1 -f Prediction_Population.txt -c False

    consensus_call_AmesDP -i OutFile.vcf -o Prediction_Out1 -f Training_Population.txt -c False

	
3) Zip the files to save space:

	gzip Training_Out1.vcf

	gzip Prediction_Out1.vcf

4) Use VCFtools v 0.1.15 to split into each chromosome and filter the SNPs for the training population (used a for loop (for chr in 1..10) with variable 'chr'):

	vcftools --gzvcf Training_Out1.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.20 --maf -0.01 --chr "$chr" --recode --recode-INFO-all --out chr"$chr"_Training_Out1 -c | gzip -c > chr"$chr"_Training_Out1.vcf.gz
	
5) Use VCFtools v 0.1.15 to split into each chromosome and filter the SNPs for the training population (used a for loop (for chr in 1..10) with variable 'chr')

	vcftools --gzvcf Prediction_Out1.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.20 --chr "$chr" --recode --recode-INFO-all --out chr"$chr"_Prediction_Out1 | gzip -c > chr"$chr"_Prediction_Out1.vcf.gz
	
6) Use VCFtools v 0.1.15 to find the common SNPs between the Training and Prediction populations (used a for loop (for chr in 1..10) with variable 'chr') and then create a file of the SNPs to exclude (meaning they are NOT in common)

	vcftools --gzvcf chr"$chr"_Training_Out1.vcf.gz --gzdiff chr"$chr"_Prediction_Out1.vcf.gz --diff-site -c | grep -v 'B' | awk '{ if($2 == ".") print$1,$3; else print $1, $2;}' > chr"$chr"_exclude_snps.txt
	
7) Use VCFtools v 0.1.15 to subset the two files using the excluded SNP file (used a for loop (for chr in 1..10) with variable 'chr') and then prepare the files to join together with the prediction population using [tabix v1.6](http://www.htslib.org/doc/tabix.html)

	vcftools --gzvcf chr"$chr"_Training_Out1.vcf.gz --exclude-positions chr"$chr"_exclude_snps.txt --out chr"$chr"_Training_Out2 --recode -c | gunzip -c > chr"$chr"_Training_Out2.vcf.gz
	tabix -p vcf chr"$chr"_Training_Out2.vcf.gz

8) Do the same thing for the prediction population (used a for loop (for chr in 1..10) with variable 'chr') 

	vcftools --gzvcf chr"$chr"_Prediction_Out1.vcf.gz --exclude-positions chr"$chr"_exclude_snps.txt --out chr"$chr"_Training_Out2 --recode -c | gunzip -c > chr"$chr"_Prediction_Out2.vcf.gz
	tabix -p vcf chr"$chr"_Prediction_Out2.vcf.gz
	
9) Use VCFtools v 0.1.15 to merge the two populations used a for loop (for chr in 1..10) with variable 'chr')

	vcf-merge chr"$chr"_Training_Out2.vcf.gz chr"$chr"_Prediction_Out2.vcf.gz | bgzip -c > chr"$chr_combined.vcf.gz
	
10) Use Beagle v 4.1 to impute the missing data, used a for loop (for chr in 1..10) with variable 'chr')

	java -Xss5m -Xmx32G -jar beagle.27Jan18.7e1.jar gt=chr"$chr"_combined.vcf.gz nthreads=16 window=10000 out=chr"$chr"_combined_imputed