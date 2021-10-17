### Workflow on how the inputed genotypic data was processed for RR-BLUP in Dzievit et al.

1) Used the merged files of the process ```Prepare_Imputed_File_Pop_Structure``` and split the vcf file into training and prediction populations

	vcftools --gzvcf all_merged_imputed.vcf.gz --out training_imputed --keep training_population_genotypes.txt --recode -c | gzip -c > training_imputed.vcf.gz
	vcftools --gzvcf all_merged_imputed.vcf.gz --out training_imputed --keep prediction_population_genotypes.txt --recode -c | gzip -c > prediction_imputed.vcf.gz

2) Need to convert the genotype calls to 0, 1, and 2, where the number represents the number of non-reference alleles. Missing genotypes are represented by -1
	
	vcftools --gzvcf training_imputed.vcf.gz --012 --out training_imputed
	vcftools --gzvcf prediction_imputed.vcf.gz --012 --out prediction_imputed

3) Need to convert these 012 files into -1,0,1 format for [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf) in R
	
	cut -f 1 --complement training_imputed.012 | sed  "s/0/-1/g; s/1/0/g; s/2/1/g; s/-0/-1/g" > training_imputed.rrblup
	cut -f 1 --complement prediction_imputed.012 | sed  "s/0/-1/g; s/1/0/g; s/2/1/g; s/-0/-1/g" > prediction_imputed.rrblup



	