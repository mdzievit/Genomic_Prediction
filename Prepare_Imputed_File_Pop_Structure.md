### Workflow on how the inputed genotypic data was processed for population structure in Dzievit et al.

1) Used the imputed files (output Step 10) of the process ```Process_Genotypic_Data_Workflow``` and made sure all the columns for each chr were sorted the same way using Use VCFtools v 0.1.15:

	for chr in{2..10}; do  vcf-shuffle-cols -t chr1_combined_imputed.vcf.gz chr"$chr"_combined_imputed.vcf.gz > gzip -c | chr"$chr"_combined_imputed_sorted.vcf.gz

2) Merged all the sorted files together using VCFtools v 0.1.15:

	vcf-concat chr1_combined_imputed.vcf.gz chr2_combined_imputed_sorted.vcf.gz ....chr10_combined_imputed_sorted.vcf.gz | gzip -c > all_merged_imputed.vcf.gz

3) Used [plink v1.9](https://www.cog-genomics.org/plink2) to convert vcf file into plink.file

	plink --vcf all_merged_imputed.vcf.gz --plink --out all_merged_imputed_plink

4) Used plink v 1.9 to prune SNPs that were in LD

	plinkv1.9 --file all_merged_imputed_plink --indep-pairwise 50 10 0.2
	plinkv1.9 --file all_merged_imputed_plink --extract plink.prune.in --make-bed --out all_merged_imputed_plink_pruned

5) Ran cross-validation according to [admixture manual](http://www.genetics.ucla.edu/software/admixture/download.html) that tested k =  {1..10}:
	
	for k in {1..10}; do admixture --cv all_merged_imputed_plink_pruned -j16 $k | tee log"$k".out; done

6) ```grep``` to pull out the CV estimates for each log file
	
	grep -h "CV" log*.out > CV.txt

7) Used plink v 1.9 to run Principal Component Analysis (PCA) and output the first 10 eigenvectors

	plinkv1.9  --vcf all_merged_imputed.vcf.gz --pca 10 --out all_merged_imputed_pca

