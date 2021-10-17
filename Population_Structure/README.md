# Population Structure Analysis
Folder contains input files for a few scripts and the corresponding output files that will be used for other analyses.

**Note**: Don't change any of the output locations or file names. This will cause problems with other scripts that read those input files.


## Admixture Folder

Some of the input files for the `Population_Structure_Comparison.R` script are located within this folder.

### Input files
- `all_merged_imputed.5.Q`: Q assignment from the admixture analysis. See `Prepare_Imputed_File_Pop_Structure.md` for more information. Needs to pair up with `all_merged_imputed_pca.nosex` to match up accession with Q value assignment.
- `all_merged_imputed_pca.nosex`: this file was an output of the merged file when preparing the genotypic files for population structure analysis. See `Prepare_Imputed_File_Pop_Structure.md` for more information. It is the accession order for the results files.
- `CV_Output.txt`: CV results from the admixture analysis. See `Prepare_Imputed_File_Pop_Structure.md` for more information.
- `Maize_Association_Panel.txt`: Another version of the training population file. It contains the accessions that are in the training population. It contains the same names as the trainining population file in the main folder, except some of the names are capitalized, etc. 
- `Validation_Panel.txt`: Another version of the validation file. It contains the accessions that are in the validation population. It contains the same names as the validation population file in the main folder, except some of the names are capitalized, etc.

### Output files
- `K5_results.txt`: Output file from the `Population_Structure_Comparison.R` script. It is used by other scripts to know what subpopulation each accession was assigned.

## PCA Folder
File that was output from the PCA analysis. Plink v1.9 was used for the analysis. See `Prepare_Imputed_File_Pop_Structure.md` for more information.

### Input files
- `all_merged_imputed_pca.eigenvec`: output file from the PCA analysis, see `Prepare_Imputed_File_Pop_Structure.md` for more information. It contains the first ten eigenvectors from the analysis. 

## Scripts
**Note** The summaries below for each script are just an overview, please see the individual script for more comments on what the individual scripts are doing.

`Combined_Population_Calculations.R` - This script takes the small empirical validation phenotypes and the training phenotypes and combines them for calculating BLUPs and heritabilities.

`Empirical_Validation_BLUP_Calculations.R` - Takes the raw phenotypic data for the empirical validation populations and calculates BLUPs and heritabilities.

`Training_BLUP_Calculations.R` - Processes the raw phenotypic data for the training populations and generates 

