# RRBLUP Input Files
Contains the input files used within the various scripts. These are the output files from the `Prepare_Imputed_File_RR-BLUP.md` workflow.

## Examples Folder
Contains example input files that can be used in the `U_Calculations/U_Code_Calculations.R`  and the `Predict_Phenotypes/Gen_Pred.R` scripts. **Typically just a subset of the larger data sets.**

Contains the following files:

- `Training_Imputed.012.indv`: GE Name corresponding to the genotypic data file, 
`Training_Imputed.rrblup`
- `Training_Imputed.rrblup`: Genotpyic data file formatted for RRBLUP, see `Prepare_Imputed_File_RR-BLUP.md` workflow
- `Validation_Imputed.rrblup`- Genotypic data file for the validation GEs.

## Input Files
- `Training_Imputed.012.indv` - The genotypic order of the training population. This was output from preparing the RRBLUP files. In this case, it was the 
- `Validation_Imputed.012.indv` - The genotypic order of the prediction population. In this case this is all of the lines in the Ames Diversity Panel. We make prediction for all 2,812 lines of the population. Later on we remove all of the duplicates from the population.