# Cross Validation
The code contained within this folder assumes that you have processed all the genotypic data (See `Process_Genotypic_Data_Workflow.md` and `Prepare_Imputed_File_RR_BLUP.md`) and you have processed all the phenotypic data for the training population (See `Phenotype_Data/README.md`), Additional comments were made within the script.

The examples files are located within `RRBLUP_Input_Files/Examples`and not intended to be the product used to generate the final files used in the manuscript. Those would be too large to store on github.


## Output files

- `Combined_Cross-Validation_All.txt`: kfold results of the combined training populations. Combined the Maize Association Population + Small Empirical Validation population.
- `Empirical-Cross_validation_results_kfold.txt`: kfold results of the empirical validation population
- `Training-Cross_validation_results_kfold.txt`: results of the kfold cross validation for the Maize Association Population.



