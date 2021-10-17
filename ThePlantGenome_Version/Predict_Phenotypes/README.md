# Predict Phenotypes
Folder contains the code used to predict trait values for the AmesDP using the maize association population as the training population or the small empirical population as the training population to predict the maize association population.

**Note**: Don't change any of the output locations or file names. This will cause problems with other scripts that read those input files.

**Note #2**: The true input files used for the analysis are not included as they are too large to store here. They can be generated using the `Prepare_Imputed_File_RR-BLUP.md` work flow. 

### Input Files
- `Training_Imputed.012.indv`- Order of the genotypic data. Used in the script to match phenotype to genotype data.
- `Training_Imputed.rrblup` - Genotypic data of the training set.
- `Validation_Imputed.rrblup` - Genotypic data of the validation or the material you want to predict

### Output Files
- `Validation_Imputed_Predicted.txt` - Output file from the `Gen_Pred.R` script when you use the example files.

## Scripts
**Note** The summaries below for each script are just an overview, please see the individual script for more comments on what the individual scripts are doing.

**Also Note** The prediction scripts were use to predict different sets of material. Generally, input files would need to be changed depending on what training and prediction sets you wanted to use. Please pay attention to that.

`Gen_Pred.R` - Takes the training genotypic and phenotypic data to estimate marker effects. Takes a prediction set to predict GEBVs using only the marker data for the prediction set. All of this is done via RRBLUP. Depending on how big your data set is, it might need to be run on a server with lots of resources. 

- **Note**: Right now this script is setup to load and save files from the `Examples` folder. The 'real' data would have to be deposited into this folder and the script updated. See `Prepare_Imputed_File_RR-BLUP.md` for more information.

## Output Files
- `Predicted_BLUPs.txt`- Output file of the `Gen_Pred.R` script. It is the GEBVs for the Ames Diversity Panel that used the Maize Association Panel as the training set. Gen ID matches the order of the genotypic data. This file is used in downstream analyses.
- `Training_Predicted_BLUPs.txt`- Output file of the `Gen_Pred.R`. It is the GEBVs for the Maize Association Population that used the Small Empirical Validation Population as the training set. Gen ID matches the order of the genotypic data. This file is used in downstream analyses.

