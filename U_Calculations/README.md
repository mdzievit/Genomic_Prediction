# U_Calculations
Contains the code needed to calculate upper bound for reliability values (U-Value). Additionally, it processes these results and creates some figures.

### Input Files
- `Training_Imputed.rrblup` - Genotypic data of the training set.
- `Validation_Imputed.rrblup` - Genotypic data of the validation or the material you want to predict.

## Scripts
- `U_Code_Calculations.R`- Script that can be run in parallel to do the U_Value calculation. More notes are available within the script. It utilizes some example data sets.
- `U_Values_Pred_Traits.R`- Script that processes the output of the `U_Code_Calculations.R` script to format the U_Values into some figures.

## Output files
`U_Results_Valid.txt`- Output from the `U_Code_Calculations.R` script on the full data set. It contains all of the U_Values for the Ames Diversity Panel.