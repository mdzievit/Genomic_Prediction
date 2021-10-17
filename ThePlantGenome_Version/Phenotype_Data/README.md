# Phenotypic data analysis
Folder contains input files for a few scripts and the corresponding output files that will be used for other analyses.

**Note**: Don't change any of the output locations or file names. This will cause problems with other scripts that read those input files.


## Folders
### Combined
This contains the results of the `Combined_Population_Calculations.R` script. It is analyzing the training + small empirical validation phenotypic data. It will generate BLUPs and heritability estimates and provide those outputs.

### Empirical_Validation
The raw phenotypic data for the small empirical validation population is within this folder (`Empirical_Validation/Empirical_Traits_Raw.txt`. See the manuscript for more information, but it was generated as part of this paper.
Additionally, it contains the results of the `Empirical_Validation_BLUP_Calculations.R` script (BLUPs and heritability estimates).

### Large_Empirical\_Validation
The BLUPs for the large empirical validation panel are included in this folder, see the manuscript for more information on where this originated from. This data is **not** used by any of the scripts in this folder, however, it made sense to store all the phenotypic data here. It will be used in the `Prediction_Accuracy` folder.

The original file can be obtained from either 
[Cornell's Site](http://cbsusrv04.tc.cornell.edu/users/panzea/download.aspx?filegroupid=11 "Cornell's Site") 
or from 
[iPlantCollaborative](http://de.iplantcollaborative.org/dl/d/A1F3CF70-AE3A-4472-845F-AA8D6D7024EA/Peiffer2014Genetics_blupPhenos20150325.xlsx) and then filtered to only contain the GEs of interest (see supplemental table 1 in the manuscript for the list of GEs that were used in the large empirical validation panel.

### Training
The true raw phenotypic data for the training population is contained within `Training/Training_Phenotypes_Raw.txt`. This file was accessed via [Panzea](http://de.iplantcollaborative.org/dl/d/A8E0DB17-879D-486B-8975-5D5FD0986656/traitMatrix_maize282NAM_v15-130212.txt). See manuscript for more information on original citation.

The `Training_BLUP_Calculations.R` script processes the raw file listed previously to generate a formatted file: `Training/Training_Phenotypes_Raw_Formatted.txt` and output a BLUP file and heritability file. All three files are used in downstream scripts, so don't rename/change their location when outputting.

## Scripts
**Note** The summaries below for each script are just an overview, please see the individual script for more comments on what the individual scripts are doing.

`Combined_Population_Calculations.R` - This script takes the small empirical validation phenotypes and the training phenotypes and combines them for calculating BLUPs and heritabilities.

`Empirical_Validation_BLUP_Calculations.R` - Takes the raw phenotypic data for the empirical validation populations and calculates BLUPs and heritabilities.

`Training_BLUP_Calculations.R` - Processes the raw phenotypic data for the training populations and generates 

