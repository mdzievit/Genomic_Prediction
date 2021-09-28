# Comprehensive analytical and empirical evaluation of genomic prediction across diverse accessions in maize

##### Matthew J. Dzievit, Xianran Li, Tingting Gui, Jianming Yu

##### Recently accepted for publication in the Plant Genome (08/2021)

This houses the scripts and analysis for the empirical validation of genomic prediction for leaf angle in maize.


## Processing Genotypic Data and Imputing

1. Start with 'Process_Genotypic_Data_Workflow' to take the raw genotypic data and apply the filters and formatting we conducted for this manuscript. This will generate the imputed genotypic data that is used for population structure and genomic prediction.
2. 'Prepare_Imputed_File_Pop_Structure' sets up the imputed data into the format needed for population structure.
3. 'Prepare Imputed_File_RR-BLUP' sets up the imputated data into the format needed for genomic prediction.


## Working with Imputed Data
Once the imputed data has been generated, it can be used for genomic prediction.
Check out the `Predict_Phenotypes` folder

**Please Note**

The actual processed and imputed genotypic data will not be included in this repository. It will need to be generate via the scripts above. The data set is too large and is reproducible following the methods above.

## Consensus_Script Folder
Contains the python script and example files used to make a consensus call for multiple samples of the same GE. Check out the README file within that folder.

## Cross-Validation Folder
Script and results of the cross-validation analysis. Check out the README within for more information.

## Excluded_SNPs Folder
Excluded SNPs from the genotypic data preparation. Check out the `Process_Genotypic_Data_Workflow.md`. Check out the README within for more information.

## Figure 1 Folder
The raw SVG file for creating figure 1. It was the only figure that was created outside of R. Check out the README within for more information.

## Final_Figures Folder
Place where all of the figures from the scripts get deposited. Check out the README within for more information.

## Phenotype_Data Folder
Contains the raw phenotypic data, and the scripts that process and summarize the data. Additionally, the downloaded phenotypic data from the Maize Association Panel is also here. Check out the README within for more information.

## Population_Structure Folder
Contains the scripts and input data to summarize the population structure results. Check out the README within for more information.

## Predict_Phenotypes Folder
Contains the script used to predict GEBVs from the training and validation data. Additionally, the predicted GEBVs for the Ames Diversity Panel are stored here. Check out the README within for more information.

## Prediction_Accuracy Folder
Contains all of the scripts that process and summarize the predicted data to determine prediction accuracy. It goes through the cross-validation results, and the small and large empirical validation populations. Check out the README within for more information.

## RRBLUP_Input\_Files Folder
Contains the GE names for the training and validation (Ames Diversity Panel). Preparing the RRBLUP files, we drop the GE name from the genotypic data file. This tracks and keeps them in order. Additionally, some example files are stores here that are used in other scripts. Check out the README within for more information.

## U_Calculations Folder
Contains the Upper Bound for Reliability (U-Value) calculation script, and the script used to process the U-Value data. Check out the README within for more information.

## Files
Summarize the files stored in the main part of this folder:

- `excluded_ges.txt`: There are duplicated GEs within the Ames Diversity Panel. These list of GEs are excluded from the analysis.
- `Genomic_Prediction.Rproj`: R Project file. Can be used to load a project within RStudio.
- `Prediction_Population.txt`: All the lines within the Ames Diversity Panel.
- `Prepare_Imputed_File_Pop_Structure.md`: See above
- `Prepare_Imputed_File_RR-BLUP.md`: See above 
- `Process_Genotypic_Data_Workflow.md`: See above
- `README.md`: Describes the contents of this main folder
- `Session_Information.pdf`: Breakdown of the packages used in this project, contains the package versions and the R Session information.
- `Session_Information.Rmd`: script that makes the .pdf file listed above
- `Training_Population.txt`: List of GEs within the training population (Maize Association Population)
- `Trait_Order.txt`: Order of the traits for plotting. Used to keep things consistent.
