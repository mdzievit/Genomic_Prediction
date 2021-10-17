# Prediction Accuracy
Contains scripts used to analyze the predictions and determine prediction accuracy from the observed phenotypic data.This will be done across the following datasets:

- Cross-validation within the training population
	- Maize Association Population
	- Small Empirical Validation Population
	- Combined Population (Maize Association Population + Small Empirical Validation Population)
- Small Empirical Validation Population 
- Large Empirical Validation Population


## Scripts
`Cross-Validation_Analysis`- Takes input from a variety of sources and determines prediction accuracy within the different training populations used throughout the study. It will output figures and tables to the appropriate locations (see `Final_Tables` and `Final_Figures`)

`Small_Empirical_Validation_Prediction_Accuracy.R` - Takes input from a variety of sources to determines prediction accuracy using the Maize Association Population as the training population and the small empirical validation population as the material we validated those predictions with. It will output figures and tables to the appropriate locations (see `Final_Tables` and `Final_Figures`)

`Large_Empirical_Validation_Prediction_Accuracy.R` - Takes input from a variety of sources to determines prediction accuracy using the Maize Association Population as the training population and the large empirical validation population as the material we validated those predictions with. It will output figures and tables to the appropriate locations (see `Final_Tables` and `Final_Figures`)
