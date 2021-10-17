#### Notes ####
#### The purpose of this package is to use the Maize Association Population to predict the Ames Diversity Panel
#### The inputs are different, otherwise the code should be the same.

#### If you want to use the Small Empirical Validation Population to predict phenotypes for the Maize Association Population
#### then just swap around the input files.

#### Package Loading ####
##These will load the packages if you are running this on a server:
# library(tibble, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/")
# library(tidyr, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/")
# library(dplyr, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/")
# library(readr, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/")
# library(rrBLUP, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/")
# library(data.table, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/")

##Load packages if running on a standard PC
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(rrBLUP)
library(data.table)

#### Training Input Data ####
#### Look at the examples to run locally. This is just a smaller data set to test the script
gen_train <- fread("Predict_Phenotypes/Examples/Training_Imputed.rrblup")

#### This would be the full imputed file to run, it is not present as the file would be too large to store here
# gen_train <- fread("RRBLUP_Input_Files/Training_Imputed.rrblup")

gen_train <- as.matrix(gen_train)

#### Load up the names of the individuals, match up phenotype + genotype data
indv_file <- "RRBLUP_Input_Files/Examples/Training_Imputed.012.indv"

#### Real file location
#indv_file <- "Predict_Phenotypes/Training_Imputed.012.indv"
genID <- read_tsv(indv_file,
                  col_names = FALSE) %>% 
  rename(Genotype = 'X1') %>% 
  mutate(Gen_Ord = row_number())

#### Load in phenotypes analyzed in the phenotypic director, see Phenotype Data folder
#### Join these with the genotypic data, then sort them so they are in the same order
phenotypes <- read_tsv("Phenotype_Data/Training/Training_Phenotypes_BLUPs.txt") %>% 
  mutate(Genotype = toupper(Genotype)) %>% 
  mutate(Genotype = case_when(
    Genotype == "W22_R-RSTD" ~ "W22R-R-STD_CS-2909-1",
    TRUE ~ Genotype
  )) %>% 
  left_join(genID,
            by = "Genotype") %>% 
  arrange(Gen_Ord) %>% 
  gather(Trait,Value,-Genotype,-Gen_Ord)

#### Pulling list of traits
traits <- phenotypes %>% 
  select(Trait) %>% 
  unique() %>% 
  pull(Trait)

### Setting up the matrices we'll be storing the data to
marker_effects <- data.frame(Marker = 1:ncol(gen_train))
marker_betas <- data.frame(Marker = 1:ncol(gen_train))

#### Estimating marker effects ####
#### Progress bar, keep track of progress
pb <- txtProgressBar(min = 0, max = length(traits), style = 3)

#### Loop to run throuh each trait and figure out marker effects
for (i in 1:length(traits)) {
  #### Filter phenotypic data to trait of interest
  y <- phenotypes %>% 
    filter(Trait == traits[i]) %>% 
    unique()
  
  #### Pull out training phenotypic data
  train_y <-y %>% 
    pull(Value)
  
  #### Extract genotypic data and make sure it is in the same order as phenotypic data
  train_Z <- gen_train[y %>% 
                 pull(Gen_Ord),]
  
  #### Run the solver from RRBLUP to estimate marker effects
  ans <- mixed.solve(y = train_y,
                     Z = train_Z)
  
  #### Pull marker effects and add marker info to it
  marker_effects <- cbind(marker_effects,data.frame(Effects = ans$u))
  names(marker_effects)[i + 1] <- traits[i]
  
  #### Add in beta, instead of resulting GEBVs being just summation of marker effects, it is that plus beta, so it 
  #### resembles an actual phenotype 
  marker_betas <- cbind(marker_betas,data.frame(Betas = ans$beta))
  names(marker_betas)[i + 1] <- traits[i]
  setTxtProgressBar(pb, i)
  gc()
}

rm(list = setdiff(ls(),c("marker_effects", "marker_betas","traits")))
gc()

marker_effects <- marker_effects %>% 
  select(-Marker)
marker_betas <- marker_betas %>% 
  select(-Marker) %>% 
  unique()

#### Loading GEs to predict ####
#### Example data set, not full one
gen_val <- fread("RRBLUP_Input_Files/Examples/Validation_Imputed.rrblup")

#### Full data set
#gen_val <- fread("RRBLUP_Input_Files/Validation_Imputed.rrblup")

gen_val <- as.matrix(gen_val)

#### Setting up data to predict
valid_set <- data.frame(Gen = 1:nrow(gen_val))

#### Generating GEBVs ####
pb <- txtProgressBar(min = 0, max = length(traits), style = 3)
for (i in 1:length(traits)) {
  
  #### Setting up GEs to predict, then summation of marker effects
  valid_set <- cbind(valid_set,data.frame(Pred = gen_val %*% as.matrix(marker_effects[,i]) + as.vector(marker_betas[,i])))
  names(valid_set)[i + 1] <- traits[i]
  
  setTxtProgressBar(pb, i)
}

#### Output results ####
write_tsv(path = "Predict_Phenotypes/Example_Validation_Imputed_Predicted.txt",
          x = valid_set)
# write_tsv(path = "Predict_Phenotypes/Validation_Imputed_Predicted.txt",
#           x = valid_set)