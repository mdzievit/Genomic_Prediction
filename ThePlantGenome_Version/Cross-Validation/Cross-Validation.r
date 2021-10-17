library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(rrBLUP)
library(data.table)
library(doSNOW)

#### Genotypic Input Data ####
#### This loads the training genotypes that were generated in the Prepare_Imputed_File-RR-BLUP
#### Loads the genotype data it is a indv x SNP matrix
#### Example data set
gen_train <- fread("RRBLUP_Input_Files/Examples/Training_Imputed.rrblup")

#### Where the larger data set would be stored
# gen_train <- fread("RRBLUP_Input_Files/Training_Imputed.rrblup")

##Converts it to matrix
gen_train <- as.matrix(gen_train)

##Loads the genotype IDs so we can match the phenotype with the genotype
genID <- read_tsv("RRBLUP_Input_Files/Examples/Training_Imputed.012.indv",
                  col_names = FALSE) %>% 
  rename(Genotype = 'X1') %>% 
  mutate(Gen_Ord = row_number())


#### Phenotypic Input ####
#### This loads the phenotype data and makes sure that a source id matched the genotype phenotype data
#### This is a indx x phenotype matrix - it is the output from 'Training_BLUP_Calculations.R'
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

#### Cross-Validation Function ####
#### This is a cross-validation function that takes an input of the Z matrix (genotype) data, a DF of phenotypes,
#### the current trait, and a fold
#### It will then split the genotype/phenotype data into approximately equal crossFolds, then perform
#### RR-BLUP on each fold iterating throuh each one. For example, if crossFold = 2, it will split data into
#### 2 groups, (A and B) A will be used to predict B, stores prediction and obs in a DF, then use B to predict A,
#### add that data to the DF. Then it will calculate the correlation and return that value.

cross_val <- function(Z,y,crossFold,curTrait) {
  #### Formats the phenotypic data and pulls current trait
  y_form <- y %>%
    filter(Trait == curTrait) %>% 
    na.omit()
  
  #### Formats and sets up the random generator for each grou
  ran <- y_form %>% 
    select(Gen_Ord) %>% 
    mutate(GenNum = sample(1:nrow(y_form),replace = FALSE),
           Group = as.numeric(cut(GenNum,crossFold)))
    
  #### Empty Validation Set
  valid_set <- NULL
  
  #### Run through each crossfold
  for (i in 1:crossFold) {
    #### Sets up the Z matrix
    train_Z <- Z[ran %>% 
                      filter(Group != i) %>% 
                      pull(Gen_Ord),]
    
    #### Pulls out the validation lines in the current fold and their phenotypic data
    train_y <- ran %>%
      filter(Group != i) %>% 
      select(Gen_Ord) %>%
      arrange(Gen_Ord) %>% 
      left_join(y_form,
                by = "Gen_Ord") %>% 
      pull(Value)
    
    #### Pulls out the genotypic data for validation
    val_Z <- Z[ran %>% 
                      filter(Group == i) %>% 
                      pull(Gen_Ord),]
    #### Pulls out the phenotypic data for the validation lines
    val_y <- ran %>%
      filter(Group == i) %>% 
      select(Gen_Ord) %>%
      arrange(Gen_Ord) %>% 
      left_join(y_form,
                by = "Gen_Ord") %>% 
      pull(Value)
    
    #### Solves for the marker effectes 
    ans <- mixed.solve(y = train_y,
                       Z = train_Z)
    
    #### Solves the validation lines and generates the GEBVs.
    #### Adds a column for the true phenotypes to make correlations easier
    valid_set <- rbind(valid_set,data.frame(Pred = val_Z %*% as.matrix(ans$u) + as.vector(ans$beta),
                                            True = val_y)) 
  }
  #### Returns just the correlation as that is what we are interested in obtaining!
  return(cor(valid_set$Pred,
             valid_set$True))
  gc()
}

#### KFold testing setup ####
#### This sets up the matrix for all of our testing, k-folds, traits, and reps
kfold <- NULL 
folds <- c(2,5,10)
reps <- 50
traits <- phenotypes %>% 
  select(Trait) %>%
  unique() %>%
  pull(Trait)

#### This builds our matrix that feeds the cross fold function
for (i in 1:length(traits)) {
  for (j in 1:length(folds)) {
    kfold <- rbind(kfold,
                   cbind(Rep = seq(1:reps),
                         Fold = rep(folds[j],reps),
                         Trait = rep(traits[i],reps),
                         Cor = NA)) 
  }
}

#### Running KFold ####
#### Converts this matrix to a tibble so we have proper numeric values
kfold <- as_tibble(kfold) %>% 
  mutate(Rep = as.numeric(Rep),
         Fold = as.numeric(Fold))

#### This sets up a progress bar for the multi-processing. This way you know how far the script is and that
#### it is actually running
pb <- txtProgressBar(max = nrow(kfold), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#### This sets the number of clusters to run. This is a very memory intensive script. For example, i could 
#### maybe run 6-7 clusters with 32GB of RAM. If you have lots of ram, then you can handle many cores.
cl <- makeSOCKcluster(2)
registerDoSNOW(cl)

#### This runs the multi-processes. The results will be stored in this variable. each trait, kfold,rep can
#### be sent to a core and multi-processed. 
kfold_out <- foreach(i = 1:nrow(kfold),
                     .combine = 'rbind',
                     .packages = c('tibble','dplyr','tidyr','rrBLUP'),
                     .options.snow = opts) %dopar% {
                       out <- c(unlist(kfold[i,1:3]), cross_val(Z = gen_train,
                                                     y = phenotypes,
                                                     curTrait = kfold$Trait[i],
                                                     crossFold = kfold$Fold[i]))
                       return(out)
                     }

#### Need this to stop the clusters
stopCluster(cl)

#### Output data ####
#### Output the results of the cross-validation
write_tsv(path = "Training-Cross_validation_results_kfold.txt",
          x = as.tibble(kfold_out))

