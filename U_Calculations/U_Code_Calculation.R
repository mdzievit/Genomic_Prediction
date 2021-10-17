#### This script was setup to run on a server, you'll need to configure it for that accordingly.
#lib_loc <- "~/R/x86_64-pc-linux-gnu-library/3.4/"

##Use this to run on desktop R studio
#lib_loc <- NULL

# library(MASS, lib.loc = (lib_loc))
# library(data.table, lib.loc = (lib_loc))
# library(doSNOW, lib.loc = (lib_loc))
# library(readr, lib.loc = (lib_loc))

library(MASS)
library(data.table)
library(doSNOW)
library(readr)

#### Input data ####
#### These are gen x markers 
#### In 1 and -1 format
#### Example Data Sets
train <- as.matrix(fread("RRBLUP_Input_Files/Examples/Training_Imputed.rrblup"))
vali <- as.matrix(fread("RRBLUP_Input_Files/Examples/Validation_Imputed.rrblup"))

#### Location of the full data if it was available. Check `Prepare_Imputed_File-RR-BLUP.md`
# train <- as.matrix(fread("RRBLUP_Input_Files/Training_Imputed.rrblup"))
# vali <- as.matrix(fread("RRBLUP_Input_Files/Validation_Imputed.rrblup"))

#### Parallelized U-Value Function ####
par_U <- function(train,
                  vali,
                  cores) {
  ##Sets the variables
  train <- train
  vali <- vali
  cores <- cores
  
  ##Determines the number if individuals in the validation population
  nind <- nrow(vali)
  
  ##Calculates the Inverse matrix of the training population first using an identity matrix
  invTTp <- ginv(tcrossprod(train))
  
  ##Registers the number of cores specified by the user
  cl <- makeSOCKcluster(cores)
  registerDoSNOW(cl)
  
  ##Sets up a progress bar to display
  pb <- txtProgressBar(max = nind, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ##Parallelizes each individual's U calculation and then rbinds them back together
  ##Keeps track of each individual by 'i'
  U_Results <- foreach(i = 1:nind,
                       .combine = 'rbind',
                       .packages = c('MASS'),
                       .options.snow = opts) %dopar% {
                         k <- as.matrix(vali[i,])
                         kHat <- crossprod(train,(invTTp %*% (train %*% k)))
                         U <- crossprod(kHat)/(crossprod(k))
                         out <- c("Indv" = i, "Result" = U[1,1])
                         return(out)
                       }
  stopCluster(cl)
  return(U_Results)
}

#### Running the function ####
U_Results <- par_U(train = train,
                   vali = vali,
                   cores = 2)

#### Output ####
#### Example output
write_tsv(x = as.data.frame(U_Results),
          path = "U_Calculations/Example_U_Results_Valid.txt")
#### Full dataset output
# write_tsv(x = as.data.frame(U_Results),
#           path = "U_Calculations/U_Results_Valid.txt")