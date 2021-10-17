library(tidyverse)
library(lme4)
library(data.table)
library(Hmisc)
library(broom)
library(corrplot)

#### Import Training phenotype data ####
#### This was accessed from panzea 11/2016

p_data <- as.matrix(read.table("Phenotype_Data/Training/Training_Phenotypes_Raw.txt",
                               sep = "\t",
                               stringsAsFactors = FALSE,
                               na.strings = "NaN")) %>% 
  t() %>% 
  as.tibble()

colnames(p_data) <- p_data[1,]
p_data <- p_data[-1,]

#### Correcting genotype names  + Formatting ####
##Two of the genotype names are slightly different from the genotypic data, so 
##we need to correct that to match
p_data_form <- p_data %>% 
  gather(Genotype,Value,-Trait,-Env) %>% 
  mutate(Value = as.numeric(Value)) %>% 
  mutate(Genotype = case_when(
    Genotype == "B2" ~ "B2-GOOD",
    Genotype == "W22_R-rstd_CS-2909-1" ~ "W22_R-RSTD",
    TRUE ~ Genotype
  )) %>%
   spread(Trait,Value) %>% 
  #### Calculating an extra trait
  mutate(TotalKernelWeight = ifelse((is.na(EarWeight) | is.na(CobWeight)),NA,(EarWeight - CobWeight)),
         TotalKernelWeight = ifelse(TotalKernelWeight < 0, NA,TotalKernelWeight)) %>% 
  gather(Trait,Value,-Genotype,-Env)

#### Collecting num envs and harmonic means #### 
##This outputs the unique number of environments this data was evaluated in:
p_data_form %>% 
  select(Env) %>% 
  unique() %>% 
  summarise(Unique_Env = n())

##Calculates the harmonic means for each trait according to Holland et al., 2010
harmonic_mean <- p_data_form %>%
  na.omit() %>% 
  group_by(Trait,Genotype) %>% 
  summarise(Num_Env = n()) %>% 
  ungroup() %>%
  group_by(Trait,Num_Env) %>% 
  summarise(n = n()) %>% 
  mutate(add = (1/Num_Env)*n) %>% 
  summarise(e = sum(n)/sum(add))

#### Trait filtering ####
#### This filters the traits based on involving a border plant and 
#### having harmonic means greater than or equal to 1.5
p_data_sel <- p_data_form %>% 
  filter(!str_detect(Trait,"Border|RowQuality")) %>% 
  mutate(Env_Code = as.numeric(as.factor(Env))) %>% 
  na.omit() %>%
  left_join(harmonic_mean) %>% 
  filter(e >= 1.5)

#### Output this filtered formatted dataset to be used in other applications
# write_tsv(x = p_data_sel,
#           path = "Phenotype_Data/Training/Training_Phenotypes_Raw_Formatted.txt")

#### Mixed model function ####
#### Custom wrapper function that has the mixed model,
#### it can calculate, BLUEs, BLUPs, or variance components
#### This assumes a P = Gen + Env + e model, where the Env is always random, and the genotype is 
#### either fixed (BLUEs) or random (BLUPs)
#### Changing the return type allows you to return var components or effects

mixed_model <- function(data, return_type,genotype_term) {
  if (genotype_term == "fixed") {
    mod <- lmer(Value ~ 1 + (1 | Env) + Genotype,
                data = data %>% 
                  mutate(Env = factor(Env)),
                REML = TRUE)
  } else if (genotype_term == "random") {
    orig_mod <- lmer(Value ~ 1 + (1 | Env) + (1 | Genotype),
                     data = data %>% 
                       mutate(Genotype = factor(Genotype),
                              Env = factor(Env)),
                     REML = TRUE)
    mod <- coef(orig_mod)$Genotype %>%
      rownames_to_column(var = "Genotype") %>%
      select(Genotype,'(Intercept)') %>%
      rename(BLUP = '(Intercept)')
  }
  if (return_type == "effects") {
    return(mod)
  } else if (return_type == "var") {
    return(print(VarCorr(orig_mod),
                 comp ="Variance"))
  }
}

#### Calculating BLUPs ####
####This calculates BLUP values for each trait and combines them into a DF
model_run <- p_data_sel %>%
  group_by(Trait) %>% 
  do((mixed_model(data = .,
                       return_type = "effects",
                       genotype_term = "random"))) %>% 
  ungroup()

#### Outputting the BLUPs ####
# write_tsv(path = "Phenotype_Data/Training/Training_Phenotypes_BLUPs.txt",
#           x = model_run %>%
#             spread(Trait,BLUP))

#### Variance components ####
#### This calculates the variance componets for each trait
model_run_var <- p_data_sel %>% 
  group_by(Trait) %>% 
  do(as.tibble((mixed_model(data = .,
                            return_type = "var",
                            genotype_term = "random")))) %>% 
  ungroup()

#### Heritability ####
#### Combines the variance component data with the harmonic means to calculate heritability
heritability <- model_run_var %>%
  select(Trait,grp,vcov) %>% 
  spread(grp, vcov) %>% 
  select(-Env) %>% 
  left_join(harmonic_mean) %>% 
  mutate(Heritability = (Genotype/(Genotype + (Residual/e))))

#### Output of varianec + heritability ####
# write_tsv(path = "Phenotype_Data/Training/Training_Heritability.txt",
#           x = heritability)