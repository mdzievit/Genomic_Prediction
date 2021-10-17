library(tidyverse)
library(broom)
library(lme4)

#### Loading raw data ####
#### Load the raw phenotypic data for the empirical population
emp_raw <- read_tsv("Phenotype_Data/Empirical_Validation/Empirical_Traits_Raw.txt") %>%
  mutate(Population = "Empirical") %>% 
  select(-Pedigree)
  
#### Load the raw phenotypic data for the training population
train_raw <- read_tsv("Phenotype_Data/Training/Training_Phenotypes_Raw_Formatted.txt") %>% 
  mutate(Population = "Training") %>% 
  select(-Env_Code,-e)

#### Pulling all of the traits of interest
traits <- emp_raw %>% 
  select(Trait) %>% 
  unique() %>% 
  pull()

#### Combining the two data sets and filtering for only traits of interest
combined_raw <- train_raw %>% 
  filter(Trait %in% traits) %>% 
  bind_rows(emp_raw)

#### Mixed Model Function ####
mixed_model <- function(data, return_type,genotype_term) {
  if (genotype_term == "fixed") {
    orig_mod <- lmer(Value ~ (1 | Env) + Genotype,
                data = data %>%
                  mutate(Env = factor(Env)),
                REML = TRUE)
    orig_geno <- data %>% 
      select(Genotype) %>% 
      unique() %>% 
      arrange()
    
    mean <- unname(fixef(orig_mod)[1])
    mod <- as.tibble(fixef(orig_mod)) %>%
      rownames_to_column(var = "Genotype") %>%
      filter(Genotype != "(Intercept)") %>% 
      separate(Genotype,into = c("Delete","Genotype"),sep = "Genotype") %>%
      select(-Delete) %>% 
      mutate(mu = mean,
             BLUE = mu + value) %>% 
      select(-mu,-value) %>% 
      bind_rows(c(orig_geno[1,1],"BLUE" = mean)) %>% 
      arrange(Genotype)
      
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
                 comp = "Variance"))
  }
}

##Running the BLUPs code
blups <- combined_raw %>%
  group_by(Trait) %>%
  do((mixed_model(data = .,
                  return_type = "effects",
                  genotype_term = "random"))) %>%
  ungroup()

#### Output of the BLUPs ####
# write_tsv("Phenotype_Data/Combined/Combined_Phenotypes_BLUPs.txt",
#           x = blups %>% 
#             spread(Trait,BLUP))


#### Calculate harmonic means, variance components, and heritability ####
var_comp <- combined_raw %>%
  group_by(Trait) %>%
  do(as_tibble(mixed_model(data = .,
                           return_type = "var",
                           genotype_term = "random"))) %>%
  ungroup()

#### Calculate harmonic means and output to the screen
(harmonic_mean <- combined_raw %>%
    na.omit() %>%
    group_by(Trait,Genotype) %>%
    summarise(Num_Env = n()) %>%
    ungroup() %>%
    group_by(Trait,Num_Env) %>%
    summarise(n = n()) %>%
    mutate(add = (1/Num_Env)*n) %>%
    summarise(e = sum(n)/sum(add)))

#### Calculate heritability and output to the screen
(herit <- var_comp %>%
  select(Trait,grp,vcov) %>%
  spread(grp,vcov) %>%
  left_join(harmonic_mean) %>%
  mutate(H = (Genotype/(Genotype + (Residual/e)))))

#### Output the heritability variance components and harmonic means ####
# write_tsv("Phenotype_Data/Combined/Combined_Heritabilities.txt",
#           x = herit)