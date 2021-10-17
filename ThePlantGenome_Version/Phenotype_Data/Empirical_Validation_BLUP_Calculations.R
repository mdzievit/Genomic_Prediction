library(tidyverse)
library(broom)
library(lme4)
library(ggpubr)
library(patchwork)


## Load the raw phenotypic data for the empirical population
emp_raw <- read_tsv("Phenotype_Data/Empirical_Validation/Empirical_Traits_Raw.txt")

#### Correlation Plot ####
#### Function that builds each row of the raw data correlation plot seen in Sup Fig 4
cor_plot <- function(data,
                     cur_trait,title = FALSE,not_bottom = TRUE) {
  if(title) {
    x_add <- 0.125
  } else {
    x_add <- 0.1
  }
  temp <- data %>% 
    filter(Trait == cur_trait) %>%
    pivot_wider(id_cols = c(Trait,Genotype,Pedigree),
                names_from = c(Env),
                values_from = Value)
  
  my_theme <- theme(plot.margin = unit(c(0,0.35,0,0), "lines"),
                    plot.title = element_text(hjust = 0.5,size = 10),
                    panel.grid = element_blank(),
                    strip.text = element_text(size = 8))
  
  plot1 <- temp %>% 
    ggscatter(x = "IA15", y = "IA16", add = "reg.line",size = 1) +
    stat_cor(aes(label = ..r.label..),
             size = 2.5,color = 'red',fontface = 'bold',
             cor.coef.name = 'r',method = 'pearson') +
    theme_bw() +
    my_theme
  
  plot2 <- temp %>% 
    ggscatter(x = "IA15", y = "IA17", add = "reg.line", size = 1) +
    stat_cor(aes(label = ..r.label..),
             size = 2.5,color = 'red',fontface = 'bold',
             cor.coef.name = 'r',method = 'pearson') +
    theme_bw() +
    my_theme
  
  plot3 <- temp %>% 
    ggscatter(x = "IA16", y = "IA17", add = "reg.line",size = 1) +
    stat_cor(aes(label = ..r.label..),
             size = 2.5,color = 'red',fontface = 'bold',
             cor.coef.name = 'r',method = 'pearson') +
    theme_bw() +
    geom_text(x = max(temp$IA16,na.rm = TRUE) + max(temp$IA16,na.rm = TRUE) * x_add,
              y = (min(temp$IA17,na.rm = TRUE) + max(temp$IA17, na.rm = TRUE))/2,
              label = cur_trait,angle = 270,
              size = 3.25) +
    coord_cartesian(clip = "off",
                    xlim = c(min(temp$IA16,na.rm = TRUE),max(temp$IA16,na.rm = TRUE))) +
    theme(plot.margin = unit(c(0,0.8,0,0), "lines"),
          plot.title = element_text(hjust = 0.5,size = 10),
          panel.grid = element_blank())
  
  if(title) {
    plot1 <- plot1 + ggtitle('IA15 vs IA16')
    plot2 <- plot2 + ggtitle('IA15 vs IA17')
    plot3 <- plot3 + ggtitle('IA16 vs IA17')
  }
  if (not_bottom) {
  }
  out <- plot1 + plot2 + plot3
  return(out)
}

#### builds the correlation row for each trait
ph <- cor_plot(data = emp_raw,
               cur_trait = 'PlantHeight')
eh <- cor_plot(data = emp_raw,
               cur_trait = 'EarHeight',
               not_bottom = FALSE)
mla <- cor_plot(data = emp_raw,
                cur_trait = 'MiddleLeafAngle')
ula <- cor_plot(data = emp_raw,
                cur_trait = 'UpperLeafAngle',
                title = TRUE)
#### Using patchwork to piece these all together
cor_out <- ula / mla / ph / eh

# ggsave(filename = 'Phenotype_Data/Empirical_Validation/Supplemental_Figure_S4.png',
#        plot = cor_out,
#        width = 7,
#        height = 8,
#        dpi = 1200)

# ggsave(filename = 'Final_Figures/Supplemental_Figure_S4.pdf',
#        plot = cor_out,
#        width = 180,
#        height = 203.2,
#        units = 'mm',
#        dpi = 1200)

#### Mixed model function for the BLUP analysis ####
mixed_model <- function(data, return_type,genotype_term) {
  if (genotype_term == "fixed") {
    mod <- lmer(Value ~ (1 | Env) + Genotype,
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

#### Running the mixed model and BLUPs ####
#### BLUP analysis
blups <- emp_raw %>%
  group_by(Trait) %>%
  do((mixed_model(data = .,
                  return_type = "effects",
                  genotype_term = "random"))) %>%
  ungroup()

#### Outputs the BLUPs for future use
# write_tsv("Phenotype_Data/Empirical_Validation/Empirical_Phenotypes_BLUPs.txt",
#           x = blups)

#### Plotting Sup Fig S3 ####
#### Plotting options
plot_theme <- theme(aspect.ratio = .75,
                    axis.title = element_text(size = 5,
                                              face = "bold"),
                    axis.text = element_text(size = 4,
                                             face = "bold",
                                             color = "black"),
                    plot.margin = unit(c(0,.1,0,.1),"cm"),
                    strip.background = element_blank(),
                    strip.text = element_text(size = 4,
                                              face = "bold"),
                    panel.grid = element_blank())

#### Combined histogram that includes the raw data across 3 env, and then the combined BLUPs
(emp_hist <- emp_raw %>%
    select(-Pedigree) %>% 
    bind_rows(blups %>% 
                mutate(Env = "Combined\nBLUPs") %>% 
                rename(Value = BLUP)) %>% 
    mutate(Trait = case_when(
      Trait == 'UpperLeafAngle' ~ 'UpperLeafAngle(°)',
      Trait == 'MiddleLeafAngle' ~ 'MiddleLeafAngle(°)',
      Trait == 'PlantHeight' ~ 'PlantHeight(cm)',
      Trait == 'EarHeight' ~ 'EarHeight(cm)'),
      Trait = factor(Trait,
                     levels = c("UpperLeafAngle(°)","MiddleLeafAngle(°)",
                                "PlantHeight(cm)","EarHeight(cm)")),
      Env = factor(Env,
                   levels = c("IA15","IA16","IA17","Combined\nBLUPs"))) %>% 
    ggplot(aes(Value)) +
    geom_histogram() +
    xlab("Phenotypic Value") +
    facet_grid(Env ~ Trait,
               scales = "free") +
    theme_bw() +
    plot_theme)

# ggsave(plot = emp_hist,
#        filename = "Phenotype_Data/Empirical_Validation/Supplemental_Figure_S3.png",
#        width = 3.5,
#        height = 2.65,
#        dpi = 1200)

# ggsave(plot = emp_hist,
#        filename = "Final_Figures/Supplemental_Figure_S3.pdf",
#        width = 80,
#        height = 65,
#        units = 'mm',
#        dpi = 1200)

#### Calculate harmonic means, variance components, and heritability ####
var_comp <- emp_raw %>%
  group_by(Trait) %>%
  do(as_tibble(mixed_model(data = .,
                           return_type = "var",
                           genotype_term = "random"))) %>%
  ungroup()

#### Calculate harmonic means and output to the screen
(harmonic_mean <- emp_raw %>%
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

#### Output the heritability
# write_tsv("Phenotype_Data/Empirical_Validation/Empirical_Heritabilities.txt",
#           x = herit)