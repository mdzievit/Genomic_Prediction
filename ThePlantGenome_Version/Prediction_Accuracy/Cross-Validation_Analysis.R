library(tidyverse)
library(broom)

#### Load input data for training ####
## This is purely for plotting purposes, it keeps the order of the traits I want them to appear in the figures.
trait_order <- read_tsv("Trait_Order.txt")

## Loads in the training heritabilities that were previously calculated
h <- read_tsv("Phenotype_Data/Training/Training_Heritability.txt") %>% 
  left_join(trait_order)

## Loads in the cross validation results that were generated from the training genotypic and phenotypic data set
kfold_results <- read_tsv("Cross-Validation/Training-Cross_validation_results_kfold.txt") %>% 
  rename(Cor = 'V1')

#### Summarize, combine, and organize ####
## Calculates avg heritability, sd, regroups traits into their different trait types
kfold_sum <- kfold_results %>% 
  group_by(Trait,Fold) %>% 
  summarise(Avg_Cor = mean(Cor),
            SD = sqrt(var(Cor))) %>% 
  ungroup() %>% 
  arrange(Avg_Cor) %>% 
  left_join(h %>% 
              select(Trait,Heritability,Order),by = "Trait") %>% 
  mutate(Trait = factor(Trait) %>% 
           fct_reorder(Order), 
         Trait_Type = factor(Trait) %>% 
           fct_reorder(Order) %>% 
           fct_collapse('Tassel' = c("MainSpikeLength","SecondaryBranchNumber",
                                     "Spikelets-MainSpike","Spikelets-PrimaryBranch",
                                     "TasselBranchLength","TasselLength",
                                     "TasselPrimaryBranches"),
                        'Ear' = c("CobDiameter","CobWeight",
                                  "EarDiameter","EarLength",
                                  "EarRankNumber","EarRowNumber",
                                  "EarWeight","TotalKernelWeight",
                                  "SeedSetLength"),
                        'Kernel' = c("20KernelWeight","GerminationCount",
                                     "NIROil","NIRProtein",
                                     "NIRStarch","TotalKernelVolume"),
                        'Flowering' = c("DaystoSilk","GDDDaystoSilk",
                                        "DaysToTassel","GDDDaystoTassel",
                                        "GDDAnthesis-SilkingInterval"),
                        'Plant Architecture' = c("UpperLeafAngle","MiddleLeafAngle",
                                                 "PlantHeight","EarHeight",
                                                 "LeafWidth","LeafLength"),
                        'Other' = c("SouthernLeafBlight","StandCount",
                                    "TilleringIndex")),
         Trait_Type = factor(Trait_Type,
                             levels = c("Plant Architecture",
                                        "Flowering",
                                        "Tassel",
                                        "Ear",
                                        "Kernel",
                                        "Other")))


#### Output of Sup Table S4: Training Results ####

# write_tsv(path = "Final_Tables/Training_Supplemental_Table_S4.txt",
#           x = kfold_sum)

#### Output Sup Fig S5 ####
(kfold_plot_all <- kfold_sum %>% 
    mutate(Fold = factor(Fold),
           Trait_Type = factor(Trait_Type,
                               levels = c("Plant Architecture",
                                          "Tassel",
                                          "Flowering",
                                          "Ear",
                                          "Kernel",
                                          "Other")),
           Trait = fct_rev(Trait)) %>%
    ggplot(aes(x = Trait,
               color = Fold,
               shape = Fold,
               y = Avg_Cor)) +
    geom_jitter(size = 1,
                width = 0.15) +
    facet_grid(rows = vars(Trait_Type),
               scales = "free_y",
               space = 'free',
               switch = 'y') +
    theme_bw() +
    scale_y_continuous(breaks = seq(0.0,1.0,.1),
                       limits = c(0.0,1.0)) +
    theme(panel.grid.minor = element_blank(),
          #panel.grid.major.y = element_blank(),
          #aspect.ratio = .5,
          legend.text = element_text(size = 4, 
                                     face = 'bold'),
          legend.title = element_text(size = 4,
                                      face = 'bold'),
          legend.key.size = unit(.45,"lines"),
          legend.margin = margin(0,0,0,-7),
          axis.title.x = element_text(size = 7,
                                      face = "bold"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 5,
                                     face = "bold"),
          axis.text.x = element_text(size = 5,
                                     face = "bold"),
          strip.text.y = element_text(size = 5,
                                      margin = margin(0,.5,0,0,unit = 'mm'),
                                      face = "bold"),
          strip.background = element_rect(fill = 'white'),
          panel.spacing = unit(0,"lines"),
          panel.grid.major = element_line(color = 'lightgray',size = 0.15),
          #legend.direction = 'horizontal',
          legend.position = c(.95,.15),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line()) +
    ylab("Predictive Ability") +
    guides(colour = guide_legend(override.aes = list(size = 1))) +
    coord_flip())

# ggsave("Final_Figures/Supplemental_Figure_S5.pdf",
#        plot = kfold_plot_all,
#        dpi = 1200,
#        width = 80,units = 'mm')
# ggsave("Final_Figures/Supplemental_Figure_S5.png",
#        plot = kfold_plot_all,
#        dpi = 1200,
#        width = 80,units = 'mm')

#### Output Sup Fig S6 ####
(kfold_plot_her <- kfold_sum %>%
    filter(Fold == 10) %>%
    select(-Fold,-SD,-Order) %>% 
    gather(Type,Value,-Trait,-Trait_Type) %>%
    mutate(Type = factor(ifelse(Type == "Avg_Cor","Predictive Ability",ifelse(Type == 'Heritability','Sqrt Heritability',Type)),
                         levels = c("Predictive Ability","Sqrt Heritability")),
           Trait_Type = factor(Trait_Type,
                               levels = c("Plant Architecture",
                                          "Flowering",
                                          "Tassel",
                                          "Ear",
                                          "Kernel",
                                          "Other")),
           Trait = fct_rev(Trait),
           Value = sqrt(Value)) %>% 
    ggplot(aes(x = Trait,
               y = Value,
               group = Type,
               fill = Type)) +
    geom_bar(stat = "identity",
             position = "dodge") +
    facet_grid(rows = vars(Trait_Type),
               scales = "free_y",
               space = 'free',
               switch = 'y') +
    theme_bw() +
    scale_y_continuous(breaks = seq(0.0,1.0,.1),
                       limits = c(0.0,1.0)) +
    theme(legend.text = element_text(size = 4, 
                                     face = 'bold'),
          legend.title = element_text(size = 4,
                                      face = 'bold'),
          legend.key.size = unit(.45,"lines"),
          legend.margin = margin(0,0,0,-7),
          axis.title.x = element_text(size = 7,
                                      face = "bold"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 5,
                                     face = "bold"),
          axis.text.x = element_text(size = 5,
                                     face = "bold"),
          strip.text.y = element_text(size = 4,
                                      margin = margin(0,.5,0,0,unit = 'mm'),
                                      face = "bold"),
          strip.background = element_rect(fill = 'white'),
          panel.spacing = unit(0,"lines"),
          #legend.direction = 'horizontal',
          legend.position = c(.95,.15),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(size = 0.15),
          panel.grid.major.y = element_blank()) +
    ylab("Predictive Ability") +
    guides(colour = guide_legend(override.aes = list(size = 1))) +
    coord_flip() +
    scale_fill_manual(breaks = c("Sqrt Heritability","Predictive Ability"),
                      values = c("#999999","#D55E00"),
                      labels = c(expression(sqrt(h^2)),'Predictive\nAbility')))

# ggsave("Final_Figures/Figures/Supplemental_Figure_S6.png",
#        plot = kfold_plot_her,
#        dpi = 1200,
#        width = 80,
#        units = 'mm')
# ggsave("Final_Figures/Supplemental_Figure_S6.pdf",
#        plot = kfold_plot_her,
#        dpi = 1200,
#        width = 80,
#        units = 'mm')

#### Output Fig 3 ####
#### Plot with only 10 fold, and then SD bars
(kfold_plot_k10 <- kfold_sum %>% 
    filter(Fold == 10) %>% 
    mutate(Trait_Type = factor(Trait_Type,
                               levels = c("Plant Architecture",
                                          "Flowering",
                                          "Tassel",
                                          "Ear",
                                          "Kernel",
                                          "Other")),
           Trait = fct_rev(Trait)) %>% 
    ggplot(aes(x = Trait,
               y = Avg_Cor)) +
    geom_bar(stat = "identity",
             width = .75) +
    geom_errorbar(aes(ymin = Avg_Cor - SD,
                      ymax = Avg_Cor + SD),
                  width = .5,
                  color = "red",
                  size = .25) +
    geom_text(aes(label = round(Avg_Cor,2)),
              hjust = -.75,
              size = 1.5,
              fontface = "bold") +
    facet_grid(rows = vars(Trait_Type),
               scales = "free_y",
               space = 'free',
               switch = 'y') +
    theme_bw() +
    scale_y_continuous(breaks = seq(0.0,1.0,.1),
                       limits = c(0.0,1.0)) +
    theme(legend.text = element_text(size = 4, 
                                     face = 'bold'),
          legend.title = element_text(size = 4,
                                      face = 'bold'),
          legend.key.size = unit(.45,"lines"),
          legend.margin = margin(0,0,0,-7),
          axis.title.x = element_text(size = 7,
                                      face = "bold"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 5,
                                     face = "bold"),
          axis.text.x = element_text(size = 5,
                                     face = "bold"),
          strip.text.y = element_text(size = 4,
                                      margin = margin(0,.5,0,0,unit = 'mm'),
                                      face = "bold"),
          strip.background = element_rect(fill = 'white'),
          panel.spacing = unit(0,"lines"),
          #legend.direction = 'horizontal',
          legend.position = c(.95,.15),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(size = 0.15),
          panel.grid.major.y = element_blank()) +
    ylab("Prediction accuracy") +
    coord_flip())

# ggsave("Final_Figures/Figures/Figure_3.png",
#        plot = kfold_plot_k10,
#        dpi = 1200,
#        width = 80,
#        units = 'mm')
# 
# ggsave("Final_Figures/Figure_3.pdf",
#        plot = kfold_plot_k10,
#        dpi = 1200,
#        width = 80,
#        units = 'mm')

#### Processing combined cross-validation data ####
#### This is the emprirical and training data combined as the 'training' population
#### Input the cross-validation results
large_kfold <- read_tsv("Cross-Validation/Combined_Cross-Validation_All.txt") %>% 
  rename(Cor = 'V1')

#### Input the heritabilities for the combined data
large_herit <- read_tsv("Phenotype_Data/Combined/Combined_Heritabilities.txt")

#### Processing the data
large_kfold_sum <- large_kfold %>% 
  group_by(Trait,Fold) %>% 
  summarise(Avg_Cor = mean(Cor),
            SD = sqrt(var(Cor))) %>% 
  ungroup() %>% 
  arrange(Avg_Cor) %>% 
  left_join(large_herit %>% 
              select(Trait,H),by = "Trait") %>% 
  rename(Heritability = H)

#### Output of Sup Table S4: Combined Results ####
# write_tsv(path = "Final_Tables/Combined_Supplemental_Table_S4.txt",
#           large_kfold_sum)


#### Processing the empirical cross-validation ####
#### Loading the kfold results
kfold_empir <- read_tsv("Cross-Validation/Empirical-Cross_validation_results_kfold.txt") %>% 
  rename(Cor = 'V1')

#### Loading the heritabilities of the empirical data
emp_herit <- read_tsv("Phenotype_Data/Empirical_Validation/Empirical_Heritabilities.txt")

#### Processing the data
kfold_empir_sum <- kfold_empir %>% 
  group_by(Trait,Fold) %>% 
  summarise(Avg_Cor = mean(Cor),
            SD = sqrt(var(Cor))) %>% 
  ungroup() %>% 
  arrange(Avg_Cor) %>% 
  left_join(emp_herit %>%
              select(Trait,H) %>% 
              mutate(Trait = case_when(
                Trait == "ULA" ~ "UpperLeafAngle",
                Trait == "MLA" ~ "MiddleLeafAngle",
                TRUE ~ Trait)),
            by = "Trait") %>% 
  rename(Heritability = H)

#### Output of Sup Table S4: Empirical Results ####
# write_tsv(path = "Final_Tables/Empirical_Supplemental_Table_S4.txt",
#           kfold_empir_sum)

#### Setting up a common theme for the plots used in the following code
plot_theme <- theme(aspect.ratio = .65,
                    legend.text = element_text(size = 8.5,
                                               face = "bold"),
                    legend.key.size = unit(0.25,'cm'),
                    legend.title = element_text(size = 9.5,
                                                face = "bold"),
                    legend.margin = margin(0,1,0,-10),
                    axis.title = element_text(size = 10,
                                              face = "bold"),
                    axis.text = element_text(size = 9.5,
                                             face = "bold",
                                             color = "black"),
                    plot.margin = unit(c(.05,0,0,0),"cm"),
                    panel.grid = element_blank())

#### Output of Sup Fig S8 ####
(large_kfold_sel <- kfold_sum %>%
    select(-Trait_Type) %>% 
    filter(str_detect(Trait,c("Angle|Height"))) %>% 
    mutate(Training_Size = "Maize\nAssociation") %>% 
    bind_rows(large_kfold_sum %>% 
                mutate(Training_Size = "Maize\nAssociation +\nSmall\nEmpirical")) %>%
    bind_rows(kfold_empir_sum %>% 
                mutate(Training_Size = "Small\nEmpirical")) %>% 
    filter(Fold == 10) %>%
    mutate(Trait = factor(Trait) %>% 
             fct_relevel(c("UpperLeafAngle","MiddleLeafAngle",
                           "PlantHeight","EarHeight")) %>% 
             fct_rev(),
           Training_Size = factor(Training_Size) %>% 
             fct_relevel(c("Maize\nAssociation",
                           "Small\nEmpirical",
                           "Maize\nAssociation +\nSmall\nEmpirical"))) %>% 
    ggplot(aes(x = Trait,
               y = Avg_Cor,
               fill = Training_Size,
               group = -as.numeric(Training_Size))) +
    geom_bar(stat = "identity",
             position = "dodge",
             width = .75) +
    geom_errorbar(aes(ymin = Avg_Cor - SD,
                      ymax = Avg_Cor + SD),
                  color = "black",
                  width = .5,
                  size = .4,
                  position = position_dodge(.75)) +
    geom_text(aes(label = format(round(Avg_Cor,3),nsmall = 3)),
              position = position_dodge(width = 0.75),
              hjust = -.5,
              size = 3.5,
              fontface = "bold") +
    ylab("Predictive Ability") +
    theme_bw() +
    scale_y_continuous(breaks = seq(0.0,1.0,.2),
                       limits = c(0.0,1.0)) +
    guides(fill = guide_legend(title = "Training\nPopulation\nSize")) +
    coord_flip() +
    plot_theme +
    theme(axis.title.y = element_blank()))

# ggsave("Final_Figures/Supplemental_Figure_S8.png",
#        plot = large_kfold_sel,
#        dpi = 1200,
#        width = 7)

# ggsave("Final_Figures/Supplemental_Figure_S8.pdf",
#        plot = large_kfold_sel,
#        dpi = 1200,
#        width = 180,
#        units = 'mm')

#### Cor test of sqrt hert w/ pred acc ####
#### Correlation test to compare relationship of sqrt of heritability with prediction accuracy
cor_test <- kfold_sum %>% 
  filter(Fold == 10) %>% 
  do(tidy(cor.test(.$Avg_Cor,sqrt(.$Heritability)))) %>% 
  mutate(X = 0,
         Y = 1,
         Sign = case_when(
           p.value <= 0.01 ~ "b",
           p.value <= 0.05 ~ "a",
           TRUE ~ NA_character_))

#### Formatting the data and adding physical spots for the text
cor_test_all <- kfold_sum %>% 
  filter(Fold == 10) %>% 
  group_by(Trait_Type) %>% 
  do(tidy(cor.test(.$Avg_Cor,sqrt(.$Heritability)))) %>%
  ungroup() %>% 
  mutate(X = 0,
         Y = case_when(
           Trait_Type == "Plant Architecture Traits" ~ .93,
           Trait_Type == "Flowering Traits" ~ 0.86,
           Trait_Type == "Tassel Traits" ~ 0.79,
           Trait_Type == "Ear Traits" ~ 0.72,
           Trait_Type == "Kernel Traits" ~ 0.65,
           Trait_Type == "Other Traits" ~ 0.58),
         Trait_Type2 = case_when(
           Trait_Type == "Plant Architecture Traits" ~ "Plant~Architecture~Traits",
           Trait_Type == "Flowering Traits" ~ "Flowering~Traits",
           Trait_Type == "Tassel Traits" ~ "Tassel~Traits",
           Trait_Type == "Ear Traits" ~ "Ear~Traits",
           Trait_Type == "Kernel Traits" ~ "Kernel~Traits",
           Trait_Type == "Other Traits" ~ "Other~Traits"))

#### Theme for the plot
plot_theme <- theme(aspect.ratio = .65,
                    legend.text = element_text(size = 9,
                                               face = "bold"),
                    legend.title = element_text(size = 10,
                                                face = "bold"),
                    legend.margin = margin(0,0.15,0,-10),
                    legend.key.size = unit(.75,"lines"),
                    axis.title = element_text(size = 12,
                                              face = "bold"),
                    axis.text = element_text(size = 11,
                                             face = "bold",
                                             color = "black"),
                    plot.margin = unit(c(0.15,0,0,0.15),"cm"),
                    panel.grid = element_blank())

#### Output of Sup Fig S7 ####
(cor_her <- kfold_sum %>%
    filter(Fold == 10) %>% 
    ggplot(aes(x = sqrt(Heritability),
               y = Avg_Cor)) +
    geom_point(aes(color = Trait_Type,shape = Trait_Type),
               size = 2) +
    scale_y_continuous(breaks = seq(0.0,1.0,.2),
                       limits = c(0.0,1.0)) +
    scale_x_continuous(breaks = seq(0.0,1.0,.2),
                       limits = c(0.0,1.0)) +
    xlab(expression(bold(sqrt(h^2)))) +
    ylab("Predictive Ability") +
    guides(color = guide_legend(title = "Trait Type"),
           shape = guide_legend(title = "Trait Type")) +
    geom_text(data = cor_test,
              aes(x = X,
                  y = Y,
                  label = paste("bold(",
                                paste("Overall~",'italic(r)',"~'='~",format(round(estimate,3),nsmall = 3),
                                      sep = ""),
                                ")",
                                sep = "")),
              color = "black",
              size = 4,
              hjust = 0,
              parse = TRUE,
              show.legend = FALSE) +
    theme_bw() +
    plot_theme)


# ggsave("Final_Figures/Supplemental_Figure_S7.png",
#        plot = cor_her,
#        dpi = 1200,
#        width = 7)

# ggsave("Final_Figures/Supplemental_Figure_S7.pdf",
#        plot = cor_her,
#        dpi = 1200,
#        width = 180,
#        units = 'mm')