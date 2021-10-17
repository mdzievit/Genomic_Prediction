library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

#### Excluded GEs ####
#### GEs that are also in the training population under a different name
excl_ges <- read_tsv("excluded_ges.txt")

gen_id <- read_tsv("RRBLUP_Input_Files/Validation_Imputed.012.indv",
                   col_names = FALSE) %>% 
  rename(Source = 'X1') %>% 
  mutate(Gen = row_number())

U_values <- read_tsv("U_Calculations/U_Results_Valid.txt") %>% 
  rename(Gen = Indv,
         U_Value = Result) %>% 
  left_join(gen_id,
            by = "Gen")
  
subpops <- read_tsv("Population_Structure/Admixture/K_5_results.txt") %>%
  rename(Source = Gen,
         Subpop = K_5) %>% 
  mutate(Subpop = factor(case_when(
    Subpop == "Pop1" ~ "Sweet Corn",
    Subpop == "Pop2" ~ "Stiff Stalk",
    Subpop == "Pop3" ~ "Tropical",
    Subpop == "Pop4" ~ "Popcorn",
    Subpop == "Pop5" ~ "Non-stiff Stalk",
    TRUE ~ as.character(Subpop))) %>% 
      fct_relevel("Stiff Stalk","Non-stiff Stalk",
                  "Tropical","Popcorn",
                  "Sweet Corn","Mixed"))


plot_theme <- theme(axis.text = element_text(size = 8,
                                             face = "bold",
                                             color = "black"),
                    axis.title = element_text(size = 8,face = "bold"),
                    plot.title = element_text(size = 7, hjust = 0.5,
                                              face = "bold"),
                    strip.text.x = element_text(size = 7,
                                              margin = margin(0,0,0,0, "cm"),
                                              face = "bold",
                                              vjust = .5),
                    strip.background = element_blank(),
                    panel.spacing.x = unit(.5,"lines"),
                    panel.spacing.y = unit(0,"lines"),
                    plot.margin = unit(c(0,0,0,0),"cm"),
                    aspect.ratio = .425,
                    panel.grid = element_blank())
(uvalue_hist <- U_values %>%
    left_join(subpops) %>%
    bind_rows(U_values %>% mutate(Subpop = "All of AmesDP")) %>% 
    mutate(Subpop = factor(Subpop) %>% 
             fct_relevel("All of AmesDP","Stiff Stalk","Non-stiff Stalk",
                                      "Tropical","Popcorn",
                                      "Sweet Corn","Mixed")) %>% 
    ggplot(aes(U_Value)) +
    geom_histogram(binwidth = .015) +
    facet_wrap(~Subpop,
               scales = "free_y",
               ncol = 2) +
    #scale_x_continuous(limits = c(0.5,1)) +
    theme_bw() +
    plot_theme)


# ggsave(plot = uvalue_hist,
#        filename = "Final_Figures/Supplemental_Figure_S2.png",
#        dpi = 1200,
#        width = 7)

# ggsave(plot = uvalue_hist,
#        filename = "Final_Figures/Supplemental_Figure_S2.pdf",
#        dpi = 1200,
#        width = 180,
#        height = 146.05,
#        units = 'mm')

trait_order <- read_tsv("Trait_Order.txt")

pred_results <- read_tsv("Predict_Phenotypes/Predicted_BLUPs.txt") %>%
  left_join(gen_id,
            by = "Gen") %>%
  select(-Gen) %>% 
  filter(!(Source %in% excl_ges$Genotype)) %>% 
  gather(Trait,Pred,-Source) %>% 
  left_join(trait_order) %>% 
  mutate(Trait = factor(Trait) %>% 
           fct_reorder(Order), 
         Trait_Type = factor(Trait) %>% 
           fct_reorder(Order) %>%
           fct_collapse('Tassel Traits' = c("MainSpikeLength","SecondaryBranchNumber",
                                            "Spikelets-MainSpike","Spikelets-PrimaryBranch",
                                            "TasselBranchLength","TasselLength",
                                            "TasselPrimaryBranches"),
                        'Ear Traits' = c("CobDiameter","CobWeight",
                                         "EarDiameter","EarLength",
                                         "EarRankNumber","EarRowNumber",
                                         "EarWeight","TotalKernelWeight",
                                         "SeedSetLength"),
                        'Kernel Traits' = c("20KernelWeight","GerminationCount",
                                            "NIROil","NIRProtein",
                                            "NIRStarch","TotalKernelVolume"),
                        'Flowering Traits' = c("DaystoSilk","GDDDaystoSilk",
                                               "DaysToTassel","GDDDaystoTassel",
                                               "GDDAnthesis-SilkingInterval"),
                        'Plant Architecture Traits' = c("UpperLeafAngle","MiddleLeafAngle",
                                                        "PlantHeight","EarHeight",
                                                        "LeafWidth","LeafLength"),
                        'Other Traits' = c("SouthernLeafBlight","StandCount",
                                           "TilleringIndex")),
         Trait_Type = factor(Trait_Type,
                             levels = c("Plant Architecture Traits",
                                        "Flowering Traits",
                                        "Tassel Traits",
                                        "Ear Traits",
                                        "Kernel Traits",
                                        "Other Traits"))) %>% 
  left_join(U_values,
            by = "Source") %>% 
  left_join(subpops,
            by = "Source")

plot_type <- function(data,
                      plot_theme,
                      ncols){
  
  plot <- data %>%
    ggplot(aes(x = Pred,
               y = U_Value)) +
    geom_point(aes(col = Subpop),
               size = .15,
               alpha = 0.25) +
    facet_wrap(Trait_Type ~ Trait,
               ncol = ncols,
               scales = "free_x") +
    theme_bw() +
    scale_y_continuous(limits = c(0.4,1.0)) +
    #ggtitle(trait_type) +
    ylab("Upper bound for reliability") +
    xlab("Predicted") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 1,
                                                    alpha = 1))) +
    plot_theme
  return(plot)
}

  
plot_theme <- theme(legend.margin = margin(0,0,0,0),
                    legend.box.margin = margin(0,-2.5,0,-5),
                    legend.text = element_text(size = 5,
                                               face = "bold"),
                    legend.title = element_blank(),
                    axis.title = element_text(size = 5,
                                              face = "bold"),
                    axis.text = element_text(size = 5,
                                             face = "bold",
                                             color = "black"),
                    strip.text.x = element_text(size = 4,
                                                margin = margin(c(0,0,0,0)),
                                                face = "bold"),
                    strip.background = element_blank(),
                    panel.spacing = unit(.15,"lines"),
                    aspect.ratio = 1.1)

(u_value_plot <-  plot_type(data = pred_results,
                   plot_theme = plot_theme,
                   ncols = 6))

# ggsave(plot = u_value_plot,
#        filename = "Final_Figures/Supplemental_Figure_S11.png",
#        dpi = 1200,
#        width = 7,
#        height = 8)

# ggsave(plot = u_value_plot,
#        filename = "Final_Figures/Supplemental_Figure_S11.pdf",
#        dpi = 1200,
#        width = 180,
#        height = 203.2,
#        units = 'mm')

###Generate summary statistics for the U-Values#####

U_values_subpop <- U_values %>%
  left_join(subpops) %>%
  bind_rows(U_values %>% mutate(Subpop = "All of AmesDP"))

(U_values_sum <- U_values_subpop %>% 
  group_by(Subpop) %>% 
  summarise(Avg = mean(U_Value),
            max = max(U_Value),
            min = min(U_Value)))

write_tsv(path = "Final_Tables/Supplemental_Table_S2.txt",
          U_values_sum)
