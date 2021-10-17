library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)

#### Output a cross-validation figure for admixture ####
#### Input CV data from the admixture analysis
cv <- read_tsv("Population_Structure/Admixture/CV_Output.txt") %>% 
  mutate(Delta = lag(CV_Error) - CV_Error)

#### Output for Fig S1 ####
#### Setting up the theme for the plots
plot_theme <- theme(aspect.ratio = .7,
                    axis.title = element_text(size = 7.5,
                                              face = "bold"),
                    axis.text = element_text(size = 6,
                                             face = "bold",
                                             color = "black"),
                    plot.title = element_text(size = 6, hjust = 0,
                                              face = "bold"),
                    plot.margin = unit(c(.1,0,.2,0),"cm"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())

(cv_plot1 <- cv %>% 
    ggplot(aes(x = K, y = CV_Error)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    plot_theme +
    scale_x_continuous(breaks = c(1:10),
                       name = "Number of subpopulations") +
    ylab("Cross-validation error") +
    ggtitle("Admixture cross-validation error \nfrom reference population"))

(cv_plot2 <- cv %>% 
    filter(!is.na(Delta)) %>% 
    ggplot(aes(x = K, y = Delta)) +
    geom_point() +
    geom_line() + 
    theme_bw() +
    plot_theme +
    scale_x_continuous(breaks = c(2:10),
                       name = "Number of subpopulations") +
    ylab(expression(bold(Delta~Cross-validation~error))) +
    ggtitle("Admixture change in cross-validation error \nfrom reference population"))

# ggsave(plot = cv_plot1 + cv_plot2,
#        filename = "Final_Figures/Supplemental_Figure_S1.png",
#        width = 7,
#        height = 3,
#        dpi = 1200)
# 
# ggsave(plot = cv_plot1 + cv_plot2,
#        filename = "Final_Figures/Supplemental_Figure_S1.pdf",
#        width = 180,
#        height = 76.2,
#        units = 'mm',
#        dpi = 1200)

#### Displaying admixture results for genotypes ####
#### Input data
#### This reads in the genotypic ID order for the results of the admixture file
gen_ids <- read_tsv("Population_Structure/Admixture/all_merged_imputed_pca.nosex",
                    col_names = FALSE) %>% 
  select(-'X2') %>% 
  rename(Gen = 'X1')

#### Reads in the population assignment
pop_str <- read_tsv("Population_Structure/Admixture/all_merged_imputed.5.Q",
                    col_names = FALSE) %>%
              separate(col = 'X1',
                       sep = " ",
                       into = paste("Pop",1:5,sep = "")) %>% 
              bind_cols(gen_ids) %>%
              select(Gen,everything()) %>% 
              mutate(K = 5)

#### Summarizes and formats the results
#### Pulls the group with the max amount of population
pop_str_results <- pop_str %>% 
  select(Gen,K, everything()) %>% 
  gather(Pop,Admixture,-Gen,-K) %>% 
  na.omit() %>% 
  arrange(K,Gen) %>% 
  group_by(K,Gen) %>% 
  filter(Admixture == max(Admixture)) %>% 
  mutate(Final_K = ifelse(Admixture >= .7,Pop,"Mixed")) %>% 
  select(-Pop,-Admixture) %>% 
  spread(K,Final_K,sep = "_")

#### Imports the PCA results and formats the data
all_pca <- read_tsv("Population_Structure/PCA/all_merged_imputed_pca.eigenvec",
                    col_names = FALSE) %>%
  separate(col = 'X1',
           into = c("Gen","Gen2",paste("PC",1:10,sep = "")),
           sep = " ") %>%
  select(-Gen2) %>% 
  mutate_at(vars(starts_with("PC")),as.numeric)

#### Import the genotypes ###
#### Imports the different GEs from respective populations, training, validation
training <- read_tsv("Population_Structure/Admixture/Maize_Association_Panel.txt") %>% 
  mutate(Type = "Training")

validation <- read_tsv("Population_Structure/Admixture/Validation_Panel.txt") %>% 
  mutate(Type = "Validation")

#### Formats and combines the ges with the PCA info
#### Assigns the different pops to names based on previous knowledge
all_pca_format <- all_pca %>% 
  left_join(pop_str_results %>% 
              select(Gen,K_5),
            by = 'Gen') %>%
  left_join(bind_rows(training,validation),
            by = 'Gen') %>% 
  mutate(K_5 = factor(case_when(
    K_5 == "Pop1" ~ "Sweet Corn",
    K_5 == "Pop2" ~ "Stiff Stalk",
    K_5 == "Pop3" ~ "Tropical",
    K_5 == "Pop4" ~ "Popcorn",
    K_5 == "Pop5" ~ "Non-stiff Stalk",
    TRUE ~ as.character(K_5))) %>% 
      fct_relevel("Stiff Stalk","Non-stiff Stalk",
                  "Tropical","Popcorn",
                  "Sweet Corn","Mixed"))

#### Output Figure 2 ####
#### Setting up and formatting the data for the plots
max_pc1 <- max(all_pca_format$PC1)
min_pc2 <- min(all_pca_format$PC2)
pop_summary_all <- all_pca_format %>% 
  rename(Subpop = K_5) %>% 
  group_by(Subpop) %>% 
  summarise(n = n()) %>% 
  mutate(Y = case_when(
    Subpop == "Stiff Stalk"  ~ min_pc2 + .06,
    Subpop == "Non-stiff Stalk" ~ min_pc2 + .05,
    Subpop == "Tropical" ~ min_pc2 + .04,
    Subpop == "Popcorn" ~ min_pc2 + .03,
    Subpop == "Sweet Corn" ~ min_pc2 + .02,
    Subpop == "Mixed" ~ min_pc2 + .01),
    X = max_pc1 - .045)

#### This is just setting up the labels for the plots and putting them on the plot
pop_summary_type <- all_pca_format %>% 
  rename(Subpop = K_5) %>%
  filter(Type %in% c("Validation","Training")) %>% 
  group_by(Type,Subpop) %>% 
  summarise(n = n()) %>% 
  mutate(Y = case_when(
    Subpop == "Stiff Stalk"  ~ min_pc2 + .06,
    Subpop == "Non-stiff Stalk" ~ min_pc2 + .05,
    Subpop == "Tropical" ~ min_pc2 + .04,
    Subpop == "Popcorn" ~ min_pc2 + .03,
    Subpop == "Sweet Corn" ~ min_pc2 + .02,
    Subpop == "Mixed" ~ min_pc2 + .01),
    X = max_pc1 - .045)

plot_theme <- theme(aspect.ratio = .7,
                    legend.text = element_text(size = 5,
                                               face = "bold"),
                    legend.key.size = unit(.5,"lines"),
                    legend.title = element_blank(),
                    legend.margin = margin(0,0,0,-7),
                    axis.title = element_text(size = 5,
                                              face = "bold"),
                    axis.text = element_text(size = 5,
                                             face = "bold",
                                             color = "black"),
                    plot.title = element_text(size = 6, hjust = 0,
                                              face = "bold"),
                    plot.subtitle = element_text(size = 4),
                    plot.margin = unit(c(0,0,0,0),"cm"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"))

(plot1 <-  all_pca_format %>% 
    ggplot(aes(x = PC1,
               y = PC2)) +
    geom_point(aes(color = K_5),
               size = .15) +
    scale_color_discrete(guide = 'none') +
    theme(legend.title = element_blank()) +
    ggtitle("Reference Population") +
    geom_text(data = pop_summary_all,
              aes(x = X,
                  y = Y,
                  label = paste(Subpop," (",n,")",sep = ""),
                  color = Subpop),
              fontface = "bold",
              size = 1.75,
              hjust = 0) +
    plot_theme)

(plot2 <- all_pca_format %>% 
    filter(Type == "Training") %>% 
    ggplot(aes(x = PC1,
               y = PC2)) +
    geom_point(aes(color = K_5),
               size = .15) +
    scale_color_discrete(guide = 'none') +
    theme(legend.title = element_blank()) +
    ggtitle("Training Population") +
    geom_text(data = pop_summary_type %>% 
                filter(Type == "Training"),
              aes(x = X,
                  y = Y,
                  label = paste(Subpop," (",n,")",sep = ""),
                  color = Subpop),
              fontface = "bold",
              size = 1.75,
              hjust = 0) +
    plot_theme)

(plot3 <- all_pca_format %>% 
    filter(Type == "Validation") %>% 
    ggplot(aes(x = PC1,
               y = PC2)) +
    geom_point(aes(color = K_5),
               size = .15) +
    theme(legend.title = element_blank()) +
    ggtitle(label = "Small Empirical Validation Population") +
    geom_text(data = pop_summary_type %>% 
                filter(Type == "Validation"),
              aes(x = X,
                  y = Y,
                  label = paste(Subpop," (",n,")",sep = ""),
                  color = Subpop),
              fontface = "bold",
              size = 1.75,
              hjust = 0) +
    scale_color_discrete(guide = 'none') +
    theme_bw() +
    plot_theme)

#### Using patchwork to put the plots together
pca_plot_out <- plot1 + plot2 + plot3


# ggsave(plot = pca_plot_out,
#        filename = "Final_Figures/Figure_2.png",
#        width = 7,
#        height = 1.9,
#        dpi = 1200)

# ggsave(plot = pca_plot_out,
#        filename = "Final_Figures/Figure_2.pdf",
#        width = 180,
#        height = 48.26,
#        units = 'mm',
#        dpi = 1200)

pop5 <- pop_str_results %>% 
  select(Gen,K_5)

#### Output data for sup table 1 ####
# write_tsv(x = pop5,
#           path = "Population_Structure/Admixture/K_5_results.txt")
