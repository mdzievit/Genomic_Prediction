library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

#### Excluded GEs ####
#### GEs that are also in the training population under a different name
excl_ges <- read_tsv("excluded_ges.txt")

#### Loading input data ####
#### Validation file and the GEs. Aligns with the genotypic data
gen_id <- read_tsv("RRBLUP_Input_Files/Validation_Imputed.012.indv",
                   col_names = FALSE) %>% 
  rename(Source = 'X1') %>% 
  mutate(Gen = row_number())

#### Predicted GEBVs
pred_results <- read_tsv("Predict_Phenotypes/Predicted_BLUPs.txt") %>% 
  left_join(gen_id,
            by = "Gen") %>% 
  filter(!(Source %in% excl_ges$Genotype))

#### U-Values
U_values <- read_tsv("U_Calculations/U_Results_Valid.txt") %>% 
  rename(Gen = Indv,
         U_Value = Result)

#### Subpopulation information for each GE
#### Also clasifies each sub pop into a type based on prior knowledge,
#### see text of manuscript for more information
subpops <- read_tsv("Population_Structure/Admixture/K_5_results.txt") %>%
  rename(Subpop = K_5,
         Source = Gen) %>% 
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

#### LA prediction accuracy #####
#### Loading up the smalle empirical phenotypes
emp_phenotypes <- read_tsv("Phenotype_Data/Empirical_Validation/Empirical_Phenotypes_BLUPs.txt") %>% 
  rename(Source = Genotype)

#### Running the correlations between predicted and observed, 
#### then setting up the text for plotting purposes (x and y coordinates)
comp <- emp_phenotypes %>%
  left_join(pred_results %>% 
              gather(Trait,Pred,-Source),
            by = c("Trait","Source")) %>% 
  left_join(subpops)

#### Overall correlations that will print to console
#### Also sets up text coordinates for plotting
(cor <- comp %>% 
    group_by(Trait) %>% 
    summarise(r = cor(BLUP,Pred,
                      use = 'complete.obs',
                      method = 'pearson')) %>% 
    mutate(X2 = case_when(
      Trait %in% c("UpperLeafAngle","MiddleLeafAngle") ~ 82.5,
      Trait == "PlantHeight" ~ 250,
      Trait == "EarHeight" ~ 150),
      Y2 = case_when(
        Trait %in% c("UpperLeafAngle","MiddleLeafAngle") ~ 90,
        Trait == "PlantHeight" ~ 275,
        Trait == "EarHeight" ~ 170)))

cor_grps <- comp %>% 
  group_by(Trait,Subpop) %>% 
  summarise(r = cor(BLUP,Pred,
                    use = 'complete.obs',
                    method = 'pearson'),
            n = n(),
            .groups = 'drop') %>% 
  ungroup() %>% 
  mutate(X = case_when(
    Trait == "UpperLeafAngle" ~ 15,
    Trait == "MiddleLeafAngle" ~ 30,
    Trait == "PlantHeight" ~ 50,
    Trait == "EarHeight" ~ 0),
    Y = case_when(
      Subpop == "Stiff Stalk" & Trait == "UpperLeafAngle"  ~ 90,
      Subpop == "Non-stiff Stalk" & Trait == "UpperLeafAngle" ~ 85.42,
      Subpop == "Tropical" & Trait == "UpperLeafAngle" ~ 80.85,
      Subpop == "Popcorn" & Trait == "UpperLeafAngle" ~ 76.28,
      Subpop == "Sweet Corn" & Trait == "UpperLeafAngle" ~ 71.71,
      Subpop == "Mixed" & Trait == "UpperLeafAngle" ~ 67.14,
      Subpop == "Stiff Stalk" & Trait == "MiddleLeafAngle"  ~ 90,
      Subpop == "Non-stiff Stalk" & Trait == "MiddleLeafAngle" ~ 86.57,
      Subpop == "Tropical" & Trait == "MiddleLeafAngle" ~ 83.14,
      Subpop == "Popcorn" & Trait == "MiddleLeafAngle" ~ 79.71,
      Subpop == "Sweet Corn" & Trait == "MiddleLeafAngle" ~ 76.28,
      Subpop == "Mixed" & Trait == "MiddleLeafAngle" ~ 72.85,
      Subpop == "Stiff Stalk" & Trait == "PlantHeight" ~ 275,
      Subpop == "Non-stiff Stalk" & Trait == "PlantHeight" ~ 262.14,
      Subpop == "Tropical" & Trait == "PlantHeight" ~ 249.28,
      Subpop == "Popcorn" & Trait == "PlantHeight" ~ 236.42,
      Subpop == "Sweet Corn" & Trait == "PlantHeight" ~ 223.57,
      Subpop == "Mixed" & Trait == "PlantHeight" ~ 210.71,
      Subpop == "Stiff Stalk" & Trait == "EarHeight" ~ 175,
      Subpop == "Non-stiff Stalk" & Trait == "EarHeight" ~ 165,
      Subpop == "Tropical" & Trait == "EarHeight" ~ 155,
      Subpop == "Popcorn" & Trait == "EarHeight" ~ 145,
      Subpop == "Sweet Corn" & Trait == "EarHeight" ~ 135,
      Subpop == "Mixed" & Trait == "EarHeight" ~ 125),
    Subpop2 = str_replace(Subpop," ","~"))

#### Function to do some basic rounding ####
mround <- function(x,base){ 
  base*round(x/base) 
} 

#### Function for plotting figures Num 1 ####
#### It puts together the multi-panel figures
plotter1 <- function(data,trait,title, xlim, ylim, incr, corData_Grps, corData,
                     theme) {
  plot <- data %>%
    filter(Trait == trait) %>% 
    ggplot(aes(x = Pred,
               y = BLUP)) +
    geom_point(size = .6,
               aes(color = Subpop),
               show.legend = FALSE) +
    geom_text(data = corData_Grps %>% 
                filter(Trait == trait),
              aes(x = X,
                  y = Y,
                  label = paste("bold(",
                                paste(Subpop2, "~(",n,")~","italic(r)","~'='~",sprintf('"%1.2f"',r),
                                      sep = ""),
                                ")",
                                sep = ""),
                  color = Subpop),
              parse = TRUE,
              size = 2.2,
              hjust = 0,
              show.legend = FALSE) +
    theme_bw() +
    scale_x_continuous(breaks = seq(mround(xlim[1],5),mround(xlim[2],5),incr),
                       limits = xlim,
                       name = "Predicted") +
    scale_y_continuous(breaks = seq(mround(ylim[1],5),mround(ylim[2],5),incr),
                       limits = ylim,
                       name = "Observed") +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(.85,.3)) +
    plot_theme +
    geom_text(data = corData %>% 
                filter(Trait == trait),
              aes(x = X2,
                  y = Y2,
                  label = paste("bold(",
                                paste("Overall~italic(r)~'='~",sprintf('"%1.2f"',r),
                                      sep = ""),
                                ")",
                                sep = "")),
              parse = TRUE,
              col = 'black',
              size = 2.2) +
    guides(colour = guide_legend(override.aes = list(size = 1),
                                 title = "Subpopulation"))
  return(plot)
}

#### Sets up the plots theme
plot_theme <- theme(aspect.ratio = .7,
                    legend.text = element_text(size = 5),
                    legend.title = element_blank(),
                    legend.margin = margin(0,0,0,-7),
                    axis.title = element_text(size = 9,
                                              face = "bold"),
                    axis.text = element_text(size = 8,
                                             face = "bold",
                                             color = "black"),
                    plot.title = element_text(size = 10, hjust = 0,
                                              face = "bold"),
                    plot.margin = unit(c(0.1,.1,0,.1),"cm"),
                    legend.key.size = unit(.75,"lines"))

#### Output Figure 4 ####
(plot_ULA <- plotter1(data = comp,
                      trait = "UpperLeafAngle",
                      title = "UpperLeafAngle",
                      xlim = c(15,90),
                      ylim = c(15,90),
                      incr = 15,
                      corData_Grps = cor_grps,
                      corData = cor,
                      theme = plot_theme))

(plot_MLA <- plotter1(data = comp,
                      trait = "MiddleLeafAngle",
                      title = "MiddleLeafAngle",
                      xlim = c(30,90),
                      ylim = c(30,90),
                      incr = 15,
                      corData_Grps = cor_grps,
                      corData = cor,
                      theme = plot_theme))

(plot_ph <- plotter1(data = comp,
                     trait = "PlantHeight",
                     title = "PlantHeight",
                     xlim = c(50,275),
                     ylim = c(50,275),
                     incr = 50,
                     corData_Grps = cor_grps,
                     corData = cor,
                     theme = plot_theme))

(plot_eh <- plotter1(data = comp,
                     trait = "EarHeight",
                     title = "EarHeight",
                     xlim = c(0,175),
                     ylim = c(0,175),
                     incr = 25,
                     corData_Grps = cor_grps,
                     corData = cor,
                     theme = plot_theme))

# ggsave(plot = grid.arrange((cbind(ggplotGrob(plot_ULA),ggplotGrob(plot_MLA),
#                                   size = "last")),
#                            (cbind(ggplotGrob(plot_ph),ggplotGrob(plot_eh),
#                                   size = "last"))),
#        filename = "Final_Figures/Figure_4.png",
#        width = 7,
#        height = 5,
#        dpi = 1200)

# ggsave(plot = grid.arrange((cbind(ggplotGrob(plot_ULA),ggplotGrob(plot_MLA),
#                                   size = "last")),
#                            (cbind(ggplotGrob(plot_ph),ggplotGrob(plot_eh),
#                                   size = "last"))),
#        filename = "Final_Figures/Figure_4.pdf",
#        width = 180,
#        height = 127,
#        units = 'mm',
#        dpi = 1200)

#### U_Plot comparisons ####
#### Formats the uvalue data with the combines gebvs + blups
#### Breaks them into groups based on top50/bot50/middle
U_values_format <- comp %>% 
  left_join(U_values %>% 
              left_join(gen_id)) %>%
  arrange(Trait,desc(U_Value)) %>% 
  group_by(Trait) %>% 
  mutate(Rank = row_number(),
         Group = case_when(
           Rank <=50 ~ "Top50",
           Rank >= 244 ~ "Bottom50",
           TRUE ~ "Middle"
         ),
         Group = factor(Group) %>% fct_relevel(c("Top50","Middle","Bottom50"))) %>% 
  ungroup()

#### Calcualtes correlations within these groups and
#### then adds in the coordinates
U_values_cor <- U_values_format %>% 
  group_by(Trait,Group) %>% 
  summarise(r = cor(BLUP,Pred)) %>% 
  ungroup() %>% 
  mutate(X = case_when(
    Trait == "UpperLeafAngle" & Group == "Top50" ~ 15,
    Trait == "MiddleLeafAngle" & Group == "Top50" ~ 30,
    Trait == "PlantHeight" & Group == "Top50" ~ 50,
    Trait == "EarHeight" & Group == "Top50" ~ 0,
    Trait == "UpperLeafAngle" & Group == "Bottom50" ~ 15,
    Trait == "MiddleLeafAngle" & Group == "Bottom50" ~ 30,
    Trait == "PlantHeight" & Group == "Bottom50" ~ 50,
    Trait == "EarHeight" & Group == "Bottom50" ~ 0,
    Trait == "UpperLeafAngle" & Group == "Middle" ~ 15,
    Trait == "MiddleLeafAngle" & Group == "Middle" ~ 30,
    Trait == "PlantHeight" & Group == "Middle" ~ 50,
    Trait == "EarHeight" & Group == "Middle" ~ 0),
    Y = case_when(
      Trait == "UpperLeafAngle" & Group == "Top50" ~ 90,
      Trait == "MiddleLeafAngle" & Group == "Top50" ~ 90,
      Trait == "PlantHeight" & Group == "Top50" ~ 275,
      Trait == "EarHeight" & Group == "Top50" ~ 160,
      Trait == "UpperLeafAngle" & Group == "Bottom50" ~ 90,
      Trait == "MiddleLeafAngle" & Group == "Bottom50" ~ 90,
      Trait == "PlantHeight" & Group == "Bottom50" ~ 275,
      Trait == "EarHeight" & Group == "Bottom50" ~ 160,
      Trait == "UpperLeafAngle" & Group == "Middle" ~ 90,
      Trait == "MiddleLeafAngle" & Group == "Middle" ~ 90,
      Trait == "PlantHeight" & Group == "Middle" ~ 275,
      Trait == "EarHeight" & Group == "Middle" ~ 160),
    X2 = case_when(
      Trait == "UpperLeafAngle" ~ 20,
      Trait == "MiddleLeafAngle" ~ 30,
      Trait == "PlantHeight" ~ 75,
      Trait == "EarHeight" ~ 25),
    Y2 = case_when(
      Trait == "UpperLeafAngle" & Group == "Top50" ~ 8,
      Trait == "UpperLeafAngle" & Group == "Bottom50" ~ 68,
      Trait == "UpperLeafAngle" & Group == "Middle" ~ 73,
      Trait == "MiddleLeafAngle" & Group == "Top50" ~ 80,
      Trait == "MiddleLeafAngle" & Group == "Bottom50" ~ 72,
      Trait == "MiddleLeafAngle" & Group == "Middle" ~ 76,
      Trait == "PlantHeight" & Group == "Top50" ~ 260,
      Trait == "PlantHeight" & Group == "Bottom50" ~ 230,
      Trait == "PlantHeight" & Group == "Middle" ~ 245,
      Trait == "EarHeight" & Group == "Top50" ~ 150,
      Trait == "EarHeight" & Group == "Bottom50" ~ 130,
      Trait == "EarHeight" & Group == "Middle" ~ 140),
    Trait = factor(Trait) %>% fct_relevel(c("UpperLeafAngle","MiddleLeafAngle",
                                            "PlantHeight","EarHeight")))
#### Plotter 2 for figures ####
plotter2 <- function(data,trait, ylim, xlim, strip_title,
                     u_cor,plot_theme, incr) {
  plot <- data %>%
    mutate(Trait = factor(Trait) %>% fct_relevel(c("UpperLeafAngle","MiddleLeafAngle",
                                                   "PlantHeight","EarHeight"))) %>%
    filter(str_detect(Trait,trait)) %>% 
    ggplot(aes(x = Pred,
               y = BLUP)) +
    geom_point(size = .5,aes(col = Subpop)) +
    facet_grid(Trait ~ Group, scales = "free") +
    geom_text(data = u_cor %>%
                filter(str_detect(Trait,trait)),
              aes(x = X,
                  y = Y,
                  label = paste("bold(",
                                paste("Overall~italic(r)","~'='~",sprintf('"%1.2f"',r),
                                      sep = ""),
                                ")",
                                sep = "")),
              parse = TRUE,
              size = 2.25,
              hjust = 0,
              show.legend = FALSE) +
    theme_bw() +
    scale_x_continuous(breaks = seq(xlim[1],xlim[2],incr),
                       limits = xlim,
                       name = "Predicted") +
    scale_y_continuous(breaks = seq(ylim[1],ylim[2],incr),
                       limits = ylim,
                       name = "Observed") +
    guides(color = guide_legend(override.aes = list(size = 1.25))) +
    plot_theme
  
  if (strip_title){
    plot <- plot +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank())
  } else {
    plot <- plot +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            strip.text.x = element_blank())
  }
}

#### Theme for the plots
plot_theme <- theme(aspect.ratio = .7,
                    legend.text = element_text(size = 7,
                                               face = "bold"),
                    legend.title = element_blank(),
                    axis.title = element_text(size = 8,
                                              face = "bold"),
                    axis.text = element_text(size = 7,
                                             face = "bold"),
                    plot.title = element_text(size = 8, hjust = 0),
                    plot.subtitle = element_text(size = 3),
                    strip.text = element_text(size = 8,
                                              margin = margin(.05,0,.05,0, "cm"),
                                              face = "bold"),
                    strip.background = element_blank(),
                    panel.spacing = unit(.1,"lines"))

#### Output Figure 5 #####
(u_plot_ula <- plotter2(data = U_values_format,
                        trait = "UpperLeafAngle",
                        ylim = c(15,90),
                        xlim = c(15,90),
                        incr = 15,
                        strip_title = TRUE,
                        u_cor = U_values_cor,
                        plot_theme = plot_theme))

#### Pulls out the legend to piece it back to the combined figure
legend <- gtable_filter(x = ggplotGrob(u_plot_ula),
                        pattern = "guide")

(u_plot_mla <- plotter2(data = U_values_format,
                        trait = "MiddleLeafAngle",
                        ylim = c(30,90),
                        xlim = c(30,90),
                        incr = 15,
                        strip_title = FALSE,
                        u_cor = U_values_cor,
                        plot_theme = plot_theme))

(u_plot_ph <- plotter2(data = U_values_format,
                       trait = "PlantHeight",
                       ylim = c(50,275),
                       xlim = c(50,275),
                       incr = 50,
                       strip_title = FALSE,
                       u_cor = U_values_cor,
                       plot_theme = plot_theme))

(u_plot_eh <- plotter2(data = U_values_format,
                       trait = "EarHeight",
                       ylim = c(0,175),
                       xlim = c(0,175),
                       incr = 25,
                       strip_title = FALSE,
                       u_cor = U_values_cor,
                       plot_theme = plot_theme))


combined_plot <- grid.arrange(arrangeGrob(rbind(ggplotGrob(u_plot_ula + theme(legend.position = "none",
                                                                              plot.margin = unit(c(0,0,0,.25),"cm"))),
                                                ggplotGrob(u_plot_mla + theme(legend.position = "none",
                                                                              plot.margin = unit(c(0,0,0,.25),"cm"))),
                                                ggplotGrob(u_plot_ph + theme(legend.position = "none",
                                                                             plot.margin = unit(c(0,0,0,.25),"cm"))),
                                                ggplotGrob(u_plot_eh + theme(legend.position = "none",
                                                                             plot.margin = unit(c(0,0,0,.25),"cm"))),
                                                size = "first"),
                                          left = 'Observed',
                                          bottom = 'Predicted'),
                              legend,
                              widths = unit.c(unit(1,"npc") - legend$width, legend$width),
                              nrow = 1)

# ggsave(plot = combined_plot,
#        filename = "Final_Figures/Figure_5.png",
#        width = 7,
#        height = 6,
#        dpi = 1200)

# ggsave(plot = combined_plot,
#        filename = "Final_Figures/Figure_5.pdf",
#        width = 180,
#        height = 152.4,
#        units = 'mm',
#        dpi = 1200)

#### Use empirical validation to predict training ####
#### Otherwise known as a reverse prediction
#### Input data for this part
train_pred_ids <- read_tsv("RRBLUP_Input_Files/Training_Imputed.012.indv",
                           col_names = FALSE) %>% 
  rename(Genotype = 'X1') %>% 
  mutate(Gen = row_number(),
         Genotype = toupper(Genotype))

#### Predicted training popolation GEBVs from 
#### using the empircial validation pop as training set
train_preds <- read_tsv("Predict_Phenotypes/Training_Predicted_BLUPs.txt") %>% 
  gather(Trait,Pred,-Gen) %>%
  left_join(train_pred_ids)

#### The training population's observed phenotypes
train_obs <- read_tsv("Phenotype_Data/Training/Training_Phenotypes_BLUPs.txt") %>% 
  gather(Trait,BLUP,-Genotype) %>% 
  mutate(Genotype = toupper(Genotype),
         Genotype = case_when(
           Genotype == "W22_R-RSTD" ~ "W22R-R-STD_CS-2909-1",
           TRUE ~ Genotype
         ))

#### Combinng the two data sets
train_comp <- train_obs %>% 
  filter(Trait %in% c("PlantHeight","EarHeight","MiddleLeafAngle","UpperLeafAngle")) %>% 
  left_join(train_preds) %>% 
  rename(Source = Genotype) %>% 
  left_join(subpops)

#### Doing some correlations and adding in the xy coordinates for the text
train_cor <- train_comp %>% 
  group_by(Trait) %>% 
  summarise(r = cor(Pred,BLUP,use = "complete.obs",method = "pearson"))  %>% 
  mutate(X2 = case_when(
    Trait %in% c("UpperLeafAngle","MiddleLeafAngle") ~ 82.5,
    Trait == "PlantHeight" ~ 250,
    Trait == "EarHeight" ~ 150),
    Y2 = case_when(
      Trait %in% c("UpperLeafAngle","MiddleLeafAngle") ~ 90,
      Trait == "PlantHeight" ~ 275,
      Trait == "EarHeight" ~ 175))

#### Breaking down by sub-populations and looking at correlations within subpop
#### Adds the xy coordinates for the text
train_cor_grps <- train_comp %>% 
  group_by(Trait,Subpop) %>% 
  summarise(r = cor(Pred,BLUP,use = "complete.obs",method = "pearson"),
            n = n(),
            .groups = 'drop') %>% 
  ungroup() %>% 
  mutate(X = case_when(
    Trait == "UpperLeafAngle" ~ 13,
    Trait == "MiddleLeafAngle" ~ 30,
    Trait == "PlantHeight" ~ 50,
    Trait == "EarHeight" ~ 0),
    Y = case_when(
      Subpop == "Stiff Stalk" & Trait == "UpperLeafAngle"  ~ 90,
      Subpop == "Non-stiff Stalk" & Trait == "UpperLeafAngle" ~ 85.42,
      Subpop == "Tropical" & Trait == "UpperLeafAngle" ~ 80.85,
      Subpop == "Popcorn" & Trait == "UpperLeafAngle" ~ 76.28,
      Subpop == "Sweet Corn" & Trait == "UpperLeafAngle" ~ 71.71,
      Subpop == "Mixed" & Trait == "UpperLeafAngle" ~ 67.14,
      Subpop == "Stiff Stalk" & Trait == "MiddleLeafAngle"  ~ 90,
      Subpop == "Non-stiff Stalk" & Trait == "MiddleLeafAngle" ~ 86.57,
      Subpop == "Tropical" & Trait == "MiddleLeafAngle" ~ 83.14,
      Subpop == "Popcorn" & Trait == "MiddleLeafAngle" ~ 79.71,
      Subpop == "Sweet Corn" & Trait == "MiddleLeafAngle" ~ 76.28,
      Subpop == "Mixed" & Trait == "MiddleLeafAngle" ~ 72.85,
      Subpop == "Stiff Stalk" & Trait == "PlantHeight" ~ 275,
      Subpop == "Non-stiff Stalk" & Trait == "PlantHeight" ~ 262.14,
      Subpop == "Tropical" & Trait == "PlantHeight" ~ 249.28,
      Subpop == "Popcorn" & Trait == "PlantHeight" ~ 236.42,
      Subpop == "Sweet Corn" & Trait == "PlantHeight" ~ 223.57,
      Subpop == "Mixed" & Trait == "PlantHeight" ~ 210.71,
      Subpop == "Stiff Stalk" & Trait == "EarHeight" ~ 175,
      Subpop == "Non-stiff Stalk" & Trait == "EarHeight" ~ 165,
      Subpop == "Tropical" & Trait == "EarHeight" ~ 155,
      Subpop == "Popcorn" & Trait == "EarHeight" ~ 145,
      Subpop == "Sweet Corn" & Trait == "EarHeight" ~ 135,
      Subpop == "Mixed" & Trait == "EarHeight" ~ 125),
    Subpop2 = str_replace(Subpop," ","~"))

#### Theme for the plots
plot_theme <- theme(aspect.ratio = .7,
                    legend.text = element_text(size = 5),
                    legend.title = element_blank(),
                    legend.margin = margin(0,0,0,-7),
                    axis.title = element_text(size = 9,
                                              face = "bold"),
                    axis.text = element_text(size = 8,
                                             face = "bold",
                                             color = "black"),
                    plot.title = element_text(size = 10, hjust = 0,
                                              face = "bold"),
                    plot.margin = unit(c(0.1,.1,0,.1),"cm"),
                    legend.key.size = unit(.75,"lines"))

#### Output for Sup Fig S9 ####
(plot_ULA_train <- plotter1(data = train_comp,
                            trait = "UpperLeafAngle",
                            title = "UpperLeafAngle",
                            xlim = c(13,90),
                            ylim = c(13,90),
                            incr = 15,
                            corData_Grps = train_cor_grps,
                            corData = train_cor,
                            theme = plot_theme))

(plot_MLA_train <- plotter1(data = train_comp,
                            trait = "MiddleLeafAngle",
                            title = "MiddleLeafAngle",
                            xlim = c(30,90),
                            ylim = c(30,90),
                            incr = 10,
                            corData_Grps = train_cor_grps,
                            corData = train_cor,
                            theme = plot_theme))

(plot_ph_train <- plotter1(data = train_comp,
                           trait = "PlantHeight",
                           title = "PlantHeight",
                           xlim = c(50,275),
                           ylim = c(50,275),
                           incr = 25,
                           corData_Grps = train_cor_grps,
                           corData = train_cor,
                           theme = plot_theme))

(plot_eh_train <- plotter1(data = train_comp,
                           trait = "EarHeight",
                           title = "EarHeight",
                           xlim = c(0,175),
                           ylim = c(0,175),
                           incr = 25,
                           corData_Grps = train_cor_grps,
                           corData = train_cor,
                           theme = plot_theme))

# ggsave(plot = grid.arrange((cbind(ggplotGrob(plot_ULA_train),ggplotGrob(plot_MLA_train),
#                                   size = "last")),
#                            (cbind(ggplotGrob(plot_ph_train),ggplotGrob(plot_eh_train),
#                                   size = "last"))),
#        filename = "Final_Figures/Supplemental_Figure_S9.png",
#        width = 7,
#        height = 5,
#        dpi = 1200)

# ggsave(plot = grid.arrange((cbind(ggplotGrob(plot_ULA_train),ggplotGrob(plot_MLA_train),
#                                   size = "last")),
#                            (cbind(ggplotGrob(plot_ph_train),ggplotGrob(plot_eh_train),
#                                   size = "last"))),
#        filename = "Final_Figures/Supplemental_Figure_S9.pdf",
#        width = 180,
#        height = 130,
#        units = 'mm',
#        dpi = 1200)