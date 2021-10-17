library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)
library(scales)
library(patchwork)

#### Excluded GEs ####
#### GEs that are also in the training population under a different name
excl_ges <- read_tsv("excluded_ges.txt")

#### Loading input data ####
#### Validation file and the GEs. Aligns with the genotypic data 
gen_id <- read_tsv("RRBLUP_Input_Files/Validation_Imputed.012.indv",
                   col_names = FALSE) %>% 
  rename(Genotype = 'X1') %>% 
  mutate(Gen = row_number())

#### Predicted GEBVs
pred_results <- read_tsv("Predict_Phenotypes/Predicted_BLUPs.txt") %>% 
  left_join(gen_id,
            by = "Gen") %>% 
  select(Gen,Genotype,everything()) %>% 
  filter(!(Genotype %in% excl_ges$Genotype))

#### U-Values
U_values <- read_tsv("U_Calculations/U_Results_Valid.txt") %>% 
  rename(Gen = Indv,
         U_Value = Result) %>% 
  left_join(gen_id) %>% 
  filter(!(Genotype %in% excl_ges$Genotype))

#### Subpopulation information for each GE
#### Also clasifies each sub pop into a type based on prior knowledge,
#### see text of manuscript for more information
subpops <- read_tsv("Population_Structure/Admixture/K_5_results.txt") %>%
  rename(Subpop = K_5,
         Genotype = Gen) %>% 
  mutate(Subpop = factor(case_when(
    Subpop == "Pop1" ~ "Sweet Corn",
    Subpop == "Pop2" ~ "Stiff Stalk",
    Subpop == "Pop3" ~ "Tropical",
    Subpop == "Pop4" ~ "Popcorn",
    Subpop == "Pop5" ~ "Non-stiff Stalk",
    TRUE ~ as.character(Subpop))) %>% 
      fct_relevel("Stiff Stalk","Non-stiff Stalk",
                  "Tropical","Popcorn",
                  "Sweet Corn","Mixed")) %>% 
  filter(!(Genotype %in% excl_ges$Genotype))

#### Pulling information for the large empirical population
#### See manuscript for citation 
### http://cbsusrv04.tc.cornell.edu/users/panzea/download.aspx?filegroupid=11
#### Downloaded the http://de.iplantcollaborative.org/dl/d/A1F3CF70-AE3A-4472-845F-AA8D6D7024EA/Peiffer2014Genetics_blupPhenos20150325.xlsx
#### file and then only pulled the data from the 'AMES' Panel
#### Most will be found, look at Peiffer name in file, it will match up with file listed above.
large_amesDP <- read_tsv("Phenotype_Data/Large_Empirical_Validation/BLUPs_Large_Empirical_Validation.txt") %>% 
  select(-Panzea_Name,-Peiffer_Name) %>% 
  gather(Trait,BLUP,-Genotype) %>% 
  filter(!(Genotype %in% excl_ges$Genotype)) %>% 
  ##Unit conversion, the data from Peiffer specifically Ames panel was listed in F instead of rest in C
  mutate(BLUP = case_when(
    Trait == 'GDDAnthesis-SilkingInterval' ~ BLUP * 5/9,
    Trait == 'GDDDaystoSilk' ~ BLUP * 5/9,
    Trait == 'GDDDaystoTassel' ~ BLUP * 5/9,
    TRUE ~ BLUP))

#### Data formatting for Sup Fig S10 and S11 ####
#### S10 is plotting the large empirical population out by trait and then looking
#### at prediction accuracy within and cross all subpopulations

#### S11 breaks these lines into groups based on their U_Value and then does the
#### same thing
#### Formatting the size of the u-value bins
u_size <- 400

#### Splitting the u-value results into bins based on top u_size, bottom u_size middle
large_uvalues <- large_amesDP %>% 
  select(Genotype) %>% 
  unique() %>% 
  left_join(U_values) %>% 
  arrange(desc(U_Value)) %>% 
  mutate(Rank = row_number(),
         Group = case_when(
           Rank <= u_size ~ paste("Top",u_size,sep = ""),
           Rank > (max(Rank)-u_size) ~ paste("Bottom",u_size,sep = ""),
           TRUE ~ "Middle")) %>% 
  select(-Rank) %>% 
  mutate(Group = factor(Group) %>% fct_relevel(c(paste("Top",u_size,sep = ""),
                                                 "Middle",
                                                 paste("Bottom",u_size,sep = ""))))

#### Formatting large emp pop and pulling in the predicted results
large_amesDP_format <- large_amesDP %>% 
  left_join(pred_results %>% 
              select(-Gen) %>% 
              gather(Trait,Pred,-Genotype)) %>% 
  left_join(subpops) %>% 
  left_join(U_values) %>% 
  left_join(large_uvalues)

#### Doing overall correlations of the predicted vs obs
#### Also adding manual text locations for the results, plotting purposes
overall_cor <- large_amesDP_format %>% 
  group_by(Trait) %>%
  summarise(Cor = cor(Pred,BLUP,use = "complete.obs",method = "pearson")) %>% 
  mutate(X = case_when(
    Trait == "DaystoSilk" ~ 107,
    Trait == "DaysToTassel" ~ 107,
    Trait == "EarHeight" ~ 118,
    Trait == 'GDDAnthesis-SilkingInterval' ~ 107,
    Trait %in% c("GDDDaystoSilk","GDDDaystoTassel") ~ 1325,
    Trait == "PlantHeight" ~ 205),
    Y = case_when(
      Trait %in% c("DaystoSilk","DaysToTassel") ~ 125,
      Trait == "EarHeight" ~ 150,
      Trait == 'GDDAnthesis-SilkingInterval' ~ 150,
      Trait %in% c("GDDDaystoSilk","GDDDaystoTassel") ~ 1600,
      Trait == "PlantHeight" ~ 250))

#### Calculating correlations within subpopulations
subpop_cor <- large_amesDP_format %>% 
  group_by(Trait,Subpop) %>%
  summarise(Cor = cor(Pred,BLUP,use = "complete.obs",method = "pearson"),
            n = n()) %>% 
  mutate(Subpop2 = str_replace(Subpop," ","~"),
         X = case_when(
           Trait == "DaystoSilk" ~ 50,
           Trait == "DaysToTassel" ~ 50,
           Trait == "EarHeight" ~ 10,
           Trait == 'GDDAnthesis-SilkingInterval' ~ -50,
           Trait %in% c("GDDDaystoSilk","GDDDaystoTassel") ~ 500,
           Trait == "PlantHeight" ~ 50),
         Y = case_when(
           Trait == "DaystoSilk" & Subpop == "Stiff Stalk" ~ 125,
           Trait == "DaystoSilk" & Subpop == "Non-stiff Stalk" ~ 120.5,
           Trait == "DaystoSilk" & Subpop == "Tropical" ~ 116,
           Trait == "DaystoSilk" & Subpop == "Popcorn" ~ 111.5,
           Trait == "DaystoSilk" & Subpop == "Sweet Corn" ~ 107,
           Trait == "DaystoSilk" & Subpop == "Mixed" ~ 102.5,
           Trait == "DaysToTassel" & Subpop == "Stiff Stalk" ~ 125,
           Trait == "DaysToTassel" & Subpop == "Non-stiff Stalk" ~ 120.38,
           Trait == "DaysToTassel" & Subpop == "Tropical" ~ 115.76,
           Trait == "DaysToTassel" & Subpop == "Popcorn" ~ 111.14,
           Trait == "DaysToTassel" & Subpop == "Sweet Corn" ~ 106.52,
           Trait == "DaysToTassel" & Subpop == "Mixed" ~ 101.9,
           Trait == "EarHeight" & Subpop == "Stiff Stalk" ~ 150,
           Trait == "EarHeight" & Subpop == "Non-stiff Stalk" ~ 141,
           Trait == "EarHeight" & Subpop == "Tropical" ~ 132,
           Trait == "EarHeight" & Subpop == "Popcorn" ~ 123,
           Trait == "EarHeight" & Subpop == "Sweet Corn" ~ 114,
           Trait == "EarHeight" & Subpop == "Mixed" ~ 105,
           Trait == "GDDAnthesis-SilkingInterval" & Subpop == "Stiff Stalk" ~ 150,
           Trait == "GDDAnthesis-SilkingInterval" & Subpop == "Non-stiff Stalk" ~ 136.68,
           Trait == "GDDAnthesis-SilkingInterval" & Subpop == "Tropical" ~ 123.36,
           Trait == "GDDAnthesis-SilkingInterval" & Subpop == "Popcorn" ~ 111.04,
           Trait == "GDDAnthesis-SilkingInterval" & Subpop == "Sweet Corn" ~ 98.72,
           Trait == "GDDAnthesis-SilkingInterval" & Subpop == "Mixed" ~ 86.4,
           Trait == "GDDDaystoSilk" & Subpop == "Stiff Stalk" ~ 1600,
           Trait == "GDDDaystoSilk" & Subpop == "Non-stiff Stalk" ~ 1532.24,
           Trait == "GDDDaystoSilk" & Subpop == "Tropical" ~ 1464.48,
           Trait == "GDDDaystoSilk" & Subpop == "Popcorn" ~ 1396.72,
           Trait == "GDDDaystoSilk" & Subpop == "Sweet Corn" ~ 1328.96,
           Trait == "GDDDaystoSilk" & Subpop == "Mixed" ~ 1261.2,
           Trait == "GDDDaystoTassel" & Subpop == "Stiff Stalk" ~ 1600,
           Trait == "GDDDaystoTassel" & Subpop == "Non-stiff Stalk" ~ 1532.24,
           Trait == "GDDDaystoTassel" & Subpop == "Tropical" ~ 1464.48,
           Trait == "GDDDaystoTassel" & Subpop == "Popcorn" ~ 1396.72,
           Trait == "GDDDaystoTassel" & Subpop == "Sweet Corn" ~ 1328.96,
           Trait == "GDDDaystoTassel" & Subpop == "Mixed" ~ 1261.2,
           Trait == "PlantHeight" & Subpop == "Stiff Stalk" ~ 250,
           Trait == "PlantHeight" & Subpop == "Non-stiff Stalk" ~ 238,
           Trait == "PlantHeight" & Subpop == "Tropical" ~ 226,
           Trait == "PlantHeight" & Subpop == "Popcorn" ~ 214,
           Trait == "PlantHeight" & Subpop == "Sweet Corn" ~ 202,
           Trait == "PlantHeight" & Subpop == "Mixed" ~ 190))


#### Custom plotter for plotting these results ####
trait_plot1 <- function(data, trait, xLim, yLim, 
                        OverallCor, plotTheme, SubpopCor, title) {
  title <- title
  plot <- data %>% 
    filter(Trait == trait) %>% 
    ggplot(aes(x = Pred,
               y = BLUP)) +
    geom_point(size = .3,
               aes(color = Subpop),
               alpha = 0.5,
               show.legend = FALSE) +
    geom_text(data = SubpopCor %>% 
                filter(Trait == trait),
              aes(x = X,
                  y = Y,
                  label = paste("bold(",
                                paste(Subpop2, "~(",n,")~","italic(r)","~'='~",
                                      sprintf('"%1.2f"',Cor),
                                      sep = ""),
                                ")",sep = ""),
                  color = Subpop),
              parse = TRUE,
              size = 2.45,
              hjust = 0,
              show.legend = FALSE) +
    geom_text(data = OverallCor %>%
                filter(Trait == trait),
              aes(x = X,
                  y = Y,
                  label = paste("bold(",
                                paste("Overall~italic(r)","~'='~",
                                      sprintf('"%1.2f"',Cor),
                                      sep = ""),
                                ")",sep = "")),
              parse = TRUE,
              size = 2.5,
              hjust = 0,
              show.legend = FALSE) +
    ylab("Observed") +
    xlab("Predicted") +
    ggtitle(paste(title,sep = "")) +
    theme_bw() +
    scale_y_continuous(limits = yLim,
                       breaks = pretty_breaks(8)) +
    scale_x_continuous(limits = xLim,
                       breaks = pretty_breaks(8)) +
    plotTheme +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(.85,.25)) +
    guides(colour = guide_legend(override.aes = list(size = 1,
                                                     alpha = 1),
                                 title = "Subpopulation"))
  
  return(plot)
}

#### Theme setup for the plots
plot_theme <- theme(aspect.ratio = .85,
                    legend.text = element_text(size = 6,
                                               face = "bold"),
                    legend.title = element_blank(),
                    legend.margin = margin(0,0,0,-7),
                    axis.title = element_text(size = 8,
                                              face = "bold"),
                    axis.text = element_text(size = 8,
                                             face = "bold",
                                             color = "black"),
                    plot.title = element_text(size = 10, hjust = 0,
                                              margin = unit(c(.1,0,0,0),"cm"),
                                              face = "bold"),
                    legend.key.size = unit(.75,"lines"))

#### Left and right plots have different themes, formatting purposes
#### Mostly for label spacing
plot_theme_l <- plot_theme + theme(plot.margin = unit(c(0,.1,0,0),"cm"))
plot_theme_r <- plot_theme + theme(plot.margin = unit(c(0,0,0,0.1),"cm"))

#### Output
(plot1 <- trait_plot1(data = large_amesDP_format,
                      trait = "PlantHeight",
                      title = "PlantHeight",
                      xLim = c(50,250),
                      yLim = c(50,250),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_l))

(plot2 <- trait_plot1(data = large_amesDP_format,
                      trait = "EarHeight",
                      title = "EarHeight",
                      xLim = c(10,150),
                      yLim = c(10,150),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_r))

(plot3 <- trait_plot1(data = large_amesDP_format,
                      trait = "DaystoSilk",
                      title = "DaystoSilk",
                      xLim = c(50,125),
                      yLim = c(50,125),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_l))

(plot4 <- trait_plot1(data = large_amesDP_format,
                      trait = "GDDDaystoSilk",
                      title = "GDDDaystoSilk",
                      xLim = c(500,1600),
                      yLim = c(500,1600),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_r))

(plot5 <- trait_plot1(data = large_amesDP_format,
                      trait = "DaysToTassel",
                      title = "DaysToTassel",
                      xLim = c(50,125),
                      yLim = c(50,125),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_l))

(plot6 <- trait_plot1(data = large_amesDP_format,
                      trait = "GDDDaystoTassel",
                      title = "GDDDaystoTassel",
                      xLim = c(500,1600),
                      yLim = c(500,1600),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_r))

(plot7 <- trait_plot1(data = large_amesDP_format,
                      trait = "GDDAnthesis-SilkingInterval",
                      title = "GDDAnthesis-SilkingInterval",
                      xLim = c(-50,150),
                      yLim = c(-50,150),
                      OverallCor = overall_cor,
                      SubpopCor = subpop_cor,
                      plotTheme = plot_theme_l))

#### Adds a blank plot at the end otherwise it tries to stretch this plot out
plot8 <- ggplot() + theme_bw() + theme(plot.margin = unit(c(0,0,0,0.1),"cm"),
                                       panel.border = element_blank())
#### Patchwork to put these all together
plot_out <- (plot1 + plot2) / (plot3 + plot4) / (plot5 + plot6) / (plot7 + plot8)

# ggsave(plot = plot_out,
#        filename = "Final_Figures/Supplemental_Figure_S10.png",
#        width = 7,
#        height = 12,
#        dpi = 1200)

# ggsave(plot = plot_out,
#        filename = "Final_Figures/Supplemental_Figure_S10.pdf",
#        width = 180,
#        height = 305,
#        units = 'mm',
#        dpi = 1200)

#### Data format for Sup Fig 11 ####
#### Breaking down prediction accuracy by U_Value group
group_cor <- large_amesDP_format %>% 
  group_by(Trait,Group) %>%
  summarise(Cor = cor(Pred,BLUP,use = "complete.obs",method = "pearson")) %>% 
  ungroup() %>% 
  mutate(X = case_when(
    Trait == "DaystoSilk" ~ 50,
    Trait == "DaysToTassel" ~ 25,
    Trait == "EarHeight" ~ 0,
    Trait == 'GDDAnthesis-SilkingInterval' ~ -50,
    Trait %in% c("GDDDaystoSilk","GDDDaystoTassel") ~ 500,
    Trait == "PlantHeight" ~ 50),
    Y = case_when(
      Trait %in% c("DaystoSilk","DaysToTassel") ~ 125,
      Trait == "EarHeight" ~ 150,
      Trait == 'GDDAnthesis-SilkingInterval' ~ 150,
      Trait %in% c("GDDDaystoSilk","GDDDaystoTassel") ~ 1600,
      Trait == "PlantHeight" ~ 250),
    Trait = ifelse(Trait == 'GDDAnthesis-SilkingInterval',
                    'GDDAnthesis-\nSilkingInterval\n',Trait))

#### New plotting function, custom setup for the results
trait_plot2 <- function(data,trait, xLim,yLim, GroupCor, plotTheme, x_strip_text) {
  plot <- data %>% 
    filter(Trait == trait) %>% 
    ggplot(aes(x = Pred,
               y = BLUP)) +
    geom_point(size = .25,aes(col = Subpop),
               alpha = 0.25) +
    facet_grid(Trait ~ Group, scales = "free") +
    geom_text(data = GroupCor %>%
                filter(Trait == trait),
              aes(x = X,
                  y = Y,
                  label = paste("bold(",
                                paste("Overall~italic(r)","~'='~",sprintf('"%1.2f"',Cor),
                                      sep = ""),
                                ")",
                                sep = "")),
              parse = TRUE,
              size = 2.25,
              hjust = 0,
              show.legend = FALSE) +
    theme_bw() +
    scale_y_continuous(limits = yLim,
                       breaks = pretty_breaks(6)) +
    scale_x_continuous(limits = xLim,
                       breaks = pretty_breaks(6)) +
    plotTheme +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 1,
                                                    alpha = 1)))
  
  if(!x_strip_text) {
    plot <- plot + theme(strip.text.x = element_blank())
  }
  return(plot)
}

#### Setting up the plot's theme
plot_theme <- theme(aspect.ratio = .7,
                    legend.text = element_text(size = 7,
                                               face = "bold"),
                    legend.title = element_blank(),
                    axis.title = element_text(size = 8,
                                              face = "bold"),
                    axis.text = element_text(size = 7,
                                             face  = "bold",
                                             color = "black"),
                    strip.text = element_text(size = 8,
                                              margin =margin(.05,0,.05,0, "cm"),
                                              face = "bold"),
                    strip.background = element_blank(),
                    panel.spacing = unit(.1,"lines"),
                    legend.margin = margin(0,0,0,1))

(plot8 <- trait_plot2(data = large_amesDP_format,
                      trait = "PlantHeight",
                      xLim = c(50,250),
                      yLim = c(50,250),
                      GroupCor = group_cor,
                      plotTheme = plot_theme,
                      x_strip_text = TRUE))

(plot9 <- trait_plot2(data = large_amesDP_format,
                      trait = "EarHeight",
                      xLim = c(0,125),
                      yLim = c(0,150),
                      GroupCor = group_cor,
                      plotTheme = plot_theme,
                      x_strip_text = FALSE))
(plot10 <- trait_plot2(data = large_amesDP_format,
                       trait = "DaystoSilk",
                       xLim = c(50,125),
                       yLim = c(50,125),
                       GroupCor = group_cor,
                       plotTheme = plot_theme,
                       x_strip_text = FALSE))
(plot11 <- trait_plot2(data = large_amesDP_format,
                       trait = "GDDDaystoSilk",
                       xLim = c(500,1600),
                       yLim = c(500,1600),
                       GroupCor = group_cor,
                       plotTheme = plot_theme,
                       x_strip_text = FALSE))
(plot12 <- trait_plot2(data = large_amesDP_format,
                       trait = "DaysToTassel",
                       xLim = c(25,125),
                       yLim = c(25,125),
                       GroupCor = group_cor,
                       plotTheme = plot_theme,
                       x_strip_text = FALSE))
(plot13 <- trait_plot2(data = large_amesDP_format,
                       trait = "GDDDaystoTassel",
                       xLim = c(500,1600),
                       yLim = c(500,1600),
                       GroupCor = group_cor,
                       plotTheme = plot_theme,
                       x_strip_text = FALSE))

(plot14 <- trait_plot2(data = large_amesDP_format %>% mutate(Trait = ifelse(Trait == "GDDAnthesis-SilkingInterval",
                                                                             "GDDAnthesis-\nSilkingInterval\n",Trait)),
                       trait = "GDDAnthesis-\nSilkingInterval\n",
                       xLim = c(-50,150),
                       yLim = c(-50,150),
                       GroupCor = group_cor,
                       plotTheme = plot_theme,
                       x_strip_text = FALSE))

#### Pulling out the legend to add it back in later
legend <- gtable_filter(x = ggplotGrob(plot8),
                        pattern = "guide")

#### Fixing for spacing of plots
left <- .4
right <- .1

combined_plot <- grid.arrange(arrangeGrob(rbind(ggplotGrob(plot8 + theme(legend.position = "none",
                                                                         plot.margin = unit(c(0,right,0,left),"cm"))),
                                                ggplotGrob(plot9 + theme(legend.position = "none",
                                                                         plot.margin = unit(c(0,right,0,left),"cm"))),
                                                ggplotGrob(plot10 + theme(legend.position = "none",
                                                                          plot.margin = unit(c(0,right,0,left),"cm"))),
                                                ggplotGrob(plot11 + theme(legend.position = "none",
                                                                          plot.margin = unit(c(0,right,0,left),"cm"))),
                                                ggplotGrob(plot12 + theme(legend.position = "none",
                                                                          plot.margin = unit(c(0,right,0,left),"cm"))),
                                                ggplotGrob(plot13 + theme(legend.position = "none",
                                                                          plot.margin = unit(c(0,right,0,left),"cm"))),
                                                ggplotGrob(plot14 + theme(legend.position = "none",
                                                                          plot.margin = unit(c(0,right,0,left),"cm"))),
                                                size = "first"),
                                          left = "Observed",
                                          bottom = "Predicted"),
                              legend,
                              widths = unit.c(unit(1,"npc") - (1.25*legend$width), legend$width),
                              nrow = 1)

# ggsave(plot = combined_plot,
#        filename = "Final_Figures/Supplemental_Figure_S11.png",
#        width = 7,
#        height = 10,
#        dpi = 1200)

# ggsave(plot = combined_plot,
#        filename = "Final_Figures/Supplemental_Figure_S11.pdf",
#        width = 7,
#        height = 10,
#        dpi = 1200)
