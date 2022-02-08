library(tidyverse)
library(readxl)
library(broom)
library(cowplot)
library(ggrepel)
library(here)

# Data from  Nov - 2022 -------------------------------------------------------


# opening data
ms_data = read_xlsx(here('data', 'for_r_analysis_GNPS_results.xlsx')) %>%
  select(-Checked, -Tags, -Annot_Source_Predicted_Compositions, 
         -Annot_Source_mzCloud_Search,
         -Annot_Source_Metabolika_Search,
         -Annot_Source_ChemSpider_Search, 
         -Annot_DeltaMass_ppm, 
         -ChemSpider_Results, 
         -mzCloud_Results, 
         -mzCloud_Best_Match, 
         -mzCloud_Best_Match_Confidence, 
         -mzVault_Best_Match, 
         -MS2) %>% 
  rename(KO = `Group_CV_ KO`, WT = Group_CV_WT)


# closer look and normalising for OD
norm_ms_data = ms_data %>%
  # drop_na(Name) %>%  ### this is also an option for NA removal
  filter(!is.na(Name)) %>% # Removes unannotated metabolites
  rename(Metabolite = Name) %>%
  select(Metabolite, Formula, Annot_Source_mzVault_Search, Area_KO, 
         Area_WT) %>%
  mutate(Area_KO = Area_KO/2.228,
         Area_WT = Area_WT/2.139)

# repeated elements, could be molecules with different groups/charges? 
norm_ms_data %>% 
  count(Metabolite) %>% 
  arrange(desc(n)) %>% view

# PCA fit
pca_fit = norm_ms_data %>% 
  select(Area_WT, Area_KO) %>% 
  prcomp(scale = TRUE)

pca_fit %>%
  tidy(matrix = "rotation")


thres_y=0.3
thres_x=3
pca_fit %>%
  augment(norm_ms_data) %>% # add original dataset back in
  mutate(new_name = case_when(.fittedPC1 > thres_x | 
                                abs(.fittedPC2) > thres_y ~ Metabolite)) %>% 
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  geom_label_repel(aes(label = new_name)) +
  theme_half_open(12) + 
  background_grid()


pca_fit %>%
  augment(norm_ms_data) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 2, alpha = .3) +
  theme_half_open(12) + 
  background_grid()




norm_ms_data %>%
  ggplot(aes(x = Area_WT, y = Area_KO)) +
  geom_point() +
  geom_segment(x = 0, y = 0, yend = 3955652465, xend = 3955652465) +
  theme_half_open(13) +
  background_grid()



norm_ms_data %>% 
  mutate(prop = Area_KO/Area_WT) %>% 
  ggplot(aes(y = prop, x = fct_reorder(Metabolite, prop, .desc=TRUE))) +
  geom_point() +
  theme_cowplot(12) +
  theme(axis.text.x = element_blank())




# Data from Jen 08-02-2022 ------------------------------------------------



# opening data
metabs =  
  read_xlsx( here('data', 
                          'Annotated_compounds_pseudomonas_test.xlsx')) %>%
  select(-Checked, -Tags, -`Annot. Source: Predicted Compositions`:-`Annot. DeltaMass [ppm]`, 
         -`# ChemSpider Results`:-`# mzVault Results`,
         -`mzCloud Best Match`:-`Reference Ion`,
         -`RT [min]`, -`Area (Max.)`,
         -`Area: LB_WT_50.raw (F1)`:-`Area: LB_WT_0.raw (F16)`) %>% 
  rename(
    calc_MW = `Calc. MW`,
    LB_WT_50 = `Group Area: LB_WT_50`,
    NGM_WT_0 = `Group Area: WT_0`,
    NGM_WT_50 = `Group Area: WT_50`,
    NGM_BioF_0 = `Group Area: BioF_0`,
    NGM_BioF_50 = `Group Area: BioF_50`,
    NGM_GacA_0 = `Group Area: GacA_0`,
    NGM_GacA_50 = `Group Area: GacA_50`,
    LB_GacA_0 = `Group Area: LB_GacA_0`,
    LB_GacA_50 = `Group Area: LB_GacA_50`,
    LB_WT_0 = `Group Area: LB_WT_0`)


# which of the metabolites are repeated
metabolites = metabs %>% drop_na(Name) %>%  pull(Name)

metabolites[duplicated(metabolites)]

# distinguish 
metabs_nondup = metabs %>% 
  drop_na(Name) %>% 
  group_by(Name) %>% 
  mutate(rep  = seq(1,n())) %>% 
  mutate(Name_unique = 
           if(n() > 1) 
             {paste(Name, rep, sep = '_')}
          else 
            {paste0(Name)},
         .before = Formula
         ) %>% 
  select(-rep)



# get the table in shape for metaboanalyst

first_step = metabs_nondup %>% 
  ungroup %>% 
  select(Name_unique, LB_WT_50:LB_WT_0)  %>% 
  data.frame()

unique_metabs = first_step$Name_unique

first_step[,1] = NULL

second_step = t(first_step)


colnames(second_step) = unique_metabs

second_step %>% 
  as_tibble(rownames = 'Sample') %>% 
  # separate(Sample, into = c('Media', 'Strain', 'Metformin_mM'), sep = '_') %>% 
  # mutate(Media = factor(Media, levels = c('NGM', 'LB')),
  #        Strain = factor(Strain, levels = c('WT', 'BioF', 'GacA')),
  #        Metformin_mM = factor(Metformin_mM, levels = c('0', '50'))) %>%
  write_csv(here('summary', 'metaboanalyst_metabs.csv'))

# convert it to a long table 
metabs_long = metabs_nondup %>% 
  pivot_longer(cols = LB_WT_50:LB_WT_0, 
               names_to = 'Sample',
               values_to = 'area') %>% 
  separate(Sample, into = c('Media', 'Strain', 'Metformin_mM'), sep = '_') %>% 
  mutate(Media = factor(Media, levels = c('NGM', 'LB')),
         Strain = factor(Strain, levels = c('WT', 'BioF', 'GacA')),
         Metformin_mM = factor(Metformin_mM, levels = c('0', '50'))) %>% 
  select(-`m/z`,-`InChI Key`)

# check integrity of data table
metabs_long










