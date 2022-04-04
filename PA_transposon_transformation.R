
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(openxlsx)

lib = read_excel("NRSetFile_v5_061004.xls", 
                 sheet = "Sheet1")


names(lib)

lib = lib %>% 
  select(Plate = `PA14 NR Set Plate (PO_PlateLabel96)`,
         Well = `PA14 NR Set Well (PO_PlateLocation96ADD)`,
         PA_ortholog_active = `PAO1 ortholog of \"Active\" gene`,
         gene_name = `\"Active\" Gene Name`,
         gene_description = `\"Active\" Gene Description`,
         PA_ortholog = `PAO1 Orthologs`,
         notes = `PA14 NR Set Comments`, 
         everything())



# are there 96 wells per plate? 
lib %>% 
  group_by(Plate) %>% count() %>% arrange((n)) %>% 
  ggplot(aes(n)) +
  geom_histogram(stat = 'count')



# how many plates
length(unique(lib$Plate))


# create new IDs for the NA --------------
lib_fix = lib %>% 
  #create shorter ids for the plates
  mutate(plate_id = paste0(str_sub(Plate, end=4), str_sub(Plate, start = -5)),
         .before=Well) %>% 
  # fix the wells with gene names and PA orthologs, and Empty wells
  mutate(new_id = case_when(!is.na(gene_name) ~ gene_name,
                            str_detect(notes, 'Empty') ~ 'Empty',
                            is.na(gene_name) ~ PA_ortholog),
         .before = gene_description) %>% 
  mutate(new_id = case_when(is.na(new_id) & str_detect(gene_description, '\\wypothetical') ~ 
                              paste0('hy_',seq(1,n()), '_', plate_id),
                            is.na(new_id) & str_detect(gene_description, '\\wutative') ~
                              paste0('pu_',seq(1,n()), '_', plate_id),
                            # weird gene that escaped somehow
                            is.na(new_id) & gene_description == 'cupD5' ~ 'cupD5',
                            is.na(new_id) & gene_description == 'transposase' ~ 'transposase',
                            is.na(new_id) & gene_description == 'pyocin protein' ~ 'pyocin protein',
                            is.na(new_id) & str_detect(gene_description, '\\wrobable') ~
                              paste0('prob_',seq(1,n()), '_', plate_id),
                            is.na(new_id) & !is.na(gene_description) ~ 
                              paste0(gene_description, '_', plate_id),
                            is.na(new_id) & is.na(gene_description) ~ paste0('Weird_',plate_id),
                            TRUE ~ new_id)
         ) 



# 96 well format hell -----------------------------------------------------


# initialise variables
plates = unique(lib$Plate)
row_names = LETTERS[1:8]
col_names = seq(1,12)

plates_list = list()
for (plate in plates) {
  
  plate_name = as.character(plate)
  
  temp_ids = lib_fix %>% 
    mutate(Column = str_extract(Well, '\\w'),
           Row = str_extract(Well, '\\d{1,}'),
           .before=PA_ortholog_active) %>% 
    filter(Plate == plate) %>% 
    mutate(Row = as.numeric(Row)) %>% 
    arrange(Column, Row) %>% 
    pull(new_id) 
  
  temp_matrix = matrix(data = temp_ids, 
                       nrow = 8, ncol = 12, byrow = T)
  
  rownames(temp_matrix) = row_names
  colnames(temp_matrix) = col_names
  
  temp_list = list(plate_name = temp_matrix)
  
  plates_list = append(plates_list, temp_list)
}

# fix naming of plates in list
names(plates_list) = plates

#check that it's ok
plates_list


sink("PA_transposon_library_onepage.csv", type="output")
invisible(lapply(names(plates_list), 
                 function(x) { print(x)
                   dput(write.csv(plates_list[[x]])) } ))
sink()


write.xlsx(plates_list, 'PA_transposon_library_multipage.xlsx',
           overwrite = T,
           row.names = T)
