
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




# B. subtilis library ----------------------------------------------------

# 05/10/2022
# Jen asked me to do the same with this B. subtilis library
# put it in 384 and 96 plate format 

# read the table
bsub = read_excel("B.subtilis library.xlsx", 
                  sheet = "BKK_2017_for distribution") %>% 
  rename(gene_name = `Gene name`, 
         strain_name = `Strain name`,
         locus_tag = `Locus tag`,
         supplement = `supplement for better growth`)

# are there any NAs in the gene names? 
bsub %>% 
  filter(is.na(gene_name))

# fill NA with well_plate info
bsub = bsub %>% 
  mutate(gene_name = case_when(is.na(gene_name) ~ paste(Plate, well, sep = '_'),
                               TRUE ~ gene_name))


## 96 well format -------
# initialise variables
plates = unique(bsub$Plate)
row_names = LETTERS[1:8]
col_names = seq(1,12)

plates_list = list()
for (plate in plates) {
  
  plate_name = as.character(plate)
  
  temp_ids = bsub %>% 
    mutate(Column = str_extract(well, '\\w'),
           Row = str_extract(well, '\\d{1,}'),
           .before=strain_name) %>% 
    filter(Plate == plate) %>% 
    mutate(Row = as.numeric(Row)) %>% 
    arrange(Column, Row) %>% 
    pull(gene_name) 
  
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

### save files ####
sink("Bsubtilis_96format_onepage.csv", type="output")
invisible(lapply(names(plates_list), 
                 function(x) { print(x)
                   dput(write.csv(plates_list[[x]])) } ))
sink()


write.xlsx(plates_list, 'Bsubtilis_96format_multipage.xlsx',
           overwrite = T,
           row.names = T)



## 384 well format -------

# import the info of the new plates

bsub_384 = read_excel("B.subtilis library.xlsx", 
                       sheet = "384 format") %>% 
  drop_na() %>% 
  rename(Plate_384 = `384 Plate name`,
         Plate = `...3`)

# as we are transferring things from a 96 to a 384 format, positions will be 
# displaced from the origin. The first top half of the 384 plate will consist on 
# plates 1 and 2 of the 96, and the second bottom half will be from plates 
# 3 and 4. 
# The positions of plate A1 will be: A1, A3, A5...
# The positions of plate A2 will be: A2, A4, A6...
# The positions of plate B1 will be: I1, I3, I5...
# The positions of plate B2 will be: I2, I4, I6...

# the strategy is to create vectors with the new positions, join them to the
# original DF and then change the coordinates for the new ones. Then
# repeat the previous loop
# create vectors with new positions

### ORIGIN ------
row_names = LETTERS[1:8]
col_names = seq(1,12, 1) 
origin_wells = expand_grid(row_names, col_names) %>% 
  # this line is because the format of wells is A01 instead of A1
  mutate(col_names = str_pad(col_names, pad=0, width = 2, side ='left')) %>% 
  unite(well,  row_names:col_names, sep = '')

### plate A1 ------
row_names = LETTERS[1:8]
col_names = seq(1,24, 2) # odd numbers
A1_wells = expand_grid(row_names, col_names) %>% 
  mutate(col_names = str_pad(col_names, pad=0, width = 2, side ='left')) %>% 
  unite(well_A1,  row_names:col_names, sep = '')

### plate A2 ------
row_names = LETTERS[1:8]
col_names = seq(2,24, 2) # even numbers
A2_wells = expand_grid(row_names, col_names) %>% 
  mutate(col_names = str_pad(col_names, pad=0, width = 2, side ='left')) %>% 
  unite(well_A2,  row_names:col_names, sep = '')


### plate B1 ------
row_names = LETTERS[9:16] # second subset of letters
col_names = seq(1,24, 2) # odd numbers
B1_wells = expand_grid(row_names, col_names) %>% 
  mutate(col_names = str_pad(col_names, pad=0, width = 2, side ='left')) %>% 
  unite(well_B1,  row_names:col_names, sep = '')

### plate B2 ------
row_names = LETTERS[9:16]
col_names = seq(2,24, 2) # even numbers
B2_wells = expand_grid(row_names, col_names) %>% 
  mutate(col_names = str_pad(col_names, pad=0, width = 2, side ='left')) %>% 
  unite(well_B2,  row_names:col_names, sep = '')


positions = bind_cols(origin_wells, 
                      A1_wells,
                      A2_wells,
                      B1_wells,
                      B2_wells)



# modify the well info 
bsub = bsub %>% 
  left_join(bsub_384) %>% 
  left_join(positions) %>% 
  mutate(well = case_when(Arrangement == 'Position A1' ~ well_A1,
                          Arrangement == 'Position A2' ~ well_A2,
                          Arrangement == 'Position B1' ~ well_B1,
                          Arrangement == 'Position B2' ~ well_B2))


bsub = bsub %>% 
  select(Plate, Plate_384, well, gene_name)

# initialise variables
plates = unique(bsub$Plate_384)
row_names = LETTERS[1:16]
col_names = seq(1,24)

plates_list = list()
for (plate in plates) {
  
  plate_name = as.character(plate)
  
  temp_ids = bsub %>% 
    mutate(Column = str_extract(well, '\\w'),
           Row = str_extract(well, '\\d{1,}'),
           .before=well) %>% 
    filter(Plate_384 == plate) %>% 
    mutate(Row = as.numeric(Row)) %>% 
    arrange(Column, Row) %>% 
    pull(gene_name) 
  
  temp_matrix = matrix(data = temp_ids, 
                       nrow = 16, ncol = 24, byrow = T)
  
  rownames(temp_matrix) = row_names
  colnames(temp_matrix) = col_names
  
  temp_list = list(plate_name = temp_matrix)
  
  plates_list = append(plates_list, temp_list)
}

# fix naming of plates in list
names(plates_list) = plates

#check that it's ok
plates_list


### save files ####
sink("BSubtilis/Bsubtilis_384format_onepage.csv", type="output")
invisible(lapply(names(plates_list), 
                 function(x) { print(x)
                   dput(write.csv(plates_list[[x]])) } ))
sink()


write.xlsx(plates_list, 'BSubtilis/Bsubtilis_384format_multipage.xlsx',
           overwrite = T,
           row.names = T)


