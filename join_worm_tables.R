library(tidyverse)
library(readr)
library(openxlsx)
library(glue)

# read the list of files and create a tibble
files = list.files('inputs/', pattern = '.csv') %>% as_tibble() %>% 
  rename(file = value)

# specify new columns
new_col_names = c(
  'file',
  'replicate',
  'date',
  'plate_type',
  'plate_number',
  'well',
  'zstack',
  'channel'
  )

# split the filenames and unnest them
file_metadata = files %>% 
  mutate(col_info = str_split(file, pattern = '_')) %>% 
  unnest_wider(col_info) %>% 
  select(-`...8`)
# change column names for the new ones
colnames(file_metadata)  = new_col_names


# loop to read each file and add the file name, it will
# help to join the metadata
merge_file = tibble()
for (file in file_metadata$file) {
  temp = read_csv(glue('inputs/{file}')) %>% 
    mutate(file = file,
           `User Label` = as.character(`User Label`))
  
  merge_file = bind_rows(merge_file, temp)
}

# chech everything is ok
file_metadata %>% 
  left_join(merge_file) %>% view

merge_file = file_metadata %>% 
  left_join(merge_file) 

# fix the plate number for ExMr 15-3
merge_file = merge_file %>% 
  mutate(plate_number = case_when(plate_number == '15-3' ~ '15_3',
                                  TRUE ~ plate_number),
         plate_type = case_when(plate_type == 'PAMR' ~ 'PAMr',
                                TRUE ~ plate_type))

# save the table
merge_file %>% 
  write_csv('merge_tables.csv')




# merge with gene info ----------------------------------------------------


library(readxl)
pa14 = read_excel("PA14_TP_annotation.xls")

controls = read_excel("PATP_Control_plate.xlsx", 
                      sheet = "control_list") %>% 
  mutate(plate_number = as.character(plate_number))

pa_metadata = pa14 %>% 
  distinct(PA14_plate_nr) %>% 
  mutate(Plate_type = str_split(PA14_plate_nr, n = 4, pattern = '_')) %>% 
  unnest_wider(Plate_type)

pa_metadata = pa_metadata %>% 
  select(-`...3`, -`...2`) %>% 
  rename(plate_type = `...1`, plate_number = `...4`)

pa14 = pa14 %>% 
  left_join(pa_metadata) %>% 
  bind_rows(controls)


pa14 %>% write_csv('PA14_gene_list_extended.csv')




# join the tables ---------------------------------------------------------



merge_file %>% 
  left_join(pa14) %>% 
  select(plate_type, well, plate_number, everything()) %>% 
  write_csv('worms_metadata_181022.csv')

