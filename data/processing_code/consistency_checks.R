source("data/processing_code/data_cleaning.R")

c("GD","KK") == unique(raw_count_data$background_library)

#Unsure of the following table because construct names in cross don't match line ID
raw_count_data %>% 
  mutate(r = row_number()) %>% 
  filter(driver == "Nos-GAL4") %>% 
  group_by(cross_name, background_library, 
           driver, construct_line,
           cross_line_name) %>% 
  summarise(replicates = n()) %>% 
  ungroup() %>%   
  group_by(construct_line) %>%
  mutate(o=n()) %>%  
  ungroup() %>% 
  filter(o > 2 & !is.na(construct_line)) %>% 
  arrange(construct_line) %>% 
  select(-o) %>% 
  write_csv("./data/data_cleaned/verify_linenumber_table.csv")

#one mismatch identified - now checks out
driver_lede_key <- raw_count_data %>% 
  group_by(driver) %>% 
  summarise(o=n()) %>% 
  filter(!is.na(driver)) %>% 
  select(-o) %>% 
  mutate(driver_lede = substr(driver,1,3))

raw_count_data %>% mutate(r=row_number()) %>% 
  filter(offspring_male + offspring_female != offspring_total)

#Recreate experiment key-----
experiment_list <- gene_driver_pair %>% 
  select(construct_line, gene_name,ends_with("library")) %>% 
  gather(-construct_line, -gene_name, key = "driver_lede", value = "libraries_used") %>% 
  mutate(driver_lede = substr(driver_lede, 1,3),
         GD = grepl("GD", libraries_used),
         KK = grepl("KK", libraries_used)) %>% 
  select(-libraries_used) %>% 
  gather(-construct_line, -gene_name, -driver_lede, key = "background_library", value = "used") %>% 
  filter(used) %>% 
  left_join(driver_lede_key,
            by = "driver_lede") %>% 
  select(-used, -driver_lede)
  
observed_experiments <- cleaned_count_data %>% 
  group_by(construct_line, gene_name, background_library, driver) %>% 
  summarise(obs= n(),
            ) %>% 
  ungroup() 

observed_experiments %>% 
  arrange(construct_line, background_library, driver)

experiment_list %>% 
  arrange(construct_line, background_library, driver)
