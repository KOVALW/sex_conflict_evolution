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
raw_count_data %>% mutate(r=row_number()) %>% 
  filter(offspring_male + offspring_female != offspring_total)
