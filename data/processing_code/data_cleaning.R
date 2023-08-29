library(tidyverse)
library(lubridate)

raw_count_data <- read_csv("./data/data_raw/final_count_data.csv", guess_max = 5e4)
gene_key <- read_csv("./data/data_raw/gene_key.csv")
gene_driver_pair <- read_csv("data/data_raw/experiment_key.csv")

#Cleaning library and line data values-----
#should coerce some values to NA
libraries_incl <- c("GD","KK")

raw_count_data <- raw_count_data %>%
  mutate(background_line = as.numeric(background_line)) %>% 
  mutate(driver_present = !is.na(driver)) %>% 
  mutate(background_library = ifelse(!background_library %in% libraries_incl , 
                                     ifelse(background_line == 60100, "KK",
                                            ifelse(background_line == 808, "GD",
                                                   NA)),
                                     background_library))

raw_count_data <- raw_count_data %>% 
  mutate(construct_line = ifelse(!background_line %in% c(60100, 808),
                                 background_line,
                                 NA)) %>% 
  mutate(background_line = ifelse(background_library == "GD",
                                  808, 60100))

raw_count_data <- raw_count_data %>%
  mutate(rnai_construct_present = !is.na(construct_line)) %>%
  mutate(cross_line_name = ifelse(
    !is.na(gene_name),
    gene_name,
    ifelse(driver == "Nos-GAL4" & rnai_construct_present,
           ifelse(background_library == "GD",
                  substr(cross_name, 14, nchar(cross_name) - 2),
                  substr(cross_name, 16, nchar(cross_name) - 2)
    ),
      as.character(construct_line)
    )
  ))

raw_count_data <- raw_count_data %>%
  mutate(cross_line_name = ifelse(
    substr(
      cross_line_name,
      nchar(cross_line_name) - 1,
      nchar(cross_line_name)
    ) == background_library,
    substr(cross_line_name,
           1, nchar(cross_line_name) - 2),
    cross_line_name
  ))

#Pairing to gene name key-----
gene_key <- gene_key %>% 
  mutate(background_library = substr(vdrc_id, 
                                     nchar(vdrc_id) - 1, 
                                     nchar(vdrc_id))) %>% 
  mutate(background_library = ifelse(!background_library %in% libraries_incl,
                                      "KK", background_library))


cleaned_count_data <- raw_count_data %>% 
  select(background_library, background_line, driver, construct_line, sex_test, count_date, bio_replicate, starts_with("offspring"), food_abnormality, parental_mortality, driver_present, rnai_construct_present) %>% 
  left_join(gene_key,
            by = c("background_library", "construct_line"))
  

#Computer-friendly variables-----
cleaned_count_data <- cleaned_count_data %>% 
  mutate(count_date = as.Date(count_date, format = "%m/%d/%Y")) %>% 
  mutate(count_month = month(count_date)) %>% 
  mutate(sex_test = ifelse(sex_test %in% c("808", "60100"), 
                           NA, tolower(sex_test))) %>%
  mutate(parental_mortality = !is.na(parental_mortality),
         food_abnormality = !is.na(food_abnormality),
         month_id = as.numeric(as.factor(as.character(count_month))),
         driver_id = as.numeric(as.factor(driver)),
         sex_id = as.numeric(as.factor(sex_test)),
         background_id = as.numeric(as.factor(background_library)),
         uniq_construct_id = as.numeric(as.factor(as.character(construct_line)))) %>% 
  mutate(driver_id = ifelse(is.na(driver_id) & rnai_construct_present, max(driver_id,na.rm=T)+1, driver_id))
  

#Final input data -----
processed_data <- cleaned_count_data %>% 
  filter(!is.na(offspring_total)) %>%
  select(background_id, driver_id, uniq_construct_id,
         sex_id, parental_mortality, month_id, 
         driver_present, rnai_construct_present,
         offspring_total)

#summary table to write-----
cleaned_count_data %>%
  filter(!is.na(offspring_total)) %>%
  group_by(background_library, background_line,
           driver, construct_line, gene_name, vdrc_id, 
           sex_test, driver_present, 
           rnai_construct_present) %>% 
  summarise(replicates = n(), net_offspring_produced = sum(offspring_total)) %>% 
  ungroup() %>% 
  write_csv("data/data_cleaned/offspring_summary.csv")
