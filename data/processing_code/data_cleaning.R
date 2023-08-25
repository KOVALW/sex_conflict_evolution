library(tidyverse)
library(lubridate)

raw_count_data <- read_csv("./data/data_raw/final_count_data.csv", guess_max = 5e4)

#Making library and line consistent-----
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

#Computer-friendly variables-----
raw_count_data <- raw_count_data %>% 
  mutate(count_date = as.Date(count_date, format = "%m/%d/%Y")) %>% 
  mutate(count_month = month(count_date)) %>% 
  mutate(sex_test = ifelse(sex_test %in% c("808", "60100"), 
                           NA, tolower(sex_test))) %>%
  mutate(parental_mortality = !is.na(parental_mortality),
         food_abnormality = !is.na(food_abnormality),
         driver_id = as.numeric(as.factor(driver)),
         sex_id = as.numeric(as.factor(sex_test)),
         background_id = as.numeric(as.factor(background_library)))
  
  