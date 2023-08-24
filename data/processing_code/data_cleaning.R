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
