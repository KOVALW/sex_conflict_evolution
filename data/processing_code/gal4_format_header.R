source("data/processing_code/data_cleaning.R")

nosgal4_data <- cleaned_count_data %>% 
  filter(!is.na(offspring_total) & !is.na(count_date))  %>%
  mutate(time_id = as.numeric(as.factor(as.character(count_date))),
         month_id = as.numeric(as.factor(month(count_date))),
         driver_id = as.numeric(as.factor(driver)),
         sex_id = as.numeric(as.factor(sex_test)),
         background_id = as.numeric(as.factor(background_library)),
         uniq_construct_id = as.numeric(as.factor(as.character(construct_line)))) 

somatic_data <- cleaned_count_data %>% 
  filter(!is.na(offspring_total) & is.na(count_date))  %>%
  mutate(driver_id = as.numeric(as.factor(driver)),
         sex_id = as.numeric(as.factor(sex_test)),
         background_id = as.numeric(as.factor(background_library)),
         uniq_construct_id = as.numeric(as.factor(as.character(construct_line)))) %>% 
  select(-month_id)

#Making stan fit data set and counting parameters------

if (germline){
data_set <- nosgal4_data %>% 
  select(background_id, driver_id, uniq_construct_id,
         sex_id, parental_mortality, time_id, month_id,
         driver_present, rnai_construct_present,
         offspring_total)
sexes <- length(unique(data_set$sex_id)) - 1
} else {
data_set <- somatic_data %>% 
  select(background_id, driver_id, uniq_construct_id,
         sex_id, parental_mortality,
         driver_present, rnai_construct_present,
         offspring_total)
sexes <- length(unique(data_set$sex_id))
}
lines <- length(unique(data_set$background_id))
drivers <- length(unique(data_set$driver_id)) - 1 #NA drivers should not be counted as a "driver" when construct is present because we consider this a separate additive effect

genes <- length(unique(data_set$uniq_construct_id)) - 1 #to exclude NA
timepoints <- length(unique(data_set$time_id))
month_count <- length(unique(data_set$month_id))

datapoints <- nrow(data_set)

construct_only <- data_set %>% 
  filter(!driver_present & rnai_construct_present) %>% 
  group_by(uniq_construct_id) %>% 
  summarise(o=n()) %>% 
  ungroup() %>% 
  select(uniq_construct_id) %>% 
  as.matrix() %>% 
  as.numeric()

#Intercept values for each treatment effect-----
eta_fx <- data_set %>% 
  mutate(i = ifelse(rnai_construct_present == 1 & driver_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (background_id - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i)))
