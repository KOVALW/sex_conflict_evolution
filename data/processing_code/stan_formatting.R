source("data/processing_code/data_cleaning.R")

#Script to prepare stan inputs from cleaned data


#Prepping beta key-----
lines <- length(unique(processed_data$background_id))
drivers <- length(unique(processed_data$driver_id)) #NA drivers are counted as a "driver"
sexes <- length(unique(processed_data$sex_id))
genes <- length(unique(processed_data$uniq_construct_id)) - 1 #to exclude NA
timepoints <- length(unique(processed_data$month_id))
datapoints <- nrow(processed_data)

tmt_spec_fx <- matrix(0, 
                      nrow = datapoints, 
                      ncol = lines + timepoints + (drivers*sexes*lines))

for(i in 1:datapoints){
  s <- processed_data[i,]$sex_id 
  d <- processed_data[i,]$driver_id
  ln <- processed_data[i,]$background_id
  m <- processed_data[i,]$month_id
  
  tmt_spec_fx[i,m] <- 1
  
  tmt_spec_fx[i,timepoints + ln] <- 1
  
  if(processed_data$driver_present[i])
    tmt_spec_fx[i, timepoints + lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes] <- 1
}

beta_value_key <- tmt_spec_fx[,(timepoints+1):ncol(tmt_spec_fx)] %>% 
  as.data.frame()

tmt_spec_fx <- tmt_spec_fx[,which(colSums(tmt_spec_fx)>0)]

eta_fx <- processed_data %>% 
  mutate(i = ifelse(rnai_construct_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (ln - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i))) 

input_count_data <- list(
  N = length(which(processed_data$offspring_total > 0)),
  NZ = length(which(processed_data$offspring_total == 0)),
  b_coefs = ncol(tmt_spec_fx),
  offspring = processed_data[which(processed_data$offspring_total > 0),]$offspring_total,
  eta_id = eta_fx[which(eta_fx$offspring_total > 0),]$eta_fctr,
  tmt_spec = tmt_spec_fx[which(processed_data$offspring_total > 0),]
)

#effect key----
#This doesn;t match up with current data set!!!
beta_value_key <- beta_value_key %>% 
  bind_cols(processed_data %>% 
              select(background_num,driver_present, sex_id,driver_id)) %>% 
  mutate(sex_num = ifelse(is.na(sex_num),0,sex_num)) %>% 
  gather(-background_num,
         -driver_present, 
         -sex_num,
         -driver_num,
         key = "beta_id",
         value = "beta_presence") %>% 
  filter(beta_presence == 1) %>% 
  mutate(tmt_score = background_num + 
           driver_present + 
           sex_num + 
           driver_num) %>% 
  group_by(beta_id) %>% 
  filter(tmt_score == min(tmt_score)) %>% 
  mutate(obs = n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(sex_num = ifelse(driver_present, sex_num, NA),
         driver_num = ifelse(driver_present, driver_num, NA)) %>% 
  select(-beta_presence, -tmt_score) %>% 
  mutate(beta_val_num = as.numeric(substr(beta_id,2,nchar(beta_id)))) %>% 
  arrange(beta_val_num) %>% 
  bind_rows(data.frame(time_mo_stanfacing = 1:timepoints)) %>% 
  mutate(beta_val_reported = row_number())
