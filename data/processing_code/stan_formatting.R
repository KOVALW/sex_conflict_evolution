source("data/processing_code/data_cleaning.R")

#Script to prepare stan inputs from cleaned data
#Switch between prcessed data for full model fits and testing data for seeing model functionality
data_set <- nosgal4_data

#Prepping beta key-----
lines <- length(unique(data_set$background_id))
drivers <- length(unique(data_set$driver_id)) - 1 #NA drivers are counted as a "driver" only when construct is present but both should be removed in beta_fx
sexes <- length(unique(data_set$sex_id)) - 1
genes <- length(unique(data_set$uniq_construct_id)) - 1 #to exclude NA
timepoints <- length(unique(data_set$month_id)) - 1
datapoints <- nrow(data_set)

tmt_spec_fx <- matrix(0, 
                      nrow = datapoints, 
                      ncol = lines + timepoints + (drivers*sexes*lines))

for(i in 1:datapoints){
  s <- data_set[i,]$sex_id 
  d <- data_set[i,]$driver_id
  ln <- data_set[i,]$background_id
  m <- data_set[i,]$month_id
  
  tmt_spec_fx[i,m] <- 1
  
  tmt_spec_fx[i,timepoints + ln] <- 1
  exp_column <- timepoints + lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes
  submit_fx <- data_set$driver_present[i]
  if(!is.na(submit_fx + exp_column) & submit_fx)
    tmt_spec_fx[i, exp_column] <- 1
}

tmt_spec_fx <- tmt_spec_fx[,which(colSums(tmt_spec_fx)>0)]

#effect key----
beta_fx_key <- data.frame()
for(ln in 1:lines){
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(background_id = ln,
                                      sex_id = NA,
                                      driver_id = NA,
                                      beta_id = timepoints + ln))
  for(d in 1:(drivers-1)) {
    for (s in 1:sexes) {
      
      beta_fx_key <-bind_rows(beta_fx_key, 
                              data.frame(background_id = ln,
                                         sex_id = s,
                                         driver_id = d,
                                         beta_id = timepoints + lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes))
      
    }
  }
}

eta_fx <- data_set %>% 
  mutate(i = ifelse(rnai_construct_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (background_id - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i))) 

#Formatted stan data----
input_count_data <- list(
  N = length(which(data_set$offspring_total > 0)),
  NZ = length(which(data_set$offspring_total == 0)),
  b_coefs = ncol(tmt_spec_fx),
  eta_coefs = length(unique(eta_fx$eta_fctr)) - 1,
  offspring = data_set[which(data_set$offspring_total > 0),]$offspring_total,
  eta_id = eta_fx[which(eta_fx$offspring_total > 0),]$eta_fctr,
  tmt_spec = tmt_spec_fx[which(data_set$offspring_total > 0),]
)

theta_set <- data_set %>% 
  group_by(sex_id, driver_id, driver_present, rnai_construct_present) %>% 
  summarise(o=n()) %>% 
  ungroup() %>% 
  filter(rnai_construct_present & driver_present) %>% 
  mutate(theta_id = row_number() + 1) %>% 
  select(-o) %>% 
  full_join(data_set %>% 
              mutate(r = row_number())) %>% 
  mutate(theta_id = ifelse(is.na(theta_id),1,theta_id)) %>% 
  arrange(r)

input_count_data_thetafx <- list(
  N = length(which(data_set$offspring_total > 0)),
  NZincl = nrow(data_set),
  Ktheta = length(unique(theta_set$theta_id)) - 1,
  b_coefs = ncol(tmt_spec_fx),
  eta_coefs = length(unique(eta_fx$eta_fctr)) - 1,
  offspring = data_set[which(data_set$offspring_total > 0),]$offspring_total,
  eta_id = eta_fx[which(eta_fx$offspring_total > 0),]$eta_fctr,
  theta_id = theta_set$theta_id,
  mortality = as.integer(data_set$parental_mortality),
  zero_bool = as.integer(data_set$offspring_total == 0),
  tmt_spec = tmt_spec_fx[which(data_set$offspring_total > 0),]
)
