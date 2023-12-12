source("data/processing_code/data_cleaning.R")

nosgal4_data <- cleaned_count_data %>% 
  filter(!is.na(offspring_total) & !is.na(count_date))  %>%
  mutate(time_id = as.numeric(as.factor(as.character(count_date))),
         driver_id = as.numeric(as.factor(driver)),
         sex_id = as.numeric(as.factor(sex_test)),
         background_id = as.numeric(as.factor(background_library)),
         uniq_construct_id = as.numeric(as.factor(as.character(construct_line)))) 

#Making stan fit data set and counting parameters------

data_set <- nosgal4_data %>% 
  select(background_id, driver_id, uniq_construct_id,
         sex_id, parental_mortality, time_id, 
         driver_present, rnai_construct_present,
         offspring_total)

lines <- length(unique(data_set$background_id))
drivers <- length(unique(data_set$driver_id)) - 1 #NA drivers should not be counted as a "driver" when construct is present because we consider this a separate additive effect
sexes <- length(unique(data_set$sex_id)) - 1
genes <- length(unique(data_set$uniq_construct_id)) - 1 #to exclude NA
timepoints <- length(unique(data_set$time_id))

datapoints <- nrow(data_set)

construct_only <- data_set %>% 
  filter(!driver_present & rnai_construct_present) %>% 
  group_by(uniq_construct_id) %>% 
  summarise(o=n()) %>% 
  ungroup() %>% 
  select(uniq_construct_id) %>% 
  as.matrix() %>% 
  as.numeric()

#Binary effect matrix for multiplying by effect size vector in regression model-----
tmt_spec_fx <- matrix(0, 
                      nrow = datapoints, 
                      ncol = lines + (drivers*sexes*lines) + length(construct_only) + timepoints)

for(i in 1:datapoints){
  s <- data_set[i,]$sex_id 
  d <- data_set[i,]$driver_id
  ln <- data_set[i,]$background_id
  con <- data_set[i,]$uniq_construct_id
  tm <- data_set[i,]$time_id
  
  tmt_spec_fx[i, ln] <- 1
  exp_column <- lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes
  submit_fx <- data_set$driver_present[i]
  if(!is.na(submit_fx + exp_column) & submit_fx)
    tmt_spec_fx[i, exp_column] <- 1
  if(!submit_fx & con %in% construct_only)
    tmt_spec_fx[i, which(construct_only == con) + lines + (drivers*sexes*lines)] <- 1
  
  tmt_spec_fx[i, tm + lines + length(construct_only) + (drivers*sexes*lines)] <- 1
}

#Simplified effect key from matrix-----
#Previously we geenrated an empty data frame that hat timepoints included. We no longer consider time effects on offspring due to conflation with earlier reproductive controls. Therefore we just make an empty data.frame
beta_fx_key <- data.frame()

for(ln in 1:lines){
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(background_id = ln,
                                      sex_id = NA,
                                      driver_id = NA,
                                      construct_id = NA,
                                      time_id = NA,
                                      beta_id = ln))
  for(d in 1:drivers) {
    for (s in 1:sexes) {
      
      beta_fx_key <-bind_rows(beta_fx_key, 
                              data.frame(background_id = ln,
                                         sex_id = s,
                                         driver_id = d,
                                         construct_id = NA,
                                         time_id = NA,
                                         beta_id = lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes))
      
      
    }
  }
}


for(con in 1:length(construct_only)) {
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(background_id = NA,
                                      sex_id = NA,
                                      driver_id = NA,
                                      construct_id = construct_only[con],
                                      time_id = NA,
                                      beta_id = con + lines + (drivers*sexes*lines)))
}

for (tm in 1:timepoints) {
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(background_id = NA,
                                      sex_id = NA,
                                      driver_id = NA,
                                      construct_id = construct_only[con],
                                      time_id = tm,
                                      beta_id = tm + length(construct_only) + lines + (drivers*sexes*lines)))
}

beta_fx_key <- beta_fx_key %>% 
  arrange(beta_id) %>% 
  mutate(experiment_count = colSums(tmt_spec_fx)) %>% 
  filter(experiment_count > 0) %>% 
  mutate(beta_val_id = row_number())

#This should be identical to the tmt_spec_fx original values. Non-identical means some driver x background x sex effects aren't present
tmt_spec_fx <- tmt_spec_fx[,which(colSums(tmt_spec_fx)>0)]

#Intercept values for each treatment effect-----
eta_fx <- data_set %>% 
  mutate(i = ifelse(rnai_construct_present == 1 & driver_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (background_id - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i)))


#Formatted stan data----
input_count_data_stan <- list(
  N = length(data_set$offspring_total),
  b_coefs = ncol(tmt_spec_fx),
  eta_coefs = length(unique(eta_fx$eta_fctr)) - 1, #1 in this case is no tmt intercept effect, fixed 0 value in eta_coefs1[b_coefs+1] and referred to appropriately through eta_id
  offspring = data_set$offspring_total,
  eta_id = eta_fx$eta_fctr,
  mortality = as.integer(data_set$parental_mortality),
  zero_bool = as.integer(data_set$offspring_total == 0),
  tmt_spec = tmt_spec_fx
)
