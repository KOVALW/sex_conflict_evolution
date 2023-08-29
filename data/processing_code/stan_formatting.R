source("data/processing_code/data_cleaning.R")

#Script to prepare stan inputs from cleaned data


#Prepping beta key-----
lines <- length(unique(processed_data$background_id))
drivers <- length(unique(processed_data$driver_id)) - 1 #NA drivers are counted as a "driver" only when construct is present but both should be removed in beta_fx
sexes <- length(unique(processed_data$sex_id)) - 1
genes <- length(unique(processed_data$uniq_construct_id)) - 1 #to exclude NA
timepoints <- length(unique(processed_data$month_id)) - 1
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
  exp_column <- timepoints + lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes
  submit_fx <- processed_data$driver_present[i]
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

eta_fx <- processed_data %>% 
  mutate(i = ifelse(rnai_construct_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (ln - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i))) 

#Formatted sstan data----
input_count_data <- list(
  N = length(which(processed_data$offspring_total > 0)),
  NZ = length(which(processed_data$offspring_total == 0)),
  b_coefs = ncol(tmt_spec_fx),
  eta_coefs = length(unique(eta_fx$eta_fctr)) - 1,
  offspring = processed_data[which(processed_data$offspring_total > 0),]$offspring_total,
  eta_id = eta_fx[which(eta_fx$offspring_total > 0),]$eta_fctr,
  tmt_spec = tmt_spec_fx[which(processed_data$offspring_total > 0),]
)

