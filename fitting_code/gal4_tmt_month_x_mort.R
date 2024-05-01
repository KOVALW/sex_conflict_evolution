#model to account for month effects
#Binary effect matrix for multiplying by effect size vector in regression model-----
tmt_spec_fx <- matrix(0, 
                      nrow = datapoints, 
                      ncol = lines + (drivers*sexes*lines) + length(construct_only) + timepoints + month_count + 1)

for(i in 1:datapoints){
  s <- data_set[i,]$sex_id 
  d <- data_set[i,]$driver_id
  ln <- data_set[i,]$background_id
  con <- data_set[i,]$uniq_construct_id
  tm <- data_set[i,]$time_id
  mo <- data_set[i,]$month_id
  
  tmt_spec_fx[i, ln] <- 1
  exp_column <- lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes
  submit_fx <- data_set$driver_present[i]
  if(!is.na(submit_fx + exp_column) & submit_fx)
    tmt_spec_fx[i, exp_column] <- 1
  if(!submit_fx & con %in% construct_only)
    tmt_spec_fx[i, which(construct_only == con) + lines + (drivers*sexes*lines)] <- 1
  
  tmt_spec_fx[i, tm + lines + length(construct_only) + (drivers*sexes*lines)] <- 1
  tmt_spec_fx[i, mo + timepoints + lines + length(construct_only) + (drivers*sexes*lines)] <- 1
  
  if(data_set[i,]$parental_mortality){
    tmt_spec_fx[i, month_count + timepoints + lines + length(construct_only) + (drivers*sexes*lines) + 1] <- 1
  }
}

#Simplified effect key from matrix-----
#Previously we geenrated an empty data frame that hat timepoints included. We no longer consider time effects on offspring due to conflation with earlier reproductive controls. Therefore we just make an empty data.frame
beta_fx_key <- data.frame()

for(ln in 1:lines){
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(background_id = ln,
                                      beta_id = ln))
  for(d in 1:drivers) {
    for (s in 1:sexes) {
      
      beta_fx_key <-bind_rows(beta_fx_key, 
                              data.frame(background_id = ln,
                                         sex_id = s,
                                         driver_id = d,
                                         beta_id = lines + d + (s-1)*drivers + (ln - 1)*drivers*sexes))
      
      
    }
  }
}


for(con in 1:length(construct_only)) {
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(construct_id = construct_only[con],
                                      beta_id = con + lines + (drivers*sexes*lines)))
}

for (tm in 1:timepoints) {
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(time_id = tm,
                                      beta_id = tm + length(construct_only) + lines + (drivers*sexes*lines)))
}

for (mo in 1:month_count) {
  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(month_id = mo,
                                      beta_id = mo + timepoints + length(construct_only) + lines + (drivers*sexes*lines)))
}

  beta_fx_key <- bind_rows(beta_fx_key,
                           data.frame(parental_mortality = T,
                                      beta_id = 1 + month_count + timepoints + length(construct_only) + lines + (drivers*sexes*lines)))
  
  
beta_fx_key <- beta_fx_key %>% 
  arrange(beta_id) %>% 
  mutate(experiment_count = colSums(tmt_spec_fx)) %>% 
  filter(experiment_count > 0) %>% 
  mutate(beta_val_id = row_number())

#This should be identical to the tmt_spec_fx original values. Non-identical means some driver x background x sex effects aren't present
tmt_spec_fx <- tmt_spec_fx[,which(colSums(tmt_spec_fx)>0)]
