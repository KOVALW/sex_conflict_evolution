source("data/processing_code/stan_formatting.R")

require(rstan)
require(loo)

options(mc.cores = parallel::detectCores())#,
#        stanc.allow_optimizations = TRUE,
#        stanc.auto_format = TRUE) #run cores in parallel

rstan_options(auto_write = TRUE) #auto-write compiled code to hard drive

model_distro <- "poisson"

model_output <- stan(file = paste0("./models/zero_hurdle_",model_distro,".stan"), 
     data = input_count_data, 
     iter = 2000)

if(data_set == testing_data){
  fp_model <- paste0("./models/",model_distro,"_testing_fit.rds")
} else if (data_set == nosgal4_data) {
  fp_model <- paste0("./models/",model_distro,"_nosgal4_fit.rds")
} else {
  fp_model <- paste0("./models/",model_distro,"_fit.rds")
}

saveRDS(model_output, file = )
