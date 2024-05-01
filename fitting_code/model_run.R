source("data/processing_code/stan_formatting.R")

require(rstan)
require(loo)

options(mc.cores = parallel::detectCores())#,
#        stanc.allow_optimizations = TRUE,
#        stanc.auto_format = TRUE) #run cores in parallel

rstan_options(auto_write = TRUE) #auto-write compiled code to hard drive

model_distro <- "poisson"
theta_fx_bio <- T

if(!theta_fx_bio) {
  model_output <-
    stan(
      file = paste0("./models/zero_hurdle_", model_distro, "_thetafix.stan"),
      data = input_count_data_theta_onefx,
      iter = 2000
    )
  
  if (identical(data_set, testing_data)) {
    fp_model <- paste0("./models/", model_distro, "_testing_fit.rds")
  } else if (identical(data_set, nosgal4_data)) {
    fp_model <- paste0("./models/", model_distro, "_nosgal4_fit.rds")
  } else if (identical(data_set, gschwend_data)) {
    fp_model <- paste0("./models/", model_distro, "_act5tubp_fit.rds")
  } else {
    fp_model <- paste0("./models/", model_distro, "_fit.rds")
  }
} else {
  model_output <-
    stan(
      file = paste0("./models/zero_hurdle_", model_distro, "_thetabio.stan"),
      data = input_count_data_thetafx,
      iter = 2000
    )
  
  if (identical(data_set, testing_data)) {
    fp_model <-
      paste0("./models/", model_distro, "_thetabio_testing_fit.rds")
  } else if (identical(data_set, nosgal4_data)) {
    fp_model <-
      paste0("./models/", model_distro, "_thetabio_nosgal4_fit.rds")
  } else if (identical(data_set, gschwend_data)) {
    fp_model <-
      paste0("./models/", model_distro, "_thetabio_act5tubp_fit.rds")
  } else {
    fp_model <- paste0("./models/", model_distro, "_thetabio_fit.rds")
  }
}

saveRDS(model_output, file = fp_model)
