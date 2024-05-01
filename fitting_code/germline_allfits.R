germline <- T
source("./data/processing_code/gal4_format_header.R")

require(rstan)
require(loo)

options(mc.cores = parallel::detectCores())

fps <- list.files("fitting_code/")
fps <- paste0("fitting_code/", fps[which(grepl("gal4_tmt", fps) & grepl(".R", fps))])

loo_list <- list()

for (j in 1:length(fps)) {
  fp <- fps[j]
  source(fp)
  
  #Formatted stan data----
  input_count_data_stan <- list(
    N = length(data_set$offspring_total),
    b_coefs = ncol(tmt_spec_fx),
    eta_coefs = length(unique(eta_fx$eta_fctr)) - 1,
    #1 in this case is no tmt intercept effect, fixed 0 value in eta_coefs1[b_coefs+1] and referred to appropriately through eta_id
    offspring = data_set$offspring_total,
    eta_id = eta_fx$eta_fctr,
    mortality = as.integer(data_set$parental_mortality),
    zero_bool = as.integer(data_set$offspring_total == 0),
    tmt_spec = tmt_spec_fx
  )
  
  fp_output <- paste0(substr(fp,1,nchar(fp)-2), "_output.rds")
  
  write_files <- file.exists(fp_output)
  #Run or load model
  if(write_files) {
    current_fit <- readRDS(fp_output)
  } else {
  current_fit <-
    stan(
      file = paste0("./models/zh_negbinom_thetafix.stan"),
      data = input_count_data_stan,
      iter = 3000
    )
  saveRDS(current_fit, fp_output)
  }
  
  data_key_postfit <- nosgal4_data  %>% 
    mutate(i = ifelse(rnai_construct_present == 1 & driver_present == 1, 
                      uniq_construct_id + (sex_id-1)*genes + 
                        (driver_id-1) * genes * sexes + (background_id - 1) * drivers * sexes * genes, 
                      0)) %>%
    mutate(eta_fctr = as.numeric(factor(i))) %>% 
    select(-i) %>% 
    group_by(background_library, background_line, driver, 
             construct_line, sex_test, time_id, month_id,
             driver_present,rnai_construct_present, 
             driver_id, sex_id, background_id, 
             uniq_construct_id, eta_fctr) %>% 
    summarise(experiment_count = n(),
              offspring_produced = sum(offspring_total)) %>% 
    ungroup() %>% 
    left_join(gene_key %>% 
                select(-background_library))
  
  parfit_summaries <- as.data.frame(summary(current_fit)$summary)
  
  parfit_summaries <- parfit_summaries %>% 
    mutate(par_name = rownames(parfit_summaries))
  
  parfit_summaries %>% 
    filter(grepl("etas",par_name)) %>% 
    mutate(eta_fctr = as.numeric(substr(par_name, 6, nchar(par_name)-1)) + 1) %>% 
    left_join(data_key_postfit) %>% 
    mutate(model_index = j) %>% 
    write_csv("./posterior_analysis/eta_summaries_germline.csv", 
              append = file.exists("./posterior_analysis/eta_summaries_germline.csv") & j > 1)
  
  if (j > 1) {
    bind_rows(
      read_csv("./posterior_analysis/beta_summaries_germline.csv"),
      parfit_summaries %>%
        filter(grepl("vals", par_name)) %>%
        mutate(beta_id = as.numeric(substr(
          par_name, 11, nchar(par_name) - 1
        ))) %>%
        left_join(beta_fx_key) %>%
        mutate(model_index = j)
    ) %>%
      write_csv("./posterior_analysis/beta_summaries_germline.csv")
  } else {
    parfit_summaries %>%
      filter(grepl("vals", par_name)) %>%
      mutate(beta_id = as.numeric(substr(par_name, 11, nchar(par_name) - 1))) %>%
      left_join(beta_fx_key) %>%
      mutate(model_index = j) %>%
      write_csv("./posterior_analysis/beta_summaries_germline.csv")
  }
  
  parfit_summaries %>% 
    filter(par_name %in% c("base_reproduction", "mort_theta", "base_theta")) %>% 
    mutate(model_index = j) %>% 
    write_csv("./posterior_analysis/base_par_summaries_germline.csv", 
              append = file.exists("./posterior_analysis/base_par_summaries_germline.csv") & j > 1)
  
  loo_list[[j]] <- loo(current_fit, cores = 1)
  remove(current_fit)
}

looic_table <- as.data.frame(loo_compare(loo_list))

looic_table %>% 
  mutate(model_id = rownames(looic_table)) %>%
  mutate(model_name_fp = fps[as.numeric(substr(model_id, 6, nchar(model_id)))]) %>%
  write_csv("./posterior_analysis/looic_table_gal4_germline.csv")
