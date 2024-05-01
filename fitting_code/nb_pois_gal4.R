source("data/processing_code/stan_formatting.R")

require(rstan)
require(loo)

options(mc.cores = parallel::detectCores())#,
#        stanc.allow_optimizations = TRUE,
#        stanc.auto_format = TRUE) #run cores in parallel

rstan_options(auto_write = TRUE) #auto-write compiled code to hard drive

fitting <- F

if (fitting){
file.remove("./models/zh_pois_thetafix.rds")
lik_output_pois <-
  stan(
    file = "./models/zh_pois_thetafix.stan",
    data = input_count_data_stan,
    iter = 3000
  )

saveRDS(lik_output_pois, file = "./models/zh_pois_thetafix_output.rds")

file.remove("./models/zh_pois_thetafix.rds")
lik_output_negbinom <-
  stan(
    file = paste0("./models/zh_negbinom_thetafix.stan"),
    data = input_count_data_stan,
    iter = 3000
  )

saveRDS(lik_output_negbinom, file = "./models/zh_negbinom_thetafix_output.rds")

} else {
  lik_output_negbinom <- readRDS("./models/zh_negbinom_thetafix_output.rds")
  lik_output_pois <- readRDS("./models/zh_pois_thetafix_output.rds")
}

loo_nb <- loo(lik_output_negbinom, cores = 1)
loo_pois <- loo(lik_output_pois, cores = 1)

loo_compare(loo_pois, loo_nb)

#-----
current_fit <- lik_output_negbinom

data_key_postfit <- nosgal4_data  %>% 
  mutate(i = ifelse(rnai_construct_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (background_id - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i))) %>% 
  select(-i) %>% 
  group_by(background_library, background_line, driver, 
           construct_line, sex_test, time_id,
           driver_present,rnai_construct_present, driver_id, sex_id, background_id, 
           uniq_construct_id, eta_fctr) %>% 
  summarise(experiment_count = n(),
            offspring_produced = sum(offspring_total)) %>% 
  ungroup() %>% left_join(gene_key %>% select(-background_library))

parfit_summaries <- as.data.frame(summary(current_fit)$summary)

parfit_summaries <- parfit_summaries %>% 
  mutate(par_name = rownames(parfit_summaries))

parfit_summaries <- parfit_summaries %>% 
  filter(grepl("etas",par_name)) %>% 
  mutate(eta_fctr = as.numeric(substr(par_name, 6, nchar(par_name)-1)) + 1) %>% 
  left_join(data_key_postfit)

parfit_summaries %>% 
  ggplot(aes(x = `50%`, xmin = `2.5%`, xmax = `97.5%`,y = eta_fctr,col = driver)) +
  geom_pointrange() +
  theme_bw()

summary_sex_fx <- parfit_summaries %>% 
  select(sex_test, background_library, driver, construct_line, time_id, gene_name, gene_age, gene_origin, `2.5%`, `50%`, mean, `97.5%`, se_mean, sd, `25%`, `75%`) %>% 
  gather(-sex_test, -background_library, -driver, -gene_name, -time_id, -construct_line, -gene_age, -gene_origin,
         key = "statistic", value = "estimate") %>% 
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate)

summary_sex_fx %>% 
  ggplot(aes(x = exp(-`male_mean`), y = exp(-`female_mean`),
             xmax = exp(-`male_2.5%`), ymax = exp(-`female_2.5%`),
             xmin = exp(-`male_97.5%`), ymin = exp(-`female_97.5%`),
             col = background_library)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "Male germline gene fitness benefit (WT / KD)", y = "Female germline gene fitness benefit (WT / KD)") +
  scale_color_manual(values = c(1,4))

summary_sex_fx %>% 
  ggplot(aes(x = exp(-`male_mean`), y = exp(-`female_mean`),
             xmax = exp(-`male_25%`), ymax = exp(-`female_25%`),
             xmin = exp(-`male_75%`), ymin = exp(-`female_75%`),
             col = background_library)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "Male germline gene fitness benefit (WT / KD)", y = "Female germline gene fitness benefit (WT / KD)") +
  scale_color_manual(values = c(1,4))

summary_sex_fx %>% 
  ggplot(aes(x = exp(-`male_mean`), y = exp(-`female_mean`),
             xmax = exp(-(male_mean - male_se_mean)), ymax = exp(-(female_mean - female_se_mean)),
             xmin = exp(-(male_mean + male_se_mean)), ymin = exp(-(female_mean + female_se_mean)),
             col = background_library)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "Male germline gene fitness benefit (WT / KD)", y = "Female germline gene fitness benefit (WT / KD)") +
  scale_color_manual(values = c(1,4))
