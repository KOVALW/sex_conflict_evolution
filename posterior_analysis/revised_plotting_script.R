source("./data/processing_code/data_cleaning.R")

#Data read-in----------------
#Best models determined from looic tables in ./posterior_analysis
soma_best_model <- 2 #mortality affects both the mean and zero hurdle probability
nosgal4_best_model <- 3 #mortality affects both the mean and zero hurdle probability, month-level effects included

#Treatment confounding factors and means
beta_vals <-
  bind_rows(
    read_csv("./posterior_analysis/beta_summaries_germline.csv") %>%
      filter(model_index == nosgal4_best_model) %>% 
      mutate(driver = ifelse(driver_id == 1, "Nos-GAL4", NA)),
    read_csv("./posterior_analysis/beta_summaries_soma.csv") %>%
      filter(model_index == soma_best_model) %>% 
      mutate(driver = ifelse(driver_id == 1, "Act5c",
                             ifelse(driver_id == 2, "Tubp-GAL4",
                                    NA)))
  )

#Base parameters
base_summaries <- bind_rows(
  read_csv("./posterior_analysis/base_par_summaries_germline.csv") %>%
    filter(model_index == nosgal4_best_model),
  read_csv("./posterior_analysis/base_par_summaries_soma.csv") %>%
    filter(model_index == soma_best_model)
)

#Treatment level effects
eta_summaries <- bind_rows(
  read_csv("./posterior_analysis/eta_summaries_germline.csv") %>%
    filter(model_index == nosgal4_best_model),
  read_csv("./posterior_analysis/eta_summaries_soma.csv") %>%
    filter(model_index == soma_best_model)
)

#Posterior manipulation-----
#Adding treatment mean effects posteriors as columns
eta_tmt_fx <- beta_vals %>%
  filter(!is.na(background_id) & !is.na(driver_id) & !is.na(sex_id)) %>% 
  rename(cross_mean = mean,
         cross_sd = sd,
         cross_lwr = `2.5%`,
         cross_upr = `97.5%`) %>% 
  select(starts_with("cross"), driver, background_id, sex_id) %>% 
  full_join(eta_summaries, by = c("driver", "background_id", "sex_id")) %>% 
  mutate(driver_shortname = substr(driver,1,3))


eta_sex_table <- eta_tmt_fx   %>%
  group_by(background_library, sex_test, gene_name, driver) %>%
  slice(1) %>% #removes duplicate rows due to multiple experiments of same sex/background/gene
  ungroup() %>%
  select(
    sex_test,
    background_library,
    driver_shortname,
    gene_name,
    gene_age,
    gene_origin,
    chromosome,
    starts_with("cross"),
    lwr_2.5 = `2.5%`,
    median = `50%`,
    mean,
    sd,
    upr_97.5 = `97.5%`
  )  %>%
  gather(
    -sex_test,-background_library,-driver_shortname,-gene_name,-gene_age,-gene_origin,
    -chromosome,
    key = "statistic",
    value = "estimate"
  ) %>%
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate)

eta_sex_table <- eta_sex_table %>% 
  mutate(male_beneficial = male_mean + male_cross_mean - 2 * male_sd > 0 & male_mean + male_cross_mean > 0,
         female_beneficial = female_mean + female_cross_mean - 2 * female_sd > 0 & female_mean + female_cross_mean > 0,
         male_detrimental = male_mean + male_cross_mean + 2 * male_sd < 0 & male_mean + male_cross_mean < 0,
         female_detrimental = female_mean + female_cross_mean + 2 * female_sd < 0 & female_mean + female_cross_mean < 0) %>% 
  mutate(signif = male_beneficial | male_detrimental | female_beneficial | female_detrimental ) 

#Writing----
write_csv(eta_sex_table, "./posterior_analysis/effect_size_summary_alldrivers.csv")

#Plots----
eta_sex_table %>% 
  ggplot(aes(y = exp(male_mean + male_cross_mean),
             x = exp(female_mean + female_cross_mean),
             col = signif)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("light grey",1)) +
  facet_wrap(driver_shortname~.) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  theme_bw() +
  labs(y = "Male KD relative fitness", x = "Female KD relative fitness")

eta_sex_table %>%
  filter(signif & chromosome != "c") %>% 
  ggplot(aes(y = exp(male_mean + male_cross_mean),
             x = exp(female_mean + female_cross_mean))) +
  geom_point(size = 3) +
  facet_grid(driver_shortname~chromosome) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  theme_bw() +
  labs(y = "Male KD relative fitness", x = "Female KD relative fitness")

#comparison of means expected and observed in data----
tmt_fx_post_means <- base_summaries %>% 
  filter(par_name == "base_reproduction") %>% 
  select(mean, sd, model_index) %>% 
  rename(base_mean = mean,
         base_sd = sd) %>% 
  left_join(beta_vals %>% 
              select(model_index, mean, sd, beta_id,
                     background_id, sex_id, driver_id,
                     construct_id, time_id, month_id,
                     beta_val_id, parental_mortality, 
                     driver))

germline <- T
source("data/processing_code/gal4_format_header.R")


mapped_fx <- eta_fx %>% 
  mutate(model_index = nosgal4_best_model) %>% 
  filter(offspring_total > 0) %>% 
  left_join(eta_summaries %>% 
              select(model_index, rnai_fx = mean, 
                     rnai_sd = sd, 
                     eta_fctr) %>% 
              group_by(eta_fctr) %>% 
              slice(1) %>% 
              ungroup()) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(month_id)) %>% 
              select(month_id, model_index, 
                     month_fx = mean,
                     month_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(time_id)) %>% 
              select(time_id, model_index, 
                     time_fx = mean,
                     time_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(background_id) & is.na(driver_id)) %>% 
              select(background_id, model_index, 
                     background_fx = mean,
                     background_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(driver_id)) %>% 
              select(driver_id, sex_id, background_id,
                     model_index, 
                     driver_fx = mean,
                     driver_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(parental_mortality)) %>% 
              select(parental_mortality, 
                     model_index, 
                     parental_mortality_fx = mean,
                     parental_mortality_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(construct_id)) %>% 
              select(uniq_construct_id = construct_id, model_index, 
                     construct_fx = mean,
                     construct_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(model_index == nosgal4_best_model) %>% 
              select(model_index,
                     base_mean,
                     base_sd) %>% 
              slice(1))

means_plot <- mapped_fx %>% 
  mutate(observation = row_number()) %>% 
  select(-ends_with("_sd"), -i) %>% 
  gather(-background_id, -driver_id, -uniq_construct_id,
         -sex_id, -parental_mortality, -time_id,
         -month_id, -driver_present, -rnai_construct_present,
         -offspring_total, -eta_fctr, -observation,
         -model_index, key = "tmt", value = "fx") %>% 
  mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
  spread(key = tmt, value = fx) %>% 
  mutate(mean_overall = rnai_fx + month_fx + time_fx + background_fx +
                            base_mean + driver_fx + parental_mortality_fx +
                            construct_fx) %>% 
  group_by(background_id, driver_id, uniq_construct_id,
           sex_id, parental_mortality, time_id,
           month_id, driver_present, rnai_construct_present,
           eta_fctr, 
           model_index, mean_overall) %>% 
  summarise(mean_obs = mean(offspring_total),
            experiments = n()) %>% 
  ungroup() %>% 
  mutate(model_mean = exp(mean_overall)) %>% 
  ggplot(aes(y = mean_obs, x = model_mean, col = parental_mortality)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(1,2)) +
  geom_abline(slope = 1, intercept = 0, col = "white",lwd = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2,lwd = 2) +
  labs(y = "Observed experimental group mean", x = "Model predicted mean",
       col = "Parental mortality")

mapped_fx %>% 
  mutate(observation = row_number()) %>% 
  select(-ends_with("_sd"), -i) %>% 
  gather(-background_id, -driver_id, -uniq_construct_id,
         -sex_id, -parental_mortality, -time_id,
         -month_id, -driver_present, -rnai_construct_present,
         -offspring_total, -eta_fctr, -observation,
         -model_index, key = "tmt", value = "fx") %>% 
  mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
  spread(key = tmt, value = fx) %>% 
  mutate(mean_overall = rnai_fx + month_fx + time_fx + background_fx +
           base_mean + driver_fx + parental_mortality_fx +
           construct_fx) %>% 
  group_by(background_id, driver_id, uniq_construct_id,
           sex_id, parental_mortality, time_id,
           month_id, driver_present, rnai_construct_present,
           eta_fctr, rnai_fx, month_fx, time_fx, background_fx,
             base_mean, driver_fx, parental_mortality_fx,
             construct_fx,
           model_index, mean_overall) %>% 
  summarise(mean_obs = mean(offspring_total),
            experiments = n()) %>% 
  ungroup() %>% 
  mutate(model_mean = exp(mean_overall)) %>% 
  ggplot(aes(x = mean_obs, y = log(mean_obs) - mean_overall, col = parental_mortality)) + 
  geom_point() +
  theme_bw() +
  labs(x = "Observed experimental group mean", y = "Log mean residual",
       alpha = "Parental mortlaity")

sex_means_plot <- mapped_fx %>% 
  mutate(observation = row_number()) %>% 
  select(-ends_with("_sd"), -i) %>% 
  gather(-background_id, -driver_id, -uniq_construct_id,
         -sex_id, -parental_mortality, -time_id,
         -month_id, -driver_present, -rnai_construct_present,
         -offspring_total, -eta_fctr, -observation,
         -model_index, key = "tmt", value = "fx") %>% 
  mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
  spread(key = tmt, value = fx) %>% 
  mutate(mean_overall = rnai_fx + month_fx + time_fx + background_fx +
           base_mean + driver_fx + parental_mortality_fx +
           construct_fx) %>% 
  group_by(background_id, driver_id, uniq_construct_id,
           sex_id, parental_mortality, time_id,
           month_id, driver_present, rnai_construct_present,
           eta_fctr, 
           model_index, mean_overall) %>% 
  summarise(mean_obs = mean(offspring_total),
            experiments = n()) %>% 
  ungroup() %>% 
  mutate(model_mean = exp(mean_overall)) %>% 
  right_join(data.frame(sex_test = c("female", "male"), sex_id = 1:2)) %>% 
  ggplot(aes(y = mean_obs, x = model_mean, col = parental_mortality)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(1,2)) +
  geom_abline(slope = 1, intercept = 0, col = "white",lwd = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2,lwd = 2) +
  labs(y = "Observed experimental group mean", x = "Model predicted mean",
       col = "Parental mortality") + 
  facet_wrap(sex_test~.)

all_data_plot <- mapped_fx %>% 
  mutate(observation = row_number()) %>% 
  select(-ends_with("_sd"), -i) %>% 
  gather(-background_id, -driver_id, -uniq_construct_id,
         -sex_id, -parental_mortality, -time_id,
         -month_id, -driver_present, -rnai_construct_present,
         -offspring_total, -eta_fctr, -observation,
         -model_index, key = "tmt", value = "fx") %>% 
  mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
  spread(key = tmt, value = fx) %>% 
  mutate(mean_overall = rnai_fx + month_fx + time_fx + background_fx +
           base_mean + driver_fx + parental_mortality_fx +
           construct_fx) %>% 
  mutate(model_mean = exp(mean_overall)) %>% 
  ggplot(aes(x = model_mean, y = offspring_total, col = parental_mortality)) + 
  geom_point(alpha = 0.7) +
  theme_bw() +
  scale_color_manual(values = c(1,2)) +
  geom_abline(slope = 1, intercept = 0, col = "white", lwd = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 2) +
  labs(x = "Model predicted mean", y = "Observed offspring count",
       col = "Parental mortality")

gridExtra::grid.arrange(all_data_plot + theme(legend.position = "none"), 
                        means_plot,
                        widths = c(0.45,0.55))

germline <- F
source("data/processing_code/gal4_format_header.R")


mapped_fx_soma <- eta_fx %>% 
  mutate(model_index = soma_best_model) %>% 
  filter(offspring_total > 0) %>% 
  left_join(eta_summaries %>% 
              select(model_index, rnai_fx = mean, 
                     rnai_sd = sd, 
                     eta_fctr) %>% 
              group_by(eta_fctr) %>% 
              slice(1) %>% 
              ungroup()) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(background_id) & is.na(driver_id)) %>% 
              select(background_id, model_index, 
                     background_fx = mean,
                     background_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(driver_id)) %>% 
              select(driver_id, sex_id, background_id,
                     model_index, 
                     driver_fx = mean,
                     driver_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(parental_mortality)) %>% 
              select(parental_mortality, 
                     model_index, 
                     parental_mortality_fx = mean,
                     parental_mortality_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(!is.na(construct_id)) %>% 
              select(uniq_construct_id = construct_id, model_index, 
                     construct_fx = mean,
                     construct_sd = sd)) %>% 
  left_join(tmt_fx_post_means %>% 
              filter(model_index == soma_best_model) %>% 
              select(model_index,
                     base_mean,
                     base_sd) %>% 
              slice(1))

means_soma_plot <-  mapped_fx_soma %>% 
  mutate(observation = row_number()) %>% 
  select(-ends_with("_sd"), -i) %>% 
  gather(-background_id, -driver_id, -uniq_construct_id,
         -sex_id, -parental_mortality, -driver_present, -rnai_construct_present,
         -offspring_total, -eta_fctr, -observation,
         -model_index, key = "tmt", value = "fx") %>% 
  mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
  spread(key = tmt, value = fx) %>% 
  mutate(mean_overall = rnai_fx + background_fx +
           base_mean + driver_fx + parental_mortality_fx +
           construct_fx) %>% 
  group_by(background_id, driver_id, uniq_construct_id,
           sex_id, parental_mortality, driver_present, rnai_construct_present,
           eta_fctr, 
           model_index, mean_overall) %>% 
  summarise(mean_obs = mean(offspring_total),
            experiments = n()) %>% 
  ungroup() %>% 
  mutate(model_mean = exp(mean_overall)) %>% 
  ggplot(aes(y = mean_obs, x = model_mean, col = parental_mortality)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(1,2)) +
  geom_abline(slope = 1, intercept = 0, col = "white",lwd = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2,lwd = 2) +
  labs(y = "Observed experimental group mean", x = "Model predicted mean",
       col = "Parental mortality")

sex_means_soma_plot <- mapped_fx_soma %>% 
  mutate(observation = row_number()) %>% 
  select(-ends_with("_sd"), -i) %>% 
  gather(-background_id, -driver_id, -uniq_construct_id,
         -sex_id, -parental_mortality, -driver_present, -rnai_construct_present,
         -offspring_total, -eta_fctr, -observation,
         -model_index, key = "tmt", value = "fx") %>% 
  mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
  spread(key = tmt, value = fx) %>% 
  mutate(mean_overall = rnai_fx + background_fx +
           base_mean + driver_fx + parental_mortality_fx +
           construct_fx) %>% 
  group_by(background_id, driver_id, uniq_construct_id,
           sex_id, parental_mortality, driver_present, rnai_construct_present,
           eta_fctr, 
           model_index, mean_overall) %>% 
  summarise(mean_obs = mean(offspring_total),
            experiments = n()) %>% 
  ungroup() %>% 
  left_join(data.frame(sex_test = c("female", "male"),
                       sex_id = 1:2)) %>% 
  mutate(model_mean = exp(mean_overall)) %>% 
  ggplot(aes(y = mean_obs, x = model_mean, col = parental_mortality)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c(1,2)) +
  geom_abline(slope = 1, intercept = 0, col = "white",lwd = 2) +
  geom_abline(slope = 1, intercept = 0, lty = 2,lwd = 2) +
  labs(y = "Observed experimental group mean", x = "Model predicted mean",
       col = "Parental mortality") +
  facet_wrap(sex_test~.)


all_data_soma_plot <-  mapped_fx_soma %>% 
   mutate(observation = row_number()) %>% 
   select(-ends_with("_sd"), -i) %>% 
   gather(-background_id, -driver_id, -uniq_construct_id,
          -sex_id, -parental_mortality, -driver_present, -rnai_construct_present,
          -offspring_total, -eta_fctr, -observation,
          -model_index, key = "tmt", value = "fx") %>% 
   mutate(fx = ifelse(is.na(fx), 0, fx)) %>% 
   spread(key = tmt, value = fx) %>% 
   mutate(mean_overall = rnai_fx + background_fx +
            base_mean + driver_fx + parental_mortality_fx +
            construct_fx) %>% 
   mutate(model_mean = exp(mean_overall)) %>% 
   ggplot(aes(x = model_mean, y = offspring_total, col = parental_mortality)) + 
   geom_point(alpha = 0.7) +
   theme_bw() +
   scale_color_manual(values = c(1,2)) +
   geom_abline(slope = 1, intercept = 0, col = "white", lwd = 2) +
   geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 2) +
   labs(x = "Model predicted mean", y = "Observed offspring count",
        col = "Parental mortality")
 
 gridExtra::grid.arrange(all_data_plot + theme(legend.position = "none"),
                         means_plot,
                         all_data_soma_plot + theme(legend.position = "none"),
                         means_soma_plot, ncol = 2,
                         widths = c(0.45,0.55))
 
 gridExtra::grid.arrange(sex_means_plot +
                           ggtitle("Nos-GAL4"),
                         sex_means_soma_plot +
                           ggtitle("Somatic drivers Act5c, Tubp"))
 