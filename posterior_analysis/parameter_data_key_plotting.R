#fctr key
current_fit <- readRDS(fp_model)

data_key_postfit <- cleaned_count_data %>% 
  # filter(!is.na(offspring_total) & 
  #          (!rnai_construct_present | construct_line %in% nosgal4_test_lines$construct_line) & 
  #          (driver == "Nos-GAL4" | is.na(driver))) %>% 
 filter(!is.na(offspring_total) & !is.na(count_date)) %>%
  filter(!is.na(offspring_total)) %>% 
  mutate(parental_mortality = !is.na(parental_mortality),
         food_abnormality = !is.na(food_abnormality),
         month_id = as.numeric(as.factor(as.character(count_month))),
         driver_id = as.numeric(as.factor(driver)),
         sex_id = as.numeric(as.factor(sex_test)),
         background_id = as.numeric(as.factor(background_library)),
         uniq_construct_id = as.numeric(as.factor(as.character(construct_line)))) %>% 
  mutate(driver_id = ifelse(is.na(driver_id) & rnai_construct_present, max(driver_id,na.rm=T)+1, driver_id))  %>% 
  mutate(i = ifelse(rnai_construct_present == 1, 
                    uniq_construct_id + (sex_id-1)*genes + 
                      (driver_id-1) * genes * sexes + (background_id - 1) * drivers * sexes * genes, 
                    0)) %>%
  mutate(eta_fctr = as.numeric(factor(i))) %>% 
  select(-i) %>% 
  group_by(background_library, background_line, driver, 
           construct_line, sex_test, 
           driver_present,rnai_construct_present,
           gene_name, vdrc_id, gene_origin, 
           gene_age, driver_id, sex_id, background_id, 
           uniq_construct_id, eta_fctr) %>% 
  summarise(experiment_count = n(),
            offspring_produced = sum(offspring_total)) %>% 
  ungroup()

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

parfit_summaries %>% 
  select(sex_test, background_library, driver, gene_name, gene_age, gene_origin, `2.5%`, `50%`, mean, `97.5%`) %>% 
  gather(-sex_test, -background_library, -driver, -gene_name, -gene_age, -gene_origin,
         key = "statistic", value = "estimate") %>% 
  mutate(estimate = exp(-estimate)) %>% 
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) %>% 
  ggplot(aes(x = `male_mean`, y = `female_mean`,
             xmin = `male_2.5%`, ymin = `female_2.5%`,
             xmax = `male_97.5%`, ymax = `female_97.5%`,
             col = background_library)) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "Male germline gene fitness benefit (offspring / expected offspring)", y = "Female germline gene fitness benefit (offspring / expected offspring)") +
  geom_smooth(aes(group=background_library), method = "glm")  +
  scale_color_manual(values = c(1,4))+
  facet_wrap(background_library~driver,scales = "free")

parfit_summaries %>%
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    `2.5%`,
    `50%`,
    mean,
    `97.5%`
  ) %>%
  gather(
    -sex_test,
    -background_library,
    -driver,
    -gene_name,
    -gene_age,
    -gene_origin,
    key = "statistic",
    value = "estimate"
  ) %>%  
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) %>% 
  filter(!is.na(male_mean) & !is.na(female_mean) &
           `female_97.5%` < 7 & `male_97.5%` < 7) %>%
  ggplot(
    aes(
      x = `male_mean`,
      y = `female_mean`,
      xmax = `male_2.5%`,
      ymax = `female_2.5%`,
      xmin = `male_97.5%`,
      ymin = `female_97.5%`,
      col = background_library
    )
  ) +
  geom_pointrange(alpha = 0.2) +
  geom_errorbarh(alpha = 0.2) + geom_point() +
  theme_bw() + 
  # scale_x_log10() + scale_y_log10() +
  # geom_hline(yintercept = 1, lty = 2) +
  # geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Male germline gene fitness benefit (offspring / expected offspring)",
       y = "Female germline gene fitness benefit (offspring / expected offspring)") +
  scale_color_manual(values = c(1, 4)) +
  facet_wrap(background_library ~ driver,
             scales = "free",
             nrow = 2)  

parfit_summaries %>%
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    `2.5%`,
    se_mean,
    `50%`,
    mean,
    `97.5%`
  ) %>%
  gather(
    -sex_test,
    -background_library,
    -driver,
    -gene_name,
    -gene_age,
    -gene_origin,
    key = "statistic",
    value = "estimate"
  ) %>%  
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) %>% 
  # filter(!is.na(male_mean) & !is.na(female_mean) &
  #          `female_97.5%` < 7 & `male_97.5%` < 7) %>%
  ggplot(
    aes(
      x = exp(-`male_mean`),
      y = exp(-`female_mean`),
      xmax = exp(-`male_mean` + 2*`male_se_mean`),
      ymax = exp(-`female_mean` + 2*`female_se_mean`),
      xmin = exp(-`male_mean` - 2*`male_se_mean`),
      ymin = exp(-`female_mean` - 2*`female_se_mean`),
      col = background_library
    )
  ) +
  geom_pointrange(alpha = 0.2) +
  geom_errorbarh(alpha = 0.2) + geom_point() +
  theme_bw() + 
  # scale_x_log10() + scale_y_log10() +
  # geom_hline(yintercept = 1, lty = 2) +
  # geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Male fitness effect (WT / KD)",
       y = "Female fitness effect (WT / KD)") +
  scale_color_manual(values = c(1, 4)) +
  facet_wrap(background_library ~ driver,
             scales = "free",
             nrow = 2)  

parfit_summaries %>%
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    `2.5%`,
    `50%`,
    mean,
    `97.5%`
  ) %>%
  gather(
    -sex_test,
    -background_library,
    -driver,
    -gene_name,
    -gene_age,
    -gene_origin,
    key = "statistic",
    value = "estimate"
  ) %>%  
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) %>% 
  filter(!is.na(male_mean) & !is.na(female_mean) &
           `female_97.5%` < 7 & `male_97.5%` < 7) %>%
  ggplot(
    aes(
      x = `male_mean`,
      y = `female_mean`,
      xmax = `male_2.5%`,
      ymax = `female_2.5%`,
      xmin = `male_97.5%`,
      ymin = `female_97.5%`,
      col = background_library
    )
  ) +
  geom_pointrange(alpha = 0.2) +
  geom_errorbarh(alpha = 0.2) + geom_point() +
  theme_bw() + 
  # scale_x_log10() + scale_y_log10() +
  # geom_hline(yintercept = 1, lty = 2) +
  # geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Male knockdown performance",
       y = "Female knockdown performance") +
  scale_color_manual(values = c(1, 4)) +
  facet_wrap(background_library ~ driver,
             scales = "free",
             nrow = 2)  

parfit_summaries %>% filter(!is.na(driver)) %>% 
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    `2.5%`,
    `50%`,
    mean,
    `97.5%`
  ) %>%
  gather(
    -sex_test,
    -background_library,
    -driver,
    -gene_name,
    -gene_age,
    -gene_origin,
    key = "statistic",
    value = "estimate"
  ) %>%  
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) %>% 
  filter(!is.na(male_mean) & !is.na(female_mean) &
           `female_97.5%` < 7 & `male_97.5%` < 7) %>%
  ggplot(
    aes(
      x = `male_mean`,
      y = `female_mean`,
      xmax = `male_2.5%`,
      ymax = `female_2.5%`,
      xmin = `male_97.5%`,
      ymin = `female_97.5%`,
      col = background_library
    )
  ) +
  geom_pointrange(alpha = 0.2) +
  geom_errorbarh(alpha = 0.2) + geom_point() +
  theme_bw() + 
  # scale_x_log10() + scale_y_log10() +
  # geom_hline(yintercept = 1, lty = 2) +
  # geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Male knockdown performance",
       y = "Female knockdown performance") +
  scale_color_manual(values = c(1, 4))   

parfit_summaries %>%
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    `2.5%`,
    `50%`,
    mean,
    `97.5%`
  ) %>%
  gather(
    -sex_test,-background_library,-driver,-gene_name,-gene_age,-gene_origin,
    key = "statistic",
    value = "estimate"
  ) %>%
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) %>%
  filter(driver == "Nos-GAL4") %>%
  mutate(
    conflict_category = ifelse(
      `male_2.5%` > 0 & `female_97.5%` < 0,
      "female_favored_conflict",
      ifelse(
        `female_2.5%` > 0 & `male_97.5%` < 0,
        "male_favored_conflict",
        ifelse(
          `female_97.5%` < 0 &
            `male_97.5%` < 0,
          "generally_beneficial",
          ifelse(
            `male_2.5%` > 0 & `female_2.5%` > 0,
            "generally_detrimental",
            "neutral_leaning"
          )
        )
      )
    )
  ) %>% 
  rename(effect_estimate_female = female_mean,
         effect_estimate_male = male_mean) %>% 
  select(-starts_with("female"), -starts_with("male")) %>% 
  write_csv("./data/results_tables/nosgal4_conflict_results.csv")

#FINAL FIGURE SHARED W LAB----

pivoted_summaries <- parfit_summaries %>%
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    `2.5%`,
    se_mean,
    `50%`,
    mean,
    `97.5%`
  ) %>%
  gather(
    -sex_test,
    -background_library,
    -driver,
    -gene_name,
    -gene_age,
    -gene_origin,
    key = "statistic",
    value = "estimate"
  ) %>%  
  pivot_wider(names_from = c(sex_test, statistic),
              values_from = estimate) 

pivoted_summaries %>% 
  filter(is.na(driver)) %>% 
  select(background_library,
         gene_name,
         female_mean,
         male_mean) %>% 
  rename(male_con_fx = male_mean,
         female_con_fx = female_mean
         ) %>% 
  right_join(pivoted_summaries %>% 
              filter(!is.na(driver))) %>% 
  mutate(female_con_fx = ifelse(is.na(female_con_fx), 0, female_con_fx),
         male_con_fx = ifelse(is.na(male_con_fx), 0, male_con_fx)) %>% 
 filter(female_mean > -log(100) & !is.na(female_mean) & !is.na(male_mean)) %>% 
  ggplot(
    aes(
      x = exp(-`male_mean` + male_con_fx),
      y = exp(-`female_mean` + female_con_fx),
      xmax = exp(-`male_mean` + male_con_fx + 2*`male_se_mean`),
      ymax = exp(-`female_mean` + female_con_fx + 2*`female_se_mean`),
      xmin = exp(-`male_mean` + male_con_fx - 2*`male_se_mean`),
      ymin = exp(-`female_mean` + female_con_fx - 2*`female_se_mean`),
      col = background_library
    )
  ) +
  geom_pointrange(alpha = 0.2) +
  geom_errorbarh(alpha = 0.2) + geom_point() +
  theme_bw() + 
  # scale_x_log10() + scale_y_log10() +
  # geom_hline(yintercept = 1, lty = 2) +
  # geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Male fitness effect (WT / KD)",
       y = "Female fitness effect (WT / KD)") +
  scale_color_manual(values = c(1, 4)) +
  facet_grid(.~driver, scales = "free")
  # facet_wrap(background_library ~ driver,
  # nrow = 2, scales = "free")

