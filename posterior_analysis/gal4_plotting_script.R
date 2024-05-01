require(tidyverse)
require(ggExtra)

nosgal4_best_model <- 4
soma_best_model <- 2


eta_summaries <- read_csv("./posterior_analysis/eta_summaries_germline.csv") %>% 
  filter(model_index == nosgal4_best_model)

eta_summaries_soma <- read_csv("./posterior_analysis/eta_summaries_soma.csv") %>% 
  filter(model_index == soma_best_model)

eta_fx <- eta_summaries  %>%
  group_by(background_library, sex_test, gene_name) %>% 
  slice(1) %>% #removes duplicate rows due to multiple experiments of same sex/background/gene
  select(
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    lwr_2.5 = `2.5%`,
    median = `50%`,
    mean,
    sd,
    upr_97.5 = `97.5%`
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

eta_fx_soma <- eta_summaries_soma  %>%
  group_by(background_library, sex_test, gene_name, construct_line) %>% 
  select(construct_line,
    sex_test,
    background_library,
    driver,
    gene_name,
    gene_age,
    gene_origin,
    lwr_2.5 = `2.5%`,
    median = `50%`,
    mean,
    sd,
    upr_97.5 = `97.5%`
  ) %>%
  gather(-construct_line,
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

scatterplot_eta <- eta_fx %>% 
  ggplot(aes(x = exp(-male_mean), y = exp(-female_mean),
             xmax = exp(-(male_mean - male_sd)), xmin = exp(-(male_mean + male_sd)),
             ymax = exp(-(female_mean - female_sd)), ymin = exp(-(female_mean + female_sd)),
             col = gene_age)) + 
  geom_point() +
  geom_pointrange(alpha = 0.1) +
  geom_errorbarh(alpha = 0.1) +
  scale_x_log10() +
  theme_bw() +
  scale_y_log10() + 
  geom_hline(yintercept = 1, lty = 3) +
  geom_vline(xintercept = 1, lty = 3) +
  theme(legend.position = "bottom") +
  labs(x = "Male gene fitness contribution (WT / KD)", y = "Female gene fitness contribution (WT / KD)")

ggMarginal(scatterplot_eta)  


eta_fx_soma %>% 
  filter(exp(-male_mean) < 100 & exp(-female_mean) < 100) %>% 
  bind_rows(eta_fx) %>% 
  ggplot(aes(x = exp(-male_mean), y = exp(-female_mean),
             xmax = exp(-(male_mean - male_sd)), xmin = exp(-(male_mean + male_sd)),
             ymax = exp(-(female_mean - female_sd)), ymin = exp(-(female_mean + female_sd)),
             col = driver)) + 
  geom_point() +
  geom_pointrange(alpha = 0.1) +
  geom_errorbarh(alpha = 0.1) +
  scale_x_log10() +
  theme_bw() +
  scale_y_log10() + 
  geom_hline(yintercept = 1, lty = 3) +
  geom_vline(xintercept = 1, lty = 3) +
  theme(legend.position = "bottom") +
  labs(x = "Male gene fitness contribution (WT / KD)", y = "Female gene fitness contribution (WT / KD)") + 
  facet_grid(driver~.)

eta_fx_soma %>% 
  bind_rows(eta_fx) %>% 
  write_csv("./posterior_analysis/effect_size_summary_alldrivers.csv")
