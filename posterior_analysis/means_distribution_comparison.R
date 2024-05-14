require(mclust)
require(tidyverse)

eta_sex_table <- read_csv("./posterior_analysis/effect_size_summary_alldrivers.csv")

driver <- "Tub" #options are ("Nos", "Act", "Tub")

#collating posterior data-----
posterior_data <- eta_sex_table %>% 
  filter(driver_shortname == driver) %>% 
  mutate(male_effect = exp(male_cross_mean + male_mean),
         female_effect = exp(female_cross_mean + female_mean),
         sex_chr = chromosome == "X",
         old_gene = gene_age == ">35MYA",
         young_gene = gene_age == "0-3MYA") %>% 
  filter(!(is.na(male_effect) | is.na(female_effect)) &
           chromosome %in% c("X", "2", "3"))

sex_effect_cluster_test <- posterior_data %>% 
  select(gene_name, gene_age,
         chromosome, driver_shortname, 
         male_effect, female_effect) %>% 
  gather(-gene_name, -gene_age, 
         -chromosome, -driver_shortname, 
         key = "sex", value = "effect") %>% 
  mutate(sex = substr(sex,1,nchar(sex)-7))

#Chisq code for assessing enrichment (or lack of enrichment) in sex conflict for each chromosome/gene age

#BIC table-----

bic_table <- data.frame()
for (driver in unique(eta_sex_table$driver_shortname)) {
  tmp <- eta_sex_table %>% 
    filter(driver_shortname == driver) %>% 
    mutate(male_effect = exp(male_cross_mean + male_mean),
           female_effect = exp(female_cross_mean + female_mean)) %>% 
    filter(!(is.na(male_effect) | is.na(female_effect)) &
             chromosome %in% c("X", "2", "3")) %>% 
    mutate(gene_age = factor(ifelse(gene_age %in% c("0-3MYA", "3-6MYA"), "0-6MYA", gene_age),
                             levels = c("0-6MYA", "6-11MYA","11-25MYA",">35MYA"))) 
  for (chr in c(T,F)) {
    for (age in unique(tmp$gene_age)) {
      subtmp <- tmp %>% 
        filter(gene_age == age & 
                 chr == (chromosome == "X"))
      
      male_ll <- ret_ll(subtmp$male_effect)
      female_ll <- ret_ll(subtmp$female_effect)
      overall_ll <- ret_ll(c(subtmp$male_effect,subtmp$female_effect))
      
      sex_bic <- 4 * log(length(c(subtmp$male_effect,
                                  subtmp$female_effect))) - 2*(female_ll + male_ll)
      null_bic <-  2 * log(length(c(subtmp$male_effect,
                                    subtmp$female_effect))) - 2*(overall_ll)
      
      bic_table <- bind_rows(bic_table, data.frame(     
        male_mean = mean(subtmp$male_effect),
        female_mean = mean(subtmp$female_effect),
        male_sd = sd(subtmp$male_effect),
        female_sd = sd(subtmp$female_effect),
        sex_chromosome = chr,
        gene_age = factor(age, levels = c("0-6MYA", "6-11MYA","11-25MYA",">35MYA")),
        driver_shortname = driver,
        sex_bic = sex_bic,
        null_bic = null_bic,
        delta_bic = null_bic - sex_bic) #positive values indicate informative difference between sex distributions
      )
    }
  }
}
#Figures-----
posterior_data %>% 
  mutate(conflict = (male_effect < 1 & female_effect > 1) | (male_effect > 1 & female_effect < 1)) %>% 
  select(gene_age, gene_name, chromosome, driver_shortname, conflict) %>% 
  right_join(sex_effect_cluster_test) %>% 
  mutate(gene_age = factor(ifelse(gene_age %in% c("0-3MYA", "3-6MYA"), "0-6MYA", gene_age),
                           levels = c("0-6MYA", "6-11MYA","11-25MYA",">35MYA"))) %>% 
  ggplot(aes(x = sex, y = effect, alpha = conflict, group = gene_name, col = gene_age)) +
  geom_hline(yintercept = 1, lty =2) +
  scale_alpha_manual(values = c(0.25,1)) +
  geom_line() +
  geom_point() +
  facet_grid(gene_age~chromosome) +
  theme_bw() +
  labs(x = "Sex", y = "RNAi KD fitness effect",
       col = "Gene age", alpha = "Under sex conflict")

#BIC differences
bic_table %>% 
  ggplot(aes(x = gene_age, y = delta_bic, group = driver_shortname, col = driver_shortname)) +
  geom_hline(yintercept = 0, lty =2) +
  geom_line()+
  geom_point() +
  theme_bw() +
  facet_wrap(sex_chromosome~.,
             labeller = labeller(sex_chromosome = as_labeller(function(x){
               ifelse(x, "Sex chromosome", "Autosomal chromosome")}))) +
  labs(x = "Gene age", y = "Sex-specific differences information (BIC)")

#
bic_table %>% 
  ggplot(aes(x = gene_age, y = female_mean - male_mean, group = driver_shortname, col = driver_shortname)) +
  geom_hline(yintercept = 0, lty =2) +
  geom_line()+
  geom_point() +
  theme_bw() +
  facet_wrap(sex_chromosome~.)  +
  labs(x = "Gene age", y = "Male relative fitness benefit")
