# Script to generate observed vs. model-estimated attack rate scatter plots for
# Qdenga, Dengvaxia, and Butantan-DV (Fig. 1a).  Inputs pre-computed
# model posterior attack rates (output of script 1) and trial data and combines 
# them into a single faceted figure (saved as an RDS for downstream use in 
# script 3).


# ── Set up script ─────────────────────────────────────────────────────────────

library(tidyverse)
library(Hmisc)
library(patchwork)

# Source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

input  = "output/M4/" # final model
output = "output/figures"

# Colour palettes
cols = c("#8EAFBF", "#D55E00","#CC79A7","#6B0040")

theme_set(
  theme_bw() +
    theme(
      text            = element_text(size = 7),
      axis.text       = element_text(size = 7),
      axis.title      = element_text(size = 7),
      plot.title      = element_text(size = 7, hjust = 0.5),
      plot.subtitle   = element_text(size = 7),
      plot.caption    = element_text(size = 7),
      legend.text     = element_text(size = 7),
      strip.text      = element_text(size = 7),
      legend.spacing.y = unit(1, "pt"),
      legend.spacing.x = unit(1, "cm"),
      legend.margin   = margin(0, 0, 0, 0),
      legend.title    = element_blank(),
      legend.key.size = unit(6, "pt"),
      plot.margin = margin(t = 0.1, b = 0, l = 0, r =0.2, "cm")
    )
)

# ── Load data ─────────────────────────────────────────────────────────────────

# Posterior model summaries (attack rate estimates) for each vaccine
AR_Q    = readRDS(paste0(input, "AR_Q.RDS"))
AR_D    = readRDS(paste0(input, "AR_De.RDS"))
AR_B    = readRDS(paste0(input, "AR_Bu.RDS"))

# Observed trial data
VCD_Q   = read.csv("data/processed_Q/vcd_data.csv")
hosp_Q  = read.csv("data/processed_Q/hosp_data.csv")
VCD_D   = readRDS("data/processed_De/cases_stan_format.RDS")
cases_B = readRDS("data/processed_Bu/cases_stan_format.RDS")

# ── Qdenga ────────────────────────────────────────────────────────────────────

# Extract posterior summaries (mean + credible intervals) from the model object
AR_model_Q = extract_model_results(AR_Q)


# Format data and calculate exact binomial confidence intervals
VCD_Q      = factor_VCD_Q(VCD_Q)
hosp_Q     = factor_VCD_Q(hosp_Q)

data_Q = VCD_Q %>% 
  select(month, age, serostatus, serotype, trial, symp = Y, N) %>%  
  left_join(hosp_Q %>% 
              select(month, age, serostatus, serotype, trial, hosp = Y)) %>% 
  filter(! (age == "all" & serotype == "all"), serostatus != "both", trial != "both") 

data_ar_Q = data_Q %>% 
  pivot_longer(cols = c(symp, hosp), values_to = "Y", names_to = "outcome") %>% 
  filter(!is.na(Y)) %>%  # hosp not at 12 and 18 months 
  mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
         lower = binconf(Y,N, method = "exact")[,2] * 100, 
         upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
  select(- c(Y,N)) %>%  
  mutate(outcome = factor(outcome,
                              levels = c("symp", "hosp"),
                              labels = c("symptomatic", "hospitalised")))



month_lookup = c("1" = 12, "2" = 18, "3" = 24, "4" = 36, "5" = 48, "6" = 54, "all" = "all")

# Format model outputs 
model_age_Q = AR_model_Q %>%
  filter(group == "AR_BVJRD_Q") %>% 
  separate(name, into = c("serostatus", "trial", "age", "outcome", "month"))

model_serotype_Q = AR_model_Q %>%
  filter(group == "AR_BVKRD_Q") %>% 
  separate(name, into = c("serostatus", "trial", "serotype", "outcome", "month"))

model_serotype_age_Q = AR_model_Q %>%
  filter(group == "AR_BVKJR_Q") %>% 
  separate(name, into = c("serostatus", "trial", "serotype", "age", "outcome"))

# Combine all model strata, fill missing stratifiers with "all",
# apply factor labels, and join to observed data

all_Q = bind_rows(model_age_Q, model_serotype_age_Q, model_serotype_Q) %>% 
  mutate(month = ifelse(is.na(month), "all", month),
         age = ifelse(is.na(age), "all", age),
         serotype = ifelse(is.na(serotype), "all", serotype)) %>% 
  mutate(
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    trial      = factor(trial, labels = c("placebo", "vaccine")),
    serotype = factor(serotype, levels = c(1,2,3,4,"all"), labels = c("DENV1", "DENV2", "DENV3", "DENV4", "all")),
    age        = factor(age, levels = c(1, 2, 3, "all"), labels = c("4-5yrs", "6-11yrs", "12-16yrs", "all")),
    month   = month_lookup[month],
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))
  ) %>%
  rename(model_lower = lower,
         model_mean = mean,
         model_upper = upper) %>% 
  left_join(data_ar_Q) %>% 
  mutate(x = paste0(trial, " ", serostatus))

# Diagnostic scatter plot: observed vs. model attack rates for Qdenga
Q_all_plot = all_Q %>%  
  ggplot(aes(x = mean, y = model_mean, color = x)) +
  geom_point(alpha = 0.2) +
  geom_errorbar(aes(ymin = model_lower, ymax = model_upper),alpha = 0.2) +
  geom_errorbar(aes(xmin = lower, xmax = upper),alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2) +
  labs(y = "model attack rate (%)", x = "data attack rate (%)")


# ── Dengvaxia ─────────────────────────────────────────────────────────────────

# Extract posterior summaries (mean + credible intervals) from the model object

AR_model_D = extract_model_results(AR_D)

# Format data and calculate exact binomial confidence intervals

AR_data_D  = lapply(VCD_D, calc_attack_rates) %>%  bind_rows() %>% 
  mutate(time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60", "all")))

# Format model outputs 
Sy_AR_VK_De = AR_model_D %>%
  filter(group == "V_AR_VK_De") %>%
  separate(name, into = c("arm", "serotype")) %>% 
  mutate(outcome = "Sy")

Sy_AR_BVJ_De = AR_model_D %>%
  filter(group == "V_AR_BVJ_De") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(outcome = "Sy")

Ho_AR_BVJD_De = AR_model_D %>%
  filter(group == "H_AR_BVJD_De") %>%
  separate(name, into = c("serostatus", "arm", "age", "time")) %>% 
  mutate(outcome = "Ho")

Ho_AR_BVKJ_De = AR_model_D %>%
  filter(group == "H_AR_BVKJ_De") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(outcome = "Ho")


# Combine all model strata, fill missing stratifiers with "all",
# apply factor labels

model_De = bind_rows(Sy_AR_VK_De, Sy_AR_BVJ_De, Ho_AR_BVJD_De, Ho_AR_BVKJ_De) %>% 
  select(- type) %>% 
  mutate(time = case_when(
    is.na(time) & grepl("V_AR", group) ~ "1",
    is.na(time) & grepl("H_AR", group) ~ "all",
    TRUE ~ as.character(time)
  ), 
         serostatus = ifelse(is.na(serostatus), "both", serostatus),
         age = ifelse(is.na(age), "all", age), 
         serotype = ifelse(is.na(serotype), "all", serotype)) %>% 
  mutate(
    serostatus = factor(serostatus, levels = c(1, 2, "both"), labels = c("seronegative", "seropositive", "both")),
    arm        = factor(arm, labels = c("placebo", "vaccine")),
    serotype   = factor(serotype, levels = c(1:4, "all"), labels = c("DENV1", "DENV2", "DENV3", "DENV4", "all")),
    age        = factor(age,  levels = c(1:2, "all"), labels = c("2-8yrs", "9-16yrs", "all")),
    time = factor(time, levels = c(1:4, "all"), labels = c("1-13", "14-24", "25-36", "37-60", "all")),
  ) %>%
  rename(model_lower = lower,
         model_mean = mean,
         model_upper = upper) 

# Combine data and model 
all_D = model_De %>% 
left_join(AR_data_D)  %>% 
  filter(serostatus != "both") %>% 
  mutate(x = paste0(arm, " ", serostatus))


# Diagnostic scatter plot: observed vs. model attack rates for Dengvaxia

D_all_plot = all_D %>%  
  ggplot(aes(x = mean, y = model_mean, color = x)) +
  geom_point(alpha = 0.2) +
  geom_errorbar(aes(ymin = model_lower, ymax = model_upper),alpha = 0.2) +
  geom_errorbar(aes(xmin = lower, xmax = upper),alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2) +
  labs(y = "model attack rate (%)", x = "data attack rate (%)")




# ── Butantan-DV ───────────────────────────────────────────────────────────────

# Extract posterior summaries (mean + credible intervals) from the model object

AR_model_B = extract_model_results(AR_B)

# Format data and calculate exact binomial confidence intervals
AR_data_B  = lapply(cases_B, calc_attack_rates) %>%  bind_rows() 

# Format model outputs 
serotype_B = AR_model_B %>%
  filter(group == "AR_BVKD_Bu") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "time"))
  
age_B = AR_model_B %>%
  filter(group == "AR_BVJD_Bu") %>%
  separate(name, into = c("serostatus", "arm", "age", "time")) 

serotype_age_B = AR_model_B %>%
  filter(group == "AR_BVKJ_Bu") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age"))

# Combine all model strata, fill missing stratifiers with "all",
# apply factor labels

model_B = bind_rows(serotype_B, age_B, serotype_age_B) %>% 
  select(- type) %>% 
  mutate(
  time = ifelse(is.na(time), "all", time),
  age = ifelse(is.na(age), "all", age), 
  serotype = ifelse(is.na(serotype), "all", serotype)
  ) %>% 
  mutate(
    arm        = factor(arm, labels = c("placebo", "vaccine")),
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    serotype   = factor(serotype, levels = c("1", "2", "all"),  labels = c("DENV1", 
                                                                           "DENV2", "all")),
    age        = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs", "all")),
    time   = factor(time, labels = c("1-24", "25-60", "1-60"))
  ) %>%
  rename(model_lower = lower,
         model_mean = mean,
         model_upper = upper) 

# Combine data and model 
all_B = model_B %>% 
  left_join(AR_data_B) %>% 
  mutate(x = paste0(arm, " ", serostatus))

# Diagnostic scatter plot: observed vs. model attack rates for Butantan-DV
B_all_plot = all_B %>%  
  ggplot(aes(x = mean, y = model_mean, color = x)) +
  geom_point(alpha = 0.2) +
  geom_errorbar(aes(ymin = model_lower, ymax = model_upper),alpha = 0.2) +
  geom_errorbar(aes(xmin = lower, xmax = upper),alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2) +
  labs(y = "model attack rate (%)", x = "data attack rate (%)")

# ── All plots ───────────────────────────────────────────────────────────────

# Plot all vaccines together 

all_plot = 
  bind_rows(all_Q, all_D, all_B, .id = "vaccine") %>% 
  mutate(vaccine = factor(vaccine, labels = c("Qdenga", "Dengvaxia", "Butantan-DV"))) %>% 
  filter(!is.na(mean)) %>%  # remove points with model but not data output
  mutate(x = factor(x, levels = c("placebo seronegative", "placebo seropositive",
                                  "vaccine seronegative", "vaccine seropositive"))) %>% 
  ggplot(aes(x = mean, y = model_mean, color = x)) +
  geom_point(alpha = 0.4) +
  geom_errorbar(aes(ymin = model_lower, ymax = model_upper),alpha = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper),alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2) +
  labs(y = "Estimated attack rate (%)", x = "Observed attack rate (%)") +
  scale_color_manual(values =cols) +
  facet_wrap(~vaccine) +
  theme(legend.position = "top")

# Save as RDS to combine with other plots in script 3. 
saveRDS(all_plot, file = "output/figures/all_attack_rates.rds")


# Save source data 
# source_data =  bind_rows(all_Q, all_D, all_B, .id = "vaccine") 
# write.csv(source_data, "output/figures/source_data/Figure1a.csv")
