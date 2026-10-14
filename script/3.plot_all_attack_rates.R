
# Script to generate aggregated observed vs. model-estimated attack rate plots 
# by serotype, serosatus, month, trial arm for each vaccine. Inputs pre-computed
# model posteriors (output of script 1) and trial data. Combines plots with 
# "output/figures/all_attack_rates.rds" - scatter plot generated in script 2.
# Saves the plot as a PDF (Fig 1).

# ── Set up script ─────────────────────────────────────────────────────────────

library(tidyverse)
library(Hmisc)
library(patchwork)
library(grid)
library(cowplot)


# Source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

input  = "output/M4/"
output = "output/figures"

# Colour palettes
serotype_fill = c("#148695", "#C51B8A", "#BDCACE", "#016C59", "#111111")
trial_fill    = c("#6B0040", "#7791D1", "#D55E00","#CC79A7") 


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
      legend.key.size = unit(6, "pt") )
)

# ── Load data ─────────────────────────────────────────────────────────────────

AR_Q    = readRDS(paste0(input, "AR_Q.RDS"))
VCD_Q   = read.csv("data/processed_Q/vcd_data.csv")
hosp_Q  = read.csv("data/processed_Q/hosp_data.csv")

AR_D    = readRDS(paste0(input, "AR_De.RDS"))
VCD_D   = readRDS("data/processed_De/cases_stan_format.RDS")

AR_B    = readRDS(paste0(input, "AR_Bu.RDS"))
cases_B = readRDS("data/processed_Bu/cases_stan_format.RDS")

all_plot = readRDS("output/figures/all_attack_rates.rds")
# ── Qdenga ────────────────────────────────────────────────────────────────────

VCD_Q      = factor_VCD_Q(VCD_Q)
hosp_Q     = factor_VCD_Q(hosp_Q)
AR_model_Q = extract_model_results(AR_Q)

AR_data_Q  = calc_attack_rates_Q(VCD_Q)
HR_data_Q  = calc_hosp_rates_Q(VCD = VCD_Q, hosp = hosp_Q)
data_Q     = bind_rows(AR_data_Q, HR_data_Q)

# Relabel outcomes consistently
relabel_outcome = function(df) {
  mutate(df, outcome = factor(outcome,
                              levels = c("VCD", "hosp"),
                              labels = c("symptomatic", "hospitalised")
  ))
}

# Time plot
AR_VRD_Q = data_Q %>%
  filter(serostatus == "both", age == "all", trial != "both", serotype == "all") %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial),
         month = as.numeric(month)) %>%
  select(-serotype, -age, -serostatus, -X, -year) %>%
  relabel_outcome()

month_lookup = c("1" = 12, "2" = 18, "3" = 24, "4" = 36, "5" = 48, "6" = 54)

time_plot_Q = AR_model_Q %>%
  filter(group == "AR_VRD_Q") %>%
  separate(name, into = c("trial", "outcome", "month")) %>%
  mutate(
    trial   = factor(trial, labels = c("placebo", "vaccine")),
    month   = month_lookup[month],
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))
  ) %>%
  bind_rows(AR_VRD_Q) %>%
  mutate(month = as.character(month)) %>% 
  mutate(x = paste0(trial, " ", outcome)) %>% 
  mutate(x = factor(x, levels = c("placebo symptomatic", "vaccine symptomatic",
                                  "placebo hospitalised", "vaccine hospitalised"),
                    labels = c("P symp", "V symp", "P hosp", "V hosp"))) %>% 
  ggplot(aes(x = month, y = mean,  color = x)) +
  geom_point(aes(shape = type),
             position = position_dodge(width = 0.4), size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, linetype = type),
                position = position_dodge(width = 0.4), width = 0, linewidth = 0.5) +
  labs(x = "Month", y = "Qdenga attack rate (%)") +
  scale_color_manual(values = trial_fill) +
  theme(legend.position = c(0.81,0.76))

# Serotype plot
AR_BVKR_Q = data_Q %>%
  filter(serotype != "all", age == "all", serostatus != "both", month == "all") %>%
  mutate(
    trial    = ifelse(trial == "TAK-003", "vaccine", trial),
    serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4"))
  ) %>%
  select(-month, -X, -year, -age) %>%
  relabel_outcome()

serotype_plot_Q = AR_model_Q %>%
  filter(group == "AR_BVKR_Q") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome")) %>%
  mutate(
    outcome    = factor(outcome, labels = c("symptomatic", "hospitalised")),
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    trial      = factor(trial, labels = c("placebo", "vaccine")),
    serotype   = factor(serotype, labels = c("D1", "D2", "D3", "D4"))
  ) %>%
  bind_rows(AR_BVKR_Q) %>%
  mutate(
    serostatus = factor(serostatus, labels = c("SN", "SP")),
    trial      = factor(trial, labels = c("P", "V"))
  ) %>%
  unite(c(trial, serostatus), sep = " ", col = "x") %>%
  filter(outcome == "symptomatic") %>% 
  ggplot(aes(x = x, y = mean, color = serotype, shape = type, linetype = type,
             group = interaction(type, serotype))) +
  geom_point(position = position_dodge(width = 0.8), size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0, linewidth = 0.5) +
  scale_color_manual(values = serotype_fill) +
  ylab("Qdenga symptomatic\nattack rate over 54 months (%)") +
  theme(legend.position = c(0.85,0.78), axis.title.x = element_blank())

# ── Dengvaxia ──────────────────────────────────────────────────────────────────

VCD_D$Ho_VD = VCD_D$Ho_BVJD %>%
  group_by(arm, time) %>%
  summarise(Y = sum(Y), N = sum(N)) %>%
  mutate(time = factor(time, labels = c("13", "24", "36", "60")))

VCD_D$Ho_BVJ = VCD_D$Ho_BVKJ %>%
  group_by(arm, serostatus, age) %>%
  summarise(Y = sum(Y), N = mean(N))

VCD_D$Ho_BVK = VCD_D$Ho_BVKJ %>%
  group_by(arm, serostatus, serotype) %>%
  summarise(Y = sum(Y), N = sum(N))

AR_data_D  = lapply(VCD_D, calc_attack_rates)
AR_model_D = extract_model_results(AR_D)

# Serotype plot
Sy_AR_VK_De = AR_model_D %>%
  filter(group == "V_AR_VK_De") %>%
  separate(name, into = c("arm", "serotype")) %>%
  mutate(
    arm      = factor(arm, labels = c("placebo", "vaccine")),
    serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")),
    outcome  = "symptomatic"
  ) %>%
  bind_rows(
    AR_data_D$Sy_VK %>%
      mutate(serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")),
             outcome = "symptomatic")
  ) %>%
  mutate(arm = factor(arm, labels = c("P", "V")))

serotype_plot_D = Sy_AR_VK_De %>% 
  ggplot(aes(x = arm, y = mean, color = serotype, shape = type, linetype = type,
             group = interaction(type, serotype))) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.5), width = 0, linewidth = 0.5) +
  scale_color_manual(values = serotype_fill) +
  ylab("Dengvaxia symptomatic\nattack rate over 13 months (%)") +
  theme(legend.position = "none",axis.title.x = element_blank())


# Time plot
time_plot_D = AR_model_D %>%
  filter(group == "H_AR_VD_De") %>%
  separate(name, into = c("arm", "time")) %>%
  mutate(
    arm  = factor(arm, labels = c("placebo", "vaccine")),
    time = factor(time, labels = c("13", "24", "36", "60"))
  ) %>%
  bind_rows(AR_data_D$Ho_VD) %>%
  mutate(outcome = "hospitalised") %>%
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(color = arm, shape = type),
             position = position_dodge(width = 0.6), size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = arm, linetype = type),
                position = position_dodge(width = 0.6), width = 0, linewidth = 0.5) +
  labs(x = "Month", y = "Dengvaxia hospitalised\nattack rate (%)") +
  scale_color_manual(values = c( "#D55E00","#CC79A7")) +
  theme(legend.position = "none")


# ── Butantan-DV ────────────────────────────────────────────────────────────────

cases_B$Sy_VD = cases_B$Sy_BVJ %>%  
  group_by(arm, time) %>%  
  summarise(Y = sum(Y),
            N = sum(N)) %>% 
  mutate(time  = factor(time, labels = c("24", "60")))

cases_B$Sy_BVK = cases_B$Sy_BVKJ %>%  
  group_by(arm, serostatus, serotype) %>%  
  summarise(Y = sum(Y),
            N = sum(N))

AR_data_B  = lapply(cases_B, calc_attack_rates)
AR_model_B         = extract_model_results(AR_B)

# Serotye plot 
serotype_plot_B = AR_model_B %>%
  filter(group == "AR_BVK_Bu") %>%
  separate(name, into = c("serostatus", "arm", "serotype")) %>%
  mutate(
    arm        = factor(arm, labels = c("placebo", "vaccine")),
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    serotype   = factor(serotype, labels = paste0("D", 1:2)),
  ) %>%
  bind_rows(AR_data_B$Sy_BVK %>% mutate(serotype = factor(serotype, labels = c("D1", "D2")))) %>%
  mutate(
    arm        = factor(arm, labels = c("P", "V")),
    serostatus = factor(serostatus, labels = c("SN", "SP")),
    outcome    = "symptomatic"
  ) %>%
  unite(c(arm, serostatus), sep = " ", col = "x") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper,
                    group = interaction(type, serotype), linetype = type, color = serotype),
                position = position_dodge(width = 0.5), width = 0, linewidth = 0.5) +
  scale_color_manual(values = serotype_fill) +
  ylab("Butantan-DV symptomatic\nattack rate over 60 months (%)") +
  theme(legend.position = "none", axis.title.x = element_blank())

# Time plot 

time_plot_B = AR_model_B %>%
  filter(group == "AR_VD_Bu") %>%
  separate(name, into = c("arm", "time")) %>%
  mutate(
    arm        = factor(arm, labels = c("placebo", "vaccine")),
    time        = factor(time, labels = c("24", "60")),
  ) %>%
  bind_rows(AR_data_B$Sy_VD) %>%
  mutate(
    arm        = factor(arm, labels = c("P", "V")),
    outcome    = "symptomatic"
  ) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(color = arm, shape = type),
             position = position_dodge(width = 0.6), size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = arm, linetype = type),
                position = position_dodge(width = 0.6), width = 0, linewidth = 0.5) +
  labs(x = "Month", y = "Butantan-DV symptomatic\nattack rate (%)") +
  scale_color_manual(values = trial_fill) +
  theme(legend.position = "none")

# ── Assemble figure ────────────────────────────────────────────────────────────

grid_plot_B = (time_plot_Q) |
  (time_plot_D) |
  (time_plot_B)  +
  theme(plot.tag = element_text(size = 12))

grid_plot_C = (serotype_plot_Q) |
  (serotype_plot_D) |
  (serotype_plot_B)  +
  theme(plot.tag = element_text(size = 12))

out = plot_grid(
  all_plot, grid_plot_B, grid_plot_C,
  ncol = 1,
  labels = c("a", "b", "c")
)

ggsave(
  plot     = out,
  filename = paste0(output, "/Fig1.pdf"),
  height   = 16,
  width    = 18,
  units    = "cm",
  dpi      = 600
)


# Save source data 
# source_data = bind_rows(AR_model_Q, AR_model_B, AR_model_D)
# library(writexl) 
# write.csv(source_data, "output/figures/source_data/Figure1b_1c_model.csv")
# write.csv(VCD_Q, "output/figures/source_data/Figure1b_1c_Q_data_vcd.csv")
# write.csv(hosp_Q, "output/figures/source_data/Figure1b_1c_Q_data_hosp.csv")
# write_xlsx(VCD_D, "output/figures/source_data/Figure1b_1c_D_data.xlsx")
# write_xlsx(cases_B, "output/figures/source_data/Figure1b_1c_B_data.xlsx")
