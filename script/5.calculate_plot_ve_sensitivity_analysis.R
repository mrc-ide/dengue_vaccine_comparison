
# Script to plot impact of sensitivity analyses on each vaccine's efficacy 
# aggregated over time. Compares VE estimates from a baseline model (M4) against 
# a set of sensitivity analysis (SA) models that vary structural assumptions
# (e.g. age-dependent VE, serotype enhancement parameters, antibody titre inputs).
# For each vaccine, the absolute difference in age-weighted mean VE between
# each SA model and the baseline is computed and visualised as a heatmap
# stratified by serostatus, serotype, and clinical outcome.
# Inputs pre-computed model posteriors of vaccine efficacy
# from different models (output of script 1). Saves the plot as a PDF (Fig. 4)

# ── Set up script ─────────────────────────────────────────────────────────────

rm(list = ls())
library(tidyverse)
library(cowplot)
library(posterior)
library(patchwork)

source("R/format_model_outputs_for_plotting.R")

n = 100

# Age weights 
Qage_weights = tibble(
  age = c("1", "2", "3"),
  n   = c(1702, 7387, 4291)
)

Dage_weights = tibble(
  age = c("1", "2"),
  n   = c(718, 2638)
)

Bage_weights = tibble(
  age = c("1", "2", "3"),
  n   = c(4931, 4980, 5694)
)

# Labels 

serostatus_labels = c(
  "1" = "Seronegative",
  "2" = "Monotypic",
  "3" = "Multitypic"
)

mycols = c("#B39EB5", "#B2D3A8", "#89A497", "#AEC6CF", "#B8B8AC", "#FFFFFF")

my_theme = theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.subtitle = element_text(size = 7),
      plot.caption = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_blank(),
      strip.text = element_text(size = 7),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)))

# ── Model folder lists ────────────────────────────────────────────────────────

Q_folder_names = c("M4", "M3", "M2", "M9", "M10", "M13", "M15")
D_folder_names = c("M4", "M5", "M11", "M14", "M15")
B_folder_names = c("M4", "M6", "M7", "M8", "M12", "M15")

# ── Load and reshape posterior VE draws ──────────────────────────────────────

VE_Q_list = process_VE_files(Q_folder_names, "Q", n)
VE_D_list = process_VE_files(D_folder_names, "De", n)
VE_B_list = process_VE_files(B_folder_names, "Bu", n)

VE_Q = process_sens(VE_Q_list, c("serostatus", "serotype", "age", "outcome", "month")) 
VE_D = process_sens(VE_D_list, c("serostatus", "serotype", "age", "outcome", "month")) 
VE_B = process_sens(VE_B_list, c("serostatus", "serotype", "age", "month")) 

# ── Compute VE differences vs. baseline (M4) ─────────────────────────────────

Q_diff = calculate_diff(
  df = VE_Q,
  age_weights = Qage_weights,
  main_model = "M4",
  serostatus_labels = serostatus_labels
) 
09
D_diff = calculate_diff(
  VE_D,
  age_weights = Dage_weights,
  main_model = "M4",
  serostatus_labels = serostatus_labels
)

VE_B$outcome = 1 # only symp, to match D and Q 

B_diff = calculate_diff(
  VE_B,
  age_weights = Bage_weights,
  main_model = "M4",
  serostatus_labels = serostatus_labels
)

# Cache the computed differences to avoid re-running the full posterior pipeline

ve_out = list(Q_diff =Q_diff,
              D_diff = D_diff,
              B_diff = B_diff)

saveRDS(ve_out, "output/ve_sens.RDS")

ve_out = readRDS("output/ve_sens.RDS")
list2env(ve_out, envir = .GlobalEnv)

# ── Relabel and order SA model codes ─────────────────────────────────────────


Q_data_plot = Q_diff %>%  
  mutate(sa = case_when(
    sa == "M2" ~ "No age VE",
    sa == "M3" ~ "Age VE\n(all age groups)",
    sa == "M4" ~ "Baseline",
    sa == "M9" ~ "No age VE\nAge FOI (4-5yrs)",
    sa == "M10" ~ "Serotype enhancement\nparameter",
    sa == "M13" ~ "Serotype hospital\nenhancement parameter",
    sa == "M15" ~ "No serotype\nantibody titres",
     TRUE   ~ NA
  )) %>% 
  mutate(sa = factor(sa, levels = c(
    "Baseline",
    "No age VE",
    "Age VE\n(all age groups)",
    "No age VE\nAge FOI (4-5yrs)",
    "No serotype\nantibody titres",
    "Serotype enhancement\nparameter",
    "Serotype hospital\nenhancement parameter"
  )))

D_data_plot = D_diff %>%  
  mutate(sa = case_when(
    sa == "M4" ~ "Baseline",
    sa == "M5" ~ "Age VE",
    sa == "M11" ~ "Serotype enhancement\nparameter",
    sa == "M14" ~ "Serotype hospital\nenhancement parameter",
    sa == "M15" ~ "No serotype\nantibody titres",
    TRUE   ~ NA
  )) %>% 
  mutate(sa = factor(sa, levels = c(
    "Baseline",
    "Age VE",
    "No serotype\nantibody titres",
    "Serotype enhancement\nparameter",
    "Serotype hospital\nenhancement parameter"
  )))


B_data_plot = B_diff %>%  
  mutate(sa = case_when(
    sa == "M4" ~ "Baseline",
    sa == "M6" ~ "Age VE",
    sa == "M7" ~ "Age VE (2-5yrs)",
    sa == "M8" ~ "Age FOI (2-5yrs)",
    sa == "M12" ~ "Serotype enhancement\nparameter",
    sa == "M15" ~ "No serotype\nantibody titres",
    TRUE   ~ NA
  )) %>% 
  mutate(sa = factor(sa, levels = c(
    "Baseline",
    "Age VE",
    "Age VE (2-5yrs)",
    "Age FOI (2-5yrs)",
    "No serotype\nantibody titres",
    "Serotype enhancement\nparameter")))

# ── Compute shared colour scale limits ────────────────────────────────────────

abs_max = max(c(abs(filter(D_data_plot, !main)$difference), 
                abs(filter(Q_data_plot, !main)$difference),
                abs(filter(B_data_plot, !main)$difference)
                ),
              na.rm = TRUE)

# Round up to nearest 5pp for clean breaks
scale_limit = ceiling(abs_max / 5) * 5

# ── Generate heatmaps ─────────────────────────────────────────────────────────

Q_plot = plot_heatmap(Q_data_plot, scale_limit)

D_plot = plot_heatmap(D_data_plot, scale_limit)

B_plot = plot_heatmap(B_data_plot, scale_limit)


# ── Assemble and save Fig. 4 ──────────────────────────────────────────────────


all = (Q_plot / D_plot / B_plot) +
  plot_annotation(tag_levels = "a") +
  plot_layout(
    guides  = "collect",
    heights = c(1.3, 1, 1)  
  ) &
  theme(legend.position = "bottom")


ggsave(
  filename = "output/figures/Fig4.pdf",
  plot     = all,
  width    = 18,
  height   = 20,   
  units    = "cm"
)


# Save source data 
# source_data = bind_rows(ve_out)
# write.csv(source_data, "output/figures/source_data/Figure4.csv")

