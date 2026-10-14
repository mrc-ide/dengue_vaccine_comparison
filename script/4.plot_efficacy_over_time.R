
# Script to plot each vaccine's efficacy over time by clinical outcome age,
# serostatus, serotype, over time for 54 months. Inputs pre-computed model 
# posteriors of vaccine efficacy (output of script 1). Saves the plots as PDFs
# (Fig. 2 and Fig 3).

# ── Set up script ─────────────────────────────────────────────────────────────

library(tidyverse)
source("R/format_model_outputs_for_plotting.R")

# colours 
cols = c("#CC79A7", "#D55E00", "#8EAFBF", "#A9D6C6")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7, hjust = 0.5),
      plot.subtitle = element_text(size = 7),
      plot.caption = element_text(size = 7),
      legend.text = element_text(size = 7),
      strip.text      = element_text(size = 7),
      legend.title = element_blank(),
      strip.clip         = "off",   # prevents label clipping if strip is narrow
      legend.position    = "top",
      legend.margin      = margin(b = 2),
      legend.spacing.x   = unit(6, "pt")
    ))

# files 
input = "output/M4/"

# ── Load posterior VE summaries ───────────────────────────────────────────────

VE_Q = readRDS(paste0(input, "VE_Q.RDS")) 
VE_D = readRDS(paste0(input, "VE_De.RDS"))
VE_B = readRDS(paste0(input, "VE_Bu.RDS"))

# ── Extract and label posterior summaries ─────────────────────────────────────

# VE 
VE_D_model = extract_model_results(VE_D)  
VE_Q_model = extract_model_results(VE_Q)  
VE_B_model = extract_model_results(VE_B)  
VE_D_model$vaccine = "Dengvaxia"
VE_Q_model$vaccine = "Qdenga"
VE_B_model$vaccine = "Butantan-DV"

# ── Combine and reshape all three vaccines ────────────────────────────────────

# Qdenga and Dengvaxia share the same name structure (serostatus/serotype/age/outcome/time)
# and are joined first; Butantan-DV lacks an outcome dimension so is appended separately.

comb_VE = VE_Q_model %>%
  bind_rows(VE_D_model) %>%
  filter(group == "VE_Q" | group == "VE_De") %>%  
  separate(name, into = c("serostatus", "serotype","age", "outcome","time")) %>%  
  filter(age == 1 | age == 2 & group == "VE_Q") %>% # only Qdenga has age
  mutate(
    month = as.numeric(time),
    serostatus = factor(serostatus, levels = 1:3, labels = c("seronegative", "monotypic", "multitypic")),
    serotype = factor(serotype, levels = 1:4, labels = c(paste0("DENV", 1:4))),
    outcome = factor(outcome, levels = 1:2,  labels = c("symptomatic", "hospitalised"))
    ) %>% 
  bind_rows(VE_B_model %>% # Add Butantan-DV after as no outcome column 
              separate(name, into = c("serostatus", "serotype","age","time")) %>% 
              mutate(
                month = as.numeric(time),
                serostatus = factor(serostatus, levels = 1:3, labels = c("seronegative", "monotypic", "multitypic")),
                serotype = factor(serotype, levels = 1:4, labels = c(paste0("DENV", 1:4))),
                outcome = factor("symptomatic")
              ) 
  ) %>%  
  # Truncate to 54 months (the common follow-up horizon across trials)
  filter(month <=54) %>% 
  mutate(age = ifelse(vaccine == "Qdenga" & age == 1, "4-5yrs", 
                      ifelse(vaccine == "Qdenga" & age == 2, "6-16yrs",
                             ifelse(vaccine == "Dengvaxia", "2-16yrs", "2-59yrs"))))


# ── Fig2: VE against symptomatic disease ─────────────────────────────────────

plot_ve_symp = comb_VE %>%   
  filter(outcome == "symptomatic") %>% 
  mutate(x = paste0(vaccine, "\n", age)) %>% 
  ggplot(aes(x = month, y = mean, color = x, fill = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha   = 0.15,
              colour  = NA) + 
  geom_line(
    aes(y = mean),
    colour      = "white",
    linewidth   = 1.1,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.75) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    colour     = "grey30",
    linewidth  = 0.4
  ) +
  labs(x = "Months post final dose", 
       y = "Vaccine efficacy\nagainst symptomatic disease (%)") +
  facet_grid(serostatus  ~ serotype, scales = "free_y") + 
  theme(legend.position = "top") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,12,24,36,48)) +
  guides(
    colour = guide_legend(
      title        = NULL,
      nrow         = 1,
      override.aes = list(linewidth = 1.2, fill = NA)
    ),
    fill = "none"
  )

ggsave(
  plot =  plot_ve_symp,
  filename = "output/figures/Fig2.pdf",
  height = 14,
  width = 18,
  units = "cm",
  dpi = 300
)


# ── Fig3: VE against hospitalisation ─────────────────────────────────────────

plot_ve_hosp = comb_VE %>%   
  filter(outcome != "symptomatic") %>% 
  mutate(x = paste0(vaccine, "\n", age)) %>% 
  ggplot(aes(x = month, y = mean, color = x, fill = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha   = 0.15,
              colour  = NA) + 
  geom_line(
    aes(y = mean),
    colour      = "white",
    linewidth   = 1.1,
    show.legend = FALSE
  ) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    colour     = "grey30",
    linewidth  = 0.4
  ) +
  labs(x = "Months post final dose", 
       y = "Vaccine efficacy\nagainst hospitalisation(%)") +
  facet_grid(serostatus  ~ serotype, scales = "free_y") + 
  theme(legend.position = "top") +
  scale_color_manual(values = cols[2:4]) +
  scale_fill_manual(values = cols[2:4]) +
  scale_x_continuous(breaks = c(0,12,24,36,48)) +
  guides(
    colour = guide_legend(
      title        = NULL,
      nrow         = 1,
      override.aes = list(linewidth = 1.2, fill = NA)
    ),
    fill = "none"
  )


ggsave(
  plot =  plot_ve_hosp,
  filename = "output/figures/Fig3.pdf",
  height = 14,
  width = 18,
  units = "cm",
  dpi = 300
)


# Save source data 
# write.csv(comb_VE, "output/figures/source_data/Figure2_3.csv")

