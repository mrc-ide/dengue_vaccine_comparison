# Dengue vaccine comparison

This repository contains the code needed to reproduce the main results presented in the manuscript:

***Modelling the efficacy of Qdenga, Dengvaxia, and the Butantan-DV dengue vaccines: a comparative analysis and open questions.***

The workflow combines publicly available clinical trial datasets with a Bayesian survival model implemented in Stan. Trial datasets are merged into a unified Stan input structure before model fitting and posterior analysis. The main analysis script (1.fit_all_trials.R) takes \~25 minutes to run 100 iterations as a demo on a standard laptop. To recreate the manuscript analyses, 6000 iterations should be run which can take several hours on a standard computer.

## Installation

### Requirements

- R (version 4.0 or later recommended)
- A working C++ toolchain, which Stan needs to compile models:
  - **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version
  - **macOS:** Xcode Command Line Tools (`xcode-select --install` in Terminal)
  - **Linux:** `g++` and `make` (e.g. `sudo apt install build-essential`)

### 1. Install CRAN packages

``` r
install.packages(c(
  "tidyverse",   # includes dplyr, tidyr, ggplot2
  "Hmisc",
  "patchwork",
  "cowplot",
  "posterior",
  "bayesplot",
  "loo",
  "truncnorm",
  "logspline"
))
```

The `grid` package is part of base R and does not need to be installed.

### 2. Install cmdstanr and CmdStan

`cmdstanr` is not on CRAN. Install it from the Stan R-universe repository:

``` r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
```

Then check your toolchain and install CmdStan, the backend that compiles and runs Stan models:

``` r
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
install_cmdstan(cores = 2)
cmdstan_version()  # should print the installed CmdStan version
```

Installing CmdStan can take several minutes. See the [cmdstanr documentation](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) if you run into problems.

This analysis has been tested on R version 4.6.1 (2026-06-24) and R version 4.4.3 (2025-02-28).

## Overview of workflow

1.  **Trial-specific data formatting**\
    Each vaccine dataset is summarised by age, serostatus, serotype, trial arm, and time. Counts are reshaped into multidimensional arrays required by Stan. Vaccine-specific formatted datasets are merged into a single Stan data list.

2.  **Model calibration**\
    The Stan model is compiled and sampled using `cmdstanr`. Posterior draws, derived quantities, diagnostics, and evaluation metrics are saved.

3.  **Posterior analysis & visualisation**\
    Model outputs are converted into interpretable vaccine efficacy and attack rate summaries and plotted.

Scripts should be run in order (1–5) to reproduce manuscript figures.

## Repository structure

### Scripts

Scripts execute the full pipeline from model calibration to figure generation:

- [**1.fit_all_trials.R**](http://1.fit)\
  Runs the multi-trial Bayesian calibration using all vaccine datasets. The script is configured to run the final model (M4) first. Alternative model variants and sensitivity analyses are implemented via configuration flags passed to the Stan model.\
  → Outputs posterior generated quantities (attack rates, vaccine efficacy, antibody titres estimates for each vaccine - e.g., `AR_Bu.RDS`, `n_Bu.RDS` and `VE_Bu.RDS`), parameter posterior chains (`posterior_chains.csv`), parameter posterior summarises (`posterior.csv`), and log-likelihood and WAIC model comparison metrics (`loglik.RDS` and `WAIC.RDS`) in each model's output directory (e.g., `output/M4/`).

- **2.plot_all_attack_rates_scatter.R**\
  For each vaccine, joins posterior attack rate estimates to observed trial data and produces observed vs. model-estimated scatter plots. Saves a combined faceted ggplot object (by vaccine) to disk for use in the next script.\
  → Outputs `output/figures/all_attack_rates.rds`

- **3.plot_all_attack_rates.R**\
  Loads the scatter plot produced by script 2 and combines it with per-vaccine time-course and serotype attack rate plots into a single assembled figure.\
  → Outputs manuscript **Fig. 1** (`output/figures/Fig1.pdf`)

- **4.plot_efficacy_over_time.R**\
  Plots model-estimated vaccine efficacy trajectories up to 54 months post final dose, stratified by serostatus, serotype, and age group, for all three vaccines. Produces separate figures for symptomatic and hospitalised outcomes.\
  → Outputs manuscript **Fig. 2** (`output/figures/Fig2.pdf`) and **Fig. 3** (`output/figures/Fig3.pdf`)

- **5.calculate_plot_ve_sensitivity_analysis.R**\
  Loads posterior VE draws from the baseline model (M4) and a set of structural sensitivity analysis variants. Computes age-weighted VE differences between each alternative model and the baseline, and visualises results as diverging heatmaps stratified by serostatus, serotype, and outcome. Individual vaccine heatmaps and a combined figure are saved.\
  → Outputs manuscript **Fig. 4** (`output/figures/Fig4.pdf`) and individual heatmaps per vaccine

### Data

All datasets are preprocessed aggregations of published phase III clinical trial results.

Where possible case counts are extracted by month, serotype, infecting serotype, age group and clinical outcome (symptomatic, hospitalised, and severe VCD).

#### Qdenga

- `processed_Q/hosp_data.csv` — hospitalised VCD cases\
- `processed_Q/vcd_data.csv` — symptomatic VCD cases\
- `processed_Q/n0_new.csv` — initial fitted neutralising titres\
- `processed_Q/seropositive_by_age_baseline.csv` — age-stratified baseline serostatus counts

#### Dengvaxia

- `processed_De/cases_stan_format.RDS` — symptomatic/hospital/severe VCD cases
- `processed_De/mu.csv` — initial fitted neutralising titres\
- `processed_De/baseline_SP.csv` — baseline serostatus counts

#### Butantan-DV

- `processed_Bu/cases_stan_format.RDS` — symptomatic VCD cases\
- `processed_Bu/mu.RDS` — neutralizing titres\
- `processed_Bu/baseline_SP.RDS` — baseline serostatus counts

------------------------------------------------------------------------

### Models

Stan models implement the Bayesian survival framework used for calibration.

- **final_model.stan**\
  Core multi-trial survival model with binomial and multinomial likelihoods. Configuration flags enable all model variants tested in the manuscript through script **1.fit_all_trials.R**

- **final_model_severe.stan**\
  Extension of the main model including Dengvaxia severe outcome likelihoods.

------------------------------------------------------------------------

### R

Functions implement data preprocessing, formatting, model fitting, and plotting:

- **format_stan_data.R**\
  Aggregates case counts, reshapes multidimensional arrays, and combines priors and flags into Stan-ready inputs.

- **run_model_fitting_joint.R**\
  High-level wrapper that preprocesses data, constructs Stan inputs, runs MCMC sampling via `cmdstanr`, and saves posterior outputs and diagnostics.

- **format_model_outputs_for_plotting.R**\
  Formats model outputs (attack rates or vaccine efficacy) ready for plotting. Calculates absolute difference in vaccine efficacy estimates for sensitivity analysis.

- **format_data_for_plotting.R**\
  Formats raw case data into observed attack rates ready for plotting.

- **calculate_bayes.R** Computes Savage-Dickey Bayes factors for the waning parameter (L) across all three vaccines. Estimates prior and posterior densities at zero using a truncated normal prior and logspline density estimation, and writes the results to BF.csv.

------------------------------------------------------------------------

## Output

Model fitting generates:

- posterior chains and summaries\
- vaccine efficacy, titres, and attack rate posterior estimates\
- WAIC/log-likelihood diagnostics
