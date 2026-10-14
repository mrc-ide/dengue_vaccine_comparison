# Script to calibrate the stan model to the publicly available data of all three 
# vaccines. The script is set up to run the final model (M4) using 100 
# iterations as a demo, with 6000 iterations used for the model fits presented 
# in the manuscript. All other model variants presented in the manuscript can be 
# run below. All model variants are run using the same stan model 
# (final_model.stan) with different parameters turned on and off, using the flags. 

# ── Set up script ─────────────────────────────────────────────────────────────

# source functions 
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

n_it = 6000 ## NOTE: this could take several hours to run depending on the machine.

dir.create("output", showWarnings = F)

# DEMO RUN ON 100 ITERATIONS ───────────────────────────────────────────────────
start <- Sys.time()
dir.create("output/DEMO", showWarnings = F)

run_model_TR (
  BF = T,
  include_beta = c(3,0,0),
  delta_KJ = c(1,2,1),
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 0,
  uniform = 1,
  tau_K = c(0,0,0),
  L_K = c(0,0,0),
  folder = "DEMO",
  n_it = 100,  # 6000 to recreate manuscript results 
  adapt_delta = 0.88,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

end <- Sys.time()
end - start   # prints e.g. "Time difference of 2.4 mins"

# # ── Run main model variants ───────────────────────────────────────────────────
# 
# 
# # M1  -----------------------------------------
# run_model_TR (
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   delta_KJ = c(1,1,1),
#   L_K = c(0,0,0),
#   include_beta = c(0,0,0),
#   folder = "M1",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M2 - M1 but delta J for Dengvaxia  -----------------------------------------
# run_model_TR (
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   include_beta = c(0,0,0),
#   folder = "M2",
#   n_it = n_it,
#   adapt_delta = 0.83,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M3 - M2 but age VE for Qdenga  -----------------------------------------------
# run_model_TR (
#   include_beta = c(1,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   folder = "M3",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
#   
# # M4 - M2 but age VE for Qdenga only youngest  ---------------------------------
# run_model_TR (
#   BF = T, 
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   folder = "M4",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
#  
# # M5 - M2 but age VE for Dengvaxia  --------------------------------------------
# run_model_TR (
#   include_beta = c(0,1,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   folder = "M5",
#   n_it = n_it,
#   adapt_delta = 0.87,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# # M6 - M2 but age VE for Butantan-DV -------------------------------------------
# run_model_TR (
#   include_beta = c(0,0,1),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   folder = "M6",
#   n_it = n_it,
#   adapt_delta = 0.89,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# # M7 - M2 but single Butantan age VE ------------------------------------
# run_model_TR (
#   include_beta = c(0,0,2),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   folder = "M7",
#   n_it = n_it,
#   adapt_delta = 0.82,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# # M8 - M2 but inc_FOIJ for B -------------------------------------------------
# 
# run_model_TR (
#   inc_FOIJ = c(0,0,1),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   include_beta = c(0,0,0),
#   folder = "M8",
#   n_it = n_it,
#   adapt_delta = 0.86,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M9 - M2 but age FOI for Qdenga   ---------------------------------
# 
# run_model_TR (
#   inc_FOIJ = c(1,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   include_beta = c(0,0,0),
#   folder = "M9",
#   n_it = n_it,
#   adapt_delta = 0.84,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# # ── Run additional sensitivity analyses ───────────────────────────────────────
# 
# 
# # M10 - M4 but serotype specific L for Qdenga   ---------------------------------
# run_model_TR (
#   L_K = c(1,0,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   folder = "M10",
#   n_it = n_it,
#   adapt_delta = 0.99,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M11 - M4 but serotype specific L for Dengvaxia   ---------------------------------
# run_model_TR (
#   L_K = c(0,1,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   folder = "M11",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M12 - M4 but serotype specific L for Butantan-DV   ---------------------------------
# run_model_TR (
#   L_K = c(0,0,1),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   folder = "M12",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M13 - M4 but serotype specific tau for Q -------------------------------------
# run_model_TR (
#   tau_K = c(1,0,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   L_K = c(0,0,0),
#   folder = "M13",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# # M14 - M4 but serotype specific tau for D -------------------------------------
# run_model_TR (
#   tau_K = c(0,1,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   L_K = c(0,0,0),
#   folder = "M14",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
#  
# 
# # M15 - M4 but not serotype-specific titres ----------------------------------
# 
# run_model_TR (
#   average_mu = 1, 
#   tau_K = c(0,0,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   L_K = c(0,0,0),
#   folder = "M15",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# 
# # M16 - M15 but not serotype-specific n50 ----------------------------------
# 
# run_model_TR (
#   mono_lc_MO = 1,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   average_mu = 1, 
#   tau_K = c(0,0,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   share_n_param = 0,
#   uniform = 1,
#   L_K = c(0,0,0),
#   folder = "M16",
#   n_it = n_it,
#   adapt_delta = 0.98,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# # M17 - M16 but serotype mu ----------------------------------
# 
# run_model_TR (
#   average_mu = 0, 
#   mono_lc_MO = 1,
#   mono_lc_MU = 1,
#   mono_lc_SN = 1,
#   tau_K = c(0,0,0),
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   share_n_param = 0,
#   uniform = 1,
#   L_K = c(0,0,0),
#   folder = "M17",
#   n_it = n_it,
#   adapt_delta = 0.88,
#   stan_model = "final_model.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
# # # MS1 - M4 but severe data and age-specific psi  -----------------------------
# run_model_TR (
#   psi_J = 1,
#   BF = T, 
#   include_beta = c(3,0,0),
#   delta_KJ = c(1,2,1),
#   mono_lc_MU = 2,
#   mono_lc_SN = 2,
#   share_n_param = 0,
#   uniform = 1,
#   tau_K = c(0,0,0),
#   L_K = c(0,0,0),
#   folder = "MS1",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "final_model_severe_2.stan",
#   mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
#   baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
#   VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
#   hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
#   baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
#   cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
#   mu_De = read.csv("data/processed_De/mu.csv"),
#   cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
#   baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
#   mu_Bu = readRDS("data/processed_Bu/mu.RDS")
# )
# 
# 
