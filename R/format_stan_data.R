# ------------------------------------------------------------------------------
# format_stan_data_Q_De
#
# Construct a unified Stan data list by combining formatted datasets from
# multiple vaccine trials (Qdenga, Dengvaxia, and Butantan-DV). This function
# combines trial-specific preprocessing outputted from the below functions
# (e.g., `format_stan_data_Q()` and `format_stan_data_De()`), and assembles
# shared model parameters and flags into a single structure.
# 
# The returned list is ready for direct use in the multi-trial Stan model.
# ------------------------------------------------------------------------------

format_stan_data_Q_De = function(baseline_SP_Q,
                                 VCD_Q,
                                 hosp_Q,
                                 mu_Q,
                                 VCD_years_Q,
                                 time_Q, 
                                 baseline_SP_De, 
                                 cases_De,
                                 mu_De,
                                 cases_Bu,
                                 baseline_SP_Bu,
                                 mu_Bu,
                                 time_Bu,
                                 VCD_years_Bu, 
                                 VCD_years_De, 
                                 time_De, 
                                 C,   
                                 B,
                                 K,
                                 V,
                                 J_Q,
                                 J_De,
                                 J_Bu,
                                 R,
                                 TR,
                                 HI,
                                 uniform,
                                 include_pK3,
                                 include_eps,
                                 L_mean,
                                 L_sd,
                                 enhancement,
                                 mono_lc_MO,
                                 mono_lc_SN,
                                 mono_lc_MU,
                                 mono_lc,
                                 rho_K,
                                 w_CK,
                                 alpha_CK,
                                 tau_K,
                                 include_beta,
                                 L_K,
                                 delta_KJ,
                                 psi_J,
                                 inc_FOIJ,
                                 share_n_param,
                                 share_omega_kappa,
                                 share_alpha
) {
  
  list_Q_data = format_stan_data_Q(
    baseline_SP = baseline_SP_Q,
    VCD = VCD_Q,
    hosp = hosp_Q,
    B = B,
    R = R,
    K = K,
    V = V,
    J = J_Q,
    C = C,
    HI = HI,
    include_pK3 = include_pK3,
    include_eps = include_eps,
    include_beta = include_beta[1],
    mono_lc_MO = mono_lc_MO,
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    mono_lc = mono_lc,
    rho_K = rho_K,
    L_K = L_K[1],
    w_CK = w_CK,
    delta_KJ = delta_KJ[1],
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    L_sd = L_sd,
    L_mean = L_mean,
    enhancement = enhancement,
    VCD_years = VCD_years_Q,
    time = time_Q,
    mu = mu_Q)
  
  list_De_data = format_stan_data_De(
    baseline_SP = baseline_SP_De,
    cases = cases_De,
    VCD_years = VCD_years_De,
    time = time_De,
    mu = mu_De,
    B = B,
    R = R,
    K = K,
    V = V,
    J = J_De,
    C = C,
    HI = HI,
    include_pK2 = include_pK3,
    include_eps = include_eps,
    include_beta = include_beta[2],
    mono_lc_MO = mono_lc_MO,
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    mono_lc = mono_lc, 
    rho_K = rho_K,
    L_K = L_K[2],
    w_CK = w_CK,
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    delta_KJ = delta_KJ[2],
    L_mean = L_mean,
    L_sd = L_sd,
    enhancement = enhancement,
    psi_J = psi_J
  )
  
  # Butantan-DV data -----------------------------------------------------------
  
  D_Bu = length(VCD_years_Bu)
  
  N_VCD_BVKD_Bu = cases_Bu$Sy_BVK %>% 
    arrange(time, serostatus, arm, serotype)
  
  VCD_BVKD_Bu = array(N_VCD_BVKD_Bu$Y, dim = c(B * 2 * V, D_Bu))
  
  N_VCD_BVJD_Bu = cases_Bu$Sy_BVJ %>% 
    arrange(time, serostatus, arm, age) 
  
  VCD_BVJD_Bu = array(N_VCD_BVJD_Bu$Y, dim = c(B * J_Bu * V, D_Bu))
  
  N_VCD_BVKJ_Bu = cases_Bu$Sy_BVKJ %>% 
    arrange(serostatus, arm, serotype, age) 
  
  VCD_BVKJ_Bu = array(N_VCD_BVKJ_Bu$Y, dim = c(B * J_Bu * V *2))
  
  N_pop_BVJD_Bu = cases_Bu$Sy_BVJ %>% 
    select(time, serostatus, age, arm, N) %>% 
    arrange(time, age, arm, serostatus)
  
  m_pop_BVJD_Bu = array(N_pop_BVJD_Bu$N, dim = c(B,V,J_Bu,D_Bu))
  
  popD_Bu = N_VCD_BVKD_Bu %>% filter(serotype == "DENV1") %>%  group_by(time) %>%  summarise(N=sum(N))
  
  VCD_D_Bu = N_VCD_BVKD_Bu %>%  group_by(time) %>%  summarise(Y=sum(Y))
  
  
  SP_J_Bu = baseline_SP_Bu %>% 
    arrange(age) %>% 
    pull(SP)
  
  pop_J_Bu = baseline_SP_Bu %>% 
    arrange(age) %>% 
    pull(N)
  
  # combine MU -----------------------------------------------------------------
  mu = array(NA, dim = c(TR, B, K)) # Create an empty array of correct size
  mu[1,,] = list_Q_data$mu  # Assign Qdenga data to TR=1
  mu[2,,] = list_De_data$mu   # Assign Dengvaxia data to TR=2 
  mu[3,,] = mu_Bu   # Assign Butantan-DV data to TR=3 
  mu[3,,3:4] = 0   # Serotypes 3 and 4 are null
  
  
  # overall stan data ----------------------------------------------------------
  
  data_list = list(
    C = C,   # No. current immune status
    B = B,   # No. baseline serostatus
    K = K,   # No. serotypes
    D_Q = list_Q_data$D, # No. Qdenga VCD data
    D_De = list_De_data$D, # No. Dengvaxia VCD data
    D_Bu = D_Bu, # No. Dengvaxia VCD data
    V = V,   # No. trial arms
    J_Q = list_Q_data$J, # No. age groups (Qdenga)
    J_De = list_De_data$J, # No. age groups (Dengvaxia)
    J_Bu = J_Bu, # No. age groups (Butantan)
    T_Q = list_Q_data$T, # No. Q model time points
    T_De = list_De_data$T, # No. D model time points
    T_Bu = length(time_Bu), # No. B model time points
    R = R,   # No. outcomes
    TR = TR,   # No. vaccine trials (Q, D, B)
    time_Q = list_Q_data$time,   # Time points for Q model
    time_De = list_De_data$time,  # Time points for D model
    time_Bu = time_Bu,  # Time points for B model
    HI = HI,       # Period of heterotypic immunity
    mu = mu,       # Titres at T0 (array[TR,B,K])
    
    # Butantan-DV data
    SP_J_Bu = SP_J_Bu, # Baseline seropositive (array[J_Bu])
    pop_J_Bu = pop_J_Bu,  # Baseline population (array[J_Bu])
    VCD_D_Bu =  VCD_D_Bu$Y, # Total VCD (array[D_Bu])
    VCD_BVKD_Bu = VCD_BVKD_Bu, # VCD (array[B*V*2, D_Bu])
    VCD_BVJD_Bu = VCD_BVJD_Bu, # VCD (array[B*V*J_Bu, D_Bu])
    VCD_BVKJ_Bu = VCD_BVKJ_Bu, # VCD (array[B*V*2*J_Bu])
    pop_BVJD_Bu = m_pop_BVJD_Bu, # Population (array[B*V*J_Bu])
    popD_Bu = popD_Bu$N, # Total population
    
    
    # Dengvaxia data
    SP_J_De = list_De_data$SP_J,  # Baseline seropositive (array[J_De])
    pop_J_De = list_De_data$pop_J, # Baseline population (array[J_De])
    VCD_De = list_De_data$VCD,   # Total VCD (int)
    HOSP_D_De = list_De_data$HOSP_D, # Total hospitalizations (array[D_De])
    VCD_VK_De = list_De_data$VCD_VK, # VCD (array[V*K])
    VCD_BVJ_De = list_De_data$VCD_BVJ, # VCD (array[B*V*J_De])
    HOSP_BVJD_De = list_De_data$HOSP_BVJD, # Hospitalizations (array[B*V*J_De,D_De])
    HOSP_BVKJ_De = list_De_data$HOSP_BVKJ, # Hospitalizations (array[B*V*K*J_De])
    pop_VCD_BVJ_De = list_De_data$pop_VCD_BVJ, # Population (array[B,V,J_De])
    pop_HOSP_BVJD_De = list_De_data$pop_HOSP_BVJD, # Population (array[B,V,J_De,D_De])
    pop_VCD_De = list_De_data$pop_VCD, # Population for VCD (int)
    pop_HOSP_De = list_De_data$pop_HOSP, # Mean hospitalization population across trial (int)
    pop_HOSP_D_De = list_De_data$pop_HOSP_D, # Hospitalization population over time (array[D_De])
    
    # Severe data for sensitivity analysis 
    SE_BVKJ_De = list_De_data$SE_BVKJ,
    SE_SN_VJD_De = list_De_data$SE_SN_VJD,
    SE_SPvJ_De = list_De_data$SE_SPvJ,
    SE_SN_D_De = list_De_data$SE_SN_D, # SN severe over time 
    SE_De = list_De_data$SE, 
    
    # Qdenga data
    SP_J_Q = list_Q_data$SP_J,  # Baseline seropositive (array[J_Q])
    pop_J_Q = list_Q_data$pop_J, # Baseline population (array[J_Q])
    VCD_D_Q = list_Q_data$VCD_D, # Total VCD (array[D_Q])
    HOSP_D_Q = list_Q_data$HOSP_D, # Total hospitalizations (array[D_Q])
    VCD_BVKD_Q = list_Q_data$VCD_BVKD, # VCD (array[B*V*K,D_Q])
    VCD_BVJD_Q = list_Q_data$VCD_BVJD, # VCD (array[B*V*J_Q,D])
    VCD_KJ2_Q = list_Q_data$VCD_KJ2, # VCD at 12 and 24 months (array[K*J_Q,2])
    HOSP_BVJD_Q = list_Q_data$HOSP_BVJD, # Hospitalizations (array[B*V*J_Q,D])
    HOSP_BVK4_Q = list_Q_data$HOSP_BVK4, # Hospitalizations (array[B*V*K,4])
    HOSP_KJ2_Q = list_Q_data$HOSP_KJ2, # Hospitalizations at 12 and 24 months (array[K*J_Q,2])
    pop_Q = list_Q_data$pop, # Population (array[B,V,J_Q,D_Q])
    pop_BVD_Q = list_Q_data$pop_BVD, # Population (array[B,V,D_Q])
    VCD_BVKJ_Q = list_Q_data$VCD_BVKJ_Q,
    HOSP_BVKJ_Q = list_Q_data$HOSP_BVKJ_Q,
    
    # Model flags
    share_alpha = share_alpha, # share alpha parameters across trial (T/F)
    share_omega_kappa = share_omega_kappa, # share omega and kappa parameters across trials (T/F)
    uniform = uniform, # Use uniform priors on probabilities of reporting? (0/1)
    include_pK3 = include_pK3, # Include serotype-specific probability? (0/1)
    include_eps = include_eps, # Include enhanced secondary hospitalization? (0/1)
    L_mean = L_mean, # Mean of L prior (real)
    L_sd = L_sd, # SD of L prior (real)
    enhancement = enhancement, # Vaccine enhancement of SN? (0/1)
    mono_lc = mono_lc, # Serostatus, age, serotype specific lc? (0/1)
    mono_lc_MO = mono_lc_MO, # Serotype-specific lc MO? (0/1)
    mono_lc_SN = mono_lc_SN, # Serotype-specific lc SN? (0/1/2)
    mono_lc_MU = mono_lc_MU, # Serotype-specific lc MU? (0/1/2)
    rho_K = rho_K, # Mono/serotype-specific rho (0/1)
    w_CK = w_CK, # Mono/serostatus/serotype w (0/1/2)
    alpha_CK = alpha_CK, # Mono/serostatus/serotype alpha (0/1/2)
    tau_K = tau_K, # Mono/serotype-specific tau (0/1)
    share_n_param = share_n_param, #  // share hs, hl, ts parameters across trials (0,1)
    
    # Trial specific flags 
    inc_FOIJ = inc_FOIJ, # Different FOI in youngest age group 
    include_beta = include_beta, # Include age-specific nc50? (array[TR], 0-3)
    L_K = L_K, # Mono/serotype-specific L (array[TR], 0/1)
    delta_KJ = delta_KJ, # Single/serotype-specific/age delta (array[TR], 0-2)
    psi_J = psi_J
  )
  return(data_list)
}


# ------------------------------------------------------------------------------
# format_stan_data_Q
#
# Prepare and aggregate trial baseline serostatus, VCD, hospitalisation, and 
# titre Qdenga trial datasets into the structured arrays required by Stan. 
# The function also passes model parameter flags and prior parameters.
# 
# Outputs a single named list formatted for Stan 
# ------------------------------------------------------------------------------

format_stan_data_Q = function(baseline_SP,
                              VCD,
                              hosp,
                              mu,
                              B,
                              K,
                              R,
                              C,
                              HI,
                              V,
                              J,
                              VCD_years,
                              time,
                              include_pK3,
                              include_eps,
                              include_beta,
                              mono_lc, 
                              mono_lc_MO,
                              mono_lc_SN,
                              mono_lc_MU,
                              rho_K,
                              L_K,
                              delta_KJ,
                              w_CK,
                              alpha_CK,
                              tau_K,
                              L_sd,
                              L_mean,
                              enhancement) {
  
  T = length(time)
  D = length(VCD_years)

  # pop by age group at trial start
  pop_J = VCD %>%
    filter(age != "all", year == 1, serotype == "all") %>%
    group_by(age) %>%
    summarise(N = sum(N))
  
# calculate # VCD over time ---------------------------------------------------
  VCD_D = VCD %>%
    filter(serostatus == "both",
           age == "all",
           serotype == "all",
           year %in% VCD_years) %>%
    group_by(year) %>%
    summarise(Y = sum(Y)) %>%
    select(Y)
  
# calculate # hosp over time ---------------------------------------------------
 
  hosp_D =
    hosp %>%
    filter(age != "all", trial != "both", month != "all") %>%
    group_by(year) %>%
    summarise(Y = sum(Y)) %>%
    select(Y) 
  
# calculate # pop by trial arm, serostatus, age group, over time ---------------
  
  N_pop_BVJD = VCD %>%
    filter(serostatus != "both",
           age != "all",
           serotype == "all",
           trial != "both") %>%
    arrange(year, age, trial , serostatus)
  
  pop_BVJD = array(N_pop_BVJD$N, dim = c(B, V, J, D))
  
# calculate # pop by trial arm, serostatus, over time  -------------------------
  
  N_pop_BVD = N_pop_BVJD %>%  
    group_by(serostatus, year, trial) %>%  
    summarise(N=sum(N))%>%
    arrange(year, trial , serostatus)
  
  pop_BVD = array(N_pop_BVD$N, dim = c(B, V, D))
  
# calculate # hosp by serostatus + trial arm + serotype, over time (24-54 months)

  N_hosp_BVK4 = hosp %>%
    filter(serotype != "all", trial != "both", month != "all") %>%
    arrange(year, serostatus, trial , serotype)
  
  hosp_BVK4_m = array(N_hosp_BVK4$Y, dim = c(B * K * V, 4))  # (D=4)
  
# calculate # VCD by serostatus + trial arm + serotype, over time --------------

    # for multinomial likelihood
  N_VCD_BVKD_m = VCD %>%
    filter(serostatus != "both",
           age == "all",
           serotype != "all",
           trial != "both",
           year %in% VCD_years) %>%
    arrange(year, serostatus, trial , serotype)
  
  VCD_BVKD_m = array(N_VCD_BVKD_m$Y, dim = c(B * K * V, D))
  
# calculate # hosp by serostatus + trial arm + age-group, over time -------------

    N_hosp_BVJD_m =  hosp %>%
    filter(serotype == "all", trial != "both", age != "all") %>%
    arrange(year, serostatus, trial , age)
  
  hosp_BVJD_m = array(N_hosp_BVJD_m$Y, dim = c(B * V * J, D))
  
# calculate # VCD by serostatus + trial arm + age-group, over time -------------

  # for multinomial likelihood
  N_VCD_BVJD_m = VCD %>%
    filter(serostatus != "both",
           age != "all",
           serotype == "all",
           trial != "both") %>%
    arrange(year, serostatus, trial , age)
  
  VCD_BVJD_m = array(N_VCD_BVJD_m$Y, dim = c(B * V * J, D))
  
# VCD age and serotype year and 2 ----------------------------------------------

  N_VCD_KJ2_m = VCD %>%
    filter(age != "all", serotype != "all", month != "all") %>%
    arrange(year, serotype, age)
  
  VCD_KJ2_m = array(N_VCD_KJ2_m$Y, dim = c(K * J, 2))
  
# hosp age and serotype year and 2 ---------------------------------------------
 
  N_hosp_KJ2_m = hosp %>%
    filter(age != "all", serotype != "all", month != "all") %>%
    arrange(year, serotype, age)
  
  hosp_KJ2_m = array(N_hosp_KJ2_m$Y, dim = c(K * J, 2))

# VCD age, serotype trial, serostatus, all time --------------------------------
  N_VCD_BVKJ_Q_m = VCD %>% 
    filter(age != "all", serotype != "all", month == "all") %>%
    arrange(serostatus, trial, serotype, age)
    

  # hosp age, serotype trial, serostatus, all time -------------------------------
 
   N_HOSP_BVKJ_Q_m = hosp %>% 
    filter(age != "all", serotype != "all", month == "all") %>%
    arrange(serostatus, trial, serotype, age)

  # data -------------------------------------------------------------------------
  
  stan_data = list(
    time = time,
    T = T,
    J = J,
    K = K,
    B = B,
    V = V,
    D = D,
    C = C,
    HI = HI,
    R = R,
    include_pK3 = include_pK3,
    include_eps = include_eps,
    include_beta = include_beta,
    mono_lc = mono_lc, 
    mono_lc_MO = mono_lc_MO,
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    rho_K = rho_K,
    L_K = L_K,
    w_CK = w_CK,
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    delta_KJ = delta_KJ, 
    L_sd = L_sd,
    L_mean = L_mean,
    enhancement = enhancement, 
    SP_J = baseline_SP$Seropositive,
    pop_J = baseline_SP$Total,
    VCD_D = as.numeric(VCD_D$Y),
    pop = pop_BVJD,
    pop_BVD = pop_BVD,
    VCD_BVKD = VCD_BVKD_m,
    VCD_BVJD = VCD_BVJD_m,
    mu = mu,
    HOSP_BVJD = hosp_BVJD_m,
    HOSP_BVK4 = hosp_BVK4_m,
    HOSP_D = as.numeric(hosp_D$Y),
    HOSP_KJ2 = hosp_KJ2_m,
    VCD_KJ2 = VCD_KJ2_m,
    VCD_BVKJ_Q = N_VCD_BVKJ_Q_m$Y,
    HOSP_BVKJ_Q = N_HOSP_BVKJ_Q_m$Y
  )
  return(stan_data)
}

# ------------------------------------------------------------------------------
# format_stan_data_De
#
# Prepare and aggregate trial baseline serostatus, VCD, hospitalisation, and 
# titre Dengvaxia trial datasets into the structured arrays required by Stan. 
# The function also passes model parameter flags and prior parameters.
#
# Outputs a single named list formatted for Stan 
# ------------------------------------------------------------------------------

format_stan_data_De = function(baseline_SP,
                               cases,
                               mu,
                               B,
                               K,
                               R,
                               C,
                               V,
                               J,
                               HI,
                               VCD_years,
                               time,
                               include_pK2,
                               include_eps,
                               include_beta,
                               mono_lc, 
                               mono_lc_MO,
                               mono_lc_SN,
                               mono_lc_MU,
                               rho_K,
                               L_K,
                               w_CK,
                               alpha_CK,
                               tau_K,
                               delta_KJ,
                               psi_J,
                               L_mean,
                               L_sd,
                               enhancement) {
  
  T = length(time)
  D = length(VCD_years)
  
  # SP at baseline 
  SP_J = baseline_SP %>%  # make sure youngest first 
    mutate(age = factor(age, levels = c("2-8yrs", "9-16yrs"))) %>% 
    arrange(age)
  
  # calculate total hosp over time
  HOSP_D = cases$Ho_BVJD %>% 
    group_by(time) %>%  
    summarise(Y=sum(Y))
  
  # calculate total severe over time 
  SE_SN_D = cases$Se_VJD %>% 
    group_by(time) %>% 
    summarise(Y=sum(Y))
  
  # Hosp BVJD 
  HOSP_BVJD = cases$Ho_BVJD %>% 
    arrange(time, serostatus, arm , age)
  
  m_HOSP_BVJD = array(HOSP_BVJD$Y, dim = c(B * V * J, D))
  
  # Hosp population over time 
  Pop_HOSP_BVJD = HOSP_BVJD %>%
    arrange(time, age, arm, serostatus)
  
  m_Pop_HOSP_BVJD = array(Pop_HOSP_BVJD$N, dim = c(B, V, J, D))
  
  # Total hosp population over time 
  Pop_HOSP_D = HOSP_BVJD %>%
    group_by(time) %>%  summarise(N = sum(N))
  
  # Hosp BVJK   
  HOSP_BVKJ = cases$Ho_BVKJ %>%
    arrange(serostatus, arm , age, serotype) 
  
  m_HOSP_BVKJ = array(HOSP_BVKJ$Y, dim = c(B * V * J * K))
  
  # Symp BVJ 
  VCD_BVJ = cases$Sy_BVJ %>%  
    arrange(serostatus, arm , age)
  
  m_VCD_BVJ = array(VCD_BVJ$Y, dim = c(B * V * J))
  
  # Symp VK (CYD14 and 15 only, PP pop)  
  VCD_VK = cases$Sy_VK %>%
    arrange(arm, serotype)
  
  m_VCD_VK = array(VCD_VK$Y, dim = c(V*K))
  
  # Symp population  (CYD14 and 15 only, ITT pop, single time point )  
  Pop_VCD_BVJ = VCD_BVJ %>%
    arrange(age, arm, serostatus)
  
  m_Pop_VCD_BVJ = array(Pop_VCD_BVJ$N, dim = c(B, V, J))
  
  # Severe SN VJD 
  SE_SN_VJD = cases$Se_VJD %>%  
    arrange(time, arm, age)
  
  m_SE_SN_VJD = array(SE_SN_VJD$Y,dim = c(V * J, D))
  
  # Severe BVJK 
  SE_BVKJ = cases$Se_BVKJ %>%   # Vaccine SP is missing (final 8)
    filter(serotype!= "all") %>%  
    arrange(serostatus, arm , age, serotype) 
  
  m_SE_BVKJ = array(SE_BVKJ$Y, dim = c(B * V * J * K - 8))
  
  # Severe SPvJ (seropos vaccine by age)
  SE_SPvJ = cases$Se_BVKJ %>%  # Vaccine SP is missing (final 8)
    filter(serotype == "all") %>%  
    arrange(age) 
  
  m_SE_SPvJ = array(SE_SPvJ$Y)
  
  # Initial antibody titres by serostatus and serotype 
  mu_t = mu %>%
    select(-X) %>% 
    as.matrix()
  
  # data -------------------------------------------------------------------------
  stan_data = list(
    time = time,
    T = T,
    J = J,
    K = K,
    B = B,
    V = V,
    D = D, 
    C = C,
    HI = HI, 
    R = R, 
    include_pK2=include_pK2,
    include_eps=include_eps,
    include_beta=include_beta,
    mono_lc = mono_lc,
    mono_lc_MO = mono_lc_MO, 
    mono_lc_SN=mono_lc_SN,
    mono_lc_MU=mono_lc_MU,
    rho_K=rho_K,
    L_K=L_K,
    w_CK=w_CK,
    alpha_CK=alpha_CK,
    tau_K=tau_K,
    delta_KJ=delta_KJ,
    psi_J=psi_J,
    L_mean=L_mean,
    L_sd=L_sd,
    enhancement=enhancement,
    SP_J = SP_J$SP,
    pop_J =SP_J$N, 
    VCD = sum(cases$Sy_VK$Y),
    HOSP_D = as.numeric(HOSP_D$Y),
    SE_SN_D = as.numeric(SE_SN_D$Y), # SN severe over time 
    SE = sum(cases$Se_BVKJ$Y), # total severe 
    VCD_BVJ = m_VCD_BVJ,  
    VCD_VK = m_VCD_VK,
    HOSP_BVJD = m_HOSP_BVJD,
    HOSP_BVKJ = m_HOSP_BVKJ,
    SE_BVKJ = m_SE_BVKJ,
    SE_SN_VJD = m_SE_SN_VJD,
    SE_SPvJ = m_SE_SPvJ,
    pop_VCD_BVJ = m_Pop_VCD_BVJ,
    pop_HOSP_BVJD = m_Pop_HOSP_BVJD,
    pop_VCD = sum(m_Pop_VCD_BVJ),
    pop_HOSP_D = as.numeric(Pop_HOSP_D$N),
    pop_HOSP = mean(Pop_HOSP_D$N),
    mu = mu_t
  )
  
  return(stan_data)
}



# ------------------------------------------------------------------------------
# factor_VCD_Q
#
# Standardise categorical fields in Qdenga VCD trial data by converting age,
# serostatus, serotype, and trial arm variables into consistently ordered
# factors - applied immediately after importing raw VCD data.
# ------------------------------------------------------------------------------

factor_VCD_Q = function(VCD){
  VCD = VCD  %>%  
    mutate(age = factor(age, levels = c("4-5yrs", "6-11yrs", "12-16yrs", "all")), 
           serostatus = factor(serostatus,
                               levels = c("seronegative", "seropositive", "both")),
           serotype = factor(serotype,levels = c("DENV1", "DENV2", "DENV3", "DENV4", "all")),
           trial = factor(trial, levels = c("placebo", "vaccine", "both")))
  
  return(VCD)
}



# ------------------------------------------------------------------------------
# factor_serology
#
# Standardise categorical fields in serology datasets by converting serostatus,
# serotype, and trial arm variables into consistently labeled factors - use
# immediately after importing data. 
# ------------------------------------------------------------------------------

factor_serology= function(data){
  data = data  %>%  
    mutate(serostatus = factor(serostatus,
                               levels = c("SN", "SP"),
                               labels = c("seronegative", "seropositive")),
           serotype = factor(serotype),
           trial = factor(trial, levels = c("Placebo", "TAK"), 
                          labels = c("placebo", "vaccine")))
  
  return(data)
}


