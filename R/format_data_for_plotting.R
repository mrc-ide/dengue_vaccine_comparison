# ------------------------------------------------------------------------------
# calc_attack_rates_Q
#
# Compute observed symptomatic attack rates from Qdenga VCD data across
# multiple stratifications (serostatus, trial arm, serotype, age, and time).
# Calculates exact binomial confidence intervals. 

# Outputs attack rates for comparison with model-derived attack rates. 
# ------------------------------------------------------------------------------

calc_attack_rates_Q = function(VCD, titre = F) {
  
  # VCD by serostatus and month 
  VCD_BD = VCD %>%  
    filter(serostatus != "both", age == "all", serotype== "all") %>% 
    group_by(age, serotype, year,month,serostatus) %>%  
    summarise(Y= sum(Y),
              N = sum(N)) %>% 
    mutate(trial = "both")
  
  # VCD by trial and month 
  VCD_VD = VCD %>%  
    filter(serostatus == "both", age == "all", serotype== "all") %>% 
    group_by(age, serotype,year, month,trial) %>%  
    summarise(Y= sum(Y),
              N = sum(N)) %>% 
    mutate(serostatus = "both")
  
  # VCD by serostatus and serotype and trial  
  VCD_BVK = VCD %>%  
    filter(serostatus != "both", age == "all", serotype != "all") %>% 
    group_by(age, serotype,trial,serostatus) %>%  
    summarise(Y= sum(Y),
              N = mean(N)) %>% 
    mutate(month = "all")
  
  # VCD by serostatus and age and trial  
  VCD_BVJ = VCD %>%  
    filter(serostatus != "both", age != "all", serotype == "all") %>% 
    group_by(age, serotype,trial,serostatus) %>%  
    summarise(Y= sum(Y),
              N = mean(N)) %>% 
    mutate(month = "all")
  
  # combine all data 
  
  out = VCD %>%  
    bind_rows(VCD_BD, VCD_VD,VCD_BVK,VCD_BVJ) %>% 
    mutate(outcome = "VCD") %>% 
    filter(!is.na(N)) %>%  # remove serotype age 
    mutate(age = factor(age)) %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  return(out)
}

# ------------------------------------------------------------------------------
# calc_hosp_rates_Q
#
# Compute observed hospitalised attack rates from Qdenga hospital data across
# multiple stratifications. Calculates exact binomial confidence intervals. 
#
# Outputs attack rates for comparison with model-derived attack rates. 
# ------------------------------------------------------------------------------

calc_hosp_rates_Q = function(hosp, VCD) {
  
  BVD_pop = VCD %>%  
    filter(age == "all", serotype == "all", serostatus != "both")  %>% 
    select(-Y, - X)
  
  # Hosp by serostatus, trial and month after months 
  BV_4 = hosp %>% 
    filter(age == "all", serotype != "all", month > 36) %>% 
    group_by (serostatus, trial, month) %>% 
    summarise(Y = sum(Y))  
  
  BVD = hosp %>%
    filter(age != "all", serotype == "all") %>% 
    group_by (serostatus, trial, month) %>% 
    summarise(Y=sum(Y)) %>% # BVD hosp upto 36 months  
    bind_rows(BV_4) %>% 
    left_join(BVD_pop)
  
  BD = BVD %>% 
    group_by(serostatus, month) %>% 
    summarise(Y = sum(Y),
              N= sum(N)) %>% 
    mutate(trial = "both", serotype = "all", age = "all") 
  
  VD = BVD %>% 
    group_by(trial, month) %>% 
    summarise(Y = sum(Y),
              N= sum(N)) %>% 
    mutate(serostatus = "both", serotype = "all", age = "all")
  
  # Hosp by age, serostatus, trial 
  
  BVJ_pop = VCD %>%  
    filter(age != "all", serotype == "all", serostatus != "both", month <= 36) %>% 
    group_by(age, trial, serostatus) %>% 
    summarise(N = mean(N))
  
  BVJ = hosp %>%  
    filter(age != "all", serotype == "all") %>% 
    group_by(serostatus, trial, age) %>% 
    summarise(Y = sum(Y)) %>% 
    left_join(BVJ_pop) %>% 
    mutate(month = "all",
           serotype = "all")
  
  # Hosp by serotype, serostatus, trial 
  
  BVK_pop = VCD %>%  
    filter(age == "all", serotype == "all", serostatus != "both") %>% 
    group_by(trial, serostatus) %>% 
    summarise(N = mean(N))
  
  BVK = hosp %>%  
    filter(age == "all", serotype != "all") %>% 
    group_by(serostatus, trial, serotype) %>% 
    summarise(Y = sum(Y)) %>% 
    left_join(BVK_pop)%>% 
    mutate(month = "all", 
           age = "all")
  
  # Hosp by age, serotype, serostatus, trial 
  
  BVKJ = hosp %>%  
    filter(age != "all", serotype != "all", month == "all") %>% 
    left_join(BVJ_pop) 
  
  # combine all data 
  
  out = bind_rows(VD, BD, BVD, BVJ, BVK, BVKJ) %>%
    mutate(outcome = "hosp") %>% 
    mutate(age = factor(age)) %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  return(out)
}


# ------------------------------------------------------------------------------
# calc_attack_rates
#
# Compute observed attack rates for Dengvaxia / Butantan-DV  
# This function assumes the input data is already aggregated by relevant
# stratifications (e.g., serostatus, serotype, trial arm, age, month).
# Calculates exact binomial confidence intervals. 
#
# Outputs attack rates for comparison with model-derived attack rates. 
# ------------------------------------------------------------------------------

calc_attack_rates = function(VCD, titre = F) {
  out = VCD %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  
  return(out)
}
