library(cmdstanr)

source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")

whales_init = whales_clean %>% filter(Number == 1) %>% filter(ID %in% c(27,12,36,32,16,
                                                                      29,7,9,13,19))

whales_seq = whales_clean%>% filter(Number != 1) %>% filter(ID %in% c(27,12,36,32,16,
                                                                      29,719,13,19))

# This is the good one #### I think... need to recheck
mod1_v1.1 = cmdstan_model("mod1_v1-1.stan")

sim = mod1_v1.1$sample(data = stan_data,
                       iter_warmup = 500,
                       iter_sampling = 500,
                       #parallel_chains = 4,
                       #init = list(init_list), 
                       adapt_delta =.9, # with .95, we get 5% divergence
                       seed = 1805314595,
                       #seed = 2026228722, 
                       chains = 1
)
#adapt_delta = .95
#2026228722, 35%
#1524746552, 13% (the best so far)
#1653525333, 5% - 27 transitions (the best so far)
#1805314595, 0% (the best so far)

#adapt_delta = .8
#28508665 91%
test = sim$metadata() ### thissss
test$seed
# Problem is with the Poisson distribution
print(sim,max_rows=52) # .95= 13% transitions ended with a divergence

# July 27 2023 ####
#saveRDS(sim,"mod1_v1_output_stan_1chain.rds")

mod1.1ch = readRDS("mod1_v1_output_stan_1chain.rds")
mod1.4ch = readRDS("mod1_v1_output_stan.rds")

mod1.1ch
print(mod1.1ch,max_rows=52) # .95= 13% transitions ended with a divergence

# With one chain everything was pretty fine, with more chains it broke, bad R-hats ###



### July 21, 2023 ### adapt_delta = .8

# All 4 chains finished successfully.
# Mean chain execution time: 1163.6 seconds.
# Total execution time: 1918.7 seconds.
# 
# Warning: 726 of 2000 (36.0%) transitions ended with a divergence.
# This may indicate insufficient exploration of the posterior distribution.
# Possible remedies include: 
#   * Increasing adapt_delta closer to 1 (default is 0.8) 
# * Reparameterizing the model (e.g. using a non-centered parameterization)
# * Using informative or weakly informative prior distributions 

sim

test

sim.85 = mod1_v1.1$sample(data = stan_data,
                       iter_warmup = 500,
                       iter_sampling = 500,
                       parallel_chains = 1,
                       adapt_delta = .95
)

print(sim.85,max_rows=82)
### July 21, 2023 ### adapt_delta = .85, no set seed

# All 4 chains finished successfully.
# Mean chain execution time: 1592.9 seconds.
# Total execution time: 2557.2 seconds.
# 
# Warning: 1574 of 2000 (79.0%) transitions ended with a divergence.
# This may indicate insufficient exploration of the posterior distribution.
# Possible remedies include: 
#   * Increasing adapt_delta closer to 1 (default is 0.8) 
# * Reparameterizing the model (e.g. using a non-centered parameterization)
# * Using informative or weakly informative prior distributions

print(sim.85,max_rows=49)

saveRDS(sim.85,"mod1_v1-1_adapt_delta_85_output_stan.rds")



mod2_mle = cmdstan_model("mod2/mod2_mle.stan")

mod2_mle_result = mod2_mle$optimize(data = stan_data,
                                          init =  list(init_list_m2),
                                          seed = 7483,
                                          # seed = 3284,
                                          # algorithm = "newton",
                                          #init_alpha = .0005
)

init_list_m2 = list(
  # duration
  mu_duration = c(140.0,334.0,516.0),
  #mu_duration = c(100.0,100.0,100.0),
  #log_mu_duration = log(c(140.0,334.0,516.0)),
  log_sigma_duration = log(c(80.0,212.0,130.0)),
  
  # surface
  mu_surface = c(70.0,86.0,151.0),
  # log_mu_surface = log(c(70.0,86.0,151.0)),
  log_sigma_surface = log(c(68.0,55.0,69.0)),

  #max depth
  mu_maxDepth = c(32.0,68.0,170.0),
  # log_mu_maxDepth = log(c(32.0,68.0,170.0)),
  log_sigma_maxDepth = log(c(24.0,65.0,60.0)),

  #step length
  mu_step = c(189.0,675.0,406.0),
  # log_mu_step = log(c(189.0,675.0,406.0)),
  log_sigma_step = log(c(134.0,305.0,287.0)),

  #turning angle
  log_kappa = log(c(1.0,3.1,.8)),

  #heading variance
  log_a = log(c(1.0,.5,1.7)),
  log_b = log(c(2.1,5.4,1.6)),

  #number of lunges
  log_lambda = log(c(.7,.05,3.4)),

  #initial distribution
  init_raw = qlogis(c(1/3,1/3)),
  # init = c(1/3,1/3,1/3),
  # init_raw = matrix(qlogis(1/3),3,2),

  # pi_raw = qlogis(c(1/3,1/3)),
  # pi_raw = qlogis(c(.9999)),

  #rows tpm
  # tpm1 = matrix(0,3,2),
  # tpm2 = matrix(0,3,2),
  # tpm3 = matrix(0,3,2))
  
  tpm1 = c(1/3,1/3,1/3),
  tpm2 = c(1/3,1/3,1/3),
  tpm3 = c(1/3,1/3,1/3))
# tpm1 = matrix(c(0.0,0.0),1,2),
# tpm2 = matrix(c(0.0,0.0),1,2),
# tpm3 = matrix(c(0.0,0.0),1,2))
