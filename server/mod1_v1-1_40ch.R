setwd("../")
library(cmdstanr)

source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")

mod1_v1.1 = cmdstan_model("mod1_v1-1.stan")

sim = mod1_v1.1$sample(data = stan_data,
                       iter_warmup = 500,
                       iter_sampling = 500,
                       #parallel_chains = 40,
                       #init = list(init_list), 
                       #adapt_delta =.9, # with .95, we get 5% divergence
                       seed = 1805314595,
                       #seed = 2026228722, 
                       chains = 1
)


sim$save_object("test123.RDS")

mod1.40ch = saveRDS(sim,"mod1_v1_output_stan_40chain.rds")

### 500 of 20000 hit the maximum depth limit I think, this should be in the diagnostics

sim1.40ch = readRDS("server/mod1_v1_output_stan_40chain.rds")

print(sim1.40ch)


temps_rds_file <- tempfile(fileext = ".RDS")
sim$save_object(file = temps_rds_file)
rm(sim)

sim = readRDS(temps_rds_file)

### This is important
sim$save_output_files()
sim$save_data_file()
