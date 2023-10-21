library(Rcpp)
library(parallel)
sourceCpp("mod2/log-lik-semi.cpp")
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")

# Modification of stan_data list
mod1.nlm = readRDS("mle/nlm.output.rds")
stan_data$K = 2
stan_data$estimate = mod1.nlm$estimate 


mod2_init_list = readRDS("mod2/mod2_init_list.rds")

mod2_init_list = mod2_init_list[1]

mod2_init_list[[1]]
theta_star_2_semi = c(0.5566620,4.2603743,-4.7024912,-0.1147913,-0.9297941,3.0402428,0.7318008,-1.2604133,-0.7158760,
  0.5930074,-0.7691263,1.9838027,-0.9877951,4.1602675,-2.3714208,1.0163141,-3.2914856)

length(theta_star_2_semi)

length(theta_star_semi)

t1 = Sys.time()
res_parallel = mclapply(mod2_init_list, function(theta_star){
  aux = nlm(log_likelihood_v2_mod2_semi,p = theta_star,list_data = stan_data)
  return(aux$minimum)
})
t2 =Sys.time()
print(t2 -t1)

res_parallel
log_likelihood_v2_mod2_semi(stan_data,mod2_init_list[[100]])

saveRDS(res_parallel,"server/mod2_mle_results.rds")

