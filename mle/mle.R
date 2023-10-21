source("dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")

# MLE for model 1

# MLE using nlm (wide data format)
t1 =Sys.time()
result = nlm(f=log.likelihood,p=theta_star,obs = list_whales,hessian = T)
t2 =Sys.time()
t2-t1


# MLE using nlm (long data format)
t1 =Sys.time()
result = nlm(f=log.likelihood_v2,p=theta_star,list_data = stan_data,hessian =T)
t2 =Sys.time()
t2-t1

result


# MLE using stan ####
# This is the good model 
mod1_mle = cmdstan_model("mle/mod1_mle.stan")
mod1_mle_result = mod1_mle$optimize(data = stan_data,
                                    init =  list(init_list),
                                    seed = 1034,
                                    #algorithm = "newton",
                                    #init_alpha = .0005,
                                    iter = 200)
