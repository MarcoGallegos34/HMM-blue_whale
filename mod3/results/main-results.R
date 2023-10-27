### Generation of MLE paramters and Hessian matrix related to it for model 3 ###

library(Rcpp)
library(gtools)
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")
seed = "2123422"
stan_data$K = 4
mod1.nlm = readRDS("mle/nlm.output.rds")
# stan_data$estimate = mod1.nlm$estimate 
stan_data$estimate = c(mod1.nlm$estimate[1:6],
                       139,342,515,log(77.8),log(217),log(129),
                       71.3,84.8,150,log(68.1),log(56.5),log(69.2),
                       32.2,67.7,169,log(23.3),log(64.4),log(60.9),
                       193,693,409,log(138),log(304),log(287),
                       log(1),log(3.1),log(.8),
                       log(.94),log(.518),log(1.66),
                       log(2.05),log(6.55),log(1.58),
                       log(.656),log(.01),log(3.31))

# Generation of initial values
set.seed(seed)
print(paste0("seed selected: ",seed))
mod3_init_mat = matrix(NA,41,500)

for(i in 1:500){
  aux_inits = rdirichlet(4,alpha=rep(1,3))
  aux_theta = rdirichlet(1,alpha=rep(1,4))
  
  mod3_init_mat[,i] = c(c(rnorm(1,-5+.997,1),rnorm(1,-5+27.1,1),
                          rnorm(1,-1.25 -17,1),rnorm(1,-5-.147,1),
                          rnorm(1,-4.5 + 1.24,1),rnorm(1,-5 - .390,1)), # tpm k = 1
                        #qlogis(aux_inits[1,2:3]), # init k = 1
                        c(-6,-1.65) + rnorm(2,0,.1), # init k = 1
                        .2642 + rnorm(1,0,.1), # pi
                        # qlogis(.449) + rnorm(1,0,1), # pi
                        c(-.525,-5), # init k = 2
                        # qlogis(c(.37,.001)) + rnorm(2,0,1), # init k = 2
                        c(rnorm(1,-8+.997,1),rnorm(1,-31.5+27.1,1),
                          rnorm(1,-5-17,1),rnorm(1,-1.35-.147,1),
                          rnorm(1,-4+1.24,1),rnorm(1,-1.2-.39,1)), # tpm k = 2
                        -1.9 + rnorm(1,0,.1), # pi for k = 3
                        # qlogis(.054) + rnorm(1,0,1), # pi for k = 3
                        c(-.525,-5) + rnorm(2,0,.1), # init k = 3
                        # qlogis(c(.5,.001)) + rnorm(2,0,1), # init k = 3
                        c(rnorm(1,-2.9+.997,1),rnorm(1,-27.52 + 27.1,1),
                          rnorm(1,-5-17,1),rnorm(1,-1.95-.147,1),
                          rnorm(1,-3.88+1.24,1),rnorm(1,-2-.39,1)), # tpm k = 3
                        -.8 + rnorm(1,0,.1), # pi for k = 4
                        # qlogis(.163) + rnorm(1,0,1), # pi for k = 4
                        c(-5,-.67) + rnorm(2,0,.1), # init k = 4
                        # qlogis(c(.001,.34)) + rnorm(1,0,1), # init k = 4
                        c(rnorm(1,-1.55+.997,1),rnorm(1,-26.55+27.1,1),
                          rnorm(1,-4-17,1),rnorm(1,.42-.147,1),
                          rnorm(1,-1.74+1.24,1),rnorm(1,1.37-.39,1)), # tpm k = 4
                        c(-.997,-27.1,17,.147,-1.24,.39)) # covariates beta
}

sourceCpp("mod3/log-lik-semi-roptim.cpp")


t1 = Sys.time()
semiMles_par = semiMleExtended_bfgs_par(stan_data$K,
                                stan_data$N,
                                stan_data$n,
                                stan_data$n_ind,
                                stan_data$ID_init,
                                stan_data$ID,
                                stan_data$x_duration_init,
                                stan_data$x_surface_init,
                                stan_data$x_maxDepth_init,
                                stan_data$x_lunges_init,
                                stan_data$x_step_init,
                                stan_data$x_angle_init,
                                stan_data$x_headVar_init,
                                stan_data$x_exposure_init,
                                stan_data$x_duration,
                                stan_data$x_surface,
                                stan_data$x_maxDepth,
                                stan_data$x_lunges,
                                stan_data$x_step,
                                stan_data$x_angle,
                                stan_data$x_headVar,
                                stan_data$x_exposure,
                                stan_data$estimate,
                                theta_test)
t2 = Sys.time()
t2-t1
# Semi Mle = 25638.102874
# It takes 9.439 mins
#saveRDS(semiMles_par,"mod3/results/semi_mle_pars.RDS")

t1 = Sys.time()
semiMles_hessian = semiMleExtended_bfgs_hessian(stan_data$K,
                                stan_data$N,
                                stan_data$n,
                                stan_data$n_ind,
                                stan_data$ID_init,
                                stan_data$ID,
                                stan_data$x_duration_init,
                                stan_data$x_surface_init,
                                stan_data$x_maxDepth_init,
                                stan_data$x_lunges_init,
                                stan_data$x_step_init,
                                stan_data$x_angle_init,
                                stan_data$x_headVar_init,
                                stan_data$x_exposure_init,
                                stan_data$x_duration,
                                stan_data$x_surface,
                                stan_data$x_maxDepth,
                                stan_data$x_lunges,
                                stan_data$x_step,
                                stan_data$x_angle,
                                stan_data$x_headVar,
                                stan_data$x_exposure,
                                stan_data$estimate,
                                mod3_init_mat[,269])
t2 = Sys.time()
t2-t1

#saveRDS(semiMles_hessian,"mod3/results/semi_mle_hessian.RDS")
semiMles_hessian = readRDS("mod3/results/semi_mle_hessian.RDS")


semiMles_par = readRDS("mod3/results/semi_mle_pars.RDS")

### tpms

# k = 1
alphas_to_tpm(semiMles_par[1:6]+semiMles_par[36:41]) 

# k = 2
alphas_to_tpm(semiMles_par[12:17]+semiMles_par[36:41]) 

# k = 3
alphas_to_tpm(semiMles_par[21:26]+semiMles_par[36:41]) 

# k = 4
alphas_to_tpm(semiMles_par[30:35]+semiMles_par[36:41]) 


sourceCpp("mod3/log-lik-mod3-roptim.cpp")


semi_theta_star = c(semiMles_par[1:6],
                    stan_data$estimate[7:42],
                    semiMles_par[7:41])

final_theta = mle_bfgs_par(stan_data$K,
                           stan_data$N,
                           stan_data$n,
                           stan_data$n_ind,
                           stan_data$ID_init,
                           stan_data$ID,
                           stan_data$x_duration_init,
                           stan_data$x_surface_init,
                           stan_data$x_maxDepth_init,
                           stan_data$x_lunges_init,
                           stan_data$x_step_init,
                           stan_data$x_angle_init,
                           stan_data$x_headVar_init,
                           stan_data$x_exposure_init,
                           stan_data$x_duration,
                           stan_data$x_surface,
                           stan_data$x_maxDepth,
                           stan_data$x_lunges,
                           stan_data$x_step,
                           stan_data$x_angle,
                           stan_data$x_headVar,
                           stan_data$x_exposure,
                           semi_theta_star)

# Final MLE = 25636.877829
#saveRDS(final_theta,"mod3/results/mle_pars.RDS")
final_theta = readRDS("mod3/results/mle_pars.RDS")
final_theta2 = readRDS("mod3/results/mle_pars2.RDS")

final_theta[6:42] - stan_data$estimate[6:42]
exp(final_theta[10:12])

### tpms
aux_final_theta = final_theta[c(1:6,43:77)]

### Baseline ###

# k = 1
alphas_to_tpm(aux_final_theta[1:6]) 

# k = 2
alphas_to_tpm(aux_final_theta[12:17]) 

# k = 3
alphas_to_tpm(aux_final_theta[21:26]) 

# k = 4
alphas_to_tpm(aux_final_theta[30:35]) 

### Acoustic exposure ###

# k = 1
alphas_to_tpm(aux_final_theta[1:6]+aux_final_theta[36:41]) 

# k = 2
alphas_to_tpm(aux_final_theta[12:17]+aux_final_theta[36:41]) 

# k = 3
alphas_to_tpm(aux_final_theta[21:26]+aux_final_theta[36:41]) 

# k = 4
alphas_to_tpm(aux_final_theta[30:35]+aux_final_theta[36:41]) 

# pi (weights for each of the contexts)

pi_vector = function(alpha){
  result = c()
   result[1]= 1/(1+sum(exp(alpha)))
   result[2]= exp(alpha[1])/(1+sum(exp(alpha)))
   result[3]= exp(alpha[2])/(1+sum(exp(alpha)))
   result[4]= exp(alpha[3])/(1+sum(exp(alpha)))
   
   return(result)
}

pi_vector(aux_final_theta[c(9,18,27)])

init_vector = function(alpha){
   result = c()
   result[1]= 1/(1+sum(exp(alpha)))
   result[2]= exp(alpha[1])/(1+sum(exp(alpha)))
   result[3]= exp(alpha[2])/(1+sum(exp(alpha)))
   
   return(result)

}

# inits

# k = 1
init_vector(aux_final_theta[7:8])
# k = 2
init_vector(aux_final_theta[10:11])
# k = 3
init_vector(aux_final_theta[19:20])
# k = 4
init_vector(aux_final_theta[28:29])

sourceCpp("mod3/log-lik-semi.cpp")
stan_data$estimate = final_theta[1:42]
stan_data$K = 4
log_likelihood_mod3_semi(stan_data,aux_final_theta)

#rnorm(1,-4.5 + 1.24,1),rnorm(1,-5 - .390,1)), # tpm k = 1
#c(-6,-1.65) + rnorm(2,0,.1), # init k = 1
#.2642 + rnorm(1,0,.1), # pi = 2
#c(-.525,-5), # init k = 2
#rnorm(1,-4+1.24,1),rnorm(1,-1.2-.39,1)), # tpm k = 2
#-1.9 + rnorm(1,0,.1), # pi for k = 3
#c(-.525,-5) + rnorm(2,0,.1), # init k = 3
#rnorm(1,-3.88+1.24,1),rnorm(1,-2-.39,1)), # tpm k = 3
#-.8 + rnorm(1,0,.1), # pi for k = 4
#c(-5,-.67) + rnorm(2,0,.1), # init k = 4
#rnorm(1,-1.74+1.24,1),rnorm(1,1.37-.39,1)), # tpm k = 4
#c(-.997,-27.1,17,.147,-1.24,.39)) # covariates beta


aux_final_theta_permuted = c(aux_final_theta[1:6], # tpm k = 1
                             aux_final_theta[7:8], # init k = 1
                             aux_final_theta[27], # pi k = 2
                             aux_final_theta[28:29], # init k = 2
                             aux_final_theta[30:35], # tpm k = 2
                             aux_final_theta[9], # pi k = 3
                             aux_final_theta[10:11], # init k =3
                             aux_final_theta[12:17], # tpm k = 3
                             aux_final_theta[18], # pi k = 4
                             aux_final_theta[19:20], # init k = 4
                             aux_final_theta[21:26], # tpm k = 4 
                             aux_final_theta[36:41])

aux_final_theta_permuted = c(aux_final_theta[1:6], # tpm k = 1
                             aux_final_theta[7:8], # init k = 1
                             0, # pi k = 2
                             0,0, # init k = 2
                             rep(0,6), # tpm k = 2
                             0, # pi k = 3
                             0,0, # init k =3
                             rep(0,6), # tpm k = 3
                             0, # pi k = 4
                             0,0, # init k = 4
                             rep(0,6), # tpm k = 4 
                             aux_final_theta[36:41])

aux_final_theta_permuted = c(rep(0,6), # tpm k = 1
                             0,0, # init k = 1
                             0, # pi k = 2
                             0,0, # init k = 2
                             rep(0,6), # tpm k = 2
                             0, # pi k = 3
                             0,0, # init k =3
                             rep(0,6), # tpm k = 3
                             0, # pi k = 4
                             aux_final_theta[7:8], # init k = 4
                             aux_final_theta[1:6], # tpm k = 4 
                             aux_final_theta[36:41])

log_likelihood_mod3_semi(stan_data,aux_final_theta_permuted)

# k = 1
alphas_to_tpm(aux_final_theta[1:6]) 

# k = 2
alphas_to_tpm(aux_final_theta[12:17]) 

# k = 3
alphas_to_tpm(aux_final_theta[21:26]) 

# k = 4
alphas_to_tpm(aux_final_theta[30:35]) 

# pi (weights for each of the contexts)

pi_vector(aux_final_theta[c(9,18,27)])

# inits

# k = 1
init_vector(aux_final_theta[7:8])
# k = 2
init_vector(aux_final_theta[10:11])
# k = 3
init_vector(aux_final_theta[19:20])
# k = 4
init_vector(aux_final_theta[28:29])



