setwd("../")
library(Rcpp)
library(gtools)
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")

# Modification of stan_data list
stan_data$K = 2
mod1.nlm = readRDS("mle/nlm.output.rds")
stan_data$estimate = mod1.nlm$estimate 

# theta_star_2_semi = c(rep(0,6), # tpm k = 1
#                       # muSigma_duration[,1],log(muSigma_duration[,2]), # duration
#                       # muSigma_surface[,1],log(muSigma_surface[,2]), # surface
#                       # muSigma_maxDepth[,1],log(muSigma_maxDepth[,2]), # maxDepth
#                       # log(muSigma_step[,1]),log(muSigma_step[,2]), # step
#                       # log(muKappa_angle[,2]), # angle
#                       # log(ab_headVar[,1]),log(ab_headVar[,2]), #varHead
#                       # log(lambda_lunges), #lunges
#                       qlogis(c(1/3,1/3)), #init distribution k = 1
#                       qlogis(c(1/2)), # pi_
#                       qlogis(c(1/3,1/3)), #init distribution k = 2
#                       rep(0,6)) # tpm k = 2


sourceCpp("mod2/log-lik-semi-roptim.cpp")

# Generation of initial values
set.seed(23984)
library(gtools)
mod2_init_mat = matrix(NA,17,15000)

for(i in 1:15000){
  aux_inits = rdirichlet(2,alpha=rep(1,3))
  
  mod2_init_mat[,i] = c(rnorm(6,0,5), # tpm k = 1
                        qlogis(aux_inits[1,2:3]), # init k = 1
                        qlogis(runif(1)), # pi
                        qlogis(aux_inits[2,2:3]), # init k = 2
                        rnorm(6,0,5)) # tpm k = 2
}


# t1 = Sys.time()
# nestim = 15
# mle_vec = c()
# ### Setting range from 1:935 would compute 14,960 different mles
# for(i in 1:10){
#   Y = mod2_init_mat[,(nestim - 14):nestim]
#   t1
#   semiMles = parallel_mle_estimator(Y,
#                                 stan_data$K,
#                                 stan_data$N,
#                                 stan_data$n,
#                                 stan_data$n_ind,
#                                 stan_data$ID_init,
#                                 stan_data$ID,
#                                 stan_data$x_duration_init,
#                                 stan_data$x_surface_init,
#                                 stan_data$x_maxDepth_init,
#                                 stan_data$x_lunges_init,
#                                 stan_data$x_step_init,
#                                 stan_data$x_angle_init,
#                                 stan_data$x_headVar_init,
#                                 stan_data$x_duration,
#                                 stan_data$x_surface,
#                                 stan_data$x_maxDepth,
#                                 stan_data$x_lunges,
#                                 stan_data$x_step,
#                                 stan_data$x_angle,
#                                 stan_data$x_headVar,
#                                 stan_data$estimate,0,15)
#   mle_vec = c(mle_vec,semiMles)
#   gc()
#   nestim = nestim + 15
#   print(paste0("iteration ",i," has finished"))
# }
# t2 = Sys.time()
# t2-t1
# 
# saveRDS(semiMles,"mod2/semi-mles.RDS")

### cols 48, 59 and 55 of the samples generated provide slower mle than the one indicated in the paper!
### We select sample 48
new_theta = semiMle_bfgs(stan_data$K,
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
                         stan_data$x_duration,
                         stan_data$x_surface,
                         stan_data$x_maxDepth,
                         stan_data$x_lunges,
                         stan_data$x_step,
                         stan_data$x_angle,
                         stan_data$x_headVar,
                         stan_data$estimate,
                         mod2_init_mat[,48])

### Completition of parameter to minimize the log-likelihood ###
complete_theta = c(new_theta[1:6],mod1.nlm$estimate[7:42],new_theta[7:17])

sourceCpp("mod2/log-lik-roptim.cpp")

final_theta = mle_bfgs(stan_data$K,
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
                       stan_data$x_duration,
                       stan_data$x_surface,
                       stan_data$x_maxDepth,
                       stan_data$x_lunges,
                       stan_data$x_step,
                       stan_data$x_angle,
                       stan_data$x_headVar,
                       complete_theta)

saveRDS(final_theta,"mod2/theta-log-mle.RDS")
