library(Rcpp)
library(gtools)
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")

# mle_files =  paste0("mod3/semiMles_mod3-",c(1:30),".RDS")
# 
# result = matrix(NA,500,30)
# 
# for(i in 1:30){
#   result[,i] = readRDS(mle_files[i])
# }
# 
# for(i in 1:30){
#   print(min(result[,i]))
# }
# 
# id_min = c()
# for(i in 1:10){
#   id_min = c(id_min,which(result[,i] == min(result[,i])))
#   
# }
# 
# saveRDS(result,"mod3/results/semi_mles.RDS")
# 
# min(result)
# 
# which( result <= 25645 & result >= 25643)

#sourceCpp("mod3/log-lik-semi.cpp")
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


### K = 4 ###

# Generation of initial values
#set.seed(212341)

mod3_init_array = array(NA,dim=c(41,500,30))

for(j in 1:30){
  
  seed = paste0("21234",j)
  set.seed(seed)
  print(paste0("seed selected = ",seed))
  # mod3_init_mat = matrix(NA,41,500)
  
  for(i in 1:500){
    aux_inits = rdirichlet(4,alpha=rep(1,3))
    aux_theta = rdirichlet(1,alpha=rep(1,4))
    
    mod3_init_array[,i,j] = c(c(rnorm(1,-5+.997,1),rnorm(1,-5+27.1,1),
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
}

sourceCpp("mod3/log-lik-semi-roptim.cpp")


semiMles_par1 = semiMleExtended_bfgs_par(stan_data$K,
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
                                mod3_init_array[,332,5])

semiMles_par2 = semiMleExtended_bfgs_par(stan_data$K,
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
                                mod3_init_array[,339,11])

semiMles_par3 = semiMleExtended_bfgs_par(stan_data$K,
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
                                mod3_init_array[,21,12])


semiMles_par1- semiMles_par2
semiMles_par1- semiMles_par3
semiMles_par2- semiMles_par3

result[332,5]
result[339,11]
result[21,12]
result[23,17]


which( result <= 25643.1)

result[8023]
result[10769]

hist(result,breaks = 80)

semiMles_par1


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

