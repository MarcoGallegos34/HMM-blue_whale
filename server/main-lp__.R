#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("-i", "--init"), type="integer", default=1, 
              help="Initial value set [default= %default]", metavar="integer"),
  make_option(c("-t", "--ntemp"), type="integer", default=10, 
              help="Number of temperatures [default= %default]", metavar="integer"),
  make_option(c("-w", "--swap"), type="character", default="y", 
              help="Incorporation of Swap proposal [default= %default]", metavar="character"),
  make_option(c("-s", "--step"), type="integer", default=1, 
              help="Step size tunning [default= %default]", metavar="integer"),
  make_option(c("-p", "--hyperpars"), type="integer", default=1, 
              help="Hyperparameter setup [default= %default]", metavar="integer"),
  make_option(c("-n", "--nsim"), type="integer", default=10000, 
              help="Number of samples [default= %default]", metavar="integer"),
  make_option(c("-d", "--diagnostics"), type="character", default="n", 
              help="Version including diagnostics [default= %default]", metavar="character")
)
### Everything is set to incorporate step size tunning for sourceCpp if needed in the future 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("../")
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
source("data/dataPrep_functions.R")
if(opt$diagnostics == "y"){
  if(opt$hyperpars == 1){
    sourceCpp("mcmc/lp__v1_1-diagnostics.cpp")
    source("stan-input/stan-data_init-values.R")
  }
  if(opt$hyperpars == 2){
    sourceCpp("mcmc/lp__v1_2-diagnostics.cpp")
    source("stan-input/stan-data_init-values-3.R")
  }
}
if(opt$diagnostics == "n"){
  if(opt$swap == "y"){
    if(opt$hyperpars == 1){
      sourceCpp("mcmc/lp__v1_1.cpp")
      source("stan-input/stan-data_init-values.R")
    }
    if(opt$hyperpars == 2){
      sourceCpp("mcmc/lp__v1_2.cpp")
      source("stan-input/stan-data_init-values-3.R")
    }
  }
  if(opt$swap == "n"){
    source("stan-input/stan-data_init-values.R")
    sourceCpp("mcmc/lp__v1_1-noSwap.cpp")
  }
  
}
init = opt$init
temp = opt$ntemp
nsimm = opt$nsim
hyperpars = opt$hyperpars
if(init == 1){
  muSigma_duration = matrix(c(140,334,516,80,212,130),3,2)
  muSigma_surface = matrix(c(70,86,151,68,55,69),3,2)
  muSigma_maxDepth = matrix(c(32,68,170,24,65,60),3,2)
  muSigma_step = matrix(c(189,675,406,134,305,287),3,2)
  muKappa_angle = matrix(c(0,0,0,1,3.1,.8),3,2)
  ab_headVar = matrix(c(1,.5,1.7,2.1,5.4,1.6),3,2)
  lambda_lunges = c(.7,.05,3.4)

  theta_star_test = c(rep(1/3,6), # tpm k = 1
                     #c(120,300,516),muSigma_duration[,2],
                     muSigma_duration[,1],muSigma_duration[,2], # duration
                     muSigma_surface[,1],muSigma_surface[,2], # surface
                     muSigma_maxDepth[,1],muSigma_maxDepth[,2], # maxDepth
                     muSigma_step[,1],muSigma_step[,2], # step
                     muKappa_angle[,2], # angle
                     ab_headVar[,1],ab_headVar[,2], #varHead
                     lambda_lunges, #lunges
                     c(1/3,1/3), #init distribution k = 1
                     rep(1/3,4),
                     rep(1/2,3)) # missing entries tpm and init
  
  rm(muSigma_duration)
  rm(muSigma_surface)
  rm(muSigma_maxDepth)
  rm(muSigma_step)
  rm(muKappa_angle)
  rm(ab_headVar)
  rm(lambda_lunges)
  
}
if(init == 2){
 
    
   muSigma_duration = matrix(c(140,515,320,80,150,155),3,2)
   muSigma_surface = matrix(c(70,140,100,68,62,69),3,2)
   muSigma_maxDepth = matrix(c(32,160,100,24,60,75),3,2)
   muSigma_step = matrix(c(189,420,600,134,300,280),3,2)
   muKappa_angle = matrix(c(0,0,0,1,1.3,1.4),3,2)
   ab_headVar = matrix(c(.95,.9,1,2.1,1,8),3,2)
   lambda_lunges = c(.7,3.3,3.4)
	  
   theta_star_test = c(rep(1/3,6), # tpm k = 1
			#c(120,300,516),muSigma_duration[,2],
			muSigma_duration[,1],muSigma_duration[,2], # duration
			muSigma_surface[,1],muSigma_surface[,2], # surface
			muSigma_maxDepth[,1],muSigma_maxDepth[,2], # maxDepth
			muSigma_step[,1],muSigma_step[,2], # step
			muKappa_angle[,2], # angle
			ab_headVar[,1],ab_headVar[,2], #varHead
			lambda_lunges, #lunges
			c(1/3,1/3), #init distribution k = 1
			rep(1/3,4), # missing entries tpm and init
			# rep(1/2,3)) 
			c(.03,.005,.88)) # weights for the zero-inflated poisson distribution 

  rm(muSigma_duration)
  rm(muSigma_surface)
  rm(muSigma_maxDepth)
  rm(muSigma_step)
  rm(muKappa_angle)
  rm(ab_headVar)
  rm(lambda_lunges)
}

gc()
t1 = Sys.time()
t1
set.seed(194)
#armadillo_set_seed(10)
print(paste0("Number of temperatures: ", temp))
print(paste0("Number of simulations ",nsimm))
print(paste0("Initial values set = ",init))
print(paste0("Swap between chains = ",opt$swap))
print(paste0("Diagnostics = ",opt$diagnostics))
print(paste0("Hyperparameter set = ",hyperpars))
sim_parallel_temp = rcpp_parallel_pt_cw_M_target_posterior(stan_data,nsim = nsimm,init = theta_star_test,
                                                           temp_vector = as.numeric(2^(0:temp)),
                                                           #data for parallel computing
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
                                                           stan_data$x_headVar)
t2 = Sys.time()
(t2-t1) # For 100 it, 10 temps, this will abe aprox. 73 min, let's see if that happens
# 1.5 hrs for 150 it, 10 temps. 

if(opt$diagnostics == "y"){
    saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_",temp,"temp-init",init,"-hyperpars",hyperpars,"-diagnostics.rds"))
}
if(opt$diagnostics == "n"){
  if(opt$swap == "y"){
    saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_",temp,"temp-init",init,"-hyperpars",hyperpars,".rds"))
  }
  
  if(opt$swap == "n"){
    saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_",temp,"temp-init",init,"-noSwap.rds"))
  }
  
  sim_target = sim_parallel_temp[[1]]
  
  variable = c(rep("tpm1",2),
               rep("tpm2",2),
               rep("tpm3",2),
               rep("mu_duration",3),
               rep("sigma_duration",3),
               rep("mu_surface",3),
               rep("sigma_surface",3),
               rep("mu_maxDepth",3),
               rep("sigma_maxDepth",3),
               rep("mu_step",3),
               rep("sigma_step",3),
               rep("kappa",3),
               rep("a",3),
               rep("b",3),
               rep("lambda",3),
               rep("init",2),
               rep("tpm1",1),
               rep("tpm2",1),
               rep("tpm3",1),
               rep("init",1),
               rep("theta",3))
  
  state = c(2:3,
            1,
            3,
            1:2,
            rep(1:3,12),
            2:3,
            1:3,
            1,
            1:3)
  
  sim_pt_df = cbind(as.data.frame(sim_target),variable,state)
  
  sim_pt_df = sim_pt_df %>% pivot_longer(-c(variable,state),names_to="iteration",values_to="values") %>% 
    mutate(iteration = as.integer(gsub("^.","",iteration))) %>% 
    pivot_wider(values_from=values, names_from=variable) %>% arrange(iteration,state)
  
  if(opt$swap == "y"){
    saveRDS(sim_pt_df,paste0("mcmc/output_rcpp_parallel_pt_m_",temp,"temp_tidy-init",init,"-hyperpars",hyperpars,".rds"))
    
  }
  
  if(opt$swap == "n"){
    saveRDS(sim_pt_df,paste0("mcmc/output_rcpp_parallel_pt_m_",temp,"temp_tidy-init",init,"-noSwap.rds"))
    
  }
}

