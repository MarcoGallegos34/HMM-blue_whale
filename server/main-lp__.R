#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("-i", "--init"), type="integer", default=1, 
              help="Initial value set [default= %default]", metavar="integer"),
  make_option(c("-c", "--constant"), type="numeric", default=2, 
              help="geometric constant [default= %default]", metavar="numeric"),
  make_option(c("-t", "--ntemp"), type="integer", default=10, 
              help="Number of temperatures [default= %default]", metavar="integer"),
  make_option(c("-w", "--swap"), type="character", default="y", 
              help="Incorporation of Swap proposal [default= %default]", metavar="character"),
  make_option(c("-l", "--likelihood"), type="character", default="n", 
              help="Avoid incorporating prior information [default= %default]", metavar="character"),
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
if(opt$likelihood == "n"){
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
      if(opt$hyperpars == 1){
        source("stan-input/stan-data_init-values.R")
        sourceCpp("mcmc/lp__v1_1-noSwap.cpp")
      }
      if(opt$hyperpars == 2)
        source("stan-input/stan-data_init-values-3.R") # have in mind this may not match to hyperparameters given in C++
        sourceCpp("mcmc/lp__v1_2-noSwap.cpp")
    }
    
  }
}
if(opt$likelihood == "y"){
  source("stan-input/stan-data_init-values.R")
  sourceCpp("mcmc/lp__v1_3.cpp")
  
}
#init = opt$init
temp = opt$ntemp
nsimm = opt$nsim
hyperpars = opt$hyperpars
c = opt$constant
if(opt$init == 1){
  ### INIT 1 ###
  init = c(0.05,0.05,0.1,0.05,0.05,0.05,
            410,420,430,140,280,400,
            140,100,65,65,50,70,
            105,110,115,48,85,130,
            320,450,460,225,300,450,
            0.5,3,1,
            2.2,0.6,0.9,1.8,5,2.1,
            3.5,2.8,0.7,
            0.25,0.15,
            0.9,0.85,0.9,
            0.6,
            0.05,0.8,0.05)
  
}
if(opt$init == 2){
  ### INIT 2 ###
  init = c(0.01,0.25,0.05,0.05,0.15,0.05,
           190,190,520,140,115,125,
           65,72,155,55,65,65,
           28,55,175,20,55,60,
           320,325,450,250,280,325,
           1.8,1,1,
           0.5,1.6,1.3,3.5,3.2,1.4,
           4,1,3.5,
           0.17,0.65,
           0.74,0.9,0.8,
           0.18,
           0.95,0.05,0.12)
    
}
if(opt$init == 3){
  
  ### INIT 3 ###
  init = c(0.05,0.05,0.05,0.2,0.05,0.1,
           90,320,520,80,225,125,
           72,80,155,70,55,65,
           28,70,175,25,65,60,
           180,500,520,125,300,380,
           1,2.8,0.8,
           1.1,0.5,1.6,2.2,5,1.5,
           0.7,2,3.5,
           0.17,0.65,
           0.9,0.75,0.85,
           0.18,
           0.05,0.9,0.05)

}
if(opt$init == 4){
  
  ### INIT 4 ###
  init = c(0.1,0.05,0.1,0.05,0.1,0.05,
           90,480,490,80,140,275,
           72,145,105,65,65,60,
           28,150,155,25,60,140,
           200,350,750,140,250,280,
           1.2,0.8,3.2,
           0.8,2.8,0.5,2,2.5,3.5,
           0.5,3.5,2.5,
           0.55,0.25,
           0.85,0.85,0.85,
           0.2,
           0.05,0.05,0.75)
  
}


gc()
t1 = Sys.time()
t1
set.seed(194)
#armadillo_set_seed(10)
print(paste0("Number of temperatures requested: ", temp))
print(paste0("Number of simulations ",nsimm))
print(paste0("Initial values set = ",opt$init))
print(paste0("Swap between chains = ",opt$swap))
print(paste0("Diagnostics = ",opt$diagnostics))
print(paste0("Hyperparameter set = ",hyperpars))
print(paste0("Geometric constant = ",c))
print(paste0("Only likelihood (no priors) = ",opt$likelihood))
aux_geom_temp_vector = as.numeric(c^(0:temp))
if(opt$swap == "y"){
  geom_temp_vector = c(aux_geom_temp_vector[aux_geom_temp_vector < 100],100)
}
if(opt$swap == "n"){
  geom_temp_vector = c(aux_geom_temp_vector[aux_geom_temp_vector < 100])
}
print(paste0("Actual number of temperaturers: ",length(geom_temp_vector)))
print(geom_temp_vector)
sim_parallel_temp = rcpp_parallel_pt_cw_M_target_posterior(stan_data,nsim = nsimm,init = init,
                                                           temp_vector = geom_temp_vector,
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
    saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp-init",opt$init,"-hyperpars",hyperpars,"-diagnostics.rds"))
}
if(opt$diagnostics == "n"){
  if(opt$swap == "y"){
    if(opt$likelihood == "n"){
      saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp-init",opt$init,"-hyperpars",hyperpars,".rds"))
    }
    if(opt$likelihood == "y"){
      saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp-init",opt$init,"-noPrior",".rds"))
    }
  }
  
  if(opt$swap == "n"){
    if(opt$likelihood == "n"){
      saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp-init",opt$init,"-hyperpars",hyperpars,"-noSwap.rds"))
    }
    if(opt$likelihood == "y"){
      saveRDS(sim_parallel_temp,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp-init",opt$init,"-noSwap.rds"))
    }
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
    if(opt$likelihood == "n"){
      saveRDS(sim_pt_df,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp_tidy-init",opt$init,"-hyperpars",hyperpars,".rds"))
    }
    if(opt$likelihood == "y"){
      saveRDS(sim_pt_df,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp_tidy-init",opt$init,"-noPrior",".rds"))
    }
    
  }
  
  if(opt$swap == "n"){
    saveRDS(sim_pt_df,paste0("mcmc/output_rcpp_parallel_pt_m_","c",c,"-",length(geom_temp_vector),"temp_tidy-init",opt$init,"-hyperpars",hyperpars,"-noSwap.rds"))
    
  }
}

