#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("-s", "--seed"), type="integer", default=1, 
              help="select seed [default= %default]", metavar="integer"),
  make_option(c("-k", "--kcontext"), type="integer", default=3, 
              help="select number of contexts [default= %default]", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("../")
library(Rcpp)
library(gtools)
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")
random.seeds = c(212341,
                 212342,
                 212343,
                 212344,
                 212345,
                 212346,
                 212347,
                 212348,
                 212349,
                 2123410)
seed = random.seeds[opt$seed]

# Modification of stan_data list
stan_data$K = opt$kcontext
mod1.nlm = readRDS("mle/nlm.output.rds")
stan_data$estimate = mod1.nlm$estimate 


# Generation of initial values
set.seed(seed)
print(paste0("number of contexts selected: ",stan_data$K))
print(paste0("seed selected: ",seed))

### K = 3 ###
if(opt$kcontext == 3){
  mod2_init_mat = matrix(NA,26,1500)
  
  for(i in 1:1500){
    aux_inits = rdirichlet(3,alpha=rep(1,3))
    aux_theta = rdirichlet(1,alpha=rep(1,3))
    
    mod2_init_mat[,i] = c(rnorm(6,-3,5), # tpm k = 1
                            qlogis(aux_inits[1,2:3]), # init k = 1
                            qlogis(aux_theta[2]), # pi
                            qlogis(aux_inits[2,2:3]), # init k = 2
                            rnorm(6,-3,5), # tpm k = 2
                            qlogis(aux_theta[3]), # pi for k = 3
                            qlogis(aux_inits[3,2:3]), # init k = 3
                            rnorm(6,-3,5)) # tpm k = 3
  }
  
}

if(opt$kcontext == 4){
  
  mod2_init_mat = matrix(NA,35,1500)
  
  for(i in 1:1500){
    aux_inits = rdirichlet(4,alpha=rep(1,3))
    aux_theta = rdirichlet(1,alpha=rep(1,4))
    
    mod2_init_mat[,i] = c(rnorm(6,-3,5), # tpm k = 1
                            qlogis(aux_inits[1,2:3]), # init k = 1
                            qlogis(aux_theta[2]), # pi
                            qlogis(aux_inits[2,2:3]), # init k = 2
                            rnorm(6,-3,5), # tpm k = 2
                            qlogis(aux_theta[3]), # pi for k = 3
                            qlogis(aux_inits[3,2:3]), # init k = 3
                            rnorm(6,-3,5), # tpm k = 3
                            qlogis(aux_theta[4]), # pi for k = 4
                            qlogis(aux_inits[4,2:3]), # init k = 4
                            rnorm(6,-3,5)) # tpm k = 3
  }
  
}

if(opt$kcontext == 5){
  
  mod2_init_mat = matrix(NA,44,1500)
  
  for(i in 1:1500){
    aux_inits = rdirichlet(5,alpha=rep(1,3))
    aux_theta = rdirichlet(1,alpha=rep(1,5))
    
    mod2_init_mat[,i] = c(rnorm(6,-3,5), # tpm k = 1
                            qlogis(aux_inits[1,2:3]), # init k = 1
                            qlogis(aux_theta[2]), # pi
                            qlogis(aux_inits[2,2:3]), # init k = 2
                            rnorm(6,-3,5), # tpm k = 2
                            qlogis(aux_theta[3]), # pi for k = 3
                            qlogis(aux_inits[3,2:3]), # init k = 3
                            rnorm(6,-3,5), # tpm k = 3
                            qlogis(aux_theta[4]), # pi for k = 4
                            qlogis(aux_inits[4,2:3]), # init k = 4
                            rnorm(6,-3,5), # tpm k = 4
                            qlogis(aux_theta[5]), # pi for k = 5
                            qlogis(aux_inits[5,2:3]), # init k = 5
                            rnorm(6,-3,5)) # tpm k = 5
  }
  
}

sourceCpp("mod2/log-lik-semi-extended-roptim.cpp")

t1 = Sys.time()
mle_vec = c()
for(i in 1:500){
  
  semiMles = semiMleExtended_bfgs(stan_data$K,
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
                                  mod2_init_mat[,i])
  print(paste0("mle for sample  ",i," has finished"))
  mle_vec = c(mle_vec,semiMles)
  aux_min = min(mle_vec)
  aux_min_id = which(aux_min == mle_vec)
  print(paste0("smallest mle is ",aux_min," and correspond to column ",aux_min_id))
}
t2 = Sys.time()
t2-t1

saveRDS(mle_vec,paste0("mod2/semiMles_",opt$kcontext,"-",opt$seed,".RDS"))
