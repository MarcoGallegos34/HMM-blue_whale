#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("-s", "--seed"), type="integer", default=1, 
              help="select seed [default= %default]", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd("../")
library(Rcpp)
library(gtools)
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")
# random.seeds = c(212341,
#                  212342,
#                  212343,
#                  212344,
#                  212345,
#                  212346,
#                  212347,
#                  212348,
#                  212349,
#                  2123410)
# seed = random.seeds[opt$seed]
seed = paste0("21234",opt$seed)
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
set.seed(seed)
print(paste0("seed selected: ",seed))
mod3_init_mat = matrix(NA,41,500)

# for(i in 1:500){
#   aux_inits = rdirichlet(4,alpha=rep(1,3))
#   aux_theta = rdirichlet(1,alpha=rep(1,4))
#   
#   mod3_init_mat[,i] = c(rnorm(6,0,10), # tpm k = 1
#                           qlogis(aux_inits[1,2:3]), # init k = 1
#                           qlogis(aux_theta[2]), # pi
#                           qlogis(aux_inits[2,2:3]), # init k = 2
#                           rnorm(6,0,10), # tpm k = 2
#                           qlogis(aux_theta[3]), # pi for k = 3
#                           qlogis(aux_inits[3,2:3]), # init k = 3
#                           rnorm(6,0,10), # tpm k = 3
#                           qlogis(aux_theta[4]), # pi for k = 4
#                           qlogis(aux_inits[4,2:3]), # init k = 4
#                           rnorm(6,0,10), # tpm k = 4
#                           rnorm(6,0,12)) # covariates beta
# }

for(i in 1:500){
  aux_inits = rdirichlet(4,alpha=rep(1,3))
  aux_theta = rdirichlet(1,alpha=rep(1,4))
  
  mod3_init_mat[,i] = c(c(rnorm(1,-5+.997,1),rnorm(1,-5+27.1,1),
                          #rnorm(1,-22,3),rnorm(1,-5,2),
                          rnorm(1,-1.25 -17,1),rnorm(1,-5-.147,1),
                          #rnorm(1,-3.5,3),rnorm(1,-10,3)), # tpm k = 1
                          rnorm(1,-4.5 + 1.24,1),rnorm(1,-5 - .390,1)), # tpm k = 1
                          qlogis(aux_inits[1,2:3]), # init k = 1
                          # qlogis(c(.001,.16)) + rnorm(2,0,1), # init k = 1
                         qlogis(aux_theta[2]), # pi
                        # qlogis(.449) + rnorm(1,0,1), # pi
                         qlogis(aux_inits[2,2:3]), # init k = 2
                        # qlogis(c(.37,.001)) + rnorm(2,0,1), # init k = 2
                        # c(rnorm(1,-4.5,3),rnorm(1,-5,3),rnorm(1,-27,5),rnorm(1,-4,3),rnorm(1,-2.8,2),rnorm(1,-1.6,2)), # tpm k = 2
                        c(rnorm(1,-8+.997,1),rnorm(1,-31.5+27.1,1),
                          rnorm(1,-5-17,1),rnorm(1,-1.35-.147,1),
                          rnorm(1,-4+1.24,1),rnorm(1,-1.2-.39,1)), # tpm k = 2
                        qlogis(aux_theta[3]), # pi for k = 3
                        # qlogis(.054) + rnorm(1,0,1), # pi for k = 3
                        qlogis(aux_inits[3,2:3]), # init k = 3
                        # qlogis(c(.5,.001)) + rnorm(2,0,1), # init k = 3
                        # c(rnorm(1,-1.9,2),rnorm(1,-.45,2),rnorm(1,-10,2),rnorm(1,-2.1,1.5),rnorm(1,2,1.5),rnorm(1,-.6,1.5)), # tpm k = 3
                        c(rnorm(1,-2.9+.997,1),rnorm(1,-27.52 + 27.1,1),
                          rnorm(1,-5-17,1),rnorm(1,-1.95-.147,1),
                          rnorm(1,-3.88+1.24,1),rnorm(1,-2-.39,1)), # tpm k = 3
                        qlogis(aux_theta[4]), # pi for k = 4
                        # qlogis(.163) + rnorm(1,0,1), # pi for k = 4
                        qlogis(aux_inits[4,2:3]), # init k = 4
                        # qlogis(c(.001,.34)) + rnorm(1,0,1), # init k = 4
                        # c(rnorm(1,-.3,1.5),rnorm(1,1,1.5),rnorm(1,-10,2),rnorm(1,.2,1.5),rnorm(1,-1,1.5),rnorm(1,-1.5,1.5)), # tpm k = 4
                        c(rnorm(1,-1.55+.997,1),rnorm(1,-26.55+27.1,1),
                          rnorm(1,-4-17,1),rnorm(1,.42-.147,1),
                          rnorm(1,-1.74+1.24,1),rnorm(1,1.37-.39,1)), # tpm k = 4
                        c(-.997,-27.1,17,.147,-1.24,.39) + rnorm(6,0,.1)) # covariates beta
}

sourceCpp("mod3/log-lik-semi-roptim.cpp")


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
                                  mod3_init_mat[,i])
  print(paste0("mle for sample  ",i," has finished"))
  mle_vec = c(mle_vec,semiMles)
  aux_min = min(mle_vec)
  aux_min_id = which(aux_min == mle_vec)
  print(paste0("smallest mle is ",aux_min," and correspond to column ",aux_min_id))
}
t2 = Sys.time()
t2-t1

saveRDS(mle_vec,paste0("mod3/semiMles_mod3-",opt$seed,".RDS"))
