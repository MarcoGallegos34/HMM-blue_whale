setwd("../")
library(Rcpp)
library(ggplot2)
library(tidyr)
library(dplyr)
#library(RcppParallel)
#library(RcppArmadillo)
source("data/dataPrep_functions.R")
#source("stan-input/stan-data_init-values.R")
source("stan-input/stan-data_init-values-3.R")
sourceCpp("mcmc/lp__v1_1.cpp")
sourceCpp("mcmc/lp__v1_2.cpp")
sourceCpp("mcmc/lp__v1_3.cpp")
#sourceCpp("mcmc/lp__v1_2-debug.cpp")
sourceCpp("mcmc/lp__v1_1-diagnostics.cpp")
sourceCpp("mcmc/lp__v1_2-diagnostics.cpp")
sourceCpp("mcmc/lp__v1_2-diagnosticsv2.cpp")
sourceCpp("mcmc/lp__v1_2-within.cpp")
sourceCpp("mcmc/lp__v1_1-conditionalMarginal.cpp")
#sourceCpp("mcmc/lp__v1_1_step2.cpp")
#sourceCpp("mcmc/lp__v1_1-debug.cpp")
#sourceCpp("mcmc/lp__v2.cpp")
if(TRUE){
  muSigma_duration = matrix(c(140,334,516,80,212,130),3,2)
  muSigma_surface = matrix(c(70,86,151,68,55,69),3,2)
  muSigma_maxDepth = matrix(c(32,68,170,24,65,60),3,2)
  muSigma_step = matrix(c(189,675,406,134,305,287),3,2)
  muKappa_angle = matrix(c(0,0,0,1,3.1,.8),3,2)
  ab_headVar = matrix(c(1,.5,1.7,2.1,5.4,1.6),3,2)
  lambda_lunges = c(.7,.05,3.4)
  
  # muSigma_duration = matrix(c(160,500,420,120,140,250),3,2)
  # muSigma_surface = matrix(c(70,140,110,68,35,69),3,2)
  # muSigma_maxDepth = matrix(c(27,150,120,17.5,20,95),3,2)
  # muSigma_step = matrix(c(240,320,620,180,175,230),3,2)
  # muKappa_angle = matrix(c(0,0,0,1.3,.8,3),3,2)
  # ab_headVar = matrix(c(.5,2.1,.4,1.8,1.5,4),3,2)
  # lambda_lunges = c(.45,3.5,2.3)
  
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
                       rep(1/2,3)) #weight for zero-inflated poisson distribution
                      #c(.1,.005,.75)) 
  
  rm(muSigma_duration)
  rm(muSigma_surface)
  rm(muSigma_maxDepth)
  rm(muSigma_step)
  rm(muKappa_angle)
  rm(ab_headVar)
  rm(lambda_lunges)
  
  muSigma_duration = matrix(c(140,515,320,80,150,155),3,2)
  muSigma_surface = matrix(c(70,140,100,68,62,69),3,2)
  muSigma_maxDepth = matrix(c(32,160,100,24,60,75),3,2)
  muSigma_step = matrix(c(189,420,600,134,300,280),3,2)
  muKappa_angle = matrix(c(0,0,0,1,1.3,1.4),3,2)
  ab_headVar = matrix(c(.95,.9,1,2.1,1,8),3,2)
  lambda_lunges = c(.7,3.3,3.4)
  
  theta_star_test_tempered = c(rep(1/3,6), # tpm k = 1
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
### Quesiton: would it be possible there's some state-switching???
### result from v1_1_step2: in many parameteres, these doesn't stay in the corresponding density areas - would
### this indicate there's actually state-switching??

cbind(theta_star_test_tempered,theta_star_test_tempered)

gc()
t1 = Sys.time()
t1
set.seed(194)
#armadillo_set_seed(10)
sim_parallel_temp = rcpp_parallel_pt_cw_M_target_posterior(stan_data,nsim=1000, init= init, #init = theta_star_test, 
                                                           #init_tempered = matrix(rep(theta_star_test_tempered,10),51,10),
                                                           # temp_vector = as.numeric(c(1.8^(0:7),100)),
                                                           temp_vector = geom_temp_vector,
                                                           #temp_vector = as.numeric(c(1:10)),
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
                                                           stan_data$x_headVar,
                                                           #within_temp = 20
                                                           )
t2 = Sys.time()
(t2-t1) # For 100 it, 10 temps, this will abe aprox. 73 min, let's see if that happens
# 1.5 hrs for 150 it, 10 temps. 
### inicie jobs a las 11:20 pm
swap = sim_parallel_temp$swap
swap_proposal = sim_parallel_temp$swap_proposal
path_init = sim_parallel_temp$path_init
swap_proposal
plot(path_init[,3],type="s")
mean(swap)
path_init[1:12,]

swap_proposal[which(swap == 1)]

sim_target = sim_parallel_temp$chains[[1]]
sim_target = sim[[9]]
#sim_target = sim_parallel_temp[[3]]

plot(sim$path_init[,9],type="s")

mean(sim$swap)

sim_parallel_temp$chain1
sim_parallel_temp$target_X_i
sim_parallel_temp$target_proposal

plot(sim$path_init[,9],type="s")



sim_target = sim_parallel_temp$chains[[1]]
sim_target = sim_parallel_temp$chain1
sim_target=sim$chains[[1]]

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

saveRDS(sim_pt_df,"mcmc/output_rcpp_parallel_pt_m_tidy-init2-conditional.rds")

### Using all variables

sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_c1-2temp_tidy-init1-hyperpars2-noSwap.rds")
sim_pt_df2 = readRDS("mcmc/output_rcpp_parallel_pt_m_c1-2temp_tidy-init2-hyperpars2-noSwap.rds")
sim_pt_df3 = readRDS("mcmc/output_rcpp_parallel_pt_m_c1-2temp_tidy-init3-hyperpars2-noSwap.rds")
sim_pt_df4 = readRDS("mcmc/output_rcpp_parallel_pt_m_c1-2temp_tidy-init4-hyperpars2-noSwap.rds")
#sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_c1-2temp_tidy-init4-hyperpars2-noSwap.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_c2.1-2temp_tidy-init2-hyperpars2-noSwap.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_8temp_tidy-init1-noPrior.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_3temp_tidy-init1-hyperpars2.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_7temp_tidy-init1-hyperpars2.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_8temp_tidy-init1-hyperpars2.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp_tidy-init1-hyperpars2.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp_tidy-init1.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp_tidy-init2.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp_tidy-init1-noSwap.rds")
# sim_pt_df2 = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp_tidy-init2-noSwap.rds")
# sim_pt_df = readRDS("mcmc/output_rcpp_parallel_pt_m_c100-2temp_tidy-init1-hyperpars2-noSwap.rds")

sim = readRDS("mcmc/output_rcpp_parallel_pt_m_c1.8-9temp-init3-hyperpars2-diagnostics.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_c1.5-13temp-init1-hyperpars2-diagnostics.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_c1.8-9temp-init1-hyperpars2-diagnostics.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp-init1-diagnostics.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp-init2-diagnostics.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_8temp-init1-hyperpars2-diagnostics.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_10temp-init1-hyperpars2.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_3temp-init1-hyperpars2.rds")
# sim = readRDS("mcmc/output_rcpp_parallel_pt_m_8temp-init1-hyperpars2.rds")
sim$swap[1:12]
sim$swap_proposal[1:12]

plot(1:10001,sim$path_init[1:10001,1],type="s")
plot(1:20001,sim$path_init[1:20001,10],type="s")
sim_target = sim$chains[[10]]

### Using all variables
sim_pt_df %>% filter(state == 3) %>% dplyr::select(iteration,
                                                   mu_duration,
                                                   sigma_duration,
                                                   mu_surface,
                                                   sigma_surface,
                                                   mu_maxDepth,
                                                   sigma_maxDepth,
                                                   mu_step,
                                                   sigma_step,
                                                   kappa,
                                                   a,b,
                                                   lambda,
                                                   theta,
                                                   init,
                                                   tpm1,
                                                   tpm2,
                                                   tpm3) %>% pivot_longer(-iteration,names_to="variable",
                                                                          values_to = "values") %>% 
  #filter(iteration > 5001) %>% 
  ggplot(aes(x=iteration,y=values)) + facet_wrap(~variable,scales = "free_y") + geom_line()
  #ggplot(aes(x=values)) + facet_wrap(~variable,scales="free") + geom_histogram()

sim_pt_df$chain = 1
sim_pt_df2$chain = 2
sim_pt_df3$chain = 3
sim_pt_df4$chain = 4

sim_pt_df_total = rbind(sim_pt_df,
                        sim_pt_df2,
                        sim_pt_df3,
                        sim_pt_df4) 

sim_pt_df_total %>% filter(state == 2) %>% dplyr::select(iteration,
                                                         chain,
                                                         mu_duration,
                                                         sigma_duration,
                                                         mu_surface,
                                                         sigma_surface,
                                                         mu_maxDepth,
                                                         sigma_maxDepth,
                                                         mu_step,
                                                         sigma_step,
                                                         kappa,
                                                         a,b,
                                                         lambda,
                                                         theta,
                                                         init,
                                                         tpm1,
                                                         tpm2,
                                                         tpm3) %>% pivot_longer(-c(iteration,chain),names_to="variable",
                                                                                                    values_to = "values") %>% 
  filter(iteration < 15000,iteration > 6000) %>% #filter(chain == 1) %>% 
  #ggplot(aes(x=iteration,y=values,color=factor(chain))) + facet_wrap(~variable,scales = "free_y") + geom_line()
  ggplot(aes(x=values,fill=factor(chain))) + facet_wrap(~variable,scales="free") + geom_histogram(position = "identity",alpha=0.5)



sim2 = sim_pt_df4 %>% #filter(iteration == 1) %>% 
  pivot_longer(-c(iteration,state),names_to = "parameter",values_to ="value") %>% 
  #mutate(par_state = paste0(parameter,".",state)) %>% 
  mutate(id_vector = case_when(parameter =="mu_duration" ~ state + 6,
                               parameter =="sigma_duration" ~ state + 9,
                               parameter == "mu_surface" ~ state + 12,
                               parameter == "sigma_surface" ~ state + 15,
                               parameter == "mu_maxDepth" ~ state + 18,
                               parameter == "sigma_maxDepth" ~ state + 21,
                               parameter == "mu_step" ~ state + 24,
                               parameter == "sigma_step" ~ state + 27,
                               parameter == "kappa" ~ state + 30,
                               parameter == "a" ~ state + 33,
                               parameter == "b" ~ state + 36,
                               parameter == "lambda" ~ state + 39,
                               parameter == "theta" ~ state + 48,
                               parameter == "tpm1" & state == 2 ~ 1,
                               parameter == "tpm1" & state == 3 ~ 2,
                               parameter == "tpm2" & state == 1 ~ 3,
                               parameter == "tpm2" & state == 3 ~ 4,
                               parameter == "tpm3" & state == 1 ~ 5,
                               parameter == "tpm3" & state == 2 ~ 6,
                               parameter == "tpm1" & state == 1 ~ 45,
                               parameter == "tpm2" & state == 2 ~ 46,
                               parameter == "tpm3" & state == 3 ~ 47,
                               parameter == "init" & state < 3 ~ state + 42,
                               parameter == "init" & state == 3 ~ 48,
                               TRUE ~ 0)) %>% 
  arrange(iteration,id_vector)

sim_lp__ = c()

sim2 %>% filter(iteration == 1) %>% dplyr::select(value) %>% pull()

for(i in 5002:20001){
  aux_x = sim2 %>% filter(iteration == i) %>% dplyr::select(value) %>% pull()
  sim_lp__ = c(sim_lp__, lp__stdVector(stan_data$N,
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
                                       aux_x))
}


hist(sim_lp__)
saveRDS(sim_lp__,"mcmc/sim_lp__init4.rds")

# sim_lp__1 = readRDS("mcmc/sim_lp__10temp-init1-noSwap.rds")
# sim_lp__2 = readRDS("mcmc/sim_lp__10temp-init2-noSwap.rds")
sim_lp__1 = readRDS("mcmc/sim_lp__init1.rds")
sim_lp__2 = readRDS("mcmc/sim_lp__init2.rds")
sim_lp__3 = readRDS("mcmc/sim_lp__init3.rds")
sim_lp__4 = readRDS("mcmc/sim_lp__init4.rds")

#sim_lp__ = rbind(cbind(sim_lp__1[6001:20001],1),cbind(sim_lp__2[6001:20001],2))
sim_lp__ = rbind(cbind(sim_lp__1,1),
                 cbind(sim_lp__2,2),
                 cbind(sim_lp__3,3),
                 cbind(sim_lp__4,4))
tibble(lp__ = sim_lp__[,1],chain=factor(sim_lp__[,2])) %>% ggplot(aes(x=lp__,fill=chain)) + 
  geom_histogram(binwidth = 4)

tibble(lp__ = sim_lp__[,1],chain=factor(sim_lp__[,2])) %>% ggplot(aes(x=lp__,fill=chain)) + 
  geom_histogram(binwidth = 1000,alpha=.5,position="identity")

tibble(lp__ = sim_lp__[,1],chain=factor(sim_lp__[,2])) %>% ggplot(aes(x=lp__,fill=chain)) + 
  geom_histogram()

length(sim_lp__1)

min(sim_lp__2)

### TO RUN ###
# 10k samples from v1_step1, with no swapping - in process
# 10k samples from v1_step2, with no swapping
# 10k samples from v1_step1, with no swapping and second set of initial values - in process
# 10k samples from v1_step2, with no swapping and second set of initial values



#### Drafts ###
# updated_theta_star = sim_target[,100]
# 
# for(k in 2:3){
#   t1 = Sys.time()
#   sim_parallel_temp_iter = rcpp_parallel_pt_cw_M_target_posterior(stan_data,nsim=100,init = updated_theta_star,
#                                                                   temp_vector = as.numeric(1:10),
#                                                                   #data for parallel computing
#                                                                   stan_data$N,
#                                                                   stan_data$n,
#                                                                   stan_data$n_ind,
#                                                                   stan_data$ID_init,
#                                                                   stan_data$ID,
#                                                                   stan_data$x_duration_init,
#                                                                   stan_data$x_surface_init,
#                                                                   stan_data$x_maxDepth_init,
#                                                                   stan_data$x_lunges_init,
#                                                                   stan_data$x_step_init,
#                                                                   stan_data$x_angle_init,
#                                                                   stan_data$x_headVar_init,
#                                                                   stan_data$x_duration,
#                                                                   stan_data$x_surface,
#                                                                   stan_data$x_maxDepth,
#                                                                   stan_data$x_lunges,
#                                                                   stan_data$x_step,
#                                                                   stan_data$x_angle,
#                                                                   stan_data$x_headVar)
#   
#   sim_target = cbind(sim_target[,-c((k-1)*100)],sim_parallel_temp_iter[[1]])
#   
#   updated_theta_star = sim_target[,k*100]
#   t2 = Sys.time()
#   print(t2-t1) # For 100 it, 10 temps, this will abe aprox. 73 min, let's see if that happens
# }
# 
# 
# saveRDS(sim_parallel_temp,"mcmc/output_rcpp_parallel_pt_m_10temp_v2.rds")
# 
# updated_theta_star
# sim_target = sim_parallel_temp_iter[[1]]


### This code is for testing ###
# sim_pt_df %>% pivot_longer(-c(state,iteration),names_to="parameter",values_to="value") %>% 
#   filter(value < 0)
# 
# sim_pt_df %>% dplyr::select(iteration,state,tpm3) %>% group_by(iteration) %>% 
#   summarise(sum = sum(tpm3))
# 
# sim_pt_df %>% dplyr::select(iteration,state,tpm2) %>% group_by(iteration) %>% 
#   summarise(sum = sum(tpm2)) #%>% dplyr::select(sum) %>% 
#   unique() %>% pull()

### K = 4 ###

# seed 1 -> 1378 (25671.0409)
# seed 2 -> 1322 (25667.8512)
# seed 3 -> 982 (25670.1908)
# seed 4 -> 396 (25661.425546) ---
# seed 5 -> 753 (25668.8220)
# seed 6 -> 754 (25670.8587)
# seed 7 -> 1408 (25669.3832)
# seed 8 -> 1058 (25671.8905)
# seed 9 -> 1195 (25664.1685)
# seed 10 -> 1041 (25672.7788)

### K = 5 ###

# seed 1 -> 1016 (25658.1807) ---
# seed 2 -> 846 (25660.2671)
# seed 3 -> 805 (25666.9617)
# seed 4 -> 650 (25661.2545)
# seed 5 -> 653 (25667.4498)
# seed 6 -> 722 (25667.8627)
# seed 7 -> 1320 (25669.8284)
# seed 8 -> 83 (25667.6982)
# seed 9 -> 341 (25666.2894)
# seed 10 -> 396 (25670.8235)

### MODEL 3 ###

# seed 1 -> 1337 (25676.8376)
# seed 2 -> 1382 (25680.5594)
# seed 3 -> 1243 (25682.4854)
# seed 4 -> 1396 (25682.3883)
# seed 5 -> 101 (25703.5097)
# seed 6 -> 745 (25665.3314) ---
# seed 7 -> 1226 (25670.3122)
# seed 8 -> 672 (25699.8331)
# seed 9 -> 731 (25671.1354)
# seed 10 -> 719 (25705.2286)

# sim_pt_df %>% pivot_longer(-c(state,iteration),names_to="parameter",values_to="value") %>% 
#   filter(value < 0)
# 
# sim_pt_df %>% filter(lambda < 0 | a < 0 | b < 0) %>% 
#   dplyr::select(iteration,state,lambda,a,b)
# 
# 
