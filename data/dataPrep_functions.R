library(dplyr)
library(tidyr)
library(CircStats)

whales = read.csv("data/bwdives_clean.csv")

whales = whales %>% filter(ID != "bw13_193a")

# Data prep for first version of likelihood

# toy_ids = unique(whales$ID)
# list_whales = list()
# 
# for(i in 1:length(toy_ids)){
#   print(id)
#   aux_whale = whales %>% filter(ID == toy_ids[i]) %>% dplyr::select(Number,
#                                                                     DIVE.TIME,
#                                                                     SURFACE.TIME,
#                                                                     MAX.DEPTH,
#                                                                     steps,
#                                                                     turns,
#                                                                     var.head,
#                                                                     LUNGES)
#   
#   list_whales[[i]] = aux_whale
#   
# }
# 
# rm(toy_ids)

#### Preparation for Stan ####

# whales_clean = whales %>% mutate(EXPOSURE = case_when(EXPOSURE %in% c("Silent") ~ 0,
whales_clean = whales %>% mutate(EXPOSURE = case_when(EXPOSURE %in% c("None", "Silent") ~ 0,
                                                      TRUE ~ 1 ))

whales_clean$ID = as.integer(as.vector(factor(whales$ID,labels = c(3,1:2,4:37))))
whales_clean[is.na(whales)] = 999

whales_init = whales_clean %>% filter(Number == 1)
whales_seq = whales_clean%>% filter(Number != 1)

# Initial distribution
#delta = c(.2,.1,.7)

# # Transition probability matrix
Gamma.matrix = matrix(c(.4,.2,.1,.3,.7,.05,.3,.1,.85),3,3)
muSigma_duration = matrix(c(140,334,516,80,212,130),3,2)
muSigma_surface = matrix(c(70,86,151,68,55,69),3,2)
muSigma_maxDepth = matrix(c(32,68,170,24,65,60),3,2)
muSigma_step = matrix(c(189,675,406,134,305,287),3,2)
muKappa_angle = matrix(c(0,0,0,1,3.1,.8),3,2)
ab_headVar = matrix(c(1,.5,1.7,2.1,5.4,1.6),3,2)
lambda_lunges = c(.7,.05,3.4)

log.likelihood = function(obs,
                          theta_star){
  
  N = 3
  
  # Creation tpm 
  
  exp_alpha1 = exp(theta_star[1:2])
  exp_alpha2 = exp(theta_star[3:4])
  exp_alpha3 = exp(theta_star[5:6])
  
  tpm_row1 = c(exp_alpha1[1],exp_alpha1[2])/(1 + sum(exp_alpha1))
  final_tpm_row1 = c(1-sum(tpm_row1),tpm_row1[1],tpm_row1[2])
  tpm_row2 = c(exp_alpha2[1],exp_alpha2[2])/(1+ sum(exp_alpha2))
  final_tpm_row2 = c(tpm_row2[1],1-sum(tpm_row2),tpm_row2[2])
  tpm_row3 = c(exp_alpha3[1],exp_alpha3[2])/(1+ sum(exp_alpha3))
  final_tpm_row3 = c(tpm_row3[1],tpm_row3[2],1-sum(tpm_row3))
  
  tpm = matrix(c(final_tpm_row1, # first row tpm
                 final_tpm_row2, # second row tpm
                 final_tpm_row3), # third row tpm
               N,N,byrow = T)
  
  # Creation initial distribution
  #init = solve(t(diag(N)-tpm + matrix(1,N,N)),rep(1,N))
  
  # Creation parameters duration distribution
  parameters_duration = matrix(c(theta_star[7:9]^2/exp(theta_star[10:12])^2, # mu parameters
                                 theta_star[7:9]/exp(theta_star[10:12])^2), # sigma parameters
                               3,2)
  
  # Creation parameters surface distribution
  parameters_surface = matrix(c(theta_star[13:15]^2/exp(theta_star[16:18])^2, # mu parameters
                                theta_star[13:15]/exp(theta_star[16:18])^2), # sigma parameters
                              3,2)
  
  # Creation parameters surface distribution
  parameters_maxDepth = matrix(c(theta_star[19:21]^2/exp(theta_star[22:24])^2, # mu parameters
                                 theta_star[19:21]/exp(theta_star[22:24])^2), # sigma parameters
                               3,2)
  
  # Creation parameters step distribution
  parameters_step = matrix(c(theta_star[25:27]^2/exp(theta_star[28:30])^2, # mu parameters
                             theta_star[25:27]/exp(theta_star[28:30])^2), # sigma parameters
                           3,2)
  
  # Creation parameters angle distribution
  parameters_angle = matrix(c(rep(0,N), # mu parameters
                              exp(theta_star[31:33])), # sigma parameters
                            3,2)
  
  # Creation parameters head variation distribution
  parameters_headVar = matrix(c(exp(theta_star[34:36]), # a parameter
                                exp(theta_star[37:39])), # b parameter
                              3,2)
  
  # Creation parameters lunges distribution
  parameters_lunges = c(exp(theta_star[40:42])) # lambda parameters
  
  # Creation parameters initial distribution
  parameters_init = c(plogis(theta_star[43:44])) # pi parameters
  
  # Creation initial distribution
  init = c(1-sum(parameters_init),parameters_init) # pi parameters
  
  
  #print(parameters_lunges)
  nwhales = length(obs)
  
  total_log_lik = c()
  # log_alpha_id = matrix(NA,N,nwhales)
  
  for(whale_i in 1:nwhales){
    
    x_duration = obs[[whale_i]]$DIVE.TIME
    x_surface = obs[[whale_i]]$SURFACE.TIME
    x_maxDepth = obs[[whale_i]]$MAX.DEPTH
    x_step = obs[[whale_i]]$steps
    x_angle = obs[[whale_i]]$turns
    x_headVar = obs[[whale_i]]$var.head
    x_lunges = obs[[whale_i]]$LUNGES
    
    n = length(x_duration)
    
    log_alpha = matrix(NA,N,n)
    
    
    log_alpha[,1] = log(init)
    
    
    if(!is.na(x_duration[1])){
      log_alpha[,1] = log_alpha[,1] + dgamma(x_duration[1],
                                             shape = parameters_duration[,1],
                                             rate = parameters_duration[,2],log = T)
    }
    
    
    if(!is.na(x_surface[1])){
      log_alpha[,1] = log_alpha[,1] + dgamma(x_surface[1],
                                             parameters_surface[,1],
                                             parameters_surface[,2],log = T)
    }
    
    
    if(!is.na(x_maxDepth[1])){
      log_alpha[,1] = log_alpha[,1] + dgamma(x_maxDepth[1],
                                             parameters_maxDepth[,1],
                                             parameters_maxDepth[,2],log = T)
    }
    
    
    if(!is.na(x_headVar[1])){
      log_alpha[,1] = log_alpha[,1] + dbeta(x_headVar[1],
                                            shape1 = parameters_headVar[,1],
                                            shape2 = parameters_headVar[,2],log = T)
    }  
    
    
    if(!is.na(x_lunges[1])){
      log_alpha[,1] = log_alpha[,1] + dpois(x_lunges[1],
                                            parameters_lunges,log = T)
    }  
    
    
    if(!is.na(x_angle[1])){
      log_alpha[,1] = log_alpha[,1] + log(dvm(x_angle[1],
                                              mu = parameters_angle[,1],
                                              kappa = parameters_angle[,2])) # Von mises
    }
    
    
    if(!is.na(x_step[1])){
      log_alpha[,1] = log_alpha[,1] + dgamma(x_step[1],
                                             parameters_step[,1],
                                             parameters_step[,2],log = T)
    }
    
    
    
    for(t in 2:n){
      for(i in 1:N){
        aux = log_alpha[,t-1] + log(tpm[,i])
        
        if(!is.na(x_duration[t])){
          aux = aux + dgamma(x_duration[t],
                             parameters_duration[i,1],
                             parameters_duration[i,2],log = T)
        }
        
        
        if(!is.na(x_surface[t])){
          aux = aux + dgamma(x_surface[t],
                             parameters_surface[i,1],
                             parameters_surface[i,2],log = T)
        }
        
        
        if(!is.na(x_maxDepth[t])){
          aux = aux + dgamma(x_maxDepth[t],
                             parameters_maxDepth[i,1],
                             parameters_maxDepth[i,2],log = T)
        }
        
        
        if(!is.na(x_headVar[t])){
          aux = aux + dbeta(x_headVar[t],
                            shape1 = parameters_headVar[i,1],
                            shape2 = parameters_headVar[i,2],log = T)
        }
        
        
        if(!is.na(x_lunges[t])){
          aux = aux + dpois(x_lunges[t],
                            parameters_lunges[i],log = T)
        }
        
        if(!is.na(x_angle[t])){
          aux = aux + log(dvm(x_angle[t],
                              mu = parameters_angle[i,1],
                              kappa = parameters_angle[i,2])) # Von mises
        }
        
        if(!is.na(x_step[t])){
          aux = aux + dgamma(x_step[t],
                             parameters_step[i,1],
                             parameters_step[i,2],log = T)
        }
        
        
        max_aux = max(aux)
        log_alpha[i,t] = max_aux+log(sum(exp(aux-max_aux)))
      }
    }
    
    max_aux = max(log_alpha[,n])
    total_log_lik[whale_i] = max_aux + log(sum(exp(log_alpha[,n] - max_aux)))
  }
  
  
  return(-sum(total_log_lik))
}

log.likelihood_v2 = function(list_data,
                             theta_star){
  
  N = list_data$N
  n = list_data$n
  n_ind = list_data$n_ind
  ID_init = list_data$ID_init
  ID = list_data$ID
  x_duration_init = list_data$x_duration_init
  x_surface_init = list_data$x_surface_init
  x_maxDepth_init = list_data$x_maxDepth_init
  x_lunges_init = list_data$x_lunges_init
  x_step_init = list_data$x_step_init
  x_angle_init = list_data$x_angle_init
  x_headVar_init = list_data$x_headVar_init
  x_duration = list_data$x_duration
  x_surface = list_data$x_surface
  x_maxDepth = list_data$x_maxDepth
  x_lunges = list_data$x_lunges
  x_step = list_data$x_step
  x_angle = list_data$x_angle
  x_headVar = list_data$x_headVar
  
  # Creation parameters duration distribution
  # log_mu_duration = theta_star[7:9]
  mu_duration = theta_star[7:9]
  log_sigma_duration = theta_star[10:12]
  
  # Creation parameters surface distribution
  # log_mu_surface = theta_star[13:15]
  mu_surface = theta_star[13:15]
  log_sigma_surface = theta_star[16:18]
  
  # Creation parameters surface distribution
  # log_mu_maxDepth = theta_star[19:21]
  mu_maxDepth = theta_star[19:21]
  log_sigma_maxDepth = theta_star[22:24]
  
  # Creation parameters surface distribution
  # log_mu_step = theta_star[25:27]
  mu_step = theta_star[25:27]
  log_sigma_step = theta_star[28:30]
  
  log_kappa = theta_star[31:33]
  
  log_a = theta_star[34:36]
  log_b = theta_star[37:39]
  
  # Creation parameters lunges distribution
  log_lambda = theta_star[40:42] # lambda parameters
  
  # Creation parameters initial distribution
  init_raw = theta_star[43:44] # pi parameters
  
  
  # parameters for duration distribution
  # alpha_duration = exp(2*log_mu_duration  - 2*log_sigma_duration)
  alpha_duration = mu_duration^2 / exp(log_sigma_duration)^2
  beta_duration = mu_duration / exp(log_sigma_duration)^2
  
  # parameters for surface distribution
  alpha_surface = mu_surface^2 /exp(log_sigma_surface)^2
  beta_surface = mu_surface / exp(log_sigma_surface)^2
  
  # parameters for maxDepth distribution
  alpha_maxDepth = mu_maxDepth^2 / exp(log_sigma_maxDepth)^2
  beta_maxDepth = mu_maxDepth / exp(log_sigma_maxDepth)^2
  
  # parameters for step distribution
  alpha_step = mu_step^2 / exp(log_sigma_step)^2
  beta_step = mu_step / exp(log_sigma_step)^2
  
  # parameters for angle distribution
  a = exp(log_a)
  b = exp(log_b)
  
  # parameters for angle distribution
  kappa = exp(log_kappa)
  
  #parameters for lunges distribution
  lambda = exp(log_lambda)
  
  # tpm
  tpm = matrix(NA,N,N)
  
  # init
  init = c()
  init[2] = 1/(1+exp(-init_raw[1]))
  init[3] = 1/(1+exp(-init_raw[2]))
  init[1] = 1- sum(init[2:3])
  
  log_alpha = matrix(NA,N,n_ind)
  log_lik_total = 0.0
  
  #Creation tpm 
  tpm[1,1] = 1 - sum(exp(theta_star[1:2]))/(1 + sum(exp(theta_star[1:2])))
  tpm[1,2] = exp(theta_star[1])/(1 + sum(exp(theta_star[1:2])))
  tpm[1,3] = exp(theta_star[2])/(1 + sum(exp(theta_star[1:2])))
  
  tpm[2,1] = exp(theta_star[3])/(1 + sum(exp(theta_star[3:4])))
  tpm[2,2] = 1 - sum(exp(theta_star[3:4]))/(1 + sum(exp(theta_star[3:4]))) 
  tpm[2,3] = exp(theta_star[4])/(1 + sum(exp(theta_star[3:4])))
  
  tpm[3,1] = exp(theta_star[5])/(1 + sum(exp(theta_star[5:6])))
  tpm[3,2] = exp(theta_star[6])/(1 + sum(exp(theta_star[5:6])))
  tpm[3,3] = 1 - sum(exp(theta_star[5:6]))/(1 + sum(exp(theta_star[5:6])))
  
  
  for(i in ID_init){
    for(state in 1:N){
      
      log_alpha[state,i] = log(init[state]) +
        dgamma(x_duration_init[i] , alpha_duration[state], beta_duration[state],log = T) +
        dgamma(x_surface_init[i] , alpha_surface[state], beta_surface[state],log = T) +
        dgamma(x_maxDepth_init[i] , alpha_maxDepth[state], beta_maxDepth[state],log=T) +
        dpois(x_lunges_init[i] , lambda[state],log =T)
      
      if(x_angle_init[i] != 999.0){
        log_alpha[state,i] = log_alpha[state,i] +
          log(dvm(x_angle_init[i] , 0.0, kappa[state]))
      }
      if(x_step_init[i] != 999.0){
        log_alpha[state,i] = log_alpha[state,i] +
          dgamma(x_step_init[i] , alpha_step[state], beta_step[state],log=T)
      }
      if(x_headVar_init[i] != 999.0){
        log_alpha[state,i] = log_alpha[state,i] +
          dbeta(x_headVar_init[i] , a[state], b[state],log=T)
      }
    }
  }
  
  for(i in 1:n){
    
    aux_log_alpha = rep(NA,N)
    
    for(state in 1:N){
      aux = log_alpha[,ID[i]] ### this is the erro!!!!!
      aux = aux + log(tpm[,state])
      aux = aux +
        dgamma(x_duration[i] , alpha_duration[state], beta_duration[state],log=T) +
        dgamma(x_surface[i] , alpha_surface[state], beta_surface[state],log=T) +
        dgamma(x_maxDepth[i] , alpha_maxDepth[state], beta_maxDepth[state],log=T) +
        dpois(x_lunges[i] , lambda[state],log=T)
      
      
      if(x_angle[i] != 999.0){
        aux = aux +
          log(dvm(x_angle[i] , 0.0, kappa[state]))
      }
      if(x_step[i] != 999.0){
        aux = aux +
          dgamma(x_step[i] , alpha_step[state], beta_step[state],log=T)
      }
      
      if(x_headVar[i] != 999.0){
        aux = aux +
          dbeta(x_headVar[i] , a[state], b[state],log=T)
      }
      
      # aux_log_alpha[state] = log_sum_exp(aux);
      aux_log_alpha[state] = max(aux) + log(sum(exp(aux-max(aux))))
      
    }
    
    # if(ID[i] == 1){
    #   print(aux_log_alpha)
    # }
    log_alpha[,ID[i]] = aux_log_alpha
  }
  
  for(i in 1:n_ind){
    #log_sum_exp(log_alpha[,i]);
    log_lik_total = log_lik_total + max(log_alpha[,i]) + log(sum(exp(log_alpha[,i]-max(log_alpha[,i]))))
  }
  return(-log_lik_total)
  #return(log_alpha)
}

theta_star = c(rep(0,6), # tpm
               muSigma_duration[,1],log(muSigma_duration[,2]), # duration
               muSigma_surface[,1],log(muSigma_surface[,2]), # surface
               muSigma_maxDepth[,1],log(muSigma_maxDepth[,2]), # maxDepth
               muSigma_step[,1],log(muSigma_step[,2]), # step
               log(muKappa_angle[,2]), # angle
               log(ab_headVar[,1]),log(ab_headVar[,2]), #varHead
               log(lambda_lunges), #lunges
               qlogis(c(1/3,1/3))) #init distribution

rm(Gamma.matrix)
rm(muSigma_duration)
rm(muSigma_surface)
rm(muSigma_maxDepth)
rm(muSigma_step)
rm(muKappa_angle)
rm(ab_headVar)
rm(lambda_lunges)

alphas_to_tpm = function(alphas){
  result = matrix(NA,3,3)
  
  result[1,1] = 1/(1 + exp(alphas[1]) + exp(alphas[2]))
  result[1,2] = exp(alphas[1])/(1 + exp(alphas[1]) + exp(alphas[2]))
  result[1,3] = exp(alphas[2])/(1 + exp(alphas[1]) + exp(alphas[2]))
  
  result[2,1] = exp(alphas[3])/(1 + exp(alphas[3]) + exp(alphas[4]))
  result[2,2] = 1/(1 + exp(alphas[3]) + exp(alphas[4]))
  result[2,3] = exp(alphas[4])/(1 + exp(alphas[3]) + exp(alphas[4]))
  
  result[3,1] = exp(alphas[5])/(1 + exp(alphas[5]) + exp(alphas[6]))
  result[3,2] = exp(alphas[6])/(1 + exp(alphas[5]) + exp(alphas[6]))
  result[3,3] = 1/(1 + exp(alphas[5]) + exp(alphas[6]))
  
  return(result)
}
