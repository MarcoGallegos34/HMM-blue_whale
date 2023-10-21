log.likelihood_v2_mod2_semi_extended = function(list_data,
                                       theta_star){
  
  K = list_data$K
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
  
  # Creation parameters initial distribution
  # init_raw = theta_star[43:44] # pi parameters
  ll = 36
  # init_raw = theta_star[(43-ll):(44-ll)] # pi parameters
  
  # Creation parameters duration distribution
  # log_mu_duration = mod1.nlm$estimate[7:9]
  mu_duration = mod1.nlm$estimate[7:9]
  log_sigma_duration = mod1.nlm$estimate[10:12]
  
  # Creation parameters surface distribution
  # log_mu_surface = mod1.nlm$estimate[13:15]
  mu_surface = mod1.nlm$estimate[13:15]
  log_sigma_surface = mod1.nlm$estimate[16:18]
  
  # Creation parameters surface distribution
  # log_mu_maxDepth = mod1.nlm$estimate[19:21]
  mu_maxDepth = mod1.nlm$estimate[19:21]
  log_sigma_maxDepth = mod1.nlm$estimate[22:24]
  
  # Creation parameters surface distribution
  # log_mu_step = mod1.nlm$estimate[25:27]
  mu_step = mod1.nlm$estimate[25:27]
  log_sigma_step = mod1.nlm$estimate[28:30]
  
  log_kappa = mod1.nlm$estimate[31:33]
  
  log_a = mod1.nlm$estimate[34:36]
  log_b = mod1.nlm$estimate[37:39]
  
  # Creation parameters lunges distribution
  log_lambda = mod1.nlm$estimate[40:42] # lambda parameters
  
  # Creation parameters initial distribution
  init_raw = theta_star[(43-ll):(44-ll)] # pi parameters
  
  
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
  tpm = array(NA,c(K,N,N))
  
  # weights pi_
  if(K == 2){
    theta = theta_star[(45-ll)] ## 36
    theta[2] = 1/(1+exp(-theta[1]))
    theta[1] = 1 - sum(theta[2])
  }
  if(K == 3){
    theta_raw2 = c(theta_star[(45-ll)],theta_star[(54-ll)]) ## 36
    theta = c()
    theta[2] = 1/(1+exp(-theta_raw2[1]))
    theta[3] = 1/(1+exp(-theta_raw2[2]))
    theta[1] = 1- sum(theta[2:3])
    
  }
  
  if(K == 4){
    theta_raw2 = c(theta_star[(45-ll)],theta_star[(54-ll)],theta_star[(63-ll)]) ## 36
    theta = c()
    theta[2] = 1/(1+exp(-theta_raw2[1]))
    theta[3] = 1/(1+exp(-theta_raw2[2]))
    theta[4] = 1/(1+exp(-theta_raw2[3]))
    theta[1] = 1- sum(theta[2:4])
    
  }
  
  if(K == 5){
    theta_raw2 = c(theta_star[(45-ll)],theta_star[(54-ll)],theta_star[(63-ll)],theta_star[(72-ll)]) ## 36
    theta = c()
    theta[2] = 1/(1+exp(-theta_raw2[1]))
    theta[3] = 1/(1+exp(-theta_raw2[2]))
    theta[4] = 1/(1+exp(-theta_raw2[3]))
    theta[5] = 1/(1+exp(-theta_raw2[4]))
    theta[1] = 1- sum(theta[2:5])
    
  }
  
  # init
  init = array(NA,c(K,N))
  init[1,2] = 1/(1+exp(-init_raw[1]))
  init[1,3] = 1/(1+exp(-init_raw[2]))
  init[1,1] = 1- sum(init[1,2:3])
  
  if(K == 2){
    init_raw2 = theta_star[(46-ll):(47-ll)]
    init[2,2] = 1/(1+exp(-init_raw2[1]))
    init[2,3] = 1/(1+exp(-init_raw2[2]))
    init[2,1] = 1- sum(init[2,2:3])
    
  }
  
  if(K == 3){
    
    init_raw2 = theta_star[(46-ll):(47-ll)]
    init[2,2] = 1/(1+exp(-init_raw2[1]))
    init[2,3] = 1/(1+exp(-init_raw2[2]))
    init[2,1] = 1- sum(init[2,2:3])
    
    
    init_raw3 = theta_star[(55-ll):(56-ll)]
    init[3,2] = 1/(1+exp(-init_raw3[1]))
    init[3,3] = 1/(1+exp(-init_raw3[2]))
    init[3,1] = 1- sum(init[3,2:3])
    
  }
  
  if(K == 4){
    
    init_raw2 = theta_star[(46-ll):(47-ll)]
    init[2,2] = 1/(1+exp(-init_raw2[1]))
    init[2,3] = 1/(1+exp(-init_raw2[2]))
    init[2,1] = 1- sum(init[2,2:3])
    
    
    init_raw3 = theta_star[(55-ll):(56-ll)]
    init[3,2] = 1/(1+exp(-init_raw3[1]))
    init[3,3] = 1/(1+exp(-init_raw3[2]))
    init[3,1] = 1- sum(init[3,2:3])
    
    init_raw4 = theta_star[(64-ll):(65-ll)]
    init[4,2] = 1/(1+exp(-init_raw4[1]))
    init[4,3] = 1/(1+exp(-init_raw4[2]))
    init[4,1] = 1- sum(init[4,2:3])
    
  }
  
  if(K == 5){
    
    init_raw2 = theta_star[(46-ll):(47-ll)]
    init[2,2] = 1/(1+exp(-init_raw2[1]))
    init[2,3] = 1/(1+exp(-init_raw2[2]))
    init[2,1] = 1- sum(init[2,2:3])
    
    
    init_raw3 = theta_star[(55-ll):(56-ll)]
    init[3,2] = 1/(1+exp(-init_raw3[1]))
    init[3,3] = 1/(1+exp(-init_raw3[2]))
    init[3,1] = 1- sum(init[3,2:3])
    
    init_raw4 = theta_star[(64-ll):(65-ll)]
    init[4,2] = 1/(1+exp(-init_raw4[1]))
    init[4,3] = 1/(1+exp(-init_raw4[2]))
    init[4,1] = 1- sum(init[4,2:3])
    
    init_raw5 = theta_star[(73-ll):(74-ll)]
    init[5,2] = 1/(1+exp(-init_raw5[1]))
    init[5,3] = 1/(1+exp(-init_raw5[2]))
    init[5,1] = 1- sum(init[5,2:3])
    
  }
  
  log_lik_total = 0.0
  
  #Creation tpm k = 1
  tpm[1,1,1] = 1 - sum(exp(theta_star[1:2]))/(1 + sum(exp(theta_star[1:2])))
  tpm[1,1,2] = exp(theta_star[1])/(1 + sum(exp(theta_star[1:2])))
  tpm[1,1,3] = exp(theta_star[2])/(1 + sum(exp(theta_star[1:2])))
  
  tpm[1,2,1] = exp(theta_star[3])/(1 + sum(exp(theta_star[3:4])))
  tpm[1,2,2] = 1 - sum(exp(theta_star[3:4]))/(1 + sum(exp(theta_star[3:4]))) 
  tpm[1,2,3] = exp(theta_star[4])/(1 + sum(exp(theta_star[3:4])))
  
  tpm[1,3,1] = exp(theta_star[5])/(1 + sum(exp(theta_star[5:6])))
  tpm[1,3,2] = exp(theta_star[6])/(1 + sum(exp(theta_star[5:6])))
  tpm[1,3,3] = 1 - sum(exp(theta_star[5:6]))/(1 + sum(exp(theta_star[5:6])))
  
  if(K == 2){
    #Creation tpm k = 2
    l = 47 - ll # 
    tpm[2,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[2,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[2,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
  }
  
  if(K == 3){

    #Creation tpm k = 2
    l = 47 - ll # 
    tpm[2,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[2,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[2,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
    #Creation tpm k = 3
    l = 47 + 9 - ll # 
    tpm[3,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[3,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[3,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[3,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[3,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[3,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[3,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[3,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[3,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
  }

  if(K == 4){

    #Creation tpm k = 2
    l = 47 - ll # 
    tpm[2,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[2,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[2,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
    #Creation tpm k = 3
    l = 47 + 9 - ll # 
    tpm[3,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[3,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[3,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[3,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[3,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[3,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[3,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[3,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[3,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
    #Creation tpm k = 3
    l = 47 + 18 - ll # 
    tpm[4,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[4,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[4,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[4,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[4,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[4,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[4,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[4,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[4,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
  }

  if(K == 5){

    #Creation tpm k = 2
    l = 47 - ll # 
    tpm[2,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[2,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[2,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[2,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[2,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[2,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
    #Creation tpm k = 3
    l = 47 + 9 - ll # 
    tpm[3,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[3,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[3,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[3,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[3,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[3,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[3,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[3,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[3,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))

    #Creation tpm k = 4
    l = 47 + 18 - ll # 
    tpm[4,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[4,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[4,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[4,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[4,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[4,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[4,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[4,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[4,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
    #Creation tpm k = 5
    l = 47 + 27 - ll # 
    tpm[5,1,1] = 1 - sum(exp(theta_star[(1+l):(2+l)]))/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[5,1,2] = exp(theta_star[(1+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    tpm[5,1,3] = exp(theta_star[(2+l)])/(1 + sum(exp(theta_star[(1+l):(2+l)])))
    
    tpm[5,2,1] = exp(theta_star[(3+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[5,2,2] = 1 - sum(exp(theta_star[(3+l):(4+l)]))/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    tpm[5,2,3] = exp(theta_star[(4+l)])/(1 + sum(exp(theta_star[(3+l):(4+l)])))
    
    tpm[5,3,1] = exp(theta_star[(5+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[5,3,2] = exp(theta_star[(6+l)])/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    tpm[5,3,3] = 1 - sum(exp(theta_star[(5+l):(6+l)]))/(1 + sum(exp(theta_star[(5+l):(6+l)])))
    
  }

  log_alpha_components = matrix(NA,K,n_ind)
  
  for(k in 1:K){
    log_alpha = matrix(NA,N,n_ind)
    
    for(i in ID_init){
      for(state in 1:N){
        
        log_alpha[state,i] = log(init[k,state]) +
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
        aux = aux + log(tpm[k,,state])
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
        
        aux_log_alpha[state] = max(aux) + log(sum(exp(aux-max(aux))))
        
      }
      
      log_alpha[,ID[i]] = aux_log_alpha
    }
    
    for(i in 1:n_ind){
      log_alpha_components[k,i] = log(theta[k]) + max(log_alpha[,i]) + log(sum(exp(log_alpha[,i]-max(log_alpha[,i]))))
    }
  }
  
  for(i in 1:n_ind){
    log_lik_total = log_lik_total + max(log_alpha_components[,i]) + log(sum(exp(log_alpha_components[,i]-max(log_alpha_components[,i]))))
  }
  return(-log_lik_total)
  # return(log_alpha_components)
}

library(Rcpp)
source("data/dataPrep_functions.R")
source("stan-input/stan-data_init-values.R")
sourceCpp("mod2/log-lik-semi-extended.cpp")

# Modification of stan_data list
stan_data$K = 2
mod1.nlm = readRDS("mle/nlm.output.rds")
stan_data$estimate = mod1.nlm$estimate 


theta_star_2_semi = c(rep(0,6), # tpm k = 1
                      # muSigma_duration[,1],log(muSigma_duration[,2]), # duration
                      # muSigma_surface[,1],log(muSigma_surface[,2]), # surface
                      # muSigma_maxDepth[,1],log(muSigma_maxDepth[,2]), # maxDepth
                      # log(muSigma_step[,1]),log(muSigma_step[,2]), # step
                      # log(muKappa_angle[,2]), # angle
                      # log(ab_headVar[,1]),log(ab_headVar[,2]), #varHead
                      # log(lambda_lunges), #lunges
                      qlogis(c(1/3,1/3)), #init distribution k = 1
                      ### These are the important ones!!!
                      qlogis(c(1/2)), # pi_
                      qlogis(c(1/3,1/3)), #init distribution k = 2
                      rep(0,6)) # tpm k = 2

log.likelihood_v2_mod2_semi_extended(stan_data,theta_star_2_semi)
log_likelihood_v2_mod2_semi_extended(stan_data,theta_star_2_semi)

theta_star_3_semi = c(rep(0,6), # tpm k = 1
                      # muSigma_duration[,1],log(muSigma_duration[,2]), # duration
                      # muSigma_surface[,1],log(muSigma_surface[,2]), # surface
                      # muSigma_maxDepth[,1],log(muSigma_maxDepth[,2]), # maxDepth
                      # log(muSigma_step[,1]),log(muSigma_step[,2]), # step
                      # log(muKappa_angle[,2]), # angle
                      # log(ab_headVar[,1]),log(ab_headVar[,2]), #varHead
                      # log(lambda_lunges), #lunges
                      qlogis(c(1/3,1/3)), #init distribution k = 1
                      ### These are the important ones!!!
                      qlogis(c(1/3)), # pi_
                      qlogis(c(1/3,1/3)), #init distribution k = 2
                      rep(0,6), # tpm k = 2
                      qlogis(c(1/3)), # pi_ (k=3)
                      qlogis(c(1/3,1/3)), #init distribution k = 3
                      rep(0,6)) # tpm k = 3

stan_data$K = 3
log.likelihood_v2_mod2_semi_extended(stan_data,theta_star_3_semi)
log_likelihood_v2_mod2_semi_extended(stan_data,theta_star_3_semi)

theta_star_4_semi = c(rep(0,6), # tpm k = 1
                      # muSigma_duration[,1],log(muSigma_duration[,2]), # duration
                      # muSigma_surface[,1],log(muSigma_surface[,2]), # surface
                      # muSigma_maxDepth[,1],log(muSigma_maxDepth[,2]), # maxDepth
                      # log(muSigma_step[,1]),log(muSigma_step[,2]), # step
                      # log(muKappa_angle[,2]), # angle
                      # log(ab_headVar[,1]),log(ab_headVar[,2]), #varHead
                      # log(lambda_lunges), #lunges
                      qlogis(c(1/3,1/3)), #init distribution k = 1
                      ### These are the important ones!!!
                      qlogis(c(1/4)), # pi_
                      qlogis(c(1/3,1/3)), #init distribution k = 2
                      rep(0,6), # tpm k = 2
                      qlogis(c(1/4)), # pi_ (k=3)
                      qlogis(c(1/3,1/3)), #init distribution k = 3
                      rep(0,6), # tpm k = 3
                      qlogis(c(1/4)), # pi_ (k=4)
                      qlogis(c(1/3,1/3)), #init distribution k = 4
                      rep(0,6)) # tpm k = 4

stan_data$K = 4
log.likelihood_v2_mod2_semi_extended(stan_data,theta_star_4_semi)
log_likelihood_v2_mod2_semi_extended(stan_data,theta_star_4_semi)

theta_star_5_semi = c(rep(0,6), # tpm k = 1
                      # muSigma_duration[,1],log(muSigma_duration[,2]), # duration
                      # muSigma_surface[,1],log(muSigma_surface[,2]), # surface
                      # muSigma_maxDepth[,1],log(muSigma_maxDepth[,2]), # maxDepth
                      # log(muSigma_step[,1]),log(muSigma_step[,2]), # step
                      # log(muKappa_angle[,2]), # angle
                      # log(ab_headVar[,1]),log(ab_headVar[,2]), #varHead
                      # log(lambda_lunges), #lunges
                      qlogis(c(1/3,1/3)), #init distribution k = 1
                      ### These are the important ones!!!
                      qlogis(c(1/5)), # pi_
                      qlogis(c(1/3,1/3)), #init distribution k = 2
                      rep(0,6), # tpm k = 2
                      qlogis(c(1/5)), # pi_ (k=3)
                      qlogis(c(1/3,1/3)), #init distribution k = 3
                      rep(0,6), # tpm k = 3
                      qlogis(c(1/5)), # pi_ (k=4)
                      qlogis(c(1/3,1/3)), #init distribution k = 4
                      rep(0,6), # tpm k = 4
                      qlogis(c(1/5)), # pi_ (k=5)
                      qlogis(c(1/3,1/3)), #init distribution k = 5
                      rep(0,6)) # tpm k = 5

stan_data$K = 5
log.likelihood_v2_mod2_semi_extended(stan_data,theta_star_5_semi)
log_likelihood_v2_mod2_semi_extended(stan_data,theta_star_5_semi)

### K = 3 ###

# Generation of initial values
set.seed(23984)
library(gtools)
mod2_init_mat_3 = matrix(NA,26,150)

for(i in 1:150){
  aux_inits = rdirichlet(3,alpha=rep(1,3))
  aux_theta = rdirichlet(1,alpha=rep(1,3))
  
  mod2_init_mat_3[,i] = c(rnorm(6,0,5), # tpm k = 1
                        qlogis(aux_inits[1,2:3]), # init k = 1
                        qlogis(aux_theta[2]), # pi
                        qlogis(aux_inits[2,2:3]), # init k = 2
                        rnorm(6,0,5), # tpm k = 2
                        qlogis(aux_theta[3]), # pi for k = 3
                        qlogis(aux_inits[3,2:3]), # init k = 3
                        rnorm(6,0,5)) # tpm k = 3
}

log.likelihood_v2_mod2_semi_extended(stan_data,mod2_init_mat_3[,1])
log_likelihood_v2_mod2_semi_extended(stan_data,mod2_init_mat_3[,1])

sourceCpp("mod2/log-lik-semi-extended-roptim.cpp")

#Y = as.matrix(theta_star_3_semi,26,1)

semiMles = parallel_mle_estimator_extended(mod2_init_mat_3[,1:5],
                                  stan_data$K,
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
                                  stan_data$estimate,0,5)


t1 = Sys.time()
nestim = 15
mle_vec = c()
### Setting range from 1:935 would compute 14,960 different mles
for(i in 1:5){
  Y = mod2_init_mat_3[,(nestim - 14):nestim]
  t1
  semiMles = parallel_mle_estimator_extended(Y,
                                stan_data$K,
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
                                stan_data$estimate,0,15)
  mle_vec = c(mle_vec,semiMles)
  gc()
  nestim = nestim + 15
  print(paste0("iteration ",i," has finished"))
}
t2 = Sys.time()
t2-t1




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

