functions{
  real  log_lik(int N, int K,
  int n,
  int n_ind,
  int[] ID,
  int[] ID_init,
  int[] x_duration_init,
  int[] x_surface_init,
  int[] x_maxDepth_init,
  int[] x_lunges_init,
  real[] x_step_init,
  real[] x_angle_init,
  real[] x_headVar_init,
  int[] x_duration,
  int[] x_surface,
  int[] x_maxDepth,
  int[] x_lunges,
  real[] x_step,
  real[] x_angle,
  real[] x_headVar,
  vector init_raw1,
  vector init_raw2,
  vector tpm11,
  vector tpm12,
  vector tpm13,
  vector tpm21,
  vector tpm22,
  vector tpm23,
  vector alpha_duration,
  vector beta_duration,
  vector alpha_surface,
  vector beta_surface,
  vector alpha_maxDepth,
  vector beta_maxDepth,
  vector kappa,
  vector a,
  vector b,
  vector alpha_step,
  vector beta_step,
  vector lambda,
  vector theta_raw){
    
  // weights conditional likelihood (pi term) 
  vector[K] theta;
  
  // init
  vector[N] init1;
  vector[N] init2;

  // tpm
  matrix[N,N] tpm1;
  matrix[N,N] tpm2;
  
  //Creation tpm 
  tpm1[1,1] = 1 - sum(exp(tpm11))/(1 + sum(exp(tpm11)));
  tpm1[1,2] = exp(tpm11[1])/(1 + sum(exp(tpm11)));
  tpm1[1,3] = exp(tpm11[2])/(1 + sum(exp(tpm11)));
  
  tpm1[2,1] = exp(tpm12[1])/(1 + sum(exp(tpm12)));
  tpm1[2,2] = 1 - sum(exp(tpm12))/(1 + sum(exp(tpm12))); 
  tpm1[2,3] = exp(tpm12[2])/(1 + sum(exp(tpm12)));
  
  tpm1[3,1] =   exp(tpm13[1])/(1 + sum(exp(tpm13)));
  tpm1[3,2] = exp(tpm13[2])/(1 + sum(exp(tpm13)));
  tpm1[3,3] = 1 - sum(exp(tpm13))/(1 + sum(exp(tpm13))); 

  //Creation tpm 
  tpm2[1,1] = 1 - sum(exp(tpm21))/(1 + sum(exp(tpm21)));
  tpm2[1,2] = exp(tpm21[1])/(1 + sum(exp(tpm21)));
  tpm2[1,3] = exp(tpm21[2])/(1 + sum(exp(tpm21)));
  
  tpm2[2,1] = exp(tpm22[1])/(1 + sum(exp(tpm22)));
  tpm2[2,2] = 1 - sum(exp(tpm22))/(1 + sum(exp(tpm22))); 
  tpm2[2,3] = exp(tpm22[2])/(1 + sum(exp(tpm22)));
  
  tpm2[3,1] =   exp(tpm23[1])/(1 + sum(exp(tpm23)));
  tpm2[3,2] = exp(tpm23[2])/(1 + sum(exp(tpm23)));
  tpm2[3,3] = 1 - sum(exp(tpm23))/(1 + sum(exp(tpm23))); 


  theta[1] = 1 - 1/(1+exp(-theta_raw[1]));
  theta[2] = 1/(1+exp(-theta_raw[1]));
  
  // init
  // init1[1] = 1 - (1/(1+exp(-init_raw1[1])) + 1/(1+exp(-init_raw1[2])));
  // init1[2] = 1/(1+exp(-init_raw1[1]));
  // init1[3] = 1/(1+exp(-init_raw1[2]));

  init1[1] = 1/(1+sum(exp(init_raw1)));
  init1[2] = exp(init_raw1[1])/(1+sum(exp(init_raw1)));
  init1[3] = exp(init_raw1[2])/(1+sum(exp(init_raw1)));

  // init2[1] = 1 - (1/(1+exp(-init_raw2[1])) + 1/(1+exp(-init_raw2[2])));
  // init2[2] = 1/(1+exp(-init_raw2[1]));
  // init2[3] = 1/(1+exp(-init_raw2[2]));

  init2[1] = 1/(1+sum(exp(init_raw2)));
  init2[2] = exp(init_raw2[1])/(1+sum(exp(init_raw2)));
  init2[3] = exp(init_raw2[2])/(1+sum(exp(init_raw2)));

  
  real log_lik_total = 0.0;


  matrix[K,n_ind] log_alpha_components;
  
  for(k in 1:K){
    
    matrix[N,n_ind] log_alpha;
    
    matrix[N,N] aux_tpm;
    vector[N] aux_init;
    
    if(k == 1){
      aux_init = init1;
      aux_tpm = tpm1;
    }
    if(k == 2){
      aux_init = init2;
      aux_tpm = tpm2;
    }

    for(i in ID_init){
      for (state in 1:N){
        
        log_alpha[state,i] = log(aux_init[state]);
        log_alpha[state,i] += gamma_lpdf(x_duration_init[i] | alpha_duration[state], beta_duration[state]);

        log_alpha[state,i] += gamma_lpdf(x_surface_init[i] | alpha_surface[state], beta_surface[state]);

        log_alpha[state,i] += gamma_lpdf(x_maxDepth_init[i] | alpha_maxDepth[state], beta_maxDepth[state]);

        log_alpha[state,i] += poisson_lpmf(x_lunges_init[i] | lambda[state]);
        
        if(x_angle_init[i] != 999.0){
          log_alpha[state,i] += von_mises_lpdf(x_angle_init[i] | 0.0, kappa[state]);

        }

        if(x_step_init[i] != 999.0){
          log_alpha[state,i] += gamma_lpdf(x_step_init[i] | alpha_step[state], beta_step[state]);
        }
        if(x_headVar_init[i] != 999.0){
          log_alpha[state,i] += beta_lpdf(x_headVar_init[i] | a[state], b[state]);
        }
      }
    }
    
    
    for(i in 1:n){
      
      vector[N] aux_log_alpha;
      
      for(state in 1:N){
        
        vector[N] aux;
        
        aux = log_alpha[,ID[i]];

        for(previous_state in 1:N){
          
          aux[previous_state] += log(aux_tpm[previous_state,state]); // if I use the same matrix, I get same results, something is weird...
          
          aux[previous_state] += gamma_lpdf(x_duration[i] | alpha_duration[state], beta_duration[state]);

          aux[previous_state] += gamma_lpdf(x_surface[i] | alpha_surface[state], beta_surface[state]);

          aux[previous_state] += gamma_lpdf(x_maxDepth[i] | alpha_maxDepth[state], beta_maxDepth[state]);

          aux[previous_state] += poisson_lpmf(x_lunges[i] | lambda[state]);

        }
        
        
        if(x_angle[i] != 999.0){
          for(previous_state in 1:N){
            aux[previous_state] += von_mises_lpdf(x_angle[i] | 0.0, kappa[state]);

          }
        }
        if(x_step[i] != 999.0){
          for(previous_state in 1:N){
            aux[previous_state] += gamma_lpdf(x_step[i] | alpha_step[state], beta_step[state]);
          }
        }
        
        if(x_headVar[i] != 999.0){
          for(previous_state in 1:N){
            aux[previous_state] += beta_lpdf(x_headVar[i] | a[state], b[state]);
          }
        }
        
        // aux_log_alpha[state] = max(aux) + log(sum(exp(aux-max(aux))));
        aux_log_alpha[state] = log_sum_exp(aux);
        
      }
      
      log_alpha[,ID[i]] = aux_log_alpha;
    }
    
    for(i in 1:n_ind){
      log_alpha_components[k,i] = log(theta[k]) + log_sum_exp(log_alpha[,i]); // max(log_alpha(_,i)) + log(sum(exp(log_alpha(_,i)-max(log_alpha(_,i)))));
    }
  }
  
  for(i in 1:n_ind){
    log_lik_total += log_sum_exp(log_alpha_components[,i]); //max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i)))));
  }
  
  return -log_lik_total;
}

  
}

data{
  
  int K;
  int N;
  int n;
  int n_ind;
  
  // First we have the first observation for each time series of each individual
  int<lower=0> ID_init[n_ind];
  int<lower=0> ID[n];
  
  int<lower=0> x_duration_init[n_ind];
  int<lower=0> x_surface_init[n_ind];
  int x_maxDepth_init[n_ind];
  int<lower=0> x_lunges_init[n_ind];
  real x_step_init[n_ind];
  real x_angle_init[n_ind];
  real x_headVar_init[n_ind];
  
  
  int<lower=0> x_duration[n];
  int<lower=0> x_surface[n];
  int x_maxDepth[n];
  int<lower=0> x_lunges[n];
  real x_step[n];
  real x_angle[n];
  real x_headVar[n];
    
  vector[44] estimate;

}

transformed data{
  
  vector[N] mu_duration;
  vector[N] log_sigma_duration;
  vector[N] mu_surface;
  vector[N] log_sigma_surface;
  vector[N] mu_maxDepth;
  vector[N] log_sigma_maxDepth;
  vector[N] mu_step;
  vector[N] log_sigma_step;
  vector[N] log_kappa;
  vector[N] log_a;
  vector[N] log_b;
  vector[N] log_lambda;

  
  // parameters for duration distribution
  vector[N] alpha_duration;
  vector[N] beta_duration;
  // vector[N] theta_duration;

  // parameters for surface distribution
  vector[N] alpha_surface;
  vector[N] beta_surface;
  // vector[N] theta_surface;
  
  // parameters for maxDepth distribution
  vector[N] alpha_maxDepth;
  vector[N] beta_maxDepth;
  // vector[N] theta_maxDepth;

  // parameters for step distribution
  vector[N] alpha_step;
  vector[N] beta_step;
  // vector[N] theta_step;
  
  // parameters for angle distribution
  vector[N] a; 
  vector[N] b; 
  
  // parameters for angle distribution
  vector[N] kappa; 
  
  // parameters for lunges distribution
  vector[N] lambda; 

  for(i in 1:3){
    
    mu_duration[i] = estimate[i+6];
    log_sigma_duration[i] = estimate[i+9];

    mu_surface[i] = estimate[i+12];
    log_sigma_surface[i] = estimate[i+15];

    mu_maxDepth[i] = estimate[i+18];
    log_sigma_maxDepth[i] = estimate[i+21];

    mu_step[i] = estimate[i+24];
    log_sigma_step[i] = estimate[i+27];

    mu_step[i] = estimate[i+24];
    log_sigma_step[i] = estimate[i+27];

    log_kappa[i] = estimate[i+30];

    log_a[i] = estimate[i+33];
    log_b[i] = estimate[i+36];

    log_lambda[i] = estimate[i + 39];

    }


  for (i in 1:N){
    alpha_duration[i] = pow(mu_duration[i],2) / exp(2*log_sigma_duration[i]);
    beta_duration[i] = mu_duration[i] / exp(2*log_sigma_duration[i]);
    // theta_duration[i] =  exp(2*log_sigma_duration[i])/mu_duration[i];

    alpha_surface[i] = pow(mu_surface[i],2) / exp(2*log_sigma_surface[i]);
    beta_surface[i] = mu_surface[i] / exp(2*log_sigma_surface[i]);
    // theta_surface[i] =  exp(2*log_sigma_surface[i])/mu_surface[i];

    alpha_maxDepth[i] = pow(mu_maxDepth[i],2) / exp(2*log_sigma_maxDepth[i]);
    beta_maxDepth[i] = mu_maxDepth[i] / exp(2*log_sigma_maxDepth[i]);
    // theta_maxDepth[i] =  exp(2*log_sigma_maxDepth[i])/mu_maxDepth[i];
    
    alpha_step[i] = pow(mu_step[i],2) / exp(2*log_sigma_step[i]);
    beta_step[i] = mu_step[i] / exp(2*log_sigma_step[i]);
    // theta_step[i] =  exp(2*log_sigma_step[i])/mu_step[i];
    
    a[i] = exp(log_a[i]);
    b[i] = exp(log_b[i]);

    kappa[i] = exp(log_kappa[i]);
    
    lambda[i] = exp(log_lambda[i]);
  }

}

parameters{
  
  // weights conditional likelihood (pi term) 
  vector[K-1] theta_raw;
  
  // initial distribution k = 1
  vector[N-1] init_raw1;
  
   // rows tpm k = 1
   vector[N-1] tpm11;
   vector[N-1] tpm12;
   vector[N-1] tpm13;
  
   // initial distribution k = 2
   vector[N-1] init_raw2;
  
   // rows tpm k = 2
   vector[N-1] tpm21;
   vector[N-1] tpm22;
   vector[N-1] tpm23;
  
}

transformed parameters{

  // // weights conditional likelihood (pi term) 
  // vector[K] theta;
  // 
  // // init
  // vector[N] init1;
  // vector[N] init2;
  // 
  // // tpm
  // matrix[N,N] tpm1;
  // matrix[N,N] tpm2;
  // 
  // //Creation tpm 
  // tpm1[1,1] = 1 - sum(exp(tpm11))/(1 + sum(exp(tpm11)));
  // tpm1[1,2] = exp(tpm11[1])/(1 + sum(exp(tpm11)));
  // tpm1[1,3] = exp(tpm11[2])/(1 + sum(exp(tpm11)));
  // 
  // tpm1[2,1] = exp(tpm12[1])/(1 + sum(exp(tpm12)));
  // tpm1[2,2] = 1 - sum(exp(tpm12))/(1 + sum(exp(tpm12))); 
  // tpm1[2,3] = exp(tpm12[2])/(1 + sum(exp(tpm12)));
  // 
  // tpm1[3,1] =   exp(tpm13[1])/(1 + sum(exp(tpm13)));
  // tpm1[3,2] = exp(tpm13[2])/(1 + sum(exp(tpm13)));
  // tpm1[3,3] = 1 - sum(exp(tpm13))/(1 + sum(exp(tpm13))); 
  // 
  // //Creation tpm 
  // tpm2[1,1] = 1 - sum(exp(tpm21))/(1 + sum(exp(tpm21)));
  // tpm2[1,2] = exp(tpm21[1])/(1 + sum(exp(tpm21)));
  // tpm2[1,3] = exp(tpm21[2])/(1 + sum(exp(tpm21)));
  // 
  // tpm2[2,1] = exp(tpm22[1])/(1 + sum(exp(tpm22)));
  // tpm2[2,2] = 1 - sum(exp(tpm22))/(1 + sum(exp(tpm22))); 
  // tpm2[2,3] = exp(tpm22[2])/(1 + sum(exp(tpm22)));
  // 
  // tpm2[3,1] =   exp(tpm23[1])/(1 + sum(exp(tpm23)));
  // tpm2[3,2] = exp(tpm23[2])/(1 + sum(exp(tpm23)));
  // tpm2[3,3] = 1 - sum(exp(tpm23))/(1 + sum(exp(tpm23))); 
  // 
  // 
  // theta[1] = 1 - 1/(1+exp(-theta_raw[1]));
  // theta[2] = 1/(1+exp(-theta_raw[1]));
  // 
  // // init
  // // init1[1] = 1 - (1/(1+exp(-init_raw1[1])) + 1/(1+exp(-init_raw1[2])));
  // // init1[2] = 1/(1+exp(-init_raw1[1]));
  // // init1[3] = 1/(1+exp(-init_raw1[2]));
  // 
  // init1[1] = 1/(1+sum(exp(init_raw1)));
  // init1[2] = exp(init_raw1[1])/(1+sum(exp(init_raw1)));
  // init1[3] = exp(init_raw1[2])/(1+sum(exp(init_raw1)));
  // 
  // // init2[1] = 1 - (1/(1+exp(-init_raw2[1])) + 1/(1+exp(-init_raw2[2])));
  // // init2[2] = 1/(1+exp(-init_raw2[1]));
  // // init2[3] = 1/(1+exp(-init_raw2[2]));
  // 
  // init2[1] = 1/(1+sum(exp(init_raw2)));
  // init2[2] = exp(init_raw2[1])/(1+sum(exp(init_raw2)));
  // init2[3] = exp(init_raw2[2])/(1+sum(exp(init_raw2)));

 
}

model{
  target += log_lik(N, K,
  n,
  n_ind,
  ID,
  ID_init,
  x_duration_init,
  x_surface_init,
  x_maxDepth_init,
  x_lunges_init,
  x_step_init,
  x_angle_init,
  x_headVar_init,
  x_duration,
  x_surface,
  x_maxDepth,
  x_lunges,
  x_step,
  x_angle,
  x_headVar,
  init_raw1,
  init_raw2,
  tpm11,
  tpm12,
  tpm13,
  tpm21,
  tpm22,
  tpm23,
  alpha_duration,
  beta_duration,
  alpha_surface,
  beta_surface,
  alpha_maxDepth,
  beta_maxDepth,
  kappa,
  a,
  b,
  alpha_step,
  beta_step,
  lambda,
  theta_raw);

}