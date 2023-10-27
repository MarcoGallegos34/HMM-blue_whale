functions {
  
  void log_lik_lp(int N,
  int n,
  int n_ind,
  int n_ind_total,
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
  vector init,
  // real[,] tpm,
  // matrix tpm,
  vector tpm1,
  vector tpm2,
  vector tpm3,
  vector mu_duration,
  vector sigma_duration,
  vector mu_surface,
  vector sigma_surface,
  vector mu_maxDepth,
  vector sigma_maxDepth,
  vector kappa,
  vector a,
  vector b,
  vector mu_step,
  vector sigma_step,
  vector lambda,
  vector theta
  ){
    

    matrix[N,n_ind_total] log_alpha;
    
    // parameters for duration distribution
    vector[N] alpha_duration = mu_duration .^2 ./ sigma_duration.^2;
    vector[N] beta_duration = mu_duration ./ sigma_duration.^2;
  
    // parameters for surface distribution
    vector[N] alpha_surface = mu_surface.^2 ./ sigma_surface.^2;
    vector[N] beta_surface = mu_surface ./ sigma_surface.^2;
  
    // parameters for maxDepth distribution
    vector[N] alpha_maxDepth = mu_maxDepth.^2 ./ sigma_maxDepth.^2;
    vector[N] beta_maxDepth = mu_maxDepth ./  sigma_maxDepth.^2;
    
    // parameters for step distribution
    vector[N] alpha_step = mu_step.^2 ./ sigma_step.^2;
    vector[N] beta_step = mu_step ./ sigma_step.^2;
    
    // tpm
    matrix[N,N] tpm;
    tpm[1,] = to_row_vector(tpm1);
    tpm[2,] = to_row_vector(tpm2);
    tpm[3,] = to_row_vector(tpm3);
    // real tpm[N,N];
    // tpm[1,] = to_array_1d(tpm1);
    // tpm[2,] = to_array_1d(tpm2);
    // tpm[3,] = to_array_1d(tpm3);
  
    
    // real log_lik_total = 0.0;
    

    for(i in 1:n_ind){
      for(state in 1:N){
        
        log_alpha[state,ID_init[i]] = log(init[state]);
        log_alpha[state,ID_init[i]] += gamma_lpdf(x_duration_init[i] | alpha_duration[state], beta_duration[state]);
        log_alpha[state,ID_init[i]] += gamma_lpdf(x_surface_init[i] | alpha_surface[state], beta_surface[state]);
        log_alpha[state,ID_init[i]] += gamma_lpdf(x_maxDepth_init[i] | alpha_maxDepth[state], beta_maxDepth[state]);
        if(x_lunges_init[i] == 0){
          log_alpha[state,ID_init[i]] += log_sum_exp(log(theta[state]) ,log1m(theta[state]) + poisson_lpmf(x_lunges_init[i] | lambda[state]));
          // log_alpha[state,ID_init[i]] += log(theta[state]);
          
        } else{
          log_alpha[state,ID_init[i]] += log1m(theta[state]) + poisson_lpmf(x_lunges_init[i] | lambda[state]);
          // log_alpha[state,ID_init[i]] += log1m(theta[state]) + poisson_lpmf(x_lunges_init[i] | lambda[state]) - log1m_exp(-lambda[state]);
        }

        if(x_angle_init[i] != 999.0){
          log_alpha[state,ID_init[i]] += von_mises_lpdf(x_angle_init[i] | 0.0, kappa[state]);
        }
        if(x_step_init[i] != 999.0){
          log_alpha[state,ID_init[i]] += gamma_lpdf(x_step_init[i] | alpha_step[state], beta_step[state]);
        }
        if(x_headVar_init[i] != 999.0){
          log_alpha[state,ID_init[i]] += beta_lpdf(x_headVar_init[i] | a[state], b[state]);
        }
      //only the variables lat,lon, TURNS, STEPS and VAR.HEAD have missing values
      }
    }
      
      for(i in 1:n){
        vector[N] aux;
        vector[N] aux_log_alpha;
        
        for(state in 1:N){
          aux = log_alpha[,ID[i]];
          aux += log(tpm[,state]);
        
        for(previous_state in 1:N){
          aux[previous_state] += gamma_lpdf(x_duration[i] | alpha_duration[state], beta_duration[state]);
          aux[previous_state] += gamma_lpdf(x_surface[i] | alpha_surface[state], beta_surface[state]);
          aux[previous_state] += gamma_lpdf(x_maxDepth[i] | alpha_maxDepth[state], beta_maxDepth[state]);
          // aux[previous_state] += poisson_lpmf(x_lunges[i] | lambda[state]);
          if(x_lunges[i] == 0){
            aux[previous_state] += log_sum_exp(log(theta[state]) ,log1m(theta[state]) + poisson_lpmf(x_lunges[i] | lambda[state]));
            // aux[previous_state] += log(theta[state]);
            
            } else{
              aux[previous_state] += log1m(theta[state]) + poisson_lpmf(x_lunges[i] | lambda[state]);
              // aux[previous_state] += log1m(theta[state]) + poisson_lpmf(x_lunges[i] | lambda[state]) - log1m_exp(-lambda[state]);

              }

            
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

        aux_log_alpha[state] = log_sum_exp(aux);
        
        }
        
        log_alpha[,ID[i]] = aux_log_alpha;
      }
      
        for(i in 1:n_ind){
          target += log_sum_exp(log_alpha[,ID_init[i]]);
        }
        // return log_lik_total;
    }    
}

data{
    int N; // Total number states in hidden process
    int n; // Total number of observations
    int n_ind; // Total number of observed individuals to be considered for fitting model
    int n_ind_total; // Total number of observed individuals
    int ID[n];
    int ID_init[n_ind];
    
    // First we have the first observation for each time series of each individual
    int<lower=0> x_duration_init[n_ind];
    int<lower=0> x_surface_init[n_ind];
    int x_maxDepth_init[n_ind];
    int<lower=0> x_lunges_init[n_ind];
    real x_step_init[n_ind];
    real x_angle_init[n_ind];
    real x_headVar_init[n_ind];
    
    // Here we have the rest of observations of the time series - splitting should help for
    // making code easier to handle
    int<lower=0> x_duration[n];
    int<lower=0> x_surface[n];
    int x_maxDepth[n];
    int<lower=0> x_lunges[n];
    real x_step[n];
    real x_angle[n];
    real x_headVar[n];
    
    // hyperparameters for the initial distributions of parameters of interest
    
    // duration
    vector<lower=0>[N] mu_duration_mean;
    vector<lower=0>[N] mu_duration_sigma;
    vector<lower=0>[N] sigma_duration_alpha;
    vector<lower=0>[N] sigma_duration_beta;
    
    // surface
    vector<lower=0>[N] mu_surface_mean;
    vector<lower=0>[N] mu_surface_sigma;
    vector<lower=0>[N] sigma_surface_alpha;
    vector<lower=0>[N] sigma_surface_beta;

    // maxDepth
    vector<lower=0>[N] mu_maxDepth_mean;
    vector<lower=0>[N] mu_maxDepth_sigma;
    vector<lower=0>[N] sigma_maxDepth_alpha;
    vector<lower=0>[N] sigma_maxDepth_beta;
    
    // step
    vector<lower=0>[N] mu_step_mean;
    vector<lower=0>[N] mu_step_sigma;
    vector<lower=0>[N] sigma_step_alpha;
    vector<lower=0>[N] sigma_step_beta;
    
    // angle
    vector<lower=0>[N] kappa_alpha;
    vector<lower=0>[N] kappa_beta;
    
    // heading variance
    vector<lower=0>[N] a_alpha;
    vector<lower=0>[N] a_beta;
    vector<lower=0>[N] b_alpha;
    vector<lower=0>[N] b_beta;
    
    // lunges
    vector<lower=0>[N] lambda_alpha;
    vector<lower=0>[N] lambda_beta;
    
  
}
parameters{
  // duration
 positive_ordered[N] mu_duration; 
 vector<lower=0>[N] sigma_duration; 
 
  // surface
 vector<lower=0>[N] mu_surface; 
 vector<lower=0>[N] sigma_surface; 
 
 // max depth
 positive_ordered[N] mu_maxDepth; 
 vector<lower=0>[N] sigma_maxDepth;
 
 // step length
 positive_ordered[N] mu_step;
 vector<lower=0>[N] sigma_step;

 // turning angle
 vector<lower=0>[N] kappa;

 // heading variance
 vector<lower=0>[N] a;
 vector<lower=0>[N] b;

 // number of lunges
 vector<lower=0>[N] lambda;
 
 // weight for zero-inflation
 vector<lower=0,upper=1>[N] theta;
 
 // initial distribution
 simplex[N] init;
 
 // rows tpm
 simplex[N] tpm1;
 simplex[N] tpm2;
 simplex[N] tpm3;
 
}
model{
  
  for(i in 1:N){
    
    // duration
    target+= normal_lpdf(mu_duration[i] | mu_duration_mean[i] , mu_duration_sigma[i]); // 100, 10 
    target+= gamma_lpdf(sigma_duration[i] | sigma_duration_alpha[i] , sigma_duration_beta[i]); // 80, .75
   
    // surface
    target+= normal_lpdf(mu_surface[i] | mu_surface_mean[i] , mu_surface_sigma[i]); // 90,10 
    target+= gamma_lpdf(sigma_surface[i] | sigma_surface_alpha[i] , sigma_surface_beta[i]); // 65,.75 
   
    // // max depth
    target+= normal_lpdf(mu_maxDepth[i] | mu_maxDepth_mean[i] , mu_maxDepth_sigma[i] ); // 70,10
    target+= gamma_lpdf(sigma_maxDepth[i] | sigma_maxDepth_alpha[i] , sigma_maxDepth_beta[i] );
   
    // // step length
    target+= normal_lpdf(mu_step[i] | mu_step_mean[i] , mu_step_sigma[i]);
    target+= gamma_lpdf(sigma_step[i] | sigma_step_alpha[i] , sigma_step_beta[i]);

    // // turning angle
    target+= gamma_lpdf(kappa[i] | kappa_alpha[i], kappa_beta[i]);

    // heading variance
    target+= gamma_lpdf(a[i] | a_alpha[i], a_beta[i]);
    target+= gamma_lpdf(b[i] | b_alpha[i], b_beta[i]);

    // number of lunges
    target+= gamma_lpdf(lambda[i] | lambda_alpha[i], lambda_beta[i]);
    

  }
  
  // initial distribution
  target += dirichlet_lpdf(init | rep_vector(1,N));
  // rows tpm
  target += dirichlet_lpdf(tpm1 | rep_vector(1,N));
  target += dirichlet_lpdf(tpm2 | rep_vector(1,N));
  target += dirichlet_lpdf(tpm3 | rep_vector(1,N));
  
  log_lik_lp(N,n,n_ind,n_ind_total,
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
  init,
  tpm1,
  tpm2,
  tpm3,
  mu_duration,
  sigma_duration,
  mu_surface,
  sigma_surface,
  mu_maxDepth,
  sigma_maxDepth,
  kappa,
  a,
  b,
  mu_step,
  sigma_step,
  lambda,
  theta
  );

}

// todo: define tpm, init
// idea: put everything together in transformed parameters (or maybe just split it)
