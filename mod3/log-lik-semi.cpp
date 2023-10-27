// #include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]] 

// The log of the Von-misses density function of x given location mu and scale kappa
double log_dvm(double x, double mu, double kappa){
    return kappa*std::cos(x-mu) -(M_LN2 + log(M_PI) + log(std::cyl_bessel_i(0,kappa)));
}


// [[Rcpp::export]]
double log_likelihood_mod3_semi(const List list_data, NumericVector theta_star){
  
  const int K = list_data["K"];
  const int N = list_data["N"];
  const int n = list_data["n"];
  const int n_ind = list_data["n_ind"];
  IntegerVector ID_init = list_data["ID_init"];
  const IntegerVector ID = list_data["ID"];
  const NumericVector x_duration_init = list_data["x_duration_init"];
  const NumericVector x_surface_init = list_data["x_surface_init"];
  const NumericVector x_maxDepth_init = list_data["x_maxDepth_init"];
  const IntegerVector x_lunges_init = list_data["x_lunges_init"];
  const NumericVector x_step_init = list_data["x_step_init"];
  const NumericVector x_angle_init = list_data["x_angle_init"];
  const NumericVector x_headVar_init = list_data["x_headVar_init"];
  const IntegerVector x_exposure_init = list_data["x_exposure_init"];
  const NumericVector x_duration = list_data["x_duration"];
  const NumericVector x_surface = list_data["x_surface"];
  const NumericVector x_maxDepth = list_data["x_maxDepth"];
  const IntegerVector x_lunges = list_data["x_lunges"];
  const NumericVector x_step = list_data["x_step"];
  const NumericVector x_angle = list_data["x_angle"];
  const NumericVector x_headVar = list_data["x_headVar"];
  const IntegerVector x_exposure = list_data["x_exposure"];
    
  const NumericVector estimate = list_data["estimate"];

  NumericVector mu_duration(N);
  NumericVector log_sigma_duration(N);
  NumericVector mu_surface(N);
  NumericVector log_sigma_surface(N);
  NumericVector mu_maxDepth(N);
  NumericVector log_sigma_maxDepth(N);
  NumericVector mu_step(N);
  NumericVector log_sigma_step(N);
  NumericVector log_kappa(N);
  NumericVector log_a(N);
  NumericVector log_b(N);
  NumericVector log_lambda(N);

  for(int i = 0; i<3; ++i){
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
  
  // # Creation parameters surface distribution
  
  // # parameters for duration distribution
  // # alpha_duration = exp(2*log_mu_duration  - 2*log_sigma_duration)
  NumericVector alpha_duration(N);
  // NumericVector beta_duration(N);
  NumericVector theta_duration(N);

  // # parameters for surface distribution
  NumericVector alpha_surface(N);
  // NumericVector beta_surface(N);
  NumericVector theta_surface(N);
  
  // # parameters for maxDepth distribution
  NumericVector alpha_maxDepth(N);
  // NumericVector beta_maxDepth(N);
  NumericVector theta_maxDepth(N);

  // # parameters for step distribution
  NumericVector alpha_step(N);
  // NumericVector beta_step(N);
  NumericVector theta_step(N);
  
  // # parameters for angle distribution
  NumericVector a(N); 
  NumericVector b(N); 
  
  // # parameters for angle distribution
  NumericVector kappa(N); 
  
  // #parameters for lunges distribution
  NumericVector lambda(N); 

  for (int i = 0; i < N; i++){
    alpha_duration[i] = pow(mu_duration[i],2) / exp(2*log_sigma_duration[i]);
    // beta_duration[i] = mu_duration[i] / exp(2*log_sigma_duration[i]);
    theta_duration[i] =  exp(2*log_sigma_duration[i])/mu_duration[i];

    alpha_surface[i] = pow(mu_surface[i],2) / exp(2*log_sigma_surface[i]);
    // beta_surface[i] = mu_surface[i] / exp(2*log_sigma_surface[i]);
    theta_surface[i] =  exp(2*log_sigma_surface[i])/mu_surface[i];

    alpha_maxDepth[i] = pow(mu_maxDepth[i],2) / exp(2*log_sigma_maxDepth[i]);
    // beta_maxDepth[i] = mu_maxDepth[i] / exp(2*log_sigma_maxDepth[i]);
    theta_maxDepth[i] =  exp(2*log_sigma_maxDepth[i])/mu_maxDepth[i];
    
    alpha_step[i] = pow(mu_step[i],2) / exp(2*log_sigma_step[i]);
    // beta_step[i] = mu_step[i] / exp(2*log_sigma_step[i]);
    theta_step[i] =  exp(2*log_sigma_step[i])/mu_step[i];
    
    a[i] = exp(log_a[i]);
    b[i] = exp(log_b[i]);

    kappa[i] = exp(log_kappa[i]);
    
    lambda[i] = exp(log_lambda[i]);
  }
  
  
  
  // from here the values will be extracted from theta_star, not estimate

    int ll = 36 + 1;
    // int ll = 1;
    // # Creation parameters initial distribution
    arma::Col<double> init_raw = {theta_star[43-ll], theta_star[44-ll]}; // pi parameters

    // # weights pi_
    arma::Col<double> theta(K); //   estimate[45-ll];  
    arma::Col<double> theta_raw2 = {theta_star[45-ll],theta_star[54-ll],theta_star[63-ll]};
    // theta[2-1] = 1/(1+exp(-theta_raw2[1-1]));
    // theta[3-1] = 1/(1+exp(-theta_raw2[2-1]));
    // theta[4-1] = 1/(1+exp(-theta_raw2[3-1]));
    // theta[1-1] = 1 - (theta[2-1] + theta[3-1] + theta[4-1]);    
    theta[1-1] = 1/(1 + exp(theta_raw2[1-1]) + exp(theta_raw2[2-1]) + exp(theta_raw2[3-1]));
    theta[2-1] = exp(theta_raw2[1-1])/(1 + exp(theta_raw2[1-1]) + exp(theta_raw2[2-1]) + exp(theta_raw2[3-1]));
    theta[3-1] = exp(theta_raw2[2-1])/(1 + exp(theta_raw2[1-1]) + exp(theta_raw2[2-1]) + exp(theta_raw2[3-1]));
    theta[4-1] = exp(theta_raw2[3-1])/(1 + exp(theta_raw2[1-1]) + exp(theta_raw2[2-1]) + exp(theta_raw2[3-1]));
    
  // return theta;
  
    // # init
    arma::mat init(K,N);
    // init(1-1,2-1) = exp(init_raw[1-1])/(1+exp(init_raw[1-1]));
    // init(1-1,3-1) = exp(init_raw[2-1])/(1+exp(init_raw[2-1]));
    // init(1-1,1-1) = 1 - (init(1-1,2-1) + init(1-1,3-1));        
    init(1-1,1-1) = 1/(1 + exp(init_raw[1-1]) + exp(init_raw[2-1]));
    init(1-1,2-1) = exp(init_raw[1-1])/(1 + exp(init_raw[1-1]) + exp(init_raw[2-1]));
    init(1-1,3-1) = exp(init_raw[2-1])/(1 + exp(init_raw[1-1]) + exp(init_raw[2-1]));
        
    arma::Col<double> init_raw2 = {theta_star[46-ll],theta_star[47-ll]};
    // init(2-1,2-1) = 1/(1+exp(-init_raw2[1-1]));
    // init(2-1,3-1) = 1/(1+exp(-init_raw2[2-1]));
    // init(2-1,1-1) = 1 - (init(2-1,2-1) + init(2-1,3-1));
    init(2-1,1-1) = 1/(1 + exp(init_raw2[1-1]) + exp(init_raw2[2-1]));
    init(2-1,2-1) = exp(init_raw2[1-1])/(1 + exp(init_raw2[1-1]) + exp(init_raw2[2-1]));
    init(2-1,3-1) = exp(init_raw2[2-1])/(1 + exp(init_raw2[1-1]) + exp(init_raw2[2-1]));
    
    
    arma::Col<double> init_raw3 = {theta_star[55-ll],theta_star[56-ll]};
    // init(3-1,2-1) = 1/(1+exp(-init_raw3[1-1]));
    // init(3-1,3-1) = 1/(1+exp(-init_raw3[2-1]));
    // init(3-1,1-1) = 1 - (init(3-1,2-1) + init(3-1,3-1));
    init(3-1,1-1) = 1/(1 + exp(init_raw3[1-1]) + exp(init_raw3[2-1]));
    init(3-1,2-1) = exp(init_raw3[1-1])/(1 + exp(init_raw3[1-1]) + exp(init_raw3[2-1]));
    init(3-1,3-1) = exp(init_raw3[2-1])/(1 + exp(init_raw3[1-1]) + exp(init_raw3[2-1]));


    arma::Col<double> init_raw4 = {theta_star[64-ll],theta_star[65-ll]};
    // init(4-1,2-1) = 1/(1+exp(-init_raw4[1-1]));
    // init(4-1,3-1) = 1/(1+exp(-init_raw4[2-1]));
    // init(4-1,1-1) = 1 - (init(4-1,2-1) + init(4-1,3-1));
    init(4-1,1-1) = 1/(1 + exp(init_raw4[1-1]) + exp(init_raw4[2-1]));
    init(4-1,2-1) = exp(init_raw4[1-1])/(1 + exp(init_raw4[1-1]) + exp(init_raw4[2-1]));
    init(4-1,3-1) = exp(init_raw4[2-1])/(1 + exp(init_raw4[1-1]) + exp(init_raw4[2-1]));
        
  
    double log_lik_total = 0.0;

    // // tpm 
    // double tpm[K][N][N];

    // // Creation tpm k = 1
    // tpm[1-1][1-1][1-1] = 1 - (exp(theta_star[1-1]) + exp(theta_star[2-1]))/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
    // tpm[1-1][1-1][2-1] = exp(theta_star[1-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
    // tpm[1-1][1-1][3-1] = exp(theta_star[2-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));

    // tpm[1-1][2-1][1-1] = exp(theta_star[3-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
    // tpm[1-1][2-1][2-1] = 1 - (exp(theta_star[3-1]) + exp(theta_star[4-1]))/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
    // tpm[1-1][2-1][3-1] = exp(theta_star[4-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));

    // tpm[1-1][3-1][1-1] = exp(theta_star[5-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
    // tpm[1-1][3-1][2-1] = exp(theta_star[6-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
    // tpm[1-1][3-1][3-1] = 1 - (exp(theta_star[5-1]) + exp(theta_star[6-1]))/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));

    // // Creation tpm k = 2
    // int l = 47 - ll; // ll = 37
    // tpm[2-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    // tpm[2-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    // tpm[2-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));

    // tpm[2-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    // tpm[2-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    // tpm[2-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));

    // tpm[2-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    // tpm[2-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    // tpm[2-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    

    // // Creation tpm k = 3
    // // l = 47 + 9 - ll 
    // l += 9; 
    // tpm[3-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    // tpm[3-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    // tpm[3-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));

    // tpm[3-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    // tpm[3-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    // tpm[3-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));

    // tpm[3-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    // tpm[3-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    // tpm[3-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    

    // // Creation tpm k = 4
    // // l = 47 + 18 - ll # 
    // l += 9;  
    // tpm[4-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    // tpm[4-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    // tpm[4-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));

    // tpm[4-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    // tpm[4-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    // tpm[4-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));

    // tpm[4-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    // tpm[4-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    // tpm[4-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    


  // Environment pkg = Environment::namespace_env("CircStats");

  // Function dvm = pkg["dvm"];
  
  NumericMatrix log_alpha_components(K,n_ind);
  
  for(int k = 0;  k < K; k++){
    
    NumericMatrix log_alpha(N,n_ind);

    for(IntegerVector::iterator i = ID_init.begin(); i != ID_init.end(); ++i){
      for (int state = 0; state < N; state++){
        
        log_alpha(state,*i-1) = log(init(k,state));
        log_alpha(state,*i-1) += R::dgamma(x_duration_init[*i-1] ,
                                            alpha_duration[state], 
                                            // 1/beta_duration[state], 
                                            theta_duration[state], 
                                            true);

        log_alpha(state,*i-1) += R::dgamma(x_surface_init[*i-1] ,
                                        alpha_surface[state], 
                                        // 1/beta_surface[state],
                                        theta_surface[state],
                                        true);

        log_alpha(state,*i-1) += R::dgamma(x_maxDepth_init[*i-1] , 
                                        alpha_maxDepth[state], 
                                        // 1/beta_maxDepth[state],
                                        theta_maxDepth[state],
                                        true);

        log_alpha(state,*i-1) += R::dpois(x_lunges_init[*i-1] , 
                                      lambda[state],
                                      true);
        
        if(x_angle_init[*i-1] != 999.0){
          // NumericVector aux_dvm = dvm(x_angle_init[*i-1] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
          // log_alpha(state,*i-1) += log(*(aux_dvm.begin()));
          log_alpha(state,*i-1) += log_dvm(x_angle_init[*i-1], 0.0, kappa[state]);

        // Rprintf("Dvm is %f", dvm(x_angle_init[*i-1] , _["mu"] = 0.0, _["kappa"] = kappa[state]));
        // Rprintf("Dvm is %f", x_lunges_init[*i-1]);
        }

        if(x_step_init[*i-1] != 999.0){
          log_alpha(state,*i-1) += R::dgamma(x_step_init[*i-1] , 
                                          alpha_step[state], 
                                          // 1/beta_step[state],
                                          theta_step[state],
                                          true);
        }
        if(x_headVar_init[*i-1] != 999.0){
          log_alpha(state,*i-1) += R::dbeta(x_headVar_init[*i-1] , 
                                      a[state], 
                                      b[state],
                                      true);
        }
      }
    }
    
    // test2 = {log_alpha(0,0),
    //         log_alpha(1,0),
    //         log_alpha(2,0)};

    for(int i = 0; i < n; i++){
      
          arma::Col<double> aux_log_alpha(N);
          arma::mat tpm(N,N);

          if(k == 0){
              // Creation tpm k = 1
              tpm(1-1,1-1) = 1/(1 + (exp(theta_star[1-1] + theta_star[72-ll]*x_exposure[i]) + exp(theta_star[2-1] + theta_star[73-ll]*x_exposure[i])));
              tpm(1-1,2-1) = exp(theta_star[1-1] + theta_star[72-ll]*x_exposure[i])/(1 + (exp(theta_star[1-1] + theta_star[72-ll]*x_exposure[i]) + exp(theta_star[2-1] + theta_star[73-ll]*x_exposure[i])));
              tpm(1-1,3-1) = exp(theta_star[2-1] + theta_star[73-ll]*x_exposure[i])/(1 + (exp(theta_star[1-1] + theta_star[72-ll]*x_exposure[i]) + exp(theta_star[2-1] + theta_star[73-ll]*x_exposure[i])));

              tpm(2-1,1-1) = exp(theta_star[3-1] + theta_star[74-ll]*x_exposure[i])/(1 + (exp(theta_star[3-1] + theta_star[74-ll]*x_exposure[i]) + exp(theta_star[4-1] + theta_star[75-ll]*x_exposure[i])));
              tpm(2-1,2-1) = 1/(1 + (exp(theta_star[3-1] + theta_star[74-ll]*x_exposure[i]) + exp(theta_star[4-1] + theta_star[75-ll]*x_exposure[i])));
              tpm(2-1,3-1) = exp(theta_star[4-1] + theta_star[75-ll]*x_exposure[i])/(1 + (exp(theta_star[3-1] + theta_star[74-ll]*x_exposure[i]) + exp(theta_star[4-1] + theta_star[75-ll]*x_exposure[i])));

              tpm(3-1,1-1) = exp(theta_star[5-1] + theta_star[76-ll]*x_exposure[i])/(1 + (exp(theta_star[5-1] + theta_star[76-ll]*x_exposure[i]) + exp(theta_star[6-1] + theta_star[77-ll]*x_exposure[i])));
              tpm(3-1,2-1) = exp(theta_star[6-1] + theta_star[77-ll]*x_exposure[i])/(1 + (exp(theta_star[5-1] + theta_star[76-ll]*x_exposure[i]) + exp(theta_star[6-1] + theta_star[77-ll]*x_exposure[i])));
              tpm(3-1,3-1) = 1/(1 + (exp(theta_star[5-1] + theta_star[76-ll]*x_exposure[i]) + exp(theta_star[6-1] + theta_star[77-ll]*x_exposure[i])));

          } else {
              int l = 47 - ll; // ll = 37
              l += 9*(k-1); 
              tpm(1-1,1-1) = 1/(1 + (exp(theta_star[1+l] + theta_star[72-ll]*x_exposure[i]) + exp(theta_star[2+l] + theta_star[73-ll]*x_exposure[i])));
              tpm(1-1,2-1) = exp(theta_star[1+l] + theta_star[72-ll]*x_exposure[i])/(1 + (exp(theta_star[1+l] + theta_star[72-ll]*x_exposure[i]) + exp(theta_star[2+l] + theta_star[73-ll]*x_exposure[i])));
              tpm(1-1,3-1) = exp(theta_star[2+l] + theta_star[73-ll]*x_exposure[i])/(1 + (exp(theta_star[1+l] + theta_star[72-ll]*x_exposure[i]) + exp(theta_star[2+l] + theta_star[73-ll]*x_exposure[i])));

              tpm(2-1,1-1) = exp(theta_star[3+l] + theta_star[74-ll]*x_exposure[i])/(1 + (exp(theta_star[3+l] + theta_star[74-ll]*x_exposure[i]) + exp(theta_star[4+l] + theta_star[75-ll]*x_exposure[i])));
              tpm(2-1,2-1) = 1/(1 + (exp(theta_star[3+l] + theta_star[74-ll]*x_exposure[i]) + exp(theta_star[4+l] + theta_star[75-ll]*x_exposure[i])));
              tpm(2-1,3-1) = exp(theta_star[4+l] + theta_star[75-ll]*x_exposure[i])/(1 + (exp(theta_star[3+l] + theta_star[74-ll]*x_exposure[i]) + exp(theta_star[4+l] + theta_star[75-ll]*x_exposure[i])));

              tpm(3-1,1-1) = exp(theta_star[5+l] + theta_star[76-ll]*x_exposure[i])/(1 + (exp(theta_star[5+l] + theta_star[76-ll]*x_exposure[i]) + exp(theta_star[6+l] + theta_star[77-ll]*x_exposure[i])));
              tpm(3-1,2-1) = exp(theta_star[6+l] + theta_star[77-ll]*x_exposure[i])/(1 + (exp(theta_star[5+l] + theta_star[76-ll]*x_exposure[i]) + exp(theta_star[6+l] + theta_star[77-ll]*x_exposure[i])));
              tpm(3-1,3-1) = 1/(1 + (exp(theta_star[5+l] + theta_star[76-ll]*x_exposure[i]) + exp(theta_star[6+l] + theta_star[77-ll]*x_exposure[i])));

          }



        // if(x_exposure[i] == 0){
        //   if(k == 0){
        //       // Creation tpm k = 1
        //       tpm(1-1,1-1) = 1 - (exp(theta_star[1-1]) + exp(theta_star[2-1]))/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
        //       tpm(1-1,2-1) = exp(theta_star[1-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
        //       tpm(1-1,3-1) = exp(theta_star[2-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));

        //       tpm(2-1,1-1) = exp(theta_star[3-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
        //       tpm(2-1,2-1) = 1 - (exp(theta_star[3-1]) + exp(theta_star[4-1]))/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
        //       tpm(2-1,3-1) = exp(theta_star[4-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));

        //       tpm(3-1,1-1) = exp(theta_star[5-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
        //       tpm(3-1,2-1) = exp(theta_star[6-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
        //       tpm(3-1,3-1) = 1 - (exp(theta_star[5-1]) + exp(theta_star[6-1]))/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));

        //   } else {
        //       int l = 47 - ll; // ll = 37
        //       l += 9*(k-1); 
        //       tpm(1-1,1-1) = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        //       tpm(1-1,2-1) = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        //       tpm(1-1,3-1) = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));

        //       tpm(2-1,1-1) = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
        //       tpm(2-1,2-1) = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
        //       tpm(2-1,3-1) = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));

        //       tpm(3-1,1-1) = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
        //       tpm(3-1,2-1) = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
        //       tpm(3-1,3-1) = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));

        //   }

        // }
        // if( x_exposure[i] == 1){
        //   if(k == 0){
        //       // Creation tpm k = 1
        //       tpm(1-1,1-1) = 1 - (exp(theta_star[1-1] + theta_star[66-1]) + exp(theta_star[2-1] + theta_star[67-1]))/(1 + (exp(theta_star[1-1] + theta_star[66-1]) + exp(theta_star[2-1] + theta_star[67-1])));
        //       tpm(1-1,2-1) = exp(theta_star[1-1] + theta_star[66-1])/(1 + (exp(theta_star[1-1] + theta_star[66-1]) + exp(theta_star[2-1] + theta_star[67-1])));
        //       tpm(1-1,3-1) = exp(theta_star[2-1] + theta_star[67-1])/(1 + (exp(theta_star[1-1] + theta_star[66-1]) + exp(theta_star[2-1] + theta_star[67-1])));

        //       tpm(2-1,1-1) = exp(theta_star[3-1] + theta_star[68-1])/(1 + (exp(theta_star[3-1] + theta_star[68-1]) + exp(theta_star[4-1] + theta_star[69-1])));
        //       tpm(2-1,2-1) = 1 - (exp(theta_star[3-1] + theta_star[68-1]) + exp(theta_star[4-1] + theta_star[69-1]))/(1 + (exp(theta_star[3-1] + theta_star[68-1]) + exp(theta_star[4-1] + theta_star[69-1])));
        //       tpm(2-1,3-1) = exp(theta_star[4-1] + theta_star[69-1])/(1 + (exp(theta_star[3-1] + theta_star[68-1]) + exp(theta_star[4-1] + theta_star[69-1])));

        //       tpm(3-1,1-1) = exp(theta_star[5-1] + theta_star[70-1])/(1 + (exp(theta_star[5-1] + theta_star[70-1]) + exp(theta_star[6-1] + theta_star[71-1])));
        //       tpm(3-1,2-1) = exp(theta_star[6-1] + theta_star[71-1])/(1 + (exp(theta_star[5-1] + theta_star[70-1]) + exp(theta_star[6-1] + theta_star[71-1])));
        //       tpm(3-1,3-1) = 1 - (exp(theta_star[5-1] + theta_star[70-1]) + exp(theta_star[6-1] + theta_star[71-1]))/(1 + (exp(theta_star[5-1] + theta_star[70-1]) + exp(theta_star[6-1] + theta_star[71-1])));

        //   } else {
        //       int l = 47 - ll; // ll = 37
        //       l += 9*(k-1); 
        //       tpm(1-1,1-1) = 1 - (exp(theta_star[1+l] + theta_star[66-1]) + exp(theta_star[2+l] + theta_star[67-1]))/(1 + (exp(theta_star[1+l] + theta_star[66-1]) + exp(theta_star[2+l] + theta_star[67-1])));
        //       tpm(1-1,2-1) = exp(theta_star[1+l] + theta_star[66-1])/(1 + (exp(theta_star[1+l] + theta_star[66-1]) + exp(theta_star[2+l] + theta_star[67-1])));
        //       tpm(1-1,3-1) = exp(theta_star[2+l] + theta_star[67-1])/(1 + (exp(theta_star[1+l] + theta_star[66-1]) + exp(theta_star[2+l] + theta_star[67-1])));

        //       tpm(2-1,1-1) = exp(theta_star[3+l] + theta_star[68-1])/(1 + (exp(theta_star[3+l] + theta_star[68-1]) + exp(theta_star[4+l] + theta_star[69-1])));
        //       tpm(2-1,2-1) = 1 - (exp(theta_star[3+l] + theta_star[68-1]) + exp(theta_star[4+l] + theta_star[69-1]))/(1 + (exp(theta_star[3+l] + theta_star[68-1]) + exp(theta_star[4+l] + theta_star[69-1])));
        //       tpm(2-1,3-1) = exp(theta_star[4+l] + theta_star[69-1])/(1 + (exp(theta_star[3+l] + theta_star[68-1]) + exp(theta_star[4+l] + theta_star[69-1])));

        //       tpm(3-1,1-1) = exp(theta_star[5+l] + theta_star[70-1])/(1 + (exp(theta_star[5+l] + theta_star[70-1]) + exp(theta_star[6+l] + theta_star[71-1])));
        //       tpm(3-1,2-1) = exp(theta_star[6+l] + theta_star[71-1])/(1 + (exp(theta_star[5+l] + theta_star[70-1]) + exp(theta_star[6+l] + theta_star[71-1])));
        //       tpm(3-1,3-1) = 1 - (exp(theta_star[5+l] + theta_star[70-1]) + exp(theta_star[6+l] + theta_star[71-1]))/(1 + (exp(theta_star[5+l] + theta_star[70-1]) + exp(theta_star[6+l] + theta_star[71-1])));

        //   }
        //   if(k == 0){
        //   //     // Creation tpm k = 1
        //   //     tpm(1-1,1-1) = 1 - (exp(theta_star[1-1]) + exp(theta_star[2-1]))/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
        //   //     tpm(1-1,2-1) = exp(theta_star[1-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
        //   //     tpm(1-1,3-1) = exp(theta_star[2-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));

        //   //     tpm(2-1,1-1) = exp(theta_star[3-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
        //   //     tpm(2-1,2-1) = 1 - (exp(theta_star[3-1]) + exp(theta_star[4-1]))/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
        //   //     tpm(2-1,3-1) = exp(theta_star[4-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));

        //   //     tpm(3-1,1-1) = exp(theta_star[5-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
        //   //     tpm(3-1,2-1) = exp(theta_star[6-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
        //   //     tpm(3-1,3-1) = 1 - (exp(theta_star[5-1]) + exp(theta_star[6-1]))/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));

        //   // } else {
        //   //     int l = 47 - ll; // ll = 37
        //   //     l += 9*(k-1); 
        //   //     tpm(1-1,1-1) = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        //   //     tpm(1-1,2-1) = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        //   //     tpm(1-1,3-1) = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));

        //   //     tpm(2-1,1-1) = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
        //   //     tpm(2-1,2-1) = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
        //   //     tpm(2-1,3-1) = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));

        //   //     tpm(3-1,1-1) = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
        //   //     tpm(3-1,2-1) = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
        //   //     tpm(3-1,3-1) = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));

        //   // }

        // }
      
      for(int state = 0; state < N; state++){
        
        NumericVector aux(N); // this is the erro!!!!!
        
        aux = log_alpha(_,ID[i]-1);
        // for(int j = 0; j < N; j++){
        //    aux[j] = log_alpha(j,ID[i]-1);
        // }
        
        for(int previous_state = 0; previous_state < N; previous_state++){

        //   aux[previous_state] += log(tpm[k][previous_state][state]); // if I use the same matrix, I get same results, something is weird...
          aux[previous_state] += log(tpm(previous_state,state)); // if I use the same matrix, I get same results, something is weird...
          
          aux[previous_state] += R::dgamma(x_duration[i] , 
                                          alpha_duration[state],
                                          // 1/beta_duration[state], 
                                          theta_duration[state], 
                                          true);

          aux[previous_state] += R::dgamma(x_surface[i] , 
                                          alpha_surface[state], 
                                          // 1/beta_surface[state], 
                                          theta_surface[state], 
                                          true);

          aux[previous_state] += R::dgamma(x_maxDepth[i] , 
                                          alpha_maxDepth[state], 
                                          // 1/beta_maxDepth[state], 
                                          theta_maxDepth[state], 
                                          true);

          aux[previous_state] += R::dpois(x_lunges[i] , 
                                          lambda[state], 
                                          true);

        }
        
        
        if(x_angle[i] != 999.0){
          // NumericVector aux_dvm = dvm(x_angle[i] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
          for(int previous_state = 0; previous_state < N; previous_state++){
            // aux[previous_state] += log(*(aux_dvm.begin()));
            aux[previous_state] += log_dvm(x_angle[i], 0.0, kappa[state]);

          }
        }
        if(x_step[i] != 999.0){
          for(int previous_state = 0; previous_state < N; previous_state++){
            aux[previous_state] += R::dgamma(x_step[i] , 
                                            alpha_step[state], 
                                            // 1/beta_step[state],
                                            theta_step[state],
                                            true);
          }
        }
        
        if(x_headVar[i] != 999.0){
          for(int previous_state = 0; previous_state < N; previous_state++){
            aux[previous_state] += R::dbeta(x_headVar[i] , a[state], b[state], true);
          }
        }
        
        // # aux_log_alpha[state] = log_sum_exp(aux);
        aux_log_alpha[state] = max(aux) + log(sum(exp(aux-max(aux))));
        
      }
      
  //     # if(ID[i] == 1){
  //     #   print(aux_log_alpha)
  //     # }
      // log_alpha(_,ID[i]-1) = aux_log_alpha;
      for(int j = 0 ; j < N; j++){
        log_alpha(j,ID[i]-1) = aux_log_alpha[j];
      }

    }
    
    for(int i=0; i < n_ind; i++){
      log_alpha_components(k,i) = log(theta[k]) + max(log_alpha(_,i)) + log(sum(exp(log_alpha(_,i)-max(log_alpha(_,i)))));
    }
    // return log_alpha;
  }
  
  for(int i=0; i < n_ind; i++){
    // #log_sum_exp(log_alpha[,i]);
    // log_lik_total[i] = max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i)))));
    log_lik_total += max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i)))));
    // Rprintf("Iteration %i : %f", i , max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i))))));
  }
  
  return -log_lik_total;
  // return {-log_lik_total,37*log(2)};
  
  // #return(log_alpha)
  // #return(log_alpha)
  // return log_alpha_components;


  // return lambda;
}
