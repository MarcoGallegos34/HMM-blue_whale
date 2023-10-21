#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]] 
// [[Rcpp::export]]
double log_likelihood(const List list_data, NumericVector theta_star){
  
//   const int K = list_data["K"];
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
  const NumericVector x_duration = list_data["x_duration"];
  const NumericVector x_surface = list_data["x_surface"];
  const NumericVector x_maxDepth = list_data["x_maxDepth"];
  const IntegerVector x_lunges = list_data["x_lunges"];
  const NumericVector x_step = list_data["x_step"];
  const NumericVector x_angle = list_data["x_angle"];
  const NumericVector x_headVar = list_data["x_headVar"];
    
//   const NumericVector estimate = list_data["estimate"];

    // idea: in case there's something not working properly, add parameters in list and extract them
  
  // # Creation parameters duration distribution
  // # log_mu_duration = theta_star[7:9]
//     const NumericVector mu_duration = NumericVector::create(140.0730, 334.1716, 516.0935);
//     const NumericVector log_sigma_duration = NumericVector::create(4.384593, 5.355711, 4.867127);
    
//   // # Creation parameters surface distribution
//   // # log_mu_surface = theta_star[13:15]
//     const NumericVector mu_surface = NumericVector::create(70.32488, 86.22668, 150.99891);
//     const NumericVector log_sigma_surface = NumericVector::create(4.216404, 4.011726, 4.236169);
//   // # 
//   // # # Creation parameters surface distribution
//   // # # log_mu_maxDepth = theta_star[19:21]
//     const NumericVector mu_maxDepth = NumericVector::create(32.37764, 68.42765, 169.74468);
//     const NumericVector log_sigma_maxDepth = NumericVector::create(3.164668, 4.173679, 4.100981);
//   // # 
//   // # # Creation parameters surface distribution
//   // # # log_mu_step = theta_star[25:27]
//     const NumericVector mu_step = NumericVector::create(187.2168, 674.9498, 406.0646);
//     const NumericVector log_sigma_step = NumericVector::create(4.884306, 5.720659, 5.658780);
//   // # 
//     const NumericVector log_kappa = NumericVector::create(0.00469386, 1.10641333, -0.20338085);
//   // # 
//     const NumericVector log_a = NumericVector::create(-0.02355296, -0.69161822, 0.51447883);
//     const NumericVector log_b = NumericVector::create(0.7448906, 1.6923719, 0.4474129);
//   // # 
//   // # # Creation parameters lunges distribution
//     const NumericVector log_lambda = NumericVector::create(-0.3998631, -3.9124458, 1.2114287); // lambda parameters
  
  // # Creation parameters initial distribution
  // # init_raw = theta_star[43:44] # pi parameters
  // # init_raw = theta_star[(43-ll):(44-ll)] # pi parameters
  
  // # Creation parameters duration distribution
  // # log_mu_duration = mod1.nlm$estimate[7:9]
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
    mu_duration[i] = theta_star[i+6];
    log_sigma_duration[i] = theta_star[i+9];

    mu_surface[i] = theta_star[i+12];
    log_sigma_surface[i] = theta_star[i+15];

    mu_maxDepth[i] = theta_star[i+18];
    log_sigma_maxDepth[i] = theta_star[i+21];

    mu_step[i] = theta_star[i+24];
    log_sigma_step[i] = theta_star[i+27];

    mu_step[i] = theta_star[i+24];
    log_sigma_step[i] = theta_star[i+27];

    log_kappa[i] = theta_star[i+30];

    log_a[i] = theta_star[i+33];
    log_b[i] = theta_star[i+36];

    log_lambda[i] = theta_star[i + 39];

    }
  
  // # Creation parameters surface distribution
  // # log_mu_surface = mod1.nlm$estimate[13:15]
  
//   // # Creation parameters surface distribution
//   // # log_mu_maxDepth = mod1.nlm$estimate[19:21]
  
//   // # Creation parameters surface distribution
//   // # log_mu_step = mod1.nlm$estimate[25:27]  
  
//   // # Creation parameters lunges distribution  
  
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
  
  
  // # tpm 
  NumericMatrix tpm(N,N);
  
  // from here the values will be extracted from theta_star, not estimate

//    int ll = 36 + 1;
  // # Creation parameters initial distribution
//   NumericVector init_raw = {theta_star[43-ll], theta_star[44-ll]}; // pi parameters
  NumericVector init_raw = {theta_star[42], theta_star[43]}; // pi parameters

  // # weights pi_
//   NumericVector theta(K); //   estimate[45-ll];
//   theta[2-1] = exp(theta_star[45-ll])/(1+exp(theta_star[45-ll])); // ll = 36 - 1
//   theta[1-1] = 1 - theta[2-1];
  
  // return theta;
  
  // # init
//   NumericMatrix init(K,N);
  NumericVector init(N);
  init[2-1] = 1/(1+exp(-init_raw[1-1]));
  init[3-1] = 1/(1+exp(-init_raw[2-1]));
  init[1-1] = 1 - (init[2-1] + init[3-1]);
  
 
//   if(K == 2){
//     NumericVector init_raw2 = {theta_star[46-ll], theta_star[47-ll]};
//     // NumericVector init_raw2(N-1); //[(46-ll):(47-ll)];
//     // init_raw2 = {theta_star[46-ll], theta_star[47-ll]};
//     init(2-1,2-1) = exp(init_raw2[1-1])/(1+exp(init_raw2[1-1]));
//     init(2-1,3-1) = exp(init_raw2[2-1])/(1+exp(init_raw2[2-1]));
//     init(2-1,1-1) = 1 - (init(2-1,2-1) + init(2-1,3-1));
    
//   }
  
  
  double log_lik_total = 0.0;
  // NumericVector log_lik_total(n_ind);
  
  // #Creation tpm k = 1
  // Here's the problem!
  tpm(1-1,1-1) = 1 - (exp(theta_star[1-1]) + exp(theta_star[2-1]))/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
  tpm(1-1,2-1) = exp(theta_star[1-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
  tpm(1-1,3-1) = exp(theta_star[2-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
  
  tpm(2-1,1-1) = exp(theta_star[3-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
  tpm(2-1,2-1) = 1 - (exp(theta_star[3-1]) + exp(theta_star[4-1]))/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
  tpm(2-1,3-1) = exp(theta_star[4-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
  
  tpm(3-1,1-1) = exp(theta_star[5-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
  tpm(3-1,2-1) = exp(theta_star[6-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
  tpm(3-1,3-1) = 1 - (exp(theta_star[5-1]) + exp(theta_star[6-1]))/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));

  // Rprintf("sum of row 1 is %f",tpm[0][0][1] + tpm[0][0][2]+tpm[0][0][3]);
  // NumericVector test = {tpm[0][1][0],
  //                       tpm[0][1][1],
  //                       tpm[0][1][2]};
  
  // tpm[1-1][1-1][1-1] = .3333333333;
  // tpm[1-1][1-1][2-1] = .3333333333;
  // tpm[1-1][1-1][3-1] = .3333333333;
  // tpm[1-1][2-1][1-1] = .3333333333;
  // tpm[1-1][2-1][2-1] = .3333333333;
  // tpm[1-1][2-1][3-1] = .3333333333;
  // tpm[1-1][3-1][1-1] = .3333333333;
  // tpm[1-1][3-1][2-1] = .3333333333;
  // tpm[1-1][3-1][3-1] = .3333333333;

//   if(K == 2){
//     // #Creation tpm k = 2
//     int l = 47 - ll; // ll = 37
//     tpm[2-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
//     tpm[2-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
//     tpm[2-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
//     tpm[2-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
//     tpm[2-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
//     tpm[2-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
//     tpm[2-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
//     tpm[2-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
//     tpm[2-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    
//     // tpm[2-1][1-1][1-1] = tpm[1-1][1-1][1-1];
//     // tpm[2-1][1-1][2-1] = tpm[1-1][1-1][2-1];
//     // tpm[2-1][1-1][3-1] = tpm[1-1][1-1][3-1];
//     // tpm[2-1][2-1][1-1] = tpm[1-1][2-1][1-1];
//     // tpm[2-1][2-1][2-1] = tpm[1-1][2-1][2-1];
//     // tpm[2-1][2-1][3-1] = tpm[1-1][2-1][3-1];
//     // tpm[2-1][3-1][1-1] = tpm[1-1][3-1][1-1];
//     // tpm[2-1][3-1][2-1] = tpm[1-1][3-1][2-1];
//     // tpm[2-1][3-1][3-1] = tpm[1-1][3-1][3-1];

//   }  

  // Rprintf("row 2 col 1 is %f \n",   tpm[1][2][0]);
  // Rprintf("row 2 col 2 is %f \n",   tpm[1][2][1]);
  // Rprintf("row 2 col 3 is %f \n",   tpm[1][2][2]);
  // Rprintf("sum of row 1 is %f \n",  tpm[1][2][0] + 
  //                                   tpm[1][2][1] +
  //                                   tpm[1][2][2]);
  // Rprintf("Difference is %f \n",  tpm[0][2][0] - tpm[1][2][0] ); 
  // Rprintf("Difference is %f \n",  tpm[0][2][1] - tpm[1][2][1]);
  // Rprintf("Difference is %f \n",  tpm[0][2][2] - tpm[1][2][2] );
  // Rprintf("sum is %f \n",tpm[1][1][_]);
  
//   NumericVector test2(N);

  Environment pkg = Environment::namespace_env("CircStats");

  Function dvm = pkg["dvm"];
  
  // NumericVector test = dvm(x_angle[0] , _["mu"] = 0.0, _["kappa"] = kappa[0]); 
  // double test = R::dgamma(x_duration[0] , alpha_duration[0],theta_duration[0] , true); 
  // return log(*(test.begin()));
  // Rprintf("dvm dist %f",*(test.begin()));
  // double  test_double = *(test.begin()); 
  // return log(*(test.begin()));
  // return test;

//   NumericMatrix log_alpha_components(K,n_ind);
    
    
    NumericMatrix log_alpha(N,n_ind);

    for(IntegerVector::iterator i = ID_init.begin(); i != ID_init.end(); ++i){
        for (int state = 0; state < N; state++){
        
        log_alpha(state,*i-1) = log(init(state));
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
            NumericVector aux_dvm = dvm(x_angle_init[*i-1] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
            log_alpha(state,*i-1) += log(*(aux_dvm.begin()));

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

    // return log_alpha;
    // test2 = {log_alpha(0,0),
    //         log_alpha(1,0),
    //         log_alpha(2,0)};

    for(int i = 0; i < n; i++){
        
        NumericVector aux_log_alpha(N);
        
        for(int state = 0; state < N; state++){
        
        NumericVector aux(N); // this is the erro!!!!!
        
        aux = log_alpha(_,ID[i]-1);
        // for(int j = 0; j < N; j++){
        //    aux[j] = log_alpha(j,ID[i]-1);
        // }
        
        for(int previous_state = 0; previous_state < N; previous_state++){
            
            aux[previous_state] += log(tpm(previous_state,state)); 
            
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
            NumericVector aux_dvm = dvm(x_angle[i] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
            for(int previous_state = 0; previous_state < N; previous_state++){
            aux[previous_state] += log(*(aux_dvm.begin()));
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
        
        log_alpha(_,ID[i]-1) = aux_log_alpha;
    }

    // for(int i=0; i < n_ind; i++){
    //     log_alpha_components(k,i) = log(theta[k]) + max(log_alpha(_,i)) + log(sum(exp(log_alpha(_,i)-max(log_alpha(_,i)))));
    // }
    // return log_alpha;

  
  for(int i=0; i < n_ind; i++){
    log_lik_total +=  max(log_alpha(_,i)) + log(sum(exp(log_alpha(_,i)-max(log_alpha(_,i)))));
    // #log_sum_exp(log_alpha[,i]);
    // log_lik_total[i] = max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i)))));
    // log_lik_total += max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i)))));
    // Rprintf("Iteration %i : %f", i , max(log_alpha_components(_,i)) + log(sum(exp(log_alpha_components(_,i)-max(log_alpha_components(_,i))))));
  }
  
  return -log_lik_total;
  // return {-log_lik_total,37*log(2)};
  
  // #return(log_alpha)
  // #return(log_alpha)
//   return log_alpha;
    // return alpha_duration;

  // return lambda;
}