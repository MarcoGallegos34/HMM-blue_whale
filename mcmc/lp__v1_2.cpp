// #include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins("cpp17")]] 

// The log of the normal density of x given mean mu and scale sd 
double normal_lpdf(double x, double mu, double sd ){
    return - .5*(M_LN2 + std::log(M_PI)) - std::log(sd) -.5*std::pow((x-mu)/sd,2);
}

// The log of the cumulative standard normal distribution of x 
double std_normal_lcdf(double x){
   return -M_LN2 + std::log(std::erfc(-x * M_SQRT1_2));
}

// The log of the gamma density of x given shape alpha and inverse scale beta
double gammapdf(double x, double alpha, double beta) {
    if(x > 0){
        return (std::pow(beta, alpha)*std::pow(x, (alpha-1))*std::pow(M_E, (-1*beta*x)))/std::tgamma(alpha);
    } else {
        return 0.0;
    }
}

// The log of the gamma density of x given shape alpha and inverse scale beta, extended to deal with big values of alpha
double log_gammapdf(double value, double alpha, double beta) {
    if(value > 0){

        if(alpha > 170){
            double aux_alpha = alpha;
            double log_gamma_alpha = 0.0;
            while(aux_alpha > 170){
                aux_alpha -= 1;
                log_gamma_alpha += std::log(aux_alpha); 
            }
            
            log_gamma_alpha += std::log(tgamma(aux_alpha));
            
            return alpha*std::log(beta) + (alpha-1)*std::log(value) - beta*value - log_gamma_alpha;
        } else{
            return alpha*std::log(beta) + (alpha-1)*std::log(value) - beta*value - std::log(tgamma(alpha));
        }
    } else {
        return 0.0;
    }
}

// The beta pdf
double betapdf(double x, double alpha, double beta) {
    if(x > 0 && x < 1){
        return (std::pow(x, (alpha-1))*std::pow(1-x, (beta-1))*std::tgamma(alpha+beta))/(std::tgamma(alpha)*std::tgamma(beta));
    } else {
        return 0.0;
    }
}

// The Zero-Inflated Poisson mass function of x given the rate lambda and probability theta (depends on Rcpp)
double dzipois(int x, double lambda, double theta, bool toLog = true){

    if(x == 0){
        NumericVector aux_vec = {log(theta),log1p(-theta) + R::dpois(x, lambda, toLog)};
        return max(aux_vec) + log(sum(exp(aux_vec - max(aux_vec))));
    } else {
        return log1p(-theta) + R::dpois(x, lambda, toLog);
    }
 
}

// The log of the Zero-Inflated Poisson mass function of x given the rate lambda and probability theta
double log_zipoispdf(int x, double lambda, double theta){

    if(x == 0){
        arma::Col<double> aux = {std::log(theta),std::log1p(-theta) +  x*std::log(lambda) -lambda - log(std::tgamma(x+1))};
        double aux_max = *std::max_element(std::begin(aux),std::end(aux));
        double sum_exp_aux_max = 0.0;
        for(int k = 0; k < 2; k++){
            aux[k] -= aux_max;
            sum_exp_aux_max += std::exp(aux[k]);
        }
         return aux_max + std::log(sum_exp_aux_max);
    } else {
        return std::log1p(-theta) + x*std::log(lambda) -lambda - log(std::tgamma(x+1));
    } 

}


// The log of the Von-misses density function of x given location mu and scale kappa
double log_dvm(double x, double mu, double kappa){
    return kappa*std::cos(x-mu) -(M_LN2 + log(M_PI) + log(std::cyl_bessel_i(0,kappa)));
}

// The log of the Dirichlet density for simplex x given vector of parameters alpha
double log_ddirichlet(const arma::Col<double> x, const arma::Col<double> alpha){

    // double sum_x = 0.0;
    // for(size_t j=0; j < x.size(); j++){
    //     sum_x += x[j];
    // }
    // if(sum_x == 1){
    //     double log_pdf = 0.0;
    //     double sum_alpha = 0.0;
    //     for(size_t i=0; i < x.size(); i++){
    //         log_pdf += (alpha[i]-1)*log(x[i]) - log(tgamma(alpha[i]));
    //         sum_alpha += alpha[i];
    //     }
        
    //     return log_pdf + log(tgamma(sum_alpha));

    // } else{
    //     return std::numeric_limits<double>::quiet_NaN();
    // }
    double log_pdf = 0.0;
    double sum_alpha = 0.0;
    for(size_t i=0; i < x.size(); i++){
        log_pdf += (alpha[i]-1)*log(x[i]) - log(tgamma(alpha[i]));
        sum_alpha += alpha[i];
    }
    
    return log_pdf + log(tgamma(sum_alpha));

}

// The log-likelihood of theta_star given the observed data (depends on Rcpp and R APIs)
double log_lik(const List list_data, NumericVector theta_star){

    // Extraction and storage of data contained in the list list_data  
    const int N = list_data["N"]; // Number of hidden states
    const int n = list_data["n"]; // Total number of observations excluding first observation of each individual
    const int n_ind = list_data["n_ind"]; // Total number of observartions considered as the first one for each individual
    // (Total number of observations = n_ind + n)
    IntegerVector ID_init = list_data["ID_init"]; // IDs corresponding to the first observation vectors  
    const IntegerVector ID = list_data["ID"]; // IDs corresponding to the observation vectors excluding first observations
    const NumericVector x_duration_init = list_data["x_duration_init"]; // first observation for the data stream "duration" of all individuals
    const NumericVector x_surface_init = list_data["x_surface_init"]; // first observation for the data stream "surface" of all individuals
    const NumericVector x_maxDepth_init = list_data["x_maxDepth_init"]; // first observation for the data stream "maximum depth" of all individuals
    const IntegerVector x_lunges_init = list_data["x_lunges_init"]; // first observation for the data stream "number of lunges" of all individuals
    const NumericVector x_step_init = list_data["x_step_init"]; // first observation for the data stream "number of steps" of all individuals
    const NumericVector x_angle_init = list_data["x_angle_init"]; // first observation for the data stream "turning angle" of all individuals
    const NumericVector x_headVar_init = list_data["x_headVar_init"]; // first observation for the data stream "heading variance" of all individuals
    const NumericVector x_duration = list_data["x_duration"]; // observations for the data stream "duration" of all individuals excluding the first observation of the sequence
    const NumericVector x_surface = list_data["x_surface"]; // observations for the data stream "surface" of all individuals excluding the first observation of the sequence
    const NumericVector x_maxDepth = list_data["x_maxDepth"]; // observations for the data stream "maximum depth" of all individuals excluding the first observation of the sequence
    const IntegerVector x_lunges = list_data["x_lunges"]; // observations for the data stream "number of lunges" of all individuals excluding the first observation of the sequence
    const NumericVector x_step = list_data["x_step"]; // observations for the data stream "number of steps" of all individuals excluding the first observation of the sequence
    const NumericVector x_angle = list_data["x_angle"]; // observations for the data stream "turning angle" of all individuals excluding the first observation of the sequence
    const NumericVector x_headVar = list_data["x_headVar"]; // observations for the data stream "heading variance" of all individuals excluding the first observation of the sequence
      
    // Initialization of vector parameters corresponding to the state-dependent distributions (the values are provided by theta_star)
    // Note: Some of the parameters has to be transformed since these doesn't match the input parameters of the density functions
    // (namely, the shape and scale parameters of the gamma density).
    NumericVector mu_duration(N);
    NumericVector sigma_duration(N);
    NumericVector mu_surface(N);
    NumericVector sigma_surface(N);
    NumericVector mu_maxDepth(N);
    NumericVector sigma_maxDepth(N);
    NumericVector mu_step(N);
    NumericVector sigma_step(N);
    NumericVector kappa(N);
    NumericVector a(N);
    NumericVector b(N);
    NumericVector lambda(N);
    NumericVector theta(N);
    
    // Setting the state-dependent distributions parameter values given the vector theta_star 
    for(int i = 0; i<3; ++i){
    mu_duration[i] = theta_star[i+6];
    sigma_duration[i] = theta_star[i+9];

    mu_surface[i] = theta_star[i+12];
    sigma_surface[i] = theta_star[i+15];

    mu_maxDepth[i] = theta_star[i+18];
    sigma_maxDepth[i] = theta_star[i+21];

    mu_step[i] = theta_star[i+24];
    sigma_step[i] = theta_star[i+27];

    kappa[i] = theta_star[i+30];

    a[i] = theta_star[i+33];
    b[i] = theta_star[i+36];

    lambda[i] = theta_star[i + 39];
    theta[i] = theta_star[i + 48];

    }

    // Initialization of the transformed vector parameters for the parameters associated to state-dependent gamma distributions
    // (data streams: duration, surface, maximum depth and number of steps)

    // parameters for duration distribution
    NumericVector alpha_duration(N);
    // NumericVector beta_duration(N);
    NumericVector theta_duration(N);

    // parameters for surface distribution
    NumericVector alpha_surface(N);
    // NumericVector beta_surface(N);
    NumericVector theta_surface(N);
    
    // parameters for maxDepth distribution
    NumericVector alpha_maxDepth(N);
    // NumericVector beta_maxDepth(N);
    NumericVector theta_maxDepth(N);

    // parameters for step distribution
    NumericVector alpha_step(N);
    // NumericVector beta_step(N);
    NumericVector theta_step(N);

    for (int i = 0; i < N; i++){
        alpha_duration[i] = pow(mu_duration[i],2) / pow(sigma_duration[i],2);
        // beta_duration[i] = mu_duration[i] / sigma_duration[i];
        theta_duration[i] =  pow(sigma_duration[i],2)/mu_duration[i];

        alpha_surface[i] = pow(mu_surface[i]/sigma_surface[i],2);
        // beta_surface[i] = mu_surface[i] / sigma_surface[i];
        theta_surface[i] =  pow(sigma_surface[i],2)/mu_surface[i];

        alpha_maxDepth[i] = pow(mu_maxDepth[i]/sigma_maxDepth[i],2);
        // beta_maxDepth[i] = mu_maxDepth[i] / sigma_maxDepth[i];
        theta_maxDepth[i] =  pow(sigma_maxDepth[i],2)/mu_maxDepth[i];
        
        alpha_step[i] = pow(mu_step[i]/sigma_step[i],2);
        // beta_step[i] = mu_step[i] / sigma_step[i];
        theta_step[i] =  pow(sigma_step[i],2)/mu_step[i];
        
    }

    // Initialization of the initial state distribution
    // initial distribution
    NumericVector init(N);
    init[2-1] = theta_star[42];
    init[3-1] = theta_star[43];

    // Initialization of the transtition probability matrix
    // tpm 
    NumericMatrix tpm(N,N);

    // creation tpm
    tpm(1-1,1-1) = theta_star[44];
    tpm(1-1,2-1) = theta_star[1-1];
    tpm(1-1,3-1) = theta_star[2-1];

    tpm(2-1,1-1) = theta_star[3-1];
    tpm(2-1,2-1) = theta_star[45];
    tpm(2-1,3-1) = theta_star[4-1];

    tpm(3-1,1-1) = theta_star[5-1];
    tpm(3-1,2-1) = theta_star[6-1];
    tpm(3-1,3-1) = theta_star[46];

    // Here we set the first entry of the initial state distribution vector
    init[1-1] = theta_star[47];
  
    // Initialization of the double variable where the log-likelihood will be stored
    double log_lik_total = 0.0;

    // Here the CircStats package is called in order to use the dvm function (von-misses density)
    // Note: the log of the Von-Misses density has been defined as a function in this file, so the code below can be omitted
    // and use our own version instead  
    Environment pkg = Environment::namespace_env("CircStats");

    Function dvm = pkg["dvm"];

    // Initialization of the matrix log_alpha. In this matrix each column i will correspond to the log-forward probabilities 
    // at time T_ind for individual i (which will end up being the log-likelihood component corresponding to individual i).
    NumericMatrix log_alpha(N,n_ind);
    
    // First we store the log-forward probabilities at time t=1 for each of the observed individuals.
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

            log_alpha(state,*i-1) += dzipois(x_lunges_init[*i-1] , 
                                            lambda[state],
                                            theta[state],
                                            true);
            
            // The data streams turning angle, number of steps and heading variance have missing values (NA's replaced with value 999.0), 
            // so an if statement is incorporated to deal with this
            if(x_angle_init[*i-1] != 999.0){
                NumericVector aux_dvm = dvm(x_angle_init[*i-1] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
                log_alpha(state,*i-1) += log(*(aux_dvm.begin()));

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
    
    // Here we compute the log-forward probablities for times t > 1 until we reach the observed time T_i corresponding to individual i.
    for(int i = 0; i < n; i++){
        
        // Auxiliary vector in which the log-forward probabilities of individual ID[i] at time t+1 (given log-forward probs in 
        // log_alpha(_,ID[i]-1) correspond to time t) will be stored
        // Note: since C++ starts indexing from 0, log_alpha(_,ID[i]-1) refers to individual ID[i] 
        NumericVector aux_log_alpha(N);
        
        for(int state = 0; state < N; state++){
        
        // Auxiliary vector which is initialized as the log-forward probabilities of individual (ID[i]-1) at time t. 
        NumericVector aux(N);
        aux = log_alpha(_,ID[i]-1);
        
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

            aux[previous_state] += dzipois(x_lunges[i],
                                            lambda[state],
                                            theta[state],
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
        
        // log-forward probabilities of individual (ID[i]-1) at time t+1
        aux_log_alpha[state] = max(aux) + log(sum(exp(aux-max(aux))));
        
        }
        
        // Here we update the log-forward probabilities for the corresponding individual
        log_alpha(_,ID[i]-1) = aux_log_alpha;
    }

    // We compute the log-likelihood using the log-sum-exp method for each one of the likelihood components
    for(int i=0; i < n_ind; i++){
    log_lik_total +=  max(log_alpha(_,i)) + log(sum(exp(log_alpha(_,i)-max(log_alpha(_,i)))));
    }

    return log_lik_total;
}

// (Unnormalilzed) log-posterior of theta_star given data (depends on Rcpp and R APIs)
double lp__(const List list_data,NumericVector theta_star){
    
    // Number of hidden states
    const int N = list_data["N"];

    NumericVector mu_duration(N);
    NumericVector sigma_duration(N);
    NumericVector mu_surface(N);
    NumericVector sigma_surface(N);
    NumericVector mu_maxDepth(N);
    NumericVector sigma_maxDepth(N);
    NumericVector mu_step(N);
    NumericVector sigma_step(N);
    NumericVector kappa(N);
    NumericVector a(N);
    NumericVector b(N);
    NumericVector lambda(N);
    NumericVector theta(N); // Since theta_i ~ U(0,1), there's no need to write down the prior for this variable

    for(int i = 0; i<N; i++){
        
        mu_duration[i] = theta_star[i+6];
        sigma_duration[i] = theta_star[i+9];

        mu_surface[i] = theta_star[i+12];
        sigma_surface[i] = theta_star[i+15];

        mu_maxDepth[i] = theta_star[i+18];
        sigma_maxDepth[i] = theta_star[i+21];

        mu_step[i] = theta_star[i+24];
        sigma_step[i] = theta_star[i+27];

        kappa[i] = theta_star[i+30];

        a[i] = theta_star[i+33];
        b[i] = theta_star[i+36];

        lambda[i] = theta_star[i + 39];
        theta[i] = theta_star[i + 48];
        
    }

    // Hyperparameters for the prior distributions

    // duration
    const NumericVector mu_duration_mean = list_data["mu_duration_mean"];
    const NumericVector mu_duration_sigma = list_data["mu_duration_sigma"];
    const NumericVector sigma_duration_alpha = list_data["sigma_duration_alpha"];
    const NumericVector sigma_duration_beta = list_data["sigma_duration_beta"];
    
    // surface
    const NumericVector mu_surface_mean = list_data["mu_surface_mean"];
    const NumericVector mu_surface_sigma = list_data["mu_surface_sigma"];
    const NumericVector sigma_surface_alpha = list_data["sigma_surface_alpha"];
    const NumericVector sigma_surface_beta = list_data["sigma_surface_beta"];

    // maxDepth
    const NumericVector mu_maxDepth_mean = list_data["mu_maxDepth_mean"];
    const NumericVector mu_maxDepth_sigma = list_data["mu_maxDepth_sigma"];
    const NumericVector sigma_maxDepth_alpha = list_data["sigma_maxDepth_alpha"];
    const NumericVector sigma_maxDepth_beta = list_data["sigma_maxDepth_beta"];
    
    // step
    const NumericVector mu_step_mean = list_data["mu_step_mean"];
    const NumericVector mu_step_sigma = list_data["mu_step_sigma"];
    const NumericVector sigma_step_alpha = list_data["sigma_step_alpha"];
    const NumericVector sigma_step_beta = list_data["sigma_step_beta"];
    
    // angle
    const NumericVector kappa_alpha = list_data["kappa_alpha"];
    const NumericVector kappa_beta = list_data["kappa_beta"];
    
    // heading variance
    const NumericVector a_alpha = list_data["a_alpha"];
    const NumericVector a_beta = list_data["a_beta"];
    const NumericVector b_alpha = list_data["b_alpha"];
    const NumericVector b_beta = list_data["b_beta"];
    
    // number of lunges
    const NumericVector lambda_alpha = list_data["lambda_alpha"];
    const NumericVector lambda_beta = list_data["lambda_beta"];    
    
    // initial distribution
    NumericVector init(N);
    init[2-1] = theta_star[42];
    init[3-1] = theta_star[43];
  
    // tpm 
    NumericMatrix tpm(N,N);
  
    // creation tpm
    tpm(1-1,1-1) = theta_star[44];
    tpm(1-1,2-1) = theta_star[1-1];
    tpm(1-1,3-1) = theta_star[2-1];
    
    tpm(2-1,1-1) = theta_star[3-1];
    tpm(2-1,2-1) = theta_star[45];
    tpm(2-1,3-1) = theta_star[4-1];    
    
    tpm(3-1,1-1) = theta_star[5-1];
    tpm(3-1,2-1) = theta_star[6-1];
    tpm(3-1,3-1) = theta_star[46];

    // Initialization of the rows of the transition probability matrix (easier to handle in this way)
    NumericVector tpm1 = tpm(0,_);
    NumericVector tpm2 = tpm(1,_);
    NumericVector tpm3 = tpm(2,_);

    init[1-1] = theta_star[47];

    // Here the gtools package is called in order to use the ddirichlet function (dirichlet density function)
    // Note: the log of the Dirichlet density has been defined as a function in this file, so the code below can be omitted
    // and use our own version instead  
    Environment pkg = Environment::namespace_env("gtools");

    Function ddirichlet = pkg["ddirichlet"];

    // Initialization of the double variable where the unnormalized log-posterior will be stored
    double target = 0.0;
    
    for(int i = 0; i < N; i++){
    
        // Prior distributions associated to the duration parameters

        // mu_duration ~ truncated from below at zero normal distribution
        target+= R::dnorm(mu_duration[i] , mu_duration_mean[i] , mu_duration_sigma[i],true) - 
                                                                            R::pnorm(-mu_duration_mean[i]/mu_duration_sigma[i],
                                                                            0,1, true,
                                                                            true);
        // sigma_duration ~ gamma distribution
        target+= R::dgamma(sigma_duration[i] , sigma_duration_alpha[i] , 1/sigma_duration_beta[i],true);
    
        // Prior distributions associated to the surface parameters

        // mu_surface ~ truncated from below at zero normal distribution
        target+= R::dnorm(mu_surface[i] , mu_surface_mean[i] , mu_surface_sigma[i], true) - 
                                                                            R::pnorm(-mu_surface_mean[i]/mu_surface_sigma[i],
                                                                            0,1, true,
                                                                            true);
        // sigma_surface ~ gamma distribution
        target+= R::dgamma(sigma_surface[i] , sigma_surface_alpha[i] , 1/sigma_surface_beta[i], true); 
    
        // Prior distributions associated to the maximum depth parameters

        // mu_maxDepth ~ truncated from below at zero normal distribution
        target+= R::dnorm(mu_maxDepth[i] , mu_maxDepth_mean[i] , mu_maxDepth_sigma[i], true) -
                                                                            R::pnorm(-mu_maxDepth_mean[i]/mu_maxDepth_sigma[i],
                                                                            0,1, true,
                                                                            true);
        // sigma_maxDepth ~ gamma distribution
        target+= R::dgamma(sigma_maxDepth[i] , sigma_maxDepth_alpha[i] , 1/sigma_maxDepth_beta[i], true );
    
        // Prior distributions associated to the number of steps parameters

        // mu_step ~ truncated from below at zero normal distribution
        target+= R::dnorm(mu_step[i] , mu_step_mean[i] , mu_step_sigma[i], true) -
                                                                            R::pnorm(-mu_step_mean[i]/mu_step_sigma[i],
                                                                            0,1, true,
                                                                            true);
        // sigma_step ~ gamma distribution
        target+= R::dgamma(sigma_step[i] , sigma_step_alpha[i] , 1/sigma_step_beta[i], true);

        // Prior distributions associated to the turning angle parameters

        // kappa ~ gamma distribution
        target+= R::dgamma(kappa[i] , kappa_alpha[i], 1/kappa_beta[i], true);

        // Prior distributions associated to the heading variance parameters

        // a ~ gamma distribution
        target+= R::dgamma(a[i] , a_alpha[i], 1/a_beta[i],true);
        // b ~ gamma distribution
        target+= R::dgamma(b[i] , b_alpha[i], 1/b_beta[i],true);

        // Prior distributions associated to the number of lunges parameters
        
        // lambda ~ gamma distribution 
        target+= R::dgamma(lambda[i] , lambda_alpha[i], 1/lambda_beta[i],true);
        // For the probability theta, we assumed theta ~ U(0,1), so there's no need to write down code to include this

    }
  
    // Prior distributions associated to the initial state distribution and transition probability matrix rows
    // initial state distribution ~ dirichlet
    NumericVector aux_init = ddirichlet(init , _["alpha"] = rep(1,N)); 
    target += log(*(aux_init.begin()));


    // rows tpm ~ dirichlet
    NumericVector aux_tpm1 = ddirichlet(tpm1 , _["alpha"] = rep(1,N)); 
    target += log(*(aux_tpm1.begin()));
    NumericVector aux_tpm2 = ddirichlet(tpm2 , _["alpha"] = rep(1,N)); 
    target += log(*(aux_tpm2.begin()));
    NumericVector aux_tpm3 = ddirichlet(tpm3 , _["alpha"] = rep(1,N)); 
    target += log(*(aux_tpm3.begin()));

    // Here we incorporate the likelihood
    target += log_lik(list_data, theta_star);
    return target;   

}

// Function that allows to sample uniformly from the simples
NumericVector next_point(const NumericVector x_old){
        
    NumericVector x_new(x_old.size());

    for(int i=0; i< x_old.size(); i++){
        x_new[i] =  -std::log(R::runif(0,1));    
    }
    double norm_const = sum(x_new);
    return x_new/norm_const;
}

// Generates a new proposal for each one of the parameters of interest and return all of them into a single vector (depends on Rcpp and R APIs)
NumericVector proposal_dist(const NumericVector theta_star, const int N = 3){
    
    // Initialization of vector parameters corresponding to the state-dependent distributions (the values are provided by theta_star)
    NumericVector mu_duration(N);
    NumericVector sigma_duration(N);
    NumericVector mu_surface(N);
    NumericVector sigma_surface(N);
    NumericVector mu_maxDepth(N);
    NumericVector sigma_maxDepth(N);
    NumericVector mu_step(N);
    NumericVector sigma_step(N);
    NumericVector kappa(N);
    NumericVector a(N);
    NumericVector b(N);
    NumericVector lambda(N);
    NumericVector theta(N);

    for(int i = 0; i<N; i++){
        
        mu_duration[i] = theta_star[i+6];
        sigma_duration[i] = theta_star[i+9];

        mu_surface[i] = theta_star[i+12];
        sigma_surface[i] = theta_star[i+15];

        mu_maxDepth[i] = theta_star[i+18];
        sigma_maxDepth[i] = theta_star[i+21];

        mu_step[i] = theta_star[i+24];
        sigma_step[i] = theta_star[i+27];

        kappa[i] = theta_star[i+30];

        a[i] = theta_star[i+33];
        b[i] = theta_star[i+36];

        lambda[i] = theta_star[i+39];
        theta[i] = theta_star[i+48];
        
    }

    // Initialization of the initial state distribution using values from theta_star
    NumericVector init(N);
    init[2-1] = theta_star[42];
    init[3-1] = theta_star[43];
  
    // Initialization of tpm 
    NumericMatrix tpm(N,N);
  
    // Defining entries of tpm from theta_star
    tpm(1-1,1-1) = theta_star[44];
    tpm(1-1,2-1) = theta_star[1-1];
    tpm(1-1,3-1) = theta_star[2-1];
    
    tpm(2-1,1-1) = theta_star[3-1];
    tpm(2-1,2-1) = theta_star[45];
    tpm(2-1,3-1) = theta_star[4-1];    
    
    tpm(3-1,1-1) = theta_star[5-1];
    tpm(3-1,2-1) = theta_star[6-1];
    tpm(3-1,3-1) = theta_star[46];

    // Initialization of the tpm rows (this is easier to handle)
    NumericVector tpm1 = tpm(0,_);
    NumericVector tpm2 = tpm(1,_);
    NumericVector tpm3 = tpm(2,_);

    init[1-1] = theta_star[47];

    // Initialization of vector where new proposals will be allocated (state-dependents parameters)
    NumericVector new_mu_duration(N);
    NumericVector new_sigma_duration(N);
    NumericVector new_mu_surface(N);
    NumericVector new_sigma_surface(N);
    NumericVector new_mu_maxDepth(N);
    NumericVector new_sigma_maxDepth(N);
    NumericVector new_mu_step(N);
    NumericVector new_sigma_step(N);
    NumericVector new_kappa(N);
    NumericVector new_a(N);
    NumericVector new_b(N);
    NumericVector new_lambda(N);
    NumericVector new_theta(N);

    // Initialization of vectors where new proposal for the initial state dist and the tpms rows will be allocated
    NumericVector new_init(N);
    NumericVector new_tpm1(N);
    NumericVector new_tpm2(N);
    NumericVector new_tpm3(N);

    // step size for the Metropolis algorithm for each one of the state-dependent distribution parameters
    // NumericVector step_size_mu_duration = {2,5,4};
    NumericVector step_size_mu_duration = {.008,.008,.008};
    // NumericVector step_size_sigma_duration = {4,6,4};
    // NumericVector step_size_sigma_duration = {.008,.008,.008};
    NumericVector step_size_sigma_duration = {.01,.01,.01};
    // NumericVector step_size_mu_surface = {2,1.5,1};
    // NumericVector step_size_mu_surface = {.008,.008,.008};
    NumericVector step_size_mu_surface = {.01,.01,.01};
    // NumericVector step_size_sigma_surface = {3,1.5,1.5};
    // NumericVector step_size_sigma_surface = {.01,.01,.01};
    NumericVector step_size_sigma_surface = {.05,.05,.05};
    // NumericVector step_size_mu_maxDepth = {.4,1,1};
    // NumericVector step_size_mu_maxDepth = {.008,.008,.008};
    NumericVector step_size_mu_maxDepth = {.01,.01,.01};
    // NumericVector step_size_sigma_maxDepth = {.5,1,1};
    NumericVector step_size_sigma_maxDepth = {.01,.01,.01};
    // NumericVector step_size_mu_step = {2,4,3};
    // NumericVector step_size_mu_step = {.008,.008,.008};
    NumericVector step_size_mu_step = {.01,.01,.01};
    // NumericVector step_size_sigma_step = {2.5,3,2};
    // NumericVector step_size_sigma_step = {.008,.008,.008};
    NumericVector step_size_sigma_step = {.01,.01,.01};
    // NumericVector step_size_kappa = {.1,.1,.2}; //first attempt, with no parallelization included
    // NumericVector step_size_kappa = {.05,.1,.05}; // second attempt, when parallelizing, but was moving toward negative values
    NumericVector step_size_kappa = {.08,.08,.08}; // third attempt, now working on the log scale
    // NumericVector step_size_a = {.05,.1,.05}; // first attempt, random walk was moving towards negative values
    NumericVector step_size_a = {.05,.05,.05}; // second attempt, now working on the log scale
    // NumericVector step_size_b = {.1,.1,.05}; // first attempt, didn't move to negative values, but just in case let's work on the log scale 
    NumericVector step_size_b = {.05,.05,.05};
    // NumericVector step_size_lambda = {.05,.1,.1}; // first attempt, random walk was moving towards negative values
    NumericVector step_size_lambda = {.05,.05,.05}; // second attempt, now working on the log scale
    // double step_size_theta = {.5,.5,.5}; // for theta, we'll take U(0,1) as the proposal distribution
    
    // double MH = 0.0;

    for(int i=0; i < N; i++){

        // new_mu_duration[i] = R::rnorm(mu_duration[i],step_size_mu_duration[i]);
        new_mu_duration[i] = std::exp(std::log(mu_duration[i]) + R::rnorm(0.0,step_size_mu_duration[i]));
        // new_sigma_duration[i] = R::rnorm(sigma_duration[i],step_size_sigma_duration[i]);
        new_sigma_duration[i] = std::exp(std::log(sigma_duration[i]) + R::rnorm(0.0,step_size_sigma_duration[i]));
        // new_mu_surface[i] = R::rnorm(mu_surface[i],step_size_mu_surface[i]);
        new_mu_surface[i] = std::exp(std::log(mu_surface[i]) + R::rnorm(0.0,step_size_mu_surface[i]));
        // new_sigma_surface[i] = R::rnorm(sigma_surface[i],step_size_sigma_surface[i]);
        new_sigma_surface[i] = std::exp(std::log(sigma_surface[i]) + R::rnorm(0.0,step_size_sigma_surface[i]));
        // new_mu_maxDepth[i] = R::rnorm(mu_maxDepth[i],step_size_mu_maxDepth[i]);
        new_mu_maxDepth[i] = std::exp(std::log(mu_maxDepth[i]) + R::rnorm(0.0,step_size_mu_maxDepth[i]));
        // new_sigma_maxDepth[i] = R::rnorm(sigma_maxDepth[i],step_size_sigma_maxDepth[i]);
        new_sigma_maxDepth[i] = std::exp(std::log(sigma_maxDepth[i]) + R::rnorm(0.0,step_size_sigma_maxDepth[i]));
        // new_mu_step[i] = R::rnorm(mu_step[i],step_size_mu_step[i]);
        new_mu_step[i] = std::exp(std::log(mu_step[i]) + R::rnorm(0.0,step_size_mu_step[i]));
        // new_sigma_step[i] = R::rnorm(sigma_step[i],step_size_sigma_step[i]);
        new_sigma_step[i] = std::exp(std::log(sigma_step[i]) + R::rnorm(0.0,step_size_sigma_step[i]));
        // new_kappa[i] = R::rnorm(kappa[i],step_size_kappa[i]); // first attempt, random walk was moving towards negative values  
        new_kappa[i] = std::exp(std::log(kappa[i]) + R::rnorm(0.0,step_size_kappa[i])); // second attempt, now working on the log scale
        // new_a[i] = R::rnorm(a[i],step_size_a[i]); // first attempt, random walk was moving towards negative values
        new_a[i] = std::exp(std::log(a[i]) + R::rnorm(0.0,step_size_a[i])); // second attempt, now working on the log scale
        // new_b[i] = R::rnorm(b[i],step_size_b[i]); // first attempt, random walk was moving towards negative values
        new_b[i] = std::exp(std::log(b[i]) + R::rnorm(0.0,step_size_b[i])); // second attempt, now working on the log scale
        // new_lambda[i] = R::rnorm(lambda[i],step_size_lambda[i]); // first attempt, random walk was moving towards negative values
        new_lambda[i] = std::exp(std::log(lambda[i]) + R::rnorm(0.0,step_size_lambda[i])); // second attempt, now working on the log scale
        new_theta[i] = R::runif(0,1);

        // // Hasting ratio
        // MH += - std::log(mu_duration[i]) - (-std::log(new_mu_duration[i]));
        // MH += - std::log(sigma_duration[i]) - (-std::log(new_sigma_duration[i]));

        // MH += - std::log(mu_surface[i]) - (-std::log(new_mu_surface[i]));
        // MH += - std::log(sigma_surface[i]) - (-std::log(new_sigma_surface[i]));
        
        // MH += - std::log(mu_maxDepth[i]) - (-std::log(new_mu_maxDepth[i]));
        // MH += - std::log(sigma_maxDepth[i]) - (-std::log(new_sigma_maxDepth[i]));
        
        // MH += - std::log(mu_step[i]) - (-std::log(new_mu_step[i]));
        // MH += - std::log(sigma_step[i]) - (-std::log(new_sigma_step[i]));

        // MH += -std::log(kappa[i]) - (-std::log(new_kappa[i]));

        // MH += -std::log(a[i]) - (-std::log(new_a[i]));
        // MH += -std::log(b[i]) - (-std::log(new_b[i]));

        // MH += -std::log(lambda[i]) - (-std::log(new_lambda[i]));
        
    }

    NumericVector ordered_new_mu_duration = new_mu_duration.sort();
    NumericVector ordered_new_mu_maxDepth = new_mu_maxDepth.sort();
    NumericVector ordered_new_mu_step = new_mu_step.sort();

    // Proposals for initial state distribution and tpm rows
    new_init = next_point(init);
    new_tpm1 = next_point(tpm1);
    new_tpm2 = next_point(tpm2);
    new_tpm3 = next_point(tpm3);

    // In this vector all the new proposals for the parameters will be stored
    NumericVector new_theta_star(theta_star.length());

    for(int i=0; i<N; i++){

        // new_theta_star[i+6] = new_mu_duration[i];
        new_theta_star[i+6] = ordered_new_mu_duration[i];
        new_theta_star[i+9] = new_sigma_duration[i];
        new_theta_star[i+12] = new_mu_surface[i];
        new_theta_star[i+15] = new_sigma_surface[i];
        // new_theta_star[i+18] = new_mu_maxDepth[i];
        new_theta_star[i+18] = ordered_new_mu_maxDepth[i];
        new_theta_star[i+21] = new_sigma_maxDepth[i];
        // new_theta_star[i+24] = new_mu_step[i];
        new_theta_star[i+24] = ordered_new_mu_step[i];
        new_theta_star[i+27] = new_sigma_step[i];
        new_theta_star[i+30] = new_kappa[i];
        new_theta_star[i+33] = new_a[i];
        new_theta_star[i+36] = new_b[i];
        new_theta_star[i+39] = new_lambda[i];
        new_theta_star[i+48] = new_theta[i];

    }

    // proposed tpm values
    new_theta_star[44] = new_tpm1[1-1];
    new_theta_star[1-1] = new_tpm1[2-1];
    new_theta_star[2-1] = new_tpm1[3-1];
    
    new_theta_star[3-1] = new_tpm2[1-1];
    new_theta_star[45] = new_tpm2[2-1];
    new_theta_star[4-1] = new_tpm2[3-1];
    
    new_theta_star[5-1] = new_tpm3[1-1];
    new_theta_star[6-1] = new_tpm3[2-1];
    new_theta_star[46] = new_tpm3[3-1];

    // proposed init values
    new_theta_star[42] = new_init[2-1];
    new_theta_star[43] = new_init[3-1];
    new_theta_star[47] = new_init[1-1];
    
    return new_theta_star;
}

// Generation of the (log and unnormalized) tempered targets (depends on Rcpp and R APIs)
double target_tau(List list_data,NumericVector theta, double thetemp , double maxtemp = 10) {
    if( (thetemp < 1) || (thetemp > maxtemp) ){
        return 0.0;
    } 
    else{
        return lp__(list_data, theta)/thetemp;
    } 
}

// Auxiliary function to create matrices (this is used to allocate matrices in a list in the function pt_cw_M_target_posterior)
NumericMatrix createMatrix(const int nrows, const int ncols){
    NumericMatrix aux_matrix(nrows,ncols);
    return aux_matrix;
}

// Parallel tempering using component-wise Metropolis algorithm (depends on Rcpp and R APIs)
List pt_cw_M_target_posterior(const List list_data, 
                        const int nsim,  
                        NumericVector init,
                        NumericVector temp_vector,
                        bool display_progress=true){
    
    // This is to display a progress bar and see how long the function it's taking to make samples (helpful for debugging and track errors)
    Progress p(nsim*temp_vector.size()*(init.size()-8), display_progress);

    // List where each element will allocate samples from the different tempered versions of the target
    List X(temp_vector.size());

    // Initial values provided for each tempered target
    for(int i = 0; i < temp_vector.size(); ++i){
        X[i] = createMatrix(init.size(),nsim+1);
        NumericMatrix aux_X = X[i];
        aux_X(_,0) = init;
    }
  
    // We will generate nsim simulations for each chain
    for(int i=0; i<nsim; ++i){
        // Part of the progress bar - it's supposed to retun -1.0 if the process is aborted
        if (Progress::check_abort() )
            return -1.0;
        
        // This loop corresponds to the component-wise step for each chain
        for(int temp = 0; temp < temp_vector.size(); ++temp){
            // For the chain temp we select the its corresponding matrix where samples are allocated 
            NumericMatrix aux_X = X[temp];
            // Generation of new proposal for the parameters
            NumericVector aux_Y_n = proposal_dist(aux_X(_,i));
            // Auxiliary vector for the component-wise step for each parameter
            NumericVector Y_n(init.length());
            
            for(int id = 0; id <  init.length(); id++){
                Y_n[id] = aux_X(id,i);
            }
            // tempered target at current iteration
            double target_at_X_i = target_tau(list_data, aux_X(_,i),temp_vector[temp],temp_vector[temp_vector.size()-1]);
            
            // component-wise steps for tpms and initial distribution
            //  component-wise step for tpm1 (tpm row 1)
            p.increment();
            double log_U_tpm1 = log(R::runif(0,1));

            Y_n[44] =  aux_Y_n[44];
            Y_n[1-1] =  aux_Y_n[1-1];
            Y_n[2-1] =  aux_Y_n[2-1];

            double log_alpha_tpm1 = 0.0;
            log_alpha_tpm1 += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            log_alpha_tpm1 -= target_at_X_i;
            double log_ratio_tpm1 = min(NumericVector::create(0,log_alpha_tpm1));
            
            if(log_U_tpm1 <= log_ratio_tpm1){
                aux_X(44,i+1) = aux_Y_n[44];
                aux_X(1-1,i+1) = aux_Y_n[1-1];
                aux_X(2-1,i+1) = aux_Y_n[2-1];
            } else {
                aux_X(44,i+1) = aux_X(44,i);
                aux_X(1-1,i+1) = aux_X(1-1,i);
                aux_X(2-1,i+1) = aux_X(2-1,i);
            }

            // return the value we had before to recycle this vector
            Y_n[44] =  aux_X(44,i);
            Y_n[1-1] =  aux_X(1-1,i);
            Y_n[2-1] =  aux_X(2-1,i);
            
            //  component-wise step for tpm2 (tpm row 2)
            p.increment();
            double log_U_tpm2 = log(R::runif(0,1));

            Y_n[3-1] =  aux_Y_n[3-1];
            Y_n[45] =  aux_Y_n[45];
            Y_n[4-1] =  aux_Y_n[4-1];

            double log_alpha_tpm2 = 0.0;
            log_alpha_tpm2 += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            log_alpha_tpm2 -= target_at_X_i;
            double log_ratio_tpm2 = min(NumericVector::create(0,log_alpha_tpm2));
            
            if(log_U_tpm2 <= log_ratio_tpm2){
                aux_X(3-1,i+1) = aux_Y_n[3-1];
                aux_X(45,i+1) = aux_Y_n[45];
                aux_X(4-1,i+1) = aux_Y_n[4-1];
            } else {
                aux_X(3-1,i+1) = aux_X(3-1,i);
                aux_X(45,i+1) = aux_X(45,i);
                aux_X(4-1,i+1) = aux_X(4-1,i);
            }

            // return the value we had before to recycle this vector
            Y_n[3-1] =  aux_X(3-1,i);
            Y_n[45] =  aux_X(45,i);
            Y_n[4-1] =  aux_X(4-1,i);

            //  component-wise step for tpm3 (tpm row 3)
            p.increment();
            double log_U_tpm3 = log(R::runif(0,1));

            Y_n[5-1] =  aux_Y_n[5-1];
            Y_n[6-1] =  aux_Y_n[6-1];
            Y_n[46] =  aux_Y_n[46];

            double log_alpha_tpm3 = 0.0;
            log_alpha_tpm3 += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            log_alpha_tpm3 -= target_at_X_i;
            double log_ratio_tpm3 = min(NumericVector::create(0,log_alpha_tpm3));
            
            if(log_U_tpm3 <= log_ratio_tpm3){
                aux_X(5-1,i+1) = aux_Y_n[5-1];
                aux_X(6-1,i+1) = aux_Y_n[6-1];
                aux_X(46,i+1) = aux_Y_n[46];
            } else {
                aux_X(5-1,i+1) = aux_X(5-1,i);
                aux_X(6-1,i+1) = aux_X(6-1,i);
                aux_X(46,i+1) = aux_X(46,i);
            }

            // return the value we had before to recycle this vector
            Y_n[5-1] =  aux_X(5-1,i);
            Y_n[6-1] =  aux_X(6-1,i);
            Y_n[46] =  aux_X(46,i);

            //  component-wise step for initial state distribution
            p.increment();
            double log_U_initd = log(R::runif(0,1));

            Y_n[42] =  aux_Y_n[42];
            Y_n[43] =  aux_Y_n[43];
            Y_n[47] = aux_Y_n[47];

            double log_alpha_initd = 0.0;
            log_alpha_initd += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            log_alpha_initd -= target_at_X_i;
            double log_ratio_initd = min(NumericVector::create(0,log_alpha_initd));
            
            if(log_U_initd <= log_ratio_initd){
                aux_X(42,i+1) = aux_Y_n[42];
                aux_X(43,i+1) = aux_Y_n[43];
                aux_X(47,i+1) = aux_Y_n[47];
            } else {
                aux_X(42,i+1) = aux_X(42,i);
                aux_X(43,i+1) = aux_X(43,i);
                aux_X(47,i+1) = aux_X(47,i);
            }

            // return the value we had before to recycle this vector
            Y_n[42] =  aux_X(42,i);
            Y_n[43] =  aux_X(43,i);
            Y_n[47] =  aux_X(47,i);

            // component-wise steps for the state-dependent parameters
            for(int j = 6; j < (init.length()-9); j++){
                p.increment();
                double log_U = log(R::runif(0,1));
                Y_n[j] =  aux_Y_n[j];
                double log_alpha = 0.0;
                log_alpha += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]); 
                log_alpha -= target_at_X_i;
                double log_ratio = min(NumericVector::create(0,log_alpha));
                
                if(log_U <= log_ratio){
                    aux_X(j,i+1) = aux_Y_n[j];
                } else {
                    aux_X(j,i+1) = aux_X(j,i);
                }

                Y_n[j] = aux_X(j,i); // return the value we had before to recycle this vector
            }

            // component-wise for theta (zero-inflated poisson)
            for(int j = init.length()-3; j < init.length(); j++){
                p.increment();
                double log_U = log(R::runif(0,1));
                Y_n[j] =  aux_Y_n[j];
                double log_alpha = 0.0;
                log_alpha += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]); 
                log_alpha -= target_at_X_i;
                double log_ratio = min(NumericVector::create(0,log_alpha));
                
                if(log_U <= log_ratio){
                    aux_X(j,i+1) = aux_Y_n[j];
                } else {
                    aux_X(j,i+1) = aux_X(j,i);
                }

                Y_n[j] = aux_X(j,i); // return the value we had before to recycle this vector
            }
        }

        // Choose randomly one of the chains
        int j = floor(R::runif(0,1)*(temp_vector.size()-1));
        int k = j + 1; // Proposed Swap;
        double log_U_swap = log(R::runif(0,1));
        // Auxiliary matrices where we will extract the values to be swapped
        NumericMatrix X_j = X[j];
        NumericMatrix X_k = X[k];
        
        double ratio = target_tau(list_data, X_j(_,i+1), k,temp_vector[temp_vector.size()-1]) + target_tau(list_data, X_k(_,i+1), j,temp_vector[temp_vector.size()-1]) - 
                        (target_tau(list_data, X_j(_,i+1), j,temp_vector[temp_vector.size()-1]) + target_tau(list_data, X_k(_,i+1), k,temp_vector[temp_vector.size()-1]));


        if(log_U_swap  < ratio){
            NumericVector tmpval = X_j(_,i);
            X_j(_,i) = X_k(_,i);
            X_k(_,i) = tmpval; // Accept Swap;
        }
    }

    return X;
} 

// NOTE: The following corresponds to the parallel version of the function above. For this the log_lik and lp__ functions
// had to be modified in order to not depend on Rcpp and R APIs (data types and functions).

// The log-likelihood of theta_star given the observed data
double log_lik_stdVector(const int N,  // Number of hidden states
                    const int n, // Total number of observations excluding first observation of each individual
                    const int n_ind, // Total number of observartions considered as the first one for each individual
                    // (Total number of observations = n_ind + n)
                    arma::Col<int> ID_init, // IDs corresponding to the first observation vectors
                    const arma::Col<int> ID, // IDs corresponding to the observation vectors excluding first observations 
                    const arma::Col<double> x_duration_init, // first observation for the data stream "duration" of all individuals
                    const arma::Col<double> x_surface_init, // first observation for the data stream "surface" of all individuals 
                    const arma::Col<double> x_maxDepth_init, // first observation for the data stream "maximum depth" of all individuals 
                    const arma::Col<int> x_lunges_init, // first observation for the data stream "number of lunges" of all individuals 
                    const arma::Col<double> x_step_init, // first observation for the data stream "number of steps" of all individuals 
                    const arma::Col<double> x_angle_init, // first observation for the data stream "turning angle" of all individuals 
                    const arma::Col<double> x_headVar_init, // first observation for the data stream "heading variance" of all individuals 
                    const arma::Col<double> x_duration, // observations for the data stream "duration" of all individuals excluding the first observation of the sequence 
                    const arma::Col<double> x_surface, // observations for the data stream "surface" of all individuals excluding the first observation of the sequence 
                    const arma::Col<double> x_maxDepth, // observations for the data stream "maximum depth" of all individuals excluding the first observation of the sequence 
                    const arma::Col<int> x_lunges, // observations for the data stream "number of lunges" of all individuals excluding the first observation of the sequence 
                    const arma::Col<double> x_step, // observations for the data stream "number of steps" of all individuals excluding the first observation of the sequence 
                    const arma::Col<double> x_angle, // observations for the data stream "turning angle" of all individuals excluding the first observation of the sequence 
                    const arma::Col<double> x_headVar, // observations for the data stream "heading variance" of all individuals excluding the first observation of the sequence 
                    arma::Col<double> theta_star){

    // Initialization of vector parameters corresponding to the state-dependent distributions (the values are provided by theta_star)
    // Note: Some of the parameters has to be transformed since these doesn't match the input parameters of the density functions
    // (namely, the shape and scale parameters of the gamma density).
    arma::Col<double> mu_duration(N);
    arma::Col<double> sigma_duration(N);
    arma::Col<double> mu_surface(N);
    arma::Col<double> sigma_surface(N);
    arma::Col<double> mu_maxDepth(N);
    arma::Col<double> sigma_maxDepth(N);
    arma::Col<double> mu_step(N);
    arma::Col<double> sigma_step(N);
    arma::Col<double> kappa(N);
    arma::Col<double> a(N);
    arma::Col<double> b(N);
    arma::Col<double> lambda(N);
    arma::Col<double> theta(N);

    // Setting the state-dependent distributions parameter values given the vector theta_star 
    for(int i = 0; i<3; ++i){
    mu_duration[i] = theta_star[i+6];
    sigma_duration[i] = theta_star[i+9];

    mu_surface[i] = theta_star[i+12];
    sigma_surface[i] = theta_star[i+15];

    mu_maxDepth[i] = theta_star[i+18];
    sigma_maxDepth[i] = theta_star[i+21];

    mu_step[i] = theta_star[i+24];
    sigma_step[i] = theta_star[i+27];

    kappa[i] = theta_star[i+30];

    a[i] = theta_star[i+33];
    b[i] = theta_star[i+36];

    lambda[i] = theta_star[i + 39];
    theta[i] = theta_star[i + 48];

    }

    // Initialization of the transformed vector parameters for the parameters associated to state-dependent gamma distributions
    // (data streams: duration, surface, maximum depth and number of steps)

    // parameters for duration distribution
    arma::Col<double> alpha_duration(N);
    arma::Col<double> beta_duration(N);
    // arma::Col<double> theta_duration(N);

    // parameters for surface distribution
    arma::Col<double> alpha_surface(N);
    arma::Col<double> beta_surface(N);
    // arma::Col<double> theta_surface(N);
    
    // parameters for maxDepth distribution
    arma::Col<double> alpha_maxDepth(N);
    arma::Col<double> beta_maxDepth(N);
    // arma::Col<double> theta_maxDepth(N);

    // parameters for step distribution
    arma::Col<double> alpha_step(N);
    arma::Col<double> beta_step(N);
    // arma::Col<double> theta_step(N);

    for (int i = 0; i < N; i++){
        alpha_duration[i] = std::pow(mu_duration[i] / sigma_duration[i],2);
        beta_duration[i] = mu_duration[i] / std::pow(sigma_duration[i],2);
        // theta_duration[i] =  std::pow(sigma_duration[i],2)/mu_duration[i];

        alpha_surface[i] = std::pow(mu_surface[i]/sigma_surface[i],2);
        beta_surface[i] = mu_surface[i] / std::pow(sigma_surface[i],2);
        // theta_surface[i] =  std::pow(sigma_surface[i],2)/mu_surface[i];

        alpha_maxDepth[i] = std::pow(mu_maxDepth[i]/sigma_maxDepth[i],2);
        beta_maxDepth[i] = mu_maxDepth[i] / std::pow(sigma_maxDepth[i],2);
        // theta_maxDepth[i] =  std::pow(sigma_maxDepth[i],2)/mu_maxDepth[i];
        
        alpha_step[i] = std::pow(mu_step[i]/sigma_step[i],2);
        beta_step[i] = mu_step[i] / std::pow(sigma_step[i],2);
        // theta_step[i] =  std::pow(sigma_step[i],2)/mu_step[i];
        
    }

    // Initialization of the initial state distribution
    // initial distribution
    arma::Col<double> init(N);
    init[2-1] = theta_star[42];
    init[3-1] = theta_star[43];

    // Initialization of the transtition probability matrix
    // tpm 
    arma::mat tpm(N,N);

    // creation tpm
    tpm(1-1,1-1) = theta_star[44];
    tpm(1-1,2-1) = theta_star[1-1];
    tpm(1-1,3-1) = theta_star[2-1];

    tpm(2-1,1-1) = theta_star[3-1];
    tpm(2-1,2-1) = theta_star[45];
    tpm(2-1,3-1) = theta_star[4-1];

    tpm(3-1,1-1) = theta_star[5-1];
    tpm(3-1,2-1) = theta_star[6-1];
    tpm(3-1,3-1) = theta_star[46];

    // Here we set the first entry of the initial state distribution vector
    init[1-1] = theta_star[47];

    // Initialization of the double variable where the log-likelihood will be stored
    double log_lik_total = 0.0;

    // Initialization of the matrix log_alpha. In this matrix each column i will correspond to the log-forward probabilities 
    // at time T_ind for individual i (which will end up being the log-likelihood component corresponding to individual i).
    arma::mat log_alpha(N,n_ind);

    // First we store the log-forward probabilities at time t=1 for each of the observed individuals.
    for(arma::Col<int>::iterator i = ID_init.begin(); i != ID_init.end(); ++i){
        for (int state = 0; state < N; state++){
        
            log_alpha(state,*i-1) = std::log(init(state));
            log_alpha(state,*i-1) += log_gammapdf(x_duration_init[*i-1] ,
                                                    alpha_duration[state], 
                                                    beta_duration[state]);

            log_alpha(state,*i-1) += log_gammapdf(x_surface_init[*i-1] ,
                                                    alpha_surface[state], 
                                                    beta_surface[state]);

            log_alpha(state,*i-1) += log_gammapdf(x_maxDepth_init[*i-1] , 
                                                    alpha_maxDepth[state], 
                                                    beta_maxDepth[state]);

            log_alpha(state,*i-1) += log_zipoispdf(x_lunges_init[*i-1] , 
                                                lambda[state],
                                                theta[state]);

            // The data streams turning angle, number of steps and heading variance have missing values (NA's replaced with value 999.0), 
            // so an if statement is incorporated to deal with this
            if(x_angle_init[*i-1] != 999.0){
                // NumericVector aux_dvm = dvm(x_angle_init[*i-1] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
                log_alpha(state,*i-1) += log_dvm(x_angle_init[*i-1], 0.0, kappa[state]);

            }

            if(x_step_init[*i-1] != 999.0){
                log_alpha(state,*i-1) += log_gammapdf(x_step_init[*i-1] , 
                                                        alpha_step[state], 
                                                        beta_step[state]);
            }
            if(x_headVar_init[*i-1] != 999.0){
                log_alpha(state,*i-1) += std::log(betapdf(x_headVar_init[*i-1] , 
                                                        a[state], 
                                                        b[state]));
            }
        }
    }

    // Here we compute the log-forward probablities for times t > 1 until we reach the observed time T_i corresponding to individual i.
    for(int i = 0; i < n; i++){

        // Auxiliary vector in which the log-forward probabilities of individual ID[i] at time t+1 (given log-forward probs in 
        // log_alpha(_,ID[i]-1) correspond to time t) will be stored
        // Note: since C++ starts indexing from 0, log_alpha(_,ID[i]-1) refers to individual ID[i] 
        std::vector<double> aux_log_alpha(N);
        
        for(int state = 0; state < N; state++){

        // Auxiliary vector which is initialized as the log-forward probabilities of individual (ID[i]-1) at time t.         
        std::vector<double> aux(N);
        
        for(int k = 0; k < N; k++){
            aux[k] = log_alpha(k,ID[i]-1);

        }
        
        for(int previous_state = 0; previous_state < N; previous_state++){
            
            aux[previous_state] += std::log(tpm(previous_state,state)); 
            
            aux[previous_state] += log_gammapdf(x_duration[i] , 
                                            alpha_duration[state],
                                            beta_duration[state]);

            aux[previous_state] += log_gammapdf(x_surface[i] , 
                                                    alpha_surface[state], 
                                                    beta_surface[state]);

            aux[previous_state] += log_gammapdf(x_maxDepth[i] , 
                                                    alpha_maxDepth[state], 
                                                    beta_maxDepth[state]);

            aux[previous_state] += log_zipoispdf(x_lunges[i],
                                                lambda[state],
                                                theta[state]);

        }
        
        
        if(x_angle[i] != 999.0){
            // NumericVector aux_dvm = dvm(x_angle[i] , _["mu"] = 0.0, _["kappa"] = kappa[state]);

            for(int previous_state = 0; previous_state < N; previous_state++){
            // aux[previous_state] += std::log(*(aux_dvm.begin()));
            aux[previous_state] += log_dvm(x_angle[i], 0.0, kappa[state]);
            }
        }
        if(x_step[i] != 999.0){
            for(int previous_state = 0; previous_state < N; previous_state++){
            aux[previous_state] += log_gammapdf(x_step[i] , 
                                                    alpha_step[state], 
                                                    beta_step[state]);
            }
        }
        
        if(x_headVar[i] != 999.0){
            for(int previous_state = 0; previous_state < N; previous_state++){
            aux[previous_state] += std::log(betapdf(x_headVar[i] , a[state], b[state]));
            }
        }
        
        double aux_max = *std::max_element(std::begin(aux),std::end(aux));
        double sum_exp_aux_max = 0.0;

        for(int k = 0; k < N; k++){
            aux[k] -= aux_max;
            sum_exp_aux_max += std::exp(aux[k]);
        }

        // log-forward probabilities of individual (ID[i]-1) at time t+1
        aux_log_alpha[state] =  aux_max + std::log(sum_exp_aux_max);
        
        }
        
        for(int k = 0; k < N; k++){
            // Here we update the log-forward probabilities for the corresponding individual
            log_alpha(k,ID[i]-1) = aux_log_alpha[k];
        }
    }

    // We compute the log-likelihood using the log-sum-exp method for each one of the likelihood components
    for(int i=0; i < n_ind; i++){
        std::vector<double> aux_i(N);
        for(int k=0; k < N; k++){
            aux_i[k] = log_alpha(k,i);
        }

        double aux_max_i = *std::max_element(std::begin(aux_i),std::end(aux_i));
        double sum_exp_aux_max_i = 0.0;

        for(int k = 0; k < N; k++){
            aux_i[k] -= aux_max_i;
            sum_exp_aux_max_i += std::exp(aux_i[k]);
        }


        log_lik_total +=  aux_max_i + std::log(sum_exp_aux_max_i);
    }

    return log_lik_total;
}

// (Unnormalilzed) log-posterior of theta_star given data
// [[Rcpp::export]]
double lp__stdVector(const int N,
                    const int n,
                    const int n_ind,
                    arma::Col<int> ID_init,
                    const arma::Col<int> ID,
                    const arma::Col<double> x_duration_init,
                    const arma::Col<double> x_surface_init,
                    const arma::Col<double> x_maxDepth_init,
                    const arma::Col<int> x_lunges_init,
                    const arma::Col<double> x_step_init,
                    const arma::Col<double> x_angle_init,
                    const arma::Col<double> x_headVar_init,
                    const arma::Col<double> x_duration,
                    const arma::Col<double> x_surface,
                    const arma::Col<double> x_maxDepth,
                    const arma::Col<int> x_lunges,
                    const arma::Col<double> x_step,
                    const arma::Col<double> x_angle,
                    const arma::Col<double> x_headVar,
                    arma::Col<double> theta_star){

    // hyperparameters for the initial distributions of parameters of interest
    // duration
    // const arma::Col<double> mu_duration_mean = {140, 334, 516};
    // const arma::Col<double> mu_duration_sigma = {10,30,30};
    const arma::Col<double> mu_duration_mean = {140,334,516};
    const arma::Col<double> mu_duration_sigma = {100,300,300};
    // const arma::Col<double> sigma_duration_alpha = {500,267,845} ;
    // const arma::Col<double> sigma_duration_beta = {6.25,1.26,6.5} ;
    const arma::Col<double> sigma_duration_alpha = {20,20,33.8};
    const arma::Col<double> sigma_duration_beta = {.25,0.09433962,.26};
    // surface
    // const arma::Col<double> mu_surface_mean = {70,86,151};
    // const arma::Col<double> mu_surface_sigma = {6,10,10};
    const arma::Col<double> mu_surface_mean = {70,86,151};
    const arma::Col<double> mu_surface_sigma = {60,50,100};
    // const arma::Col<double> sigma_surface_alpha = {361.25,151.25,661.25};
    // const arma::Col<double> sigma_surface_beta = {5.3125,2.75,9.6};
    const arma::Col<double> sigma_surface_alpha = {9.248,67.222,9.9146189};
    const arma::Col<double> sigma_surface_beta = {.136,1.222,0.1436901};
    // maxDepth
    // const arma::Col<double> mu_maxDepth_mean = {32,68,170};
    // const arma::Col<double> mu_maxDepth_sigma = {5,10,5};
    const arma::Col<double> mu_maxDepth_mean = {32,68,170};
    const arma::Col<double> mu_maxDepth_sigma = {25,50,100};
    // const arma::Col<double> sigma_maxDepth_alpha = {10,30,40};
    // const arma::Col<double> sigma_maxDepth_beta = {2,2,2};
    const arma::Col<double> sigma_maxDepth_alpha = {180,845,45};
    const arma::Col<double> sigma_maxDepth_beta = {7.5,13,.75};
    // step
    // const arma::Col<double> mu_step_mean = {189,675,406};
    // const arma::Col<double> mu_step_sigma = {15,45,30};
    const arma::Col<double> mu_step_mean = {189,675,406};
    const arma::Col<double> mu_step_sigma = {150, 450, 300};
    // const arma::Col<double> sigma_step_alpha = {134.0,305.0,287.0};
    // const arma::Col<double> sigma_step_beta = {2,2,2};
    const arma::Col<double> sigma_step_alpha = {8.978, 46.5125, 41.1845};
    const arma::Col<double> sigma_step_beta = {.067, .1525, .1435};
    // angle
    // const arma::Col<double> kappa_alpha = {125,133.47,80};
    // const arma::Col<double> kappa_beta = {125,43.05,100};
    const arma::Col<double> kappa_alpha = {.3125,1.922,.20}; // Remember--- for state 2, this values are different to those used in Stan
    const arma::Col<double> kappa_beta = {.3125,.620,.25};
    // heading variance
    // const arma::Col<double> a_alpha = {125,125,361.25};
    // const arma::Col<double> a_beta = {125,250,212.5};
    // const arma::Col<double> b_alpha = {245,5.4,50.45};
    // const arma::Col<double> b_beta = {116.6,1,9.34};
    const arma::Col<double> a_alpha = {1.25,.2,14.45};
    const arma::Col<double> a_beta = {1.25,.4,8.5};
    const arma::Col<double> b_alpha = {5.5125,145.8,12.8};
    const arma::Col<double> b_beta = {2.6250,27,8};
    // lunges
    // const arma::Col<double> lambda_alpha = {245,.0138,1445};
    // const arma::Col<double> lambda_beta = {350,2.77,425};
    const arma::Col<double> lambda_alpha = {.6125,.000125,14.45};
    const arma::Col<double> lambda_beta = {.8750,.025,4.25};


    arma::Col<double> mu_duration(N);
    arma::Col<double> sigma_duration(N);
    arma::Col<double> mu_surface(N);
    arma::Col<double> sigma_surface(N);
    arma::Col<double> mu_maxDepth(N);
    arma::Col<double> sigma_maxDepth(N);
    arma::Col<double> mu_step(N);
    arma::Col<double> sigma_step(N);
    arma::Col<double> kappa(N);
    arma::Col<double> a(N);
    arma::Col<double> b(N);
    arma::Col<double> lambda(N);
    arma::Col<double> theta(N); // Since theta_i ~ U(0,1), there's no need to write down the prior for this variable

    for(int i = 0; i<N; i++){
        
        mu_duration[i] = theta_star[i+6];
        sigma_duration[i] = theta_star[i+9];

        mu_surface[i] = theta_star[i+12];
        sigma_surface[i] = theta_star[i+15];

        mu_maxDepth[i] = theta_star[i+18];
        sigma_maxDepth[i] = theta_star[i+21];

        mu_step[i] = theta_star[i+24];
        sigma_step[i] = theta_star[i+27];

        kappa[i] = theta_star[i+30];

        a[i] = theta_star[i+33];
        b[i] = theta_star[i+36];

        lambda[i] = theta_star[i + 39];
        theta[i] = theta_star[i + 48];
        
    }

    // NumericVector init_raw = {theta_star[42], theta_star[43]}; // pi parameters
    
    
    // initial distribution
    arma::Col<double> init(N);
    init[2-1] = theta_star[42];
    init[3-1] = theta_star[43];
  
    // tpm 
    arma::mat tpm(N,N);
  
    // creation tpm
    tpm(1-1,1-1) = theta_star[44];
    tpm(1-1,2-1) = theta_star[1-1];
    tpm(1-1,3-1) = theta_star[2-1];
    
    tpm(2-1,1-1) = theta_star[3-1];
    tpm(2-1,2-1) = theta_star[45];
    tpm(2-1,3-1) = theta_star[4-1];    
    
    tpm(3-1,1-1) = theta_star[5-1];
    tpm(3-1,2-1) = theta_star[6-1];
    tpm(3-1,3-1) = theta_star[46];

    // Initialization of the rows of the transition probability matrix (easier to handle in this way)
    arma::Col<double> tpm1(N);
    arma::Col<double> tpm2(N);
    arma::Col<double> tpm3(N);

    for(int j =0; j < N; j++){
        tpm1[j] = tpm(1-1,j); 
        tpm2[j] = tpm(2-1,j); 
        tpm3[j] = tpm(3-1,j); 
    }

    init[1-1] = theta_star[47];

    // Initialization of the double variable where the unnormalized log-posterior will be stored
    double target = 0.0;
    
    for(int i = 0; i < N; i++){
        
        // Prior distributions associated to the duration parameters

        // mu_duration ~ truncated from below at zero normal distribution
        target+= normal_lpdf(mu_duration[i] , mu_duration_mean[i] , mu_duration_sigma[i]) - std_normal_lcdf(-mu_duration_mean[i]/mu_duration_sigma[i]);
        // sigma_duration ~ gamma distribution
        target+= log_gammapdf(sigma_duration[i] , sigma_duration_alpha[i], sigma_duration_beta[i]);
    
        // Prior distributions associated to the surface parameters

        // mu_surface ~ truncated from below at zero normal distribution
        target+= normal_lpdf(mu_surface[i] , mu_surface_mean[i] , mu_surface_sigma[i]) -  std_normal_lcdf(-mu_surface_mean[i]/mu_surface_sigma[i]);
        // sigma_surface ~ gamma distribution        
        target+= log_gammapdf(sigma_surface[i], sigma_surface_alpha[i], sigma_surface_beta[i]); 
    
        // Prior distributions associated to the maximum depth parameters
        
        // mu_maxDepth ~ truncated from below at zero normal distribution
        target+= normal_lpdf(mu_maxDepth[i] , mu_maxDepth_mean[i] , mu_maxDepth_sigma[i]) - std_normal_lcdf(-mu_maxDepth_mean[i]/mu_maxDepth_sigma[i]);
        // sigma_maxDepth ~ gamma distribution
        target+= log_gammapdf(sigma_maxDepth[i] , sigma_maxDepth_alpha[i], sigma_maxDepth_beta[i]);

        // Prior distributions associated to the number of steps parameters

        // mu_step ~ truncated from below at zero normal distribution
        target+= normal_lpdf(mu_step[i] , mu_step_mean[i] , mu_step_sigma[i]) - std_normal_lcdf(-mu_step_mean[i]/mu_step_sigma[i]);
        // sigma_step ~ gamma distribution        
        target+= log_gammapdf(sigma_step[i], sigma_step_alpha[i], sigma_step_beta[i]);

        // Prior distributions associated to the turning angle parameters

        // kappa ~ gamma distribution
        target+= log_gammapdf(kappa[i] , kappa_alpha[i], kappa_beta[i]);

        // Prior distributions associated to the heading variance parameters
        // a ~ gamma distribution
        target+= log_gammapdf(a[i] , a_alpha[i], a_beta[i]);
        // b ~ gamma distribution
        target+= log_gammapdf(b[i] , b_alpha[i], b_beta[i]);

        // Prior distributions associated to the number of lunges parameters

        // lambda ~ gamma distribution 
        target+= log_gammapdf(lambda[i] , lambda_alpha[i], lambda_beta[i]);

        // For the probability theta, we assumed theta ~ U(0,1), so there's no need to write down code to include this
        

    }
  
    arma::Col<double> aux_alpha_dir(N);
    for(int j = 0 ; j < N; j++){
        aux_alpha_dir[j] += 1;
    } 

    // Prior distributions associated to the initial state distribution and transition probability matrix rows
    // initial state distribution ~ dirichlet
    target += log_ddirichlet(init,aux_alpha_dir);

    // rows tpm ~ dirichlet
    target += log_ddirichlet(tpm1,aux_alpha_dir);
    target += log_ddirichlet(tpm2,aux_alpha_dir);
    target += log_ddirichlet(tpm3,aux_alpha_dir);

    // Here we incorporate the likelihood
    target += log_lik_stdVector(N,
                                n,
                                n_ind,
                                ID_init,
                                ID,
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
                                x_headVar, theta_star);
    return target;   

}


// For the implementation of RcppParallel, we need to define a Worker object
struct ComponentWiseStep : public Worker
{

    const RMatrix<double> the_Y_n_matrix; // input matrix
    
    // data
    const int N;
    const int n;
    const int n_ind;
    const RVector<int> ID_init;
    const RVector<int> ID;
    const RVector<double> x_duration_init;
    const RVector<double> x_surface_init;
    const RVector<double> x_maxDepth_init;
    const RVector<int> x_lunges_init;
    const RVector<double> x_step_init;
    const RVector<double> x_angle_init;
    const RVector<double> x_headVar_init;
    const RVector<double> x_duration;
    const RVector<double> x_surface;
    const RVector<double> x_maxDepth;
    const RVector<int> x_lunges;
    const RVector<double> x_step;
    const RVector<double> x_angle;
    const RVector<double> x_headVar;

    const RVector<double> X_i; // samples from iteration i
    const double the_target_at_X_i; // value of the tempered target at the iteration i
    const std::size_t size; // size of the parameter vector
    const double the_temp; // the temperature related to the tempered version of the target
    const double max_temp;
    const RVector<double> log_U; // log-uniform samples
    // double log_alpha;
    // double log_ratio;
    // double the_target_at_proposal; // value fo the proposal at the tempered target
    
    RVector<double> output; // output vector 

    // void  target_tau_(arma::Col<int> aux_ID_init_in,
    double  target_tau_(arma::Col<int> aux_ID_init_in,
                    const arma::Col<int> aux_ID,
                    const arma::Col<double> aux_x_duration_init,
                    const arma::Col<double> aux_x_surface_init,
                    const arma::Col<double> aux_x_maxDepth_init,
                    const arma::Col<int> aux_x_lunges_init,
                    const arma::Col<double> aux_x_step_init,
                    const arma::Col<double> aux_x_angle_init,
                    const arma::Col<double> aux_x_headVar_init,
                    const arma::Col<double> aux_x_duration,
                    const arma::Col<double> aux_x_surface,
                    const arma::Col<double> aux_x_maxDepth,
                    const arma::Col<int> aux_x_lunges,
                    const arma::Col<double> aux_x_step,
                    const arma::Col<double> aux_x_angle,
                    const arma::Col<double> aux_x_headVar,
                    //here is the col of Y matrix
                    arma::Col<double> Y_matrix);

    // convert RVector/RMatrix into arma type for target_tau_ function and the following arma data will be shared in parallel computing
    arma::mat convert(){
        RMatrix<double> tmp_the_Y_n_matrix = the_Y_n_matrix;
        arma::mat MAT(tmp_the_Y_n_matrix.begin(), size, size, false);
        return MAT;
    }
    
    arma::Col<int> convert_ID_init(){
        RVector<int> tmp_ID_init = ID_init;
        arma::Col<int> vec(tmp_ID_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<int> convert_ID(){
        RVector<int> tmp_ID = ID;
        arma::Col<int> vec(tmp_ID.begin(),n,false);
        return vec;
    }

    arma::Col<double> convert_x_duration_init(){
        RVector<double> tmp_x_duration_init = x_duration_init;
        arma::Col<double> vec(tmp_x_duration_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<double> convert_x_surface_init(){
        RVector<double> tmp_x_surface_init = x_surface_init;
        arma::Col<double> vec(tmp_x_surface_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<double> convert_x_maxDepth_init(){
        RVector<double> tmp_x_maxDepth_init = x_maxDepth_init;
        arma::Col<double> vec(tmp_x_maxDepth_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<int> convert_x_lunges_init(){
        RVector<int> tmp_x_lunges_init = x_lunges_init;
        arma::Col<int> vec(tmp_x_lunges_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<double> convert_x_step_init(){
        RVector<double> tmp_x_step_init = x_step_init;
        arma::Col<double> vec(tmp_x_step_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<double> convert_x_angle_init(){
        RVector<double> tmp_x_angle_init = x_angle_init;
        arma::Col<double> vec(tmp_x_angle_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<double> convert_x_headVar_init(){
        RVector<double> tmp_x_headVar_init = x_headVar_init;
        arma::Col<double> vec(tmp_x_headVar_init.begin(),n_ind,false);
        return vec;
    }

    arma::Col<double> convert_x_duration(){
        RVector<double> tmp_x_duration = x_duration;
        arma::Col<double> vec(tmp_x_duration.begin(),n,false);
        return vec;
    }

    arma::Col<double> convert_x_surface(){
        RVector<double> tmp_x_surface = x_surface;
        arma::Col<double> vec(tmp_x_surface.begin(),n,false);
        return vec;
    }

    arma::Col<double> convert_x_maxDepth(){
        RVector<double> tmp_x_maxDepth = x_maxDepth;
        arma::Col<double> vec(tmp_x_maxDepth.begin(),n,false);
        return vec;
    }

    arma::Col<int> convert_x_lunges(){
        RVector<int> tmp_x_lunges = x_lunges;
        arma::Col<int> vec(tmp_x_lunges.begin(),n,false);
        return vec;
    }

    arma::Col<double> convert_x_step(){
        RVector<double> tmp_x_step = x_step;
        arma::Col<double> vec(tmp_x_step.begin(),n,false);
        return vec;
    }

    arma::Col<double> convert_x_angle(){
        RVector<double> tmp_x_angle = x_angle;
        arma::Col<double> vec(tmp_x_angle.begin(),n,false);
        return vec;
    }

    arma::Col<double> convert_x_headVar(){
        RVector<double> tmp_x_headVar = x_headVar;
        arma::Col<double> vec(tmp_x_headVar.begin(),n,false);
        return vec;
    }

    // The main constructor
    ComponentWiseStep (
            const NumericMatrix the_Y_n_matrix_in,
            const int N_in,
            const int n_in,
            const int n_ind_in,
            IntegerVector ID_init_in,
            const IntegerVector ID_in,
            const NumericVector x_duration_init_in,
            const NumericVector x_surface_init_in,
            const NumericVector x_maxDepth_init_in,
            const IntegerVector x_lunges_init_in,
            const NumericVector x_step_init_in,
            const NumericVector x_angle_init_in,
            const NumericVector x_headVar_init_in,
            const NumericVector x_duration_in,
            const NumericVector x_surface_in,
            const NumericVector x_maxDepth_in,
            const IntegerVector x_lunges_in,
            const NumericVector x_step_in,
            const NumericVector x_angle_in,
            const NumericVector x_headVar_in,
            // other parameters
            const NumericVector X_i_in,
            const double the_target_at_X_i_in,
            const size_t size_in,
            const double the_temp_in,
            const double max_temp_in,
            const NumericVector log_U_in,
            NumericVector output_in) :
        the_Y_n_matrix(the_Y_n_matrix_in),
        N(N_in),
        n(n_in),
        n_ind(n_ind_in),
        ID_init(ID_init_in),
        ID(ID_in),
        x_duration_init(x_duration_init_in),
        x_surface_init(x_surface_init_in),
        x_maxDepth_init(x_maxDepth_init_in),
        x_lunges_init(x_lunges_init_in),
        x_step_init(x_step_init_in),
        x_angle_init(x_angle_init_in),
        x_headVar_init(x_headVar_init_in),
        x_duration(x_duration_in),
        x_surface(x_surface_in),
        x_maxDepth(x_maxDepth_in),
        x_lunges(x_lunges_in),
        x_step(x_step_in),
        x_angle(x_angle_in),
        x_headVar(x_headVar_in),
        X_i(X_i_in),
        the_target_at_X_i(the_target_at_X_i_in),
        size(size_in),the_temp(the_temp_in),
        max_temp(max_temp_in), log_U(log_U_in),
        // log_alpha(0), log_ratio(0), the_target_at_proposal(0),
        output(output_in)
    {
        // log_alpha -= the_target_at_X_i;

    }


    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        // rows we will operate on
        for (std::size_t i = begin; i < end; i++)
        {
            arma::Col<int> aux_ID_init = convert_ID_init();
            const arma::Col<int> aux_ID = convert_ID();
            const arma::Col<double> aux_x_duration_init = convert_x_duration_init();
            const arma::Col<double> aux_x_surface_init = convert_x_surface_init();
            const arma::Col<double> aux_x_maxDepth_init = convert_x_maxDepth_init();
            const arma::Col<int> aux_x_lunges_init = convert_x_lunges_init();
            const arma::Col<double> aux_x_step_init = convert_x_step_init();
            const arma::Col<double> aux_x_angle_init = convert_x_angle_init();
            const arma::Col<double> aux_x_headVar_init = convert_x_headVar_init();
            const arma::Col<double> aux_x_duration = convert_x_duration();
            const arma::Col<double> aux_x_surface = convert_x_surface();
            const arma::Col<double> aux_x_maxDepth = convert_x_maxDepth();
            const arma::Col<int> aux_x_lunges = convert_x_lunges();
            const arma::Col<double> aux_x_step = convert_x_step();
            const arma::Col<double> aux_x_angle = convert_x_angle();
            const arma::Col<double> aux_x_headVar = convert_x_headVar();

            arma::mat MAT = convert();

            if( i == 2){

                arma::Col<double> col_i = MAT.col(0);

                double log_alpha = 0.0;
                double log_ratio = 0.0;
                double the_target_at_proposal = target_tau_(aux_ID_init,
                            aux_ID,
                            aux_x_duration_init,
                            aux_x_surface_init,
                            aux_x_maxDepth_init,
                            aux_x_lunges_init,
                            aux_x_step_init,
                            aux_x_angle_init,
                            aux_x_headVar_init,
                            aux_x_duration,
                            aux_x_surface,
                            aux_x_maxDepth,
                            aux_x_lunges,
                            aux_x_step,
                            aux_x_angle,
                            aux_x_headVar,
                            col_i);
                
                log_alpha += the_target_at_proposal;
                log_alpha -= the_target_at_X_i;
                
                if(log_alpha < 0){
                    log_ratio += log_alpha;
                }

                if(log_U[2] <= log_ratio){
                    output[1-1] = the_Y_n_matrix(1-1,1-1);
                    output[2-1] = the_Y_n_matrix(2-1,1-1);
                    output[44] = the_Y_n_matrix(44,1-1);
                } else {
                    output[1-1] = X_i[1-1];
                    output[2-1] = X_i[2-1];
                    output[44] = X_i[44];
                }

            }  else if ( i == 3)
            {
                arma::Col<double> col_i = MAT.col(2);
                double log_alpha = 0.0;
                double log_ratio = 0.0;
                double the_target_at_proposal = target_tau_(aux_ID_init,
                            aux_ID,
                            aux_x_duration_init,
                            aux_x_surface_init,
                            aux_x_maxDepth_init,
                            aux_x_lunges_init,
                            aux_x_step_init,
                            aux_x_angle_init,
                            aux_x_headVar_init,
                            aux_x_duration,
                            aux_x_surface,
                            aux_x_maxDepth,
                            aux_x_lunges,
                            aux_x_step,
                            aux_x_angle,
                            aux_x_headVar,
                            col_i);
                
                log_alpha += the_target_at_proposal;
                log_alpha -= the_target_at_X_i;
                
                if(log_alpha < 0){
                    log_ratio += log_alpha;
                }

                if(log_U[3] <= log_ratio){
                    output[3-1] = the_Y_n_matrix(3-1,3-1);
                    output[4-1] = the_Y_n_matrix(4-1,3-1);
                    output[45] = the_Y_n_matrix(45,3-1);
                } else {
                    output[3-1] = X_i[3-1];
                    output[4-1] = X_i[4-1];
                    output[45] = X_i[45];
                }

            } else if ( i == 4)
            {
                arma::Col<double> col_i = MAT.col(4);
                double log_alpha = 0.0;
                double log_ratio = 0.0;
                double the_target_at_proposal = target_tau_(aux_ID_init,
                            aux_ID,
                            aux_x_duration_init,
                            aux_x_surface_init,
                            aux_x_maxDepth_init,
                            aux_x_lunges_init,
                            aux_x_step_init,
                            aux_x_angle_init,
                            aux_x_headVar_init,
                            aux_x_duration,
                            aux_x_surface,
                            aux_x_maxDepth,
                            aux_x_lunges,
                            aux_x_step,
                            aux_x_angle,
                            aux_x_headVar,
                            col_i);
                
                log_alpha += the_target_at_proposal;
                log_alpha -= the_target_at_X_i;
                
                if(log_alpha < 0){
                    log_ratio += log_alpha;
                }

                if(log_U[4] <= log_ratio){
                    output[5-1] = the_Y_n_matrix(5-1,5-1);
                    output[6-1] = the_Y_n_matrix(6-1,5-1);
                    output[46] = the_Y_n_matrix(46,5-1);
                } else {
                    output[5-1] = X_i[5-1];
                    output[6-1] = X_i[6-1];
                    output[46] = X_i[46];
                }

            } else if ( i == 5)
            {
                arma::Col<double> col_i = MAT.col(42);
                double log_alpha = 0.0;
                double log_ratio = 0.0;
                double the_target_at_proposal = target_tau_(aux_ID_init,
                            aux_ID,
                            aux_x_duration_init,
                            aux_x_surface_init,
                            aux_x_maxDepth_init,
                            aux_x_lunges_init,
                            aux_x_step_init,
                            aux_x_angle_init,
                            aux_x_headVar_init,
                            aux_x_duration,
                            aux_x_surface,
                            aux_x_maxDepth,
                            aux_x_lunges,
                            aux_x_step,
                            aux_x_angle,
                            aux_x_headVar,
                            col_i);
                
                log_alpha += the_target_at_proposal;
                log_alpha -= the_target_at_X_i;
                
                if(log_alpha < 0){
                    log_ratio += log_alpha;
                }

                if(log_U[5] <= log_ratio){
                    output[42] = the_Y_n_matrix(42,42);
                    output[43] = the_Y_n_matrix(43,42);
                    output[47] = the_Y_n_matrix(47,42);
                } else {
                    output[42] = X_i[42];
                    output[43] = X_i[43];
                    output[47] = X_i[47];
                }

            } else{
                arma::Col<double> col_i = MAT.col(i);

                double log_alpha = 0.0;
                double log_ratio = 0.0;

                double the_target_at_proposal = target_tau_(aux_ID_init,
                            aux_ID,
                            aux_x_duration_init,
                            aux_x_surface_init,
                            aux_x_maxDepth_init,
                            aux_x_lunges_init,
                            aux_x_step_init,
                            aux_x_angle_init,
                            aux_x_headVar_init,
                            aux_x_duration,
                            aux_x_surface,
                            aux_x_maxDepth,
                            aux_x_lunges,
                            aux_x_step,
                            aux_x_angle,
                            aux_x_headVar,
                            col_i);

                log_alpha += the_target_at_proposal;
                log_alpha -= the_target_at_X_i;
                // Hasting ratio
                // if(i != 48 && i != 49 && i !=50){
                if(i < (size - 3)){
                    log_alpha +=  std::log(the_Y_n_matrix(i,i));
                    log_alpha -=  std::log(X_i[i]);
                }

                if(log_alpha < 0){
                    log_ratio += log_alpha;
                }

                if(log_U[i] <= log_ratio){
                    output[i] = the_Y_n_matrix(i,i);
                } else {
                    output[i] = X_i[i];    
                }
            }
    
        }
    } // end parallel function operator

};

// void ComponentWiseStep::target_tau_(arma::Col<int> aux_ID_init,
double ComponentWiseStep::target_tau_(arma::Col<int> aux_ID_init,
                                    const arma::Col<int> aux_ID,
                                    const arma::Col<double> aux_x_duration_init,
                                    const arma::Col<double> aux_x_surface_init,
                                    const arma::Col<double> aux_x_maxDepth_init,
                                    const arma::Col<int> aux_x_lunges_init,
                                    const arma::Col<double> aux_x_step_init,
                                    const arma::Col<double> aux_x_angle_init,
                                    const arma::Col<double> aux_x_headVar_init,
                                    const arma::Col<double> aux_x_duration,
                                    const arma::Col<double> aux_x_surface,
                                    const arma::Col<double> aux_x_maxDepth,
                                    const arma::Col<int> aux_x_lunges,
                                    const arma::Col<double> aux_x_step,
                                    const arma::Col<double> aux_x_angle,
                                    const arma::Col<double> aux_x_headVar,
                                    //here is the col of Y matrix
                                    arma::Col<double> Y_matrix){
    
    // the_target_at_proposal = lp__stdVector(N,
    return lp__stdVector(N,
                                        n,
                                        n_ind,
                                        aux_ID_init,
                                        aux_ID,
                                        aux_x_duration_init,
                                        aux_x_surface_init,
                                        aux_x_maxDepth_init,
                                        aux_x_lunges_init,
                                        aux_x_step_init,
                                        aux_x_angle_init,
                                        aux_x_headVar_init,
                                        aux_x_duration,
                                        aux_x_surface,
                                        aux_x_maxDepth,
                                        aux_x_lunges,
                                        aux_x_step,
                                        aux_x_angle,
                                        aux_x_headVar, 
                                        Y_matrix)/the_temp;
}

// This functions calls the ComponentWiseStep Worker we defined above
void vector_aggregator(NumericMatrix Y,
                                const int N,
                                const int n,
                                const int n_ind,
                                IntegerVector ID_init,
                                const IntegerVector ID,
                                const NumericVector x_duration_init,
                                const NumericVector x_surface_init,
                                const NumericVector x_maxDepth_init,
                                const IntegerVector x_lunges_init,
                                const NumericVector x_step_init,
                                const NumericVector x_angle_init,
                                const NumericVector x_headVar_init,
                                const NumericVector x_duration,
                                const NumericVector x_surface,
                                const NumericVector x_maxDepth,
                                const IntegerVector x_lunges,
                                const NumericVector x_step,
                                const NumericVector x_angle,
                                const NumericVector x_headVar,
                                NumericVector output_parallel,
                                NumericVector X_i, NumericVector vec_log_U, double target_X_i, double temp, double maxtemp, int begin_in, int end_in){

    const std::size_t size = static_cast <size_t> (Y.nrow());
    const std::size_t begin = static_cast <size_t> (begin_in);
    const std::size_t end = static_cast <size_t> (end_in);        
    
    ComponentWiseStep componentWiseStep(Y,
                                        //data
                                        N,
                                        n,
                                        n_ind,
                                        ID_init,
                                        ID,
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
                                        //other parameters    
                                        X_i,target_X_i,size, temp, maxtemp,vec_log_U,
                                        output_parallel);

    parallelFor(begin,end,componentWiseStep); 
    // return output_parallel;

}

// Parallelized version of the function pt_cw_M_target_posterior
// [[Rcpp::export]]
List rcpp_parallel_pt_cw_M_target_posterior(const List list_data, 
                        const int nsim,  
                        NumericVector init,
                        NumericVector temp_vector,
                        const int N,
                        const int n,
                        const int n_ind,
                        IntegerVector ID_init,
                        const  IntegerVector ID,
                        const  NumericVector x_duration_init,
                        const  NumericVector x_surface_init,
                        const  NumericVector x_maxDepth_init,
                        const  IntegerVector x_lunges_init,
                        const  NumericVector x_step_init,
                        const  NumericVector x_angle_init,
                        const  NumericVector x_headVar_init,
                        const  NumericVector x_duration,
                        const  NumericVector x_surface,
                        const  NumericVector x_maxDepth,
                        const  IntegerVector x_lunges,
                        const  NumericVector x_step,
                        const  NumericVector x_angle,
                        const  NumericVector x_headVar,
                        bool display_progress=true){
    
    // Progress p(nsim*temp_vector.size()*init.size(), display_progress);

    // This is to display a progress bar and see how long the function it's taking to make samples (helpful for debugging and track errors)
    Progress p(nsim*temp_vector.size(), display_progress);

    // List where each element will allocate samples from the different tempered versions of the target
    List X(temp_vector.size());

    // Initial values provided for each tempered target
    for(int i = 0; i < temp_vector.size(); ++i){
        X[i] = createMatrix(init.size(),nsim+1);
        NumericMatrix aux_X = X[i];
        aux_X(_,0) = init;
    }

    // We will generate nsim simulations for each chain
    for(int i=0; i<nsim; ++i){
        // Part of the progress bar - it's supposed to retun -1.0 if the process is aborted
        if (Progress::check_abort() )
            return -1.0;

        // This loop corresponds to the component-wise step for each chain        
        for(int temp = 0; temp < temp_vector.size(); ++temp){
            // For the chain temp we select the its corresponding matrix where samples are allocated 
            NumericMatrix aux_X = X[temp];
            // Generation of new proposal for the parameters
            NumericVector aux_Y_n = proposal_dist(aux_X(_,i));

            // NumericVector aux_Y_n(init.size());

            // // Initialization of vectors where new proposal for the initial state dist and the tpms rows will be allocated
            // NumericVector new_init(N);
            // NumericVector new_tpm1(N);
            // NumericVector new_tpm2(N);
            // NumericVector new_tpm3(N);

            // // step size for the Metropolis algorithm for each one of the state-dependent distribution parameters
            // NumericVector step_size_mu_duration = {.008,.008,.008};
            // NumericVector step_size_sigma_duration = {.01,.01,.01};
            // NumericVector step_size_mu_surface = {.01,.01,.01};
            // NumericVector step_size_sigma_surface = {.05,.05,.05};
            // NumericVector step_size_mu_maxDepth = {.01,.01,.01};
            // NumericVector step_size_sigma_maxDepth = {.01,.01,.01};
            // NumericVector step_size_mu_step = {.01,.01,.01};
            // NumericVector step_size_sigma_step = {.01,.01,.01};
            // NumericVector step_size_kappa = {.08,.08,.08}; // third attempt, now working on the log scale
            // NumericVector step_size_a = {.05,.05,.05}; // second attempt, now working on the log scale
            // NumericVector step_size_b = {.05,.05,.05};
            // NumericVector step_size_lambda = {.05,.05,.05}; // second attempt, now working on the log scale    

            // for(int j=0; j < N; j++){

            //     aux_Y_n[j+6] = std::exp(std::log(aux_X(j+6,i)) + R::rnorm(0.0,step_size_mu_duration[j]));
            //     aux_Y_n[j+9] = std::exp(std::log(aux_X(j+9,i)) + R::rnorm(0.0,step_size_sigma_duration[j]));
            //     aux_Y_n[j+12] = std::exp(std::log(aux_X(j+12,i)) + R::rnorm(0.0,step_size_mu_surface[j]));
            //     aux_Y_n[j+15] = std::exp(std::log(aux_X(j+15,i)) + R::rnorm(0.0,step_size_sigma_surface[j]));
            //     aux_Y_n[j+18] = std::exp(std::log(aux_X(j+18,i)) + R::rnorm(0.0,step_size_mu_maxDepth[j]));
            //     aux_Y_n[j+21] = std::exp(std::log(aux_X(j+21,i)) + R::rnorm(0.0,step_size_sigma_maxDepth[j]));
            //     aux_Y_n[j+24] = std::exp(std::log(aux_X(j+24,i)) + R::rnorm(0.0,step_size_mu_step[j]));
            //     aux_Y_n[j+27] = std::exp(std::log(aux_X(j+27,i)) + R::rnorm(0.0,step_size_sigma_step[j]));
            //     aux_Y_n[j+30] = std::exp(std::log(aux_X(j+30,i)) + R::rnorm(0.0,step_size_kappa[j]));
            //     aux_Y_n[j+33] = std::exp(std::log(aux_X(j+33,i)) + R::rnorm(0.0,step_size_a[j]));
            //     aux_Y_n[j+36] = std::exp(std::log(aux_X(j+36,i)) + R::rnorm(0.0,step_size_b[j]));
            //     aux_Y_n[j+39] = std::exp(std::log(aux_X(j+39,i)) + R::rnorm(0.0,step_size_lambda[j]));
            //     aux_Y_n[j+48] = R::runif(0,1);
            // }

            // // Proposals for initial state distribution and tpm rows
            // NumericVector dummy_vec = NumericVector::create(1/3,1/3,1/3);
            // NumericVector aux_new_init(N);
            // NumericVector aux_new_tpm1(N);
            // NumericVector aux_new_tpm2(N);
            // NumericVector aux_new_tpm3(N);

            // for(int j=0; j< dummy_vec.size(); j++){
            //     aux_new_init[j] =  -std::log(R::runif(0,1));    
            //     aux_new_tpm1[j] =  -std::log(R::runif(0,1));    
            //     aux_new_tpm2[j] =  -std::log(R::runif(0,1));    
            //     aux_new_tpm3[j] =  -std::log(R::runif(0,1));
            //     // if(temp == 0){
            //     //     Rprintf("sim for tpm1, temp %i is %f \n",temp,aux_new_tpm1[j]);
            //     //     Rprintf("sim for tpm2, temp %i is %f \n",temp,aux_new_tpm2[j]);
            //     //     Rprintf("sim for tpm3, temp %i is %f \n",temp,aux_new_tpm3[j]);
            //     //     Rprintf("sim for init, temp %i is %f \n",temp,aux_new_init[j]);

            //     // }

            // }

            // double aux_sum_init = 0.0;
            // double aux_sum_tpm1 = 0.0;
            // double aux_sum_tpm2 = 0.0;
            // double aux_sum_tpm3 = 0.0;

            // for(int j=0; j < dummy_vec.size(); j++){
            //     aux_sum_init += aux_new_init[j]; 
            //     aux_sum_tpm1 += aux_new_tpm1[j]; 
            //     aux_sum_tpm2 += aux_new_tpm2[j]; 
            //     aux_sum_tpm3 += aux_new_tpm3[j]; 
            // }
            // for(int j=0; j < dummy_vec.size();j++){
            //     new_init[j] = aux_new_init[j]/aux_sum_init;
            //     new_tpm1[j] = aux_new_tpm1[j]/aux_sum_tpm1;
            //     new_tpm2[j] = aux_new_tpm2[j]/aux_sum_tpm2;
            //     new_tpm3[j] = aux_new_tpm3[j]/aux_sum_tpm3;
            //     // if(temp == 0 && i ==10){
            //     //     Rprintf("prop for tpm1, temp %i is %f \n",temp,new_tpm1[j]);
            //     //     Rprintf("prop for tpm2, temp %i is %f \n",temp,new_tpm2[j]);
            //     //     Rprintf("prop for tpm3, temp %i is %f \n",temp,new_tpm3[j]);
            //     //     Rprintf("prop for init, temp %i is %f \n",temp,new_init[j]);

            //     // }
            // }

            // // proposed tpm values
            // aux_Y_n[44] = new_tpm1[1-1];
            // aux_Y_n[1-1] = new_tpm1[2-1];
            // aux_Y_n[2-1] = new_tpm1[3-1];
            
            // aux_Y_n[3-1] = new_tpm2[1-1];
            // aux_Y_n[45] = new_tpm2[2-1];
            // aux_Y_n[4-1] = new_tpm2[3-1];
            
            // aux_Y_n[5-1] = new_tpm3[1-1];
            // aux_Y_n[6-1] = new_tpm3[2-1];
            // aux_Y_n[46] = new_tpm3[3-1];

            // // proposed init values
            // aux_Y_n[42] = new_init[2-1];
            // aux_Y_n[43] = new_init[3-1];
            // aux_Y_n[47] = new_init[1-1];



            // Auxiliary vector for the component-wise step for each parameter            
            // NumericVector Y_n(init.length());
            
            // for(int id = 0; id <  init.length(); id++){
            //     Y_n[id] = aux_X(id,i);
            // }

            // tempered target at current iteration
            // double target_at_X_i = target_tau(list_data, aux_X(_,i),temp_vector[temp],temp_vector[temp_vector.size()-1]);
            
            double target_at_X_i = lp__stdVector(N,
                    n,
                    n_ind,
                    ID_init,
                    ID,
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
                    aux_X(_,i))/temp_vector[temp];

            // component-wise for tpms and initial distribution
            //  component-wise step for tpm1 (tpm row 1)
            // p.increment();
            // double log_U_tpm1 = log(R::runif(0,1));

            // Y_n[44] =  aux_Y_n[44];
            // Y_n[1-1] =  aux_Y_n[1-1];
            // Y_n[2-1] =  aux_Y_n[2-1];

            // double log_alpha_tpm1 = 0.0;
            // log_alpha_tpm1 += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            // log_alpha_tpm1 -= target_at_X_i;
            // double log_ratio_tpm1 = min(NumericVector::create(0,log_alpha_tpm1));
            
            // if(log_U_tpm1 <= log_ratio_tpm1){
            //     aux_X(44,i+1) = aux_Y_n[44];
            //     aux_X(1-1,i+1) = aux_Y_n[1-1];
            //     aux_X(2-1,i+1) = aux_Y_n[2-1];
            // } else {
            //     aux_X(44,i+1) = aux_X(44,i);
            //     aux_X(1-1,i+1) = aux_X(1-1,i);
            //     aux_X(2-1,i+1) = aux_X(2-1,i);
            // }

            // // return the value we had before to recycle this vector
            // Y_n[44] =  aux_X(44,i);
            // Y_n[1-1] =  aux_X(1-1,i);
            // Y_n[2-1] =  aux_X(2-1,i);
            
            // //  component-wise step for tpm2 (tpm row 2)
            // // p.increment();
            // double log_U_tpm2 = log(R::runif(0,1));

            // Y_n[3-1] =  aux_Y_n[3-1];
            // Y_n[45] =  aux_Y_n[45];
            // Y_n[4-1] =  aux_Y_n[4-1];

            // double log_alpha_tpm2 = 0.0;
            // log_alpha_tpm2 += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            // log_alpha_tpm2 -= target_at_X_i;
            // double log_ratio_tpm2 = min(NumericVector::create(0,log_alpha_tpm2));
            
            // if(log_U_tpm2 <= log_ratio_tpm2){
            //     aux_X(3-1,i+1) = aux_Y_n[3-1];
            //     aux_X(45,i+1) = aux_Y_n[45];
            //     aux_X(4-1,i+1) = aux_Y_n[4-1];
            // } else {
            //     aux_X(3-1,i+1) = aux_X(3-1,i);
            //     aux_X(45,i+1) = aux_X(45,i);
            //     aux_X(4-1,i+1) = aux_X(4-1,i);
            // }

            // // return the value we had before to recycle this vector
            // Y_n[3-1] =  aux_X(3-1,i);
            // Y_n[45] =  aux_X(45,i);
            // Y_n[4-1] =  aux_X(4-1,i);

            // //  component-wise step for tpm3 (tpm row 3)
            // // p.increment();
            // double log_U_tpm3 = log(R::runif(0,1));

            // Y_n[5-1] =  aux_Y_n[5-1];
            // Y_n[6-1] =  aux_Y_n[6-1];
            // Y_n[46] =  aux_Y_n[46];

            // double log_alpha_tpm3 = 0.0;
            // log_alpha_tpm3 += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            // log_alpha_tpm3 -= target_at_X_i;
            // double log_ratio_tpm3 = min(NumericVector::create(0,log_alpha_tpm3));
            
            // if(log_U_tpm3 <= log_ratio_tpm3){
            //     aux_X(5-1,i+1) = aux_Y_n[5-1];
            //     aux_X(6-1,i+1) = aux_Y_n[6-1];
            //     aux_X(46,i+1) = aux_Y_n[46];
            // } else {
            //     aux_X(5-1,i+1) = aux_X(5-1,i);
            //     aux_X(6-1,i+1) = aux_X(6-1,i);
            //     aux_X(46,i+1) = aux_X(46,i);
            // }

            // // return the value we had before to recycle this vector
            // Y_n[5-1] =  aux_X(5-1,i);
            // Y_n[6-1] =  aux_X(6-1,i);
            // Y_n[46] =  aux_X(46,i);

            // //  component-wise step for initial state distribution
            // // p.increment();
            // double log_U_initd = log(R::runif(0,1));

            // Y_n[42] =  aux_Y_n[42];
            // Y_n[43] =  aux_Y_n[43];
            // Y_n[47] = aux_Y_n[47];

            // double log_alpha_initd = 0.0;
            // log_alpha_initd += target_tau(list_data, Y_n, temp_vector[temp],temp_vector[temp_vector.size()-1]);
            // log_alpha_initd -= target_at_X_i;
            // double log_ratio_initd = min(NumericVector::create(0,log_alpha_initd));
            
            // if(log_U_initd <= log_ratio_initd){
            //     aux_X(42,i+1) = aux_Y_n[42];
            //     aux_X(43,i+1) = aux_Y_n[43];
            //     aux_X(47,i+1) = aux_Y_n[47];
            // } else {
            //     aux_X(42,i+1) = aux_X(42,i);
            //     aux_X(43,i+1) = aux_X(43,i);
            //     aux_X(47,i+1) = aux_X(47,i);
            // }

            // // return the value we had before to recycle this vector
            // Y_n[42] =  aux_X(42,i);
            // Y_n[43] =  aux_X(43,i);
            // Y_n[47] =  aux_X(47,i);

            // component-wise for state-dependent parameters  

            // Creation of input matrix for the parallelization 
            // (the column i is the vector Y_n with the ith value being replaced by the new proposal for parameter i)
            NumericMatrix matrix_Y_n(init.size(),init.size());
            for(int row = 0; row < init.size(); row++){
                for(int col = 0; col < init.size(); col++){
                    // matrix_Y_n(row,col) = Y_n[row];
                    matrix_Y_n(row,col) = aux_X(row,i);
                }
            }
            for(int diag = 0; diag < init.size(); diag++){
                matrix_Y_n(diag,diag) = aux_Y_n[diag];
            }

            // completition tpm1 in matrix_Y_n
            matrix_Y_n(2-1,1-1) =  aux_Y_n[2-1];
            matrix_Y_n(44,1-1) =  aux_Y_n[44];

            // completition tpm2 in matrix_Y_n
            matrix_Y_n(4-1,3-1) =  aux_Y_n[4-1];
            matrix_Y_n(45,3-1) =  aux_Y_n[45];

            // completition tpm2 in matrix_Y_n
            matrix_Y_n(6-1,5-1) =  aux_Y_n[6-1];
            matrix_Y_n(46,5-1) =  aux_Y_n[46];


            // completition init in matrix_Y_n
            matrix_Y_n(43,42) = aux_Y_n[43];
            matrix_Y_n(47,42) = aux_Y_n[47];

            // Vector of log-uniform samples to be used in the parallelization
            NumericVector vec_log_U(init.size());
            for(int j = 2; j < (init.length()-9); j++){
                vec_log_U[j] = log(R::runif(0,1));
            }

            for(int j = (init.length()-3); j < init.length(); j++){
                vec_log_U[j] = log(R::runif(0,1));
            }
            
            // Rprintf("log_U for tpm1, temp %i is %f \n",temp,vec_log_U[2]);
            // Rprintf("log_U for tpm2, temp %i is %f \n",temp,vec_log_U[3]);
            // Rprintf("log_U for tpm3, temp %i is %f \n",temp,vec_log_U[4]);
            // Rprintf("log_U for init, temp %i is %f \n",temp,vec_log_U[5]);
            // Vector where the new sample generated from the parallelization will be stored
            NumericVector output_parallel(init.size());
            // Using vector_aggregator to get new samples from state-dependent distributions
            // NumericVector aux_parallel_Y_n = vector_aggregator(matrix_Y_n,
            vector_aggregator(matrix_Y_n,
                                                                N,
                                                                n,
                                                                n_ind,
                                                                ID_init,
                                                                ID,
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
                                                                output_parallel,
                                                                aux_X(_,i),vec_log_U,target_at_X_i,temp_vector[temp],temp_vector[temp_vector.size()-1],2,init.length()-9);
            // Using vector_aggregator for the theta parameter
            // NumericVector aux_parallel_Y_n_theta = vector_aggregator(matrix_Y_n,
            vector_aggregator(matrix_Y_n,
                                                    N,
                                                    n,
                                                    n_ind,
                                                    ID_init,
                                                    ID,
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
                                                    // aux_parallel_Y_n,
                                                    output_parallel,
                                                    aux_X(_,i),vec_log_U,target_at_X_i,temp_vector[temp],temp_vector[temp_vector.size()-1],init.length()-3,init.length());

            // for(int j = 6; j < (init.length()-9); j++){
            //     aux_X(j,i+1) = aux_X(j,i);
            // }                            

            // for(int j = (init.length()-3); j < init.length(); j++){
            //     aux_X(j,i+1) = aux_X(j,i);
            // }

            // for(int j = 0; j < 6; j++){
            //     aux_X(j,i+1) = output_parallel[j];
            // }

            // for(int j = (init.length()-9); j < (init.length()-3); j++){
            //     aux_X(j,i+1) = output_parallel[j];
            // }

            

            for(int j = 0; j < init.length(); j++){
                aux_X(j,i+1) = output_parallel[j];
            }

            p.increment();

        }

        // Choose randomly one of the chains
        int j = floor(R::runif(0,1)*(temp_vector.size()-1));
        int k = j + 1; // Proposed Swap;
        double log_U_swap = log(R::runif(0,1));

        // Auxiliary matrices where we will extract the values to be swapped
        NumericMatrix X_j = X[j];
        NumericMatrix X_k = X[k];
        // double ratio = target_tau(list_data, X_j(_,i+1), k,temp_vector[temp_vector.size()-1]) + target_tau(list_data, X_k(_,i+1), j,temp_vector[temp_vector.size()-1]) - 
        //                 (target_tau(list_data, X_j(_,i+1), j,temp_vector[temp_vector.size()-1]) + target_tau(list_data, X_k(_,i+1), k,temp_vector[temp_vector.size()-1]));
        double ratio = 0.0;
        ratio += lp__stdVector(N,
                    n,
                    n_ind,
                    ID_init,
                    ID,
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
                    X_j(_,i+1))/temp_vector[k];

        ratio += lp__stdVector(N,
                    n,
                    n_ind,
                    ID_init,
                    ID,
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
                    X_k(_,i+1))/temp_vector[j];

        ratio -= lp__stdVector(N,
                    n,
                    n_ind,
                    ID_init,
                    ID,
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
                    X_j(_,i+1))/temp_vector[j];

        ratio -= lp__stdVector(N,
                    n,
                    n_ind,
                    ID_init,
                    ID,
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
                    X_k(_,i+1))/temp_vector[k];


        if(log_U_swap  < ratio){
            NumericVector tmpval = X_j(_,i+1);
            X_j(_,i+1) = X_k(_,i+1);
            X_k(_,i+1) = tmpval; // Accept Swap;
        }
    }

    return X;
} 
