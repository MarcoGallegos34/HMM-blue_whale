// #include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(roptim)]]
#include <roptim.h>

using namespace roptim;

// using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]] 

// The log of the Von-misses density function of x given location mu and scale kappa
double log_dvm(double x, double mu, double kappa){
    return kappa*std::cos(x-mu) -(M_LN2 + log(M_PI) + log(std::cyl_bessel_i(0,kappa)));
}


struct Mle : public Functor {
    
    // private:
    public:
        // data
        int K;
        int N;
        int n;
        int n_ind;
        arma::Col<int> ID_init;
        arma::Col<int> ID;
        arma::Col<double> x_duration_init;
        arma::Col<double> x_surface_init;
        arma::Col<double> x_maxDepth_init;
        arma::Col<int> x_lunges_init;
        arma::Col<double> x_step_init;
        arma::Col<double> x_angle_init;
        arma::Col<double> x_headVar_init;
        arma::Col<double> x_duration;
        arma::Col<double> x_surface;
        arma::Col<double> x_maxDepth;
        arma::Col<int> x_lunges;
        arma::Col<double> x_step;
        arma::Col<double> x_angle;
        arma::Col<double> x_headVar;
        // arma::Col<double> estimate;


    Mle (
        const int K_in,
        const int N_in,
        const int n_in,
        const int n_ind_in,
        arma::Col<int> ID_init_in,
        const arma::Col<int> ID_in,
        const arma::Col<double> x_duration_init_in,
        const arma::Col<double> x_surface_init_in,
        const arma::Col<double> x_maxDepth_init_in,
        const arma::Col<int> x_lunges_init_in,
        const arma::Col<double> x_step_init_in,
        const arma::Col<double> x_angle_init_in,
        const arma::Col<double> x_headVar_init_in,
        const arma::Col<double> x_duration_in,
        const arma::Col<double> x_surface_in,
        const arma::Col<double> x_maxDepth_in,
        const arma::Col<int> x_lunges_in,
        const arma::Col<double> x_step_in,
        const arma::Col<double> x_angle_in,
        const arma::Col<double> x_headVar_in) :
        K(K_in),
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
        x_headVar(x_headVar_in)
        {

        }

  double operator()(const arma::Col<double> &theta_star) override {
    
    arma::Col<double> mu_duration(N);
    arma::Col<double> log_sigma_duration(N);
    arma::Col<double> mu_surface(N);
    arma::Col<double> log_sigma_surface(N);
    arma::Col<double> mu_maxDepth(N);
    arma::Col<double> log_sigma_maxDepth(N);
    arma::Col<double> mu_step(N);
    arma::Col<double> log_sigma_step(N);
    arma::Col<double> log_kappa(N);
    arma::Col<double> log_a(N);
    arma::Col<double> log_b(N);
    arma::Col<double> log_lambda(N);

    for(int i = 0; i<3; ++i){
        mu_duration[i] = theta_star[i+6];
        log_sigma_duration[i] = theta_star[i+9];

        mu_surface[i] = theta_star[i+12];
        log_sigma_surface[i] = theta_star[i+15];

        mu_maxDepth[i] = theta_star[i+18];
        log_sigma_maxDepth[i] = theta_star[i+21];

        mu_step[i] = theta_star[i+24];
        log_sigma_step[i] = theta_star[i+27];

        log_kappa[i] = theta_star[i+30];

        log_a[i] = theta_star[i+33];
        log_b[i] = theta_star[i+36];

        log_lambda[i] = theta_star[i + 39];

    }
    
    arma::Col<double> alpha_duration(N);
    // arma::Col<double> beta_duration(N);
    arma::Col<double> theta_duration(N);

    // # parameters for surface distribution
    arma::Col<double> alpha_surface(N);
    // arma::Col<double> beta_surface(N);
    arma::Col<double> theta_surface(N);
    
    // # parameters for maxDepth distribution
    arma::Col<double> alpha_maxDepth(N);
    // arma::Col<double> beta_maxDepth(N);
    arma::Col<double> theta_maxDepth(N);

    // # parameters for step distribution
    arma::Col<double> alpha_step(N);
    // arma::Col<double> beta_step(N);
    arma::Col<double> theta_step(N);
    
    // # parameters for angle distribution
    arma::Col<double> a(N); 
    arma::Col<double> b(N); 
    
    // # parameters for angle distribution
    arma::Col<double> kappa(N); 
    
    // #parameters for lunges distribution
    arma::Col<double> lambda(N); 

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
    double tpm[K][N][N];
    
    // from here the values will be extracted from theta_star, not estimate

    // int ll = 36 + 1;
    int ll = 1;
    // # Creation parameters initial distribution
    arma::Col<double> init_raw = {theta_star[43-ll], theta_star[44-ll]}; // pi parameters

    // # weights pi_
    arma::Col<double> theta(K); //   estimate[45-ll];
    if(K == 2){
        theta[2-1] = exp(theta_star[45-ll])/(1+exp(theta_star[45-ll])); // ll = 36 - 1
        theta[1-1] = 1 - theta[2-1];

    }

  if(K == 3){
    arma::Col<double> theta_raw2 = {theta_star[45-ll],theta_star[54-ll]};
    theta[2-1] = 1/(1+exp(-theta_raw2[1-1]));
    theta[3-1] = 1/(1+exp(-theta_raw2[2-1]));
    theta[1-1] = 1 - (theta[2-1] + theta[3-1]);
    
  }
  
  if(K == 4){
    arma::Col<double> theta_raw2 = {theta_star[45-ll],theta_star[54-ll],theta_star[63-ll]};
    theta[2-1] = 1/(1+exp(-theta_raw2[1-1]));
    theta[3-1] = 1/(1+exp(-theta_raw2[2-1]));
    theta[4-1] = 1/(1+exp(-theta_raw2[3-1]));
    theta[1-1] = 1 - (theta[2-1] + theta[3-1] + theta[4-1]);
    
  }
  
  if(K == 5){
    arma::Col<double> theta_raw2 = {theta_star[45-ll],theta_star[54-ll],theta_star[63-ll],theta_star[72-ll]};
    theta[2-1] = 1/(1+exp(-theta_raw2[1-1]));
    theta[3-1] = 1/(1+exp(-theta_raw2[2-1]));
    theta[4-1] = 1/(1+exp(-theta_raw2[3-1]));
    theta[5-1] = 1/(1+exp(-theta_raw2[4-1]));
    theta[1-1] = 1 - (theta[2-1] + theta[3-1] + theta[4-1] + theta[5-1]);
    
  }
    
    // return theta;
    
    // # init
    arma::mat init(K,N);
    init(1-1,2-1) = exp(init_raw[1-1])/(1+exp(init_raw[1-1]));
    init(1-1,3-1) = exp(init_raw[2-1])/(1+exp(init_raw[2-1]));
    init(1-1,1-1) = 1 - (init(1-1,2-1) + init(1-1,3-1));
    
    
    if(K == 2){
        arma::Col<double> init_raw2 = {theta_star[46-ll], theta_star[47-ll]};
        // NumericVector init_raw2(N-1); //[(46-ll):(47-ll)];
        // init_raw2 = {theta_star[46-ll], theta_star[47-ll]};
        init(2-1,2-1) = exp(init_raw2[1-1])/(1+exp(init_raw2[1-1]));
        init(2-1,3-1) = exp(init_raw2[2-1])/(1+exp(init_raw2[2-1]));
        init(2-1,1-1) = 1 - (init(2-1,2-1) + init(2-1,3-1));
        
    }

    if(K == 3){
        
        arma::Col<double> init_raw2 = {theta_star[46-ll],theta_star[47-ll]};
        init(2-1,2-1) = 1/(1+exp(-init_raw2[1-1]));
        init(2-1,3-1) = 1/(1+exp(-init_raw2[2-1]));
        init(2-1,1-1) = 1 - (init(2-1,2-1) + init(2-1,3-1));
        
        
        arma::Col<double> init_raw3 = {theta_star[55-ll],theta_star[56-ll]};
        init(3-1,2-1) = 1/(1+exp(-init_raw3[1-1]));
        init(3-1,3-1) = 1/(1+exp(-init_raw3[2-1]));
        init(3-1,1-1) = 1 - (init(3-1,2-1) + init(3-1,3-1));
        
    }
    
    if(K == 4){
        
        arma::Col<double> init_raw2 = {theta_star[46-ll],theta_star[47-ll]};
        init(2-1,2-1) = 1/(1+exp(-init_raw2[1-1]));
        init(2-1,3-1) = 1/(1+exp(-init_raw2[2-1]));
        init(2-1,1-1) = 1 - (init(2-1,2-1) + init(2-1,3-1));
        
        
        arma::Col<double> init_raw3 = {theta_star[55-ll],theta_star[56-ll]};
        init(3-1,2-1) = 1/(1+exp(-init_raw3[1-1]));
        init(3-1,3-1) = 1/(1+exp(-init_raw3[2-1]));
        init(3-1,1-1) = 1 - (init(3-1,2-1) + init(3-1,3-1));
        
        arma::Col<double> init_raw4 = {theta_star[64-ll],theta_star[65-ll]};
        init(4-1,2-1) = 1/(1+exp(-init_raw4[1-1]));
        init(4-1,3-1) = 1/(1+exp(-init_raw4[2-1]));
        init(4-1,1-1) = 1 - (init(4-1,2-1) + init(4-1,3-1));
        
    }
    
    if(K == 5){
        
        arma::Col<double> init_raw2 = {theta_star[46-ll],theta_star[47-ll]};
        init(2-1,2-1) = 1/(1+exp(-init_raw2[1-1]));
        init(2-1,3-1) = 1/(1+exp(-init_raw2[2-1]));
        init(2-1,1-1) = 1 - (init(2-1,2-1) + init(2-1,3-1));
        
        
        arma::Col<double> init_raw3 = {theta_star[55-ll],theta_star[56-ll]};
        init(3-1,2-1) = 1/(1+exp(-init_raw3[1-1]));
        init(3-1,3-1) = 1/(1+exp(-init_raw3[2-1]));
        init(3-1,1-1) = 1 - (init(3-1,2-1) + init(3-1,3-1));
        
        arma::Col<double> init_raw4 = {theta_star[64-ll],theta_star[65-ll]};
        init(4-1,2-1) = 1/(1+exp(-init_raw4[1-1]));
        init(4-1,3-1) = 1/(1+exp(-init_raw4[2-1]));
        init(4-1,1-1) = 1 - (init(4-1,2-1) + init(4-1,3-1));
        
        arma::Col<double> init_raw5 = {theta_star[73-ll],theta_star[74-ll]};
        init(5-1,2-1) = 1/(1+exp(-init_raw5[1-1]));
        init(5-1,3-1) = 1/(1+exp(-init_raw5[2-1]));
        init(5-1,1-1) = 1 - (init(5-1,2-1) + init(5-1,3-1));
        
    }
    
    
    double log_lik_total = 0.0;
    // NumericVector log_lik_total(n_ind);
    
    // #Creation tpm k = 1
    // Here's the problem!
    tpm[1-1][1-1][1-1] = 1 - (exp(theta_star[1-1]) + exp(theta_star[2-1]))/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
    tpm[1-1][1-1][2-1] = exp(theta_star[1-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
    tpm[1-1][1-1][3-1] = exp(theta_star[2-1])/(1 + (exp(theta_star[1-1]) + exp(theta_star[2-1])));
    
    tpm[1-1][2-1][1-1] = exp(theta_star[3-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
    tpm[1-1][2-1][2-1] = 1 - (exp(theta_star[3-1]) + exp(theta_star[4-1]))/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
    tpm[1-1][2-1][3-1] = exp(theta_star[4-1])/(1 + (exp(theta_star[3-1]) + exp(theta_star[4-1])));
    
    tpm[1-1][3-1][1-1] = exp(theta_star[5-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
    tpm[1-1][3-1][2-1] = exp(theta_star[6-1])/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));
    tpm[1-1][3-1][3-1] = 1 - (exp(theta_star[5-1]) + exp(theta_star[6-1]))/(1 + (exp(theta_star[5-1]) + exp(theta_star[6-1])));

    if(K == 2){
        // #Creation tpm k = 2
        int l = 47 - ll; // ll = 37
        tpm[2-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        tpm[2-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        tpm[2-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
        
        tpm[2-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
        tpm[2-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
        tpm[2-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
        
        tpm[2-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
        tpm[2-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
        tpm[2-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    

    }  

  if(K == 3){

    // Creation tpm k = 2
    int l = 47 - ll; // ll = 37
    tpm[2-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[2-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[2-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[2-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[2-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[2-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[2-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[2-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[2-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
    // Creation tpm k = 3
    // l = 47 + 9 - ll 
    l += 9; 
    tpm[3-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[3-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[3-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[3-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[3-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[3-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[3-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[3-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[3-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
  }

  if(K == 4){

    // Creation tpm k = 2
    int l = 47 - ll; // ll = 37
    tpm[2-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[2-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[2-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[2-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[2-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[2-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[2-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[2-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[2-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
    // Creation tpm k = 3
    // l = 47 + 9 - ll 
    l += 9; 
    tpm[3-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[3-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[3-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[3-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[3-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[3-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[3-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[3-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[3-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
    // Creation tpm k = 4
    // l = 47 + 18 - ll # 
    l += 9;  
    tpm[4-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[4-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[4-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[4-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[4-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[4-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[4-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[4-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[4-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
  }

  if(K == 5){

    // Creation tpm k = 2
    int l = 47 - ll; // ll = 37
    tpm[2-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[2-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[2-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[2-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[2-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[2-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[2-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[2-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[2-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
    // Creation tpm k = 3
    // l = 47 + 9 - ll 
    l += 9; 
    tpm[3-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[3-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[3-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[3-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[3-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[3-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[3-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[3-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[3-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
    // Creation tpm k = 4
    // l = 47 + 18 - ll # 
    l += 9;  
    tpm[4-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[4-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[4-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[4-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[4-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[4-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[4-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[4-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[4-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
    // Creation tpm k = 5
    // l = 47 + 27 - ll # 
    l += 9; 
    tpm[5-1][1-1][1-1] = 1 - (exp(theta_star[1+l]) + exp(theta_star[2+l]))/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[5-1][1-1][2-1] = exp(theta_star[1+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    tpm[5-1][1-1][3-1] = exp(theta_star[2+l])/(1 + (exp(theta_star[1+l]) + exp(theta_star[2+l])));
    
    tpm[5-1][2-1][1-1] = exp(theta_star[3+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    tpm[5-1][2-1][2-1] = 1 - (exp(theta_star[3+l]) + exp(theta_star[4+l]))/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l]))); 
    tpm[5-1][2-1][3-1] = exp(theta_star[4+l])/(1 + (exp(theta_star[3+l]) + exp(theta_star[4+l])));
    
    tpm[5-1][3-1][1-1] = exp(theta_star[5+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[5-1][3-1][2-1] = exp(theta_star[6+l])/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));
    tpm[5-1][3-1][3-1] = 1 - (exp(theta_star[5+l]) + exp(theta_star[6+l]))/(1 + (exp(theta_star[5+l]) + exp(theta_star[6+l])));    
    
  }

    // Environment pkg = Environment::namespace_env("CircStats");

    // Function dvm = pkg["dvm"];  

    arma::mat log_alpha_components(K,n_ind);
    
    for(int k = 0;  k < K; k++){
        
        arma::mat log_alpha(N,n_ind);

        for(arma::Col<int>::iterator i = ID_init.begin(); i != ID_init.end(); ++i){
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
            // arma::Col<double> aux_dvm = dvm(x_angle_init[*i-1] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
            // log_alpha(state,*i-1) += log(*(aux_dvm.begin()));
            log_alpha(state,*i-1) += log_dvm(x_angle_init[*i-1], 0.0, kappa[state]);


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
        
        for(int i = 0; i < n; i++){
            
            arma::Col<double> aux_log_alpha(N);
        
            for(int state = 0; state < N; state++){
                
                arma::Col<double> aux(N);
                
                // aux = log_alpha(_,ID[i]-1);
                for(int k = 0; k < N; k++){
                    aux[k] = log_alpha(k,ID[i]-1);

                }

                
                for(int previous_state = 0; previous_state < N; previous_state++){
                
                aux[previous_state] += log(tpm[k][previous_state][state]); // if I use the same matrix, I get same results, something is weird...
                
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
                    // arma::Col<double> aux_dvm = dvm(x_angle[i] , _["mu"] = 0.0, _["kappa"] = kappa[state]); 
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
                
                aux_log_alpha[state] = max(aux) + log(sum(exp(aux-max(aux))));
                
            }
            
            // log_alpha(_,ID[i]-1) = aux_log_alpha;
            for(int j = 0 ; j < N; j++){
                log_alpha(j,ID[i]-1) = aux_log_alpha[j];
            }
        }
        
        for(int i=0; i < n_ind; i++){
            log_alpha_components(k,i) = log(theta[k]) + max(log_alpha.col(i)) + log(sum(exp(log_alpha.col(i)-max(log_alpha.col(i)))));
        }
    }
    
    for(int i=0; i < n_ind; i++){
        log_lik_total += max(log_alpha_components.col(i)) + log(sum(exp(log_alpha_components.col(i)-max(log_alpha_components.col(i)))));
    }
    
    return -log_lik_total;
  }
};

// [[Rcpp::export]]
arma::Col<double> mle_extended_bfgs(    
    const int K,
    const int N,
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
    arma::Col<double> theta_star
)
{
  Mle sm(
    K,
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
    x_headVar);

  Roptim<Mle> opt("BFGS");
  opt.control.trace = 1;
  // opt.set_hessian(true);
  // arma::vec x = {-1.2, 1};
  opt.minimize(sm, theta_star);
  Rcpp::Rcout << "-------------------------" << std::endl;
  opt.print();
//   return opt.value();
  return opt.par();
}