// #include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(roptim)]]
#include <roptim.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>


using namespace Rcpp;
using namespace roptim;
using namespace RcppParallel;

// using namespace Rcpp;

// [[Rcpp::plugins("cpp17")]] 

// The log of the Von-misses density function of x given location mu and scale kappa
double log_dvm(double x, double mu, double kappa){
    return kappa*std::cos(x-mu) -(M_LN2 + log(M_PI) + log(std::cyl_bessel_i(0,kappa)));
}


struct SemiMle : public Functor {
// class SemiMle {
    
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
        arma::Col<double> estimate;


    SemiMle (
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
        const arma::Col<double> x_headVar_in,
        const arma::Col<double> estimate_in) :
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
        x_headVar(x_headVar_in),
        estimate(estimate_in)
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

    int ll = 36 + 1;
    // # Creation parameters initial distribution
    arma::Col<double> init_raw = {theta_star[43-ll], theta_star[44-ll]}; // pi parameters

    // # weights pi_
    arma::Col<double> theta(K); //   estimate[45-ll];
    theta[2-1] = exp(theta_star[45-ll])/(1+exp(theta_star[45-ll])); // ll = 36 - 1
    theta[1-1] = 1 - theta[2-1];
    
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
arma::Col<double> semiMle_bfgs(    
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
    const arma::Col<double> estimate,
    arma::Col<double> theta_star
)
{
  SemiMle sm(
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
    x_headVar,
    estimate);

  Roptim<SemiMle> opt("BFGS");
  opt.control.trace = 1;
  // opt.set_hessian(true);
  // arma::vec x = {-1.2, 1};
  opt.minimize(sm, theta_star);
  Rcpp::Rcout << "-------------------------" << std::endl;
  opt.print();
//   return opt.value();
  return opt.par();
}

// For the implementation of RcppParallel, we need to define a Worker object
struct MleEstim : public Worker
{
    const RMatrix<double> the_init_matrix; // input matrix
    // data
    const int K;
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
    const RVector<double> estimate;

    const std::size_t size; // size of the parameter vector
    const std::size_t nsim; // size of the parameter vector
    
    RVector<double> output; // output vector 
    double mle;

    // void semiMle_bfgs(    
    double semiMle_bfgs(    
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
        const arma::Col<double> estimate,
        arma::Col<double> theta_star);

    // convert RVector/RMatrix into arma type for target_tau_ function and the following arma data will be shared in parallel computing
    arma::mat convert(){
        RMatrix<double> tmp_the_init_matrix = the_init_matrix;
        arma::mat MAT(tmp_the_init_matrix.begin(), size, nsim, false);
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

    arma::Col<double> convert_estimate(){
        RVector<double> tmp_estimate = estimate;
        arma::Col<double> vec(tmp_estimate.begin(),n,false);
        return vec;
    }

    // The main constructor
    MleEstim (
            const NumericMatrix the_init_matrix_in,
            const int K_in,
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
            const NumericVector estimate_in,
            const size_t size_in,
            const size_t nsim_in,
            NumericVector output_in) :
        the_init_matrix(the_init_matrix_in),
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
        x_headVar(x_headVar_in),
        estimate(estimate_in),
        size(size_in),
        nsim(nsim_in),
        output(output_in),
        mle(0)
    {

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
            const arma::Col<double> aux_estimate = convert_estimate();
            arma::mat MAT = convert();

            arma::Col<double> col_i = MAT.col(i);

            output[i] = semiMle_bfgs(aux_ID_init,
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
                        aux_estimate,
                        col_i);
            
            // output[i] += mle;

        }
    }

};

// void MleEstim::semiMle_bfgs(    
double MleEstim::semiMle_bfgs(    
    arma::Col<int> aux_ID_init,
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
    const arma::Col<double> aux_estimate,
    arma::Col<double> theta_star){
        SemiMle sm(
            K,
            N,
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
            aux_estimate);

        Roptim<SemiMle> opt("BFGS");
        opt.control.trace = 1;
        // opt.set_hessian(true);
        // arma::vec x = {-1.2, 1};
        opt.minimize(sm, theta_star);
        // Rcpp::Rcout << "-------------------------" << std::endl;
        // opt.print();
        //   mle += opt.value();
        return opt.value();
    }

// [[Rcpp::export]]
NumericVector parallel_mle_estimator(NumericMatrix Y,
                                const int K,
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
                                const NumericVector estimate,
                                int begin_in, int end_in){

    const std::size_t size = static_cast <size_t> (Y.nrow());
    const std::size_t nsim = static_cast <size_t> (Y.ncol());
    const std::size_t begin = static_cast <size_t> (begin_in);
    const std::size_t end = static_cast <size_t> (end_in);        
    const NumericVector output_parallel(Y.ncol());

    MleEstim mleEstim(Y,
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
                    x_headVar,
                    estimate,
                    size,
                    nsim,
                    output_parallel);

    parallelFor(begin,end,mleEstim); 

    return output_parallel;

}
