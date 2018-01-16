// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
//#include <progress.hpp>
using namespace arma;
//using namespace Rcpp;

int N = 0;   //num individuals
int W = 0;   // num of waves
int G = 0;   // num of groups
int K = 0;   // num of questions
int P = 3;   // num of factors
int N_Iterations = 0;  // num iterations
mat IM[K, M];    //factors each question belongs to
Vec IG[R];   //group the observation belongs to
Vec IR[N];   // unique id of the respondent
Vec IW[R];   // wave the respondent belongs to
Vec IK[N];  // question being answered
Vec Y[N];    // answers as 0/1 choices
  

arma::mat mvrnormArma(int n, arma::vec& mu, arma::mat& sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

double rtcauchy(double location, double scale, double a, double b){
  return(R::qcauchy(R::runif(R::pcauchy(a, location, scale, true, false), 
                             R::pcauchy(b, location, scale, true, false)),
                             location, scale, true, false));
  
}

arma::mat rwishart(unsigned int df, const arma::mat& S){
  unsigned int m = S.n_rows;
  
  arma::mat Z(m,m);
  
  for(unsigned int i = 0; i < m; i++){
    Z(i,i) = std::sqrt(R::rchisq(df-i));
  }
  
  for(unsigned int j = 0; j < m; j++){  
    for(unsigned int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1);
    }
  }
  
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  
  return C.t()*C;
}

arma::mat riwishart(unsigned int df, const arma::mat& S){
  return rwishart(df,S.i()).i();
}

double inversegamma{
}
double inversewishart{
  
}

void run_sampler(){
  
}

void initialize_parameters(mat& Sigma, mat& Ximp, vec& gamma, 
                           double& tau, vec& beta, vec& Psi, vec& lambda, vec& nu, double& xi,
                           double& sigma2, double& phi2, double& ka, double& kb, vec& a, vec& b,
                           mat& S0, double& nu0, vec& mu, vec& omega, 
                           double& a_gamma, double& b_gamma, double& a_sigma2, double& b_sigma2,
                           double& a_phi2, double& b_phi2, double& a_kappa, double& b_kappa,
                           mat& X, mat& Imis, int& N, int& P,
                           int& K, int C, int T, vec& isControl){
  
  ka = 1/R::rgamma(a_kappa, 1/b_kappa);
  kb = 1/R::rgamma(a_kappa, 1/b_kappa);
  phi2 = 1/R::rgamma(a_phi2, 1/b_phi2);
  sigma2 = 1/R::rgamma(a_sigma2, 1/b_sigma2);
  
  a = rnorm(C, 0, ka);
  b = rnorm(T, 0, kb);
  
  Sigma = riwishart(nu0, S0);
  
  Ximp = zeros(N, P);
  for (int i=0; i<N; i++){
    for(int j=0; j<P; j++){
      if(Imis(i, j) == 1){
        Ximp(i, j) = as<double>(rnorm(1, mu[j], omega[j]));
      }
      else{
        Ximp(i, j) = X(i, j);
      }
    }
  }
  
  tau = rtcauchy(0, 1, 0, 1);
  
  for (int j=0; j<P; j++){
    if(isControl(j)==1)
      gamma(j) = rtcauchy(0, tau, 0, 1);
    else
      gamma(j) = 1/as<double>(rgamma(1, a_gamma, b_gamma));
  }
  
  for(int j=0; j<P; j++){
    beta(j) = as<double>(rnorm(1, 0, gamma(j, 0)));
  }
  
  vec Z = zeros<vec>(K);
  Psi = mvrnormArma(1, Z, Sigma).t();
  Psi(0) = 1.0;
  
  for(int i=0; i<N; i++){
    mat mu = Ximp.row(i)*beta;
    lambda(i, 0) = R::rnorm(mu(0, 0), sigma2 + ka + kb);
  }
  for(int i=0; i<P; i++)
    nu(i) = 1/R::rgamma(0.5, 1);
  xi = 1/R::rgamma(0.5, 1);
  return;
}

void save_samples(mat& Sigma, mat& Ximp, vec& gamma, 
                  double& tau, vec& beta, vec& Psi, vec& lambda, vec& nu, double& xi,
                  double& sigma2, double& phi2, double& ka, double& kb, vec& a, vec& b,
                  cube& Sigma_s, cube& Ximp_s, mat& gamma_s,
                  vec& tau_s, mat& beta_s, mat& Psi_s, mat& lambda_s, mat& nu_s, vec& xi_s,
                  vec& sigma2_s, vec& phi2_s, vec& ka_s, vec& kb_s, mat& a_s, mat& b_s, int& iter){
  
  Sigma_s.slice(iter) = Sigma;
  Ximp_s.slice(iter) = Ximp;
  gamma_s.col(iter) = gamma;
  tau_s(iter) = tau;
  beta_s.col(iter) = beta;
  Psi_s.col(iter) = Psi;
  lambda_s.col(iter) = lambda;
  nu_s.col(iter) = nu;
  xi_s(iter) = xi;
  sigma2_s(iter) = sigma2;
  phi2_s(iter) = phi2;
  ka_s(iter) = ka;
  kb_s(iter) = kb;
  a_s.col(iter) = a;
  b_s.col(iter) = b;
  
  return;
}
void update_y(modelparams, secondparams, i){}
void update_miu(modelparams, secondparams, i){}
void update_ak(modelparams, secondparams, i){}
void update_bk(modelparams, secondparams, i){}
void update_lambda(modelparams, secondparams, i){}
void update_theta(modelparams, secondparams, i) {}
void update_beta(modelparams, secondparams, i){}
void update_T(modelparams, secondparams, i){}
                                                      
                                                        
void sample_parameters(mat& Sigma, mat& Ximp, vec& gamma, 
                       double& tau, vec& beta, vec& Psi, vec& lambda, vec& nu, double& xi,
                       double& sigma2, double& phi2, double& ka, double& kb, vec& a, vec& b,
                       cube& Sigma_s, cube& Ximp_s, mat& gamma_s, 
                       vec& tau_s, mat& beta_s, mat& Psi_s, mat& lambda_s, mat& nu_s, vec& xi_s,
                       vec& sigma2_s, vec& phi2_s, vec& ka_s, vec& kb_s, mat& a_s, mat& b_s,
                       mat& S0, double& nu0, vec& mu, vec& omega, 
                       double& a_gamma, double& b_gamma, double& a_sigma2, double& b_sigma2,
                       double& a_phi2, double& b_phi2, double& a_kappa, double& b_kappa,
                       mat& X, mat& Y, mat& Imis, vec& isControl, vec& IC, vec& IT,
                       int& N, int& P, int& K, int& C, int& T, int& n_burnin, int& n_sampling, bool& verbose){
  
  //Progress prog(n_burnin + n_sampling, true);
  
  for(int iter=0; iter<n_burnin + n_sampling; iter++){
  //  if(Progress::check_abort())
     // return;
   // prog.increment();
  if(verbose) {
    Rcout<<"###################################################"<<endl;
    Rcout<<"Iteration "<<iter<<endl;
  }
    update_ak(modelparams,secondparams,i)
    update_bk(modelparams,secondparams,i)
    update_lambda(modelparams,secondparams,i)
    update_beta(modelparams, secondparams, i)
    update_theta(modelparams, secondparams, i)
    update_T(modelparams, secondparams, i)
    update_y(modelparams, secondparams, i)
    update_miu(modelparams, secondparams, i)
    
    vec a_vec(N);
    vec b_vec(N);
    for(int i=0; i<N; i++){
      a_vec(i) = a(IC(i)-1);
      b_vec(i) = b(IT(i)-1);
    }
    
    //kappa a
    ka = 1/R::rgamma(a_kappa + 0.5*C, 1/(b_kappa + 0.5 * accu(a.t()*a)));
    if (verbose) Rcout<<"Done sampling ka"<<endl;
    if(verbose) Rcout<<"ka="<<ka<<endl;
    
    //kappa b
    kb = 1/R::rgamma(a_kappa + 0.5*T, 1/(b_kappa + 0.5 * accu(b.t()*b)));
    if (verbose) Rcout<<"Done sampling kb"<<endl;
    if(verbose) Rcout<<"kb="<<kb<<endl;
    
    //phi2
    // double RSS_lambda = accu((lambda - a_vec - b_vec - Ximp * beta).t() * (lambda - a_vec - b_vec - Ximp * beta));
    // if(verbose) Rcout<<"RSS_lambda="<<RSS_lambda<<endl;
    // phi2 = 1/R::rgamma(a_phi2 + 0.5*N, 1/(b_phi2 + 0.5 * RSS_lambda));
    // if (verbose) Rcout<<"Done sampling phi2"<<endl;
    // if(verbose) Rcout<<"phi2="<<phi2<<endl;
    phi2 = 10.0;
    
    //sigma2
    double RSS_y = 0;
    for(int i=0; i<N; i++){
      RSS_y += accu((Y.row(i).t() - lambda(i) * Psi).t() * (Y.row(i).t() - lambda(i)*Psi));
    }
    if(verbose) Rcout<<"RSS_Y="<<RSS_y<<endl;
    sigma2 = 1/R::rgamma(a_sigma2 + 0.5*N, 1/(b_sigma2 + 0.5 * RSS_y));
    if (verbose) Rcout<<"Done sampling Sigma2"<<endl;
    if(verbose) Rcout<<"sigma2="<<sigma2<<endl;
    
    //a 
    double VaI = 1/(C/phi2 + 1/ka);
    vec XBb = lambda - Ximp*beta - b_vec;
    for(int c=1; c<=C; c++){
      double Ma = (1/phi2)*accu(XBb.elem(find(IC==c)));
      a(c-1) = R::rnorm(VaI*Ma, VaI);
    }
    if (verbose) Rcout<<"Done sampling a"<<endl;
    if(verbose) Rcout<<"a: "<<"mean="<<mean(a)
                     <<", var="<<var(a)
                     <<", min="<<min(a)
                     <<", max="<<max(a)<<endl;
    
    //b 
    double VbI = 1/(T/phi2 + 1/kb);
    vec XBa = lambda - Ximp * beta - a_vec;
    for(int t=1; t<=T; t++){
      double Mb = (1/phi2)*accu(XBa.elem(find(IT==t)));
      b(t-1) = R::rnorm(VbI * Mb, VbI);
    }
    if (verbose) Rcout<<"Done sampling b"<<endl;
    if(verbose) Rcout<<"b: "<<"mean="<<mean(b)
                     <<", var="<<var(b)
                     <<", min="<<min(b)
                     <<", max="<<max(b)<<endl;
    
    // missing Xs
    for(int i=0; i<N; i++){
      for(int j=0; j<P; j++){
        if (Imis(i,j)==0){
          Ximp(i,j) = X(i, j);
        }else{
          double vj = std::pow(beta(j), 2) * (1/phi2) + 1/omega(j);
          mat mj = (beta(j)/phi2) * 
            (lambda(i) - Ximp.row(i)*beta - a_vec(i) - b_vec(i) +
            Ximp(i,j)*beta(j)) + (mu[j]/omega[j]);
          mj(0,0) = mj(0,0) * 1/vj;
          Ximp(i, j) = as<double>(rnorm(1, mj(0,0), 1/vj));
        }
      }
    }
    if(verbose){
      Rcout<<"Done sampling X"<<endl;
      for(int j=1; j< P; j++){
        Rcout<<"Ximp col "<<j<<": "<<"mean="<<mean(Ximp.col(j))
             <<", var="<<var(Ximp.col(j))
             <<", min="<<min(Ximp.col(j))
             <<", max="<<max(Ximp.col(j))<<endl;
      }
    }
    
    // nu
    for(int j=0; j<P; j++){
      nu(j) = 1/R::rgamma(1, 1/(1 + 1/gamma(j)));
    }
    if(verbose) Rcout<<"Done sampling nu"<<endl; 
    if(verbose) Rcout<<"nu: "<<"mean="<<mean(nu)
                     <<", var="<<var(nu)
                     <<", min="<<min(nu)
                     <<", max="<<max(nu)<<endl;
    
    //xi
    xi = 1/R::rgamma(1, 1/(1 + 1/tau));
    if(verbose) Rcout<<"Done sampling xi"<<endl;
    if(verbose) Rcout<<"xi="<<xi<<endl;
    
    // tau
    double BLS = 0;
    for(int j=0; j<P; j++)
      BLS += std::pow(beta(j), 2) * gamma(j);
    tau = 1/R::rgamma((P+1)/2, 1/(1/xi + 0.5 * BLS));
    if(verbose) Rcout<<"Done sampling tau"<<endl;
    if(verbose) Rcout<<"tau="<<tau<<endl;
    
    // gamma
    for(int j=0; j<P; j++){
      if(isControl(j)==1){
        gamma(j) = 1/R::rgamma(1, 1/(1/nu(j) + std::pow(beta(j), 2)/(2*tau)));
      }else{
        gamma(j) = 1/R::rgamma(a_gamma + 0.5, 1/(b_gamma + 0.5 * std::pow(beta(j), 2)));
      }
    }
    if(verbose) Rcout<<"Done sampling gamma"<<endl;
    if(verbose) Rcout<<"gamma: "<<"mean="<<mean(gamma)
                     <<", var="<<var(gamma)
                     <<", min="<<min(gamma)
                     <<", max="<<max(gamma)<<endl;
    
    // Sigma
    mat S1 = Psi * Psi.t() + S0;
    Sigma = riwishart(nu0 + 1, inv(S1));
    if(verbose) Rcout<<"Done sampling Sigma"<<endl;
    if(verbose) Sigma.print();
    
    //Psi
    mat Vpsi = inv(Sigma) + (1/sigma2)*accu((lambda.t()*lambda));
    mat VpsiI = inv(Vpsi);
    if(verbose) Rcout<<"VpsiI"<<endl;
    if(verbose) VpsiI.print();
    vec mpsi = VpsiI * ((1/sigma2)*lambda.t()*Y).t();
    Psi = mvrnormArma(1, mpsi, VpsiI).t();
    
    //Set the loading on GER to 0 for identification
    Psi(0) = 0.0;
    
    if(verbose) Rcout<<"Done sampling Psi"<<endl;
    if(verbose) Psi.print();
    
    //Lambda
    mat Vlambda = (1/phi2) + (1/sigma2)*Psi.t()*Psi;
    double VlambdaI = 1/Vlambda(0,0);
    for(int i=0; i<N; i++){
      double mlambda = accu((1/phi2)*(Ximp.row(i)*beta + a_vec(i) + b_vec(i)) + 
                            (1/sigma2)*Psi.t()*Y.row(i).t());
      lambda(i) = R::rnorm(VlambdaI * mlambda, VlambdaI);
    }
    if(verbose) Rcout<<"Done sampling lambda"<<endl;
    if(verbose) Rcout<<"lambda: "<<"mean="<<mean(lambda)
                     <<", var="<<var(lambda)
                     <<", min="<<min(lambda)
                     <<", max="<<max(lambda)<<endl;
    
    //Beta
    vec gammaTau(P);
    for(int j=0; j<P; j++){
      if(isControl[j]==1)
        gammaTau(j) = gamma(j) * tau;
      else
        gammaTau(j) = gamma(j);
    }
    
    mat Gamma = diagmat(gammaTau) * eye(P, P)*0.0001;
    mat GammaI = inv(Gamma);
    mat Vbeta = GammaI + (1/phi2) * Ximp.t() * Ximp;
    mat VbetaI = inv(Vbeta);
    vec mbeta = VbetaI * ((1/phi2) * Ximp.t() * (lambda - a_vec - b_vec));
    beta = mvrnormArma(1, mbeta, VbetaI).t();
    
    if(verbose) Rcout<<"Done sampling beta"<<endl;
    if(verbose) Rcout<<"beta: "<<"mean="<<mean(beta)
                     <<", var="<<var(beta)
                     <<", min="<<min(beta)
                     <<", max="<<max(beta)<<endl;
    
    // Save samples
    if(iter >= n_burnin){
      int is = iter - n_burnin;
      save_samples(Sigma, Ximp, gamma, tau, beta, Psi, lambda, nu, xi,
                   sigma2, phi2, ka, kb, a, b, Sigma_s, Ximp_s, gamma_s, 
                   tau_s, beta_s, Psi_s, lambda_s, nu_s, xi_s,
                   sigma2_s,  phi2_s, ka_s, kb_s, a_s, b_s, is);
    }
    
  }
  return;
}


//[[Rcpp::export]]
List sample_two_step(mat X, mat Y, vec IC, vec IT, mat Imis, vec isControl, 
                     int n_burnin, int n_sampling, bool verbose, List hyperparameters){
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Y.n_cols;
  int C = max(IC);
  int T = max(IT);
  
  //Parameter storage declarations
  cube Sigma_s(K, K, n_sampling);
  cube Ximp_s(N, P, n_sampling);
  mat gamma_s(P, n_sampling);
  vec tau_s(n_sampling);
  mat beta_s(P, n_sampling);
  mat Psi_s(K, n_sampling);
  mat lambda_s(N, n_sampling);
  mat nu_s(P, n_sampling);
  vec xi_s(n_sampling);
  vec sigma2_s(n_sampling);
  vec phi2_s(n_sampling);
  vec ka_s(n_sampling);
  vec kb_s(n_sampling);
  mat a_s(C, n_sampling); 
  mat b_s(T, n_sampling);
  
  // Burnin samples
  mat Sigma(K, K);
  mat Ximp(N, P);
  vec gamma(P);
  double tau;
  vec beta(P);
  vec Psi(K);
  vec lambda(N);
  vec nu(P);
  double xi(P);
  double sigma2;
  double phi2;
  double ka;
  double kb;
  vec a(C);
  vec b(T);
  
  // get hyperparameters out of list
  mat S0 = as<mat>(hyperparameters["S0"]);
  double nu0 = hyperparameters["nu0"];
  vec mu = as<vec>(hyperparameters["mu"]);
  vec omega = as<vec>(hyperparameters["omega"]);
  double a_gamma = hyperparameters["a_gamma"];
  double b_gamma = hyperparameters["b_gamma"];
  double a_sigma2 = hyperparameters["a_sigma2"];
  double b_sigma2 = hyperparameters["b_sigma2"];
  double a_phi2 = hyperparameters["a_phi2"];
  double b_phi2 = hyperparameters["b_phi2"];
  double a_kappa = hyperparameters["a_kappa"];
  double b_kappa = hyperparameters["b_kappa"];
  
  Rcout<<"Initializing Parameters..."<<endl;
  initialize_parameters(Sigma, Ximp, gamma, tau, beta, Psi, lambda, nu, xi,
                        sigma2, phi2, ka, kb, a, b,
                        S0, nu0, mu, omega, a_gamma, b_gamma, a_sigma2, b_sigma2,
                        a_phi2, b_phi2, a_kappa, b_kappa,
                        X, Imis, N, P, K, C, T, isControl);  
  
  Rcout<<"Sampling..."<<endl;
  sample_parameters(Sigma, Ximp, gamma, tau, beta, Psi, lambda, nu, xi, 
                    sigma2, phi2, ka, kb, a, b,
                    Sigma_s, Ximp_s, gamma_s, tau_s, beta_s, Psi_s, lambda_s, nu_s, xi_s, 
                    sigma2_s, phi2_s, ka_s, kb_s, a_s, b_s,
                    S0, nu0, mu, omega, a_gamma, b_gamma, a_sigma2, b_sigma2,
                    a_phi2, b_phi2, a_kappa, b_kappa,
                    X, Y, Imis, isControl, IC, IT, N, P, K, C, T, n_burnin, n_sampling, verbose);
  
  
  return List::create(_["Sigma"]=Sigma_s, 
                      _["Ximp"]=Ximp_s, 
                      _["gamma"]=gamma_s, 
                      _["tau"]=tau_s, 
                      _["beta"]=beta_s, 
                      _["Psi"]=Psi_s, 
                      _["lambda"]=lambda_s,
                      _["sigma2"]=sigma2_s,
                      _["phi2"]=phi2_s,
                      _["ka"]=ka_s,
                      _["kb"]=kb_s);
}