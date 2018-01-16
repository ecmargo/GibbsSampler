data{
  int<lower=1> N;                               ## total number of observations
  int<lower=1> G;                               ## total number of groups
  int<lower=1> W;                               ## total number of waves/time periods
  int<lower=1> K;                               ## total number of questions
  int<lower=1> M;                               ## number of intermediate factors     
  int<lower=1> R;
  
  int<lower=0, upper=1> IM[K, M];                ## factors each question belongs to
  int<lower=1, upper=G> IG[R];                   ## group the observation belongs to
  int<lower=0, upper=R> IR[N];                   ## unique id of the respondent
  int<lower=1, upper=W> IW[R];                   ## wave the respondent belongs to
  int<lower=1, upper=K> IK[N];                   ## question being answered
  
  int<lower=0, upper=1> Y[N];                   ## answers as 0/1 choices
}

parameters{
  # Model parameters
  matrix[R, M] theta;                           ## intermediate-level factors
  vector[M] lambda;                             ## intermediate factor loadings
  row_vector[K] a;                              ## 2PL IRT parameters
  row_vector[K] b;
  matrix[G, W] T;                               ## the second level factor
  vector<lower=0>[M-1] beta;                             ## second level loading

  ## Hyperparameters
  row_vector[G] nu;
  real<lower=0> kappa;                          ## intermediate factor loading variance
  real<lower=0> tau;                            ## second level factor loading variance
  real<lower=0> phi;                            ## second level factor variance
  real<lower=0> eta;                            ## question variance
  cholesky_factor_corr[M] Psi;                  ## first level factor correlation matrix
  real<lower=0> omega;                          ## question variances
}

transformed parameters{
  matrix[N, M] theta_zero;
  cov_matrix[M] Sigma;
  vector[M] beta1;
  Sigma = quad_form_diag(multiply_lower_tri_self_transpose(Psi), rep_vector(1, M));
  for(i in 1:N){
    for(j in 1:M){
      if(IM[IK[i], j]==0)
        theta_zero[i, j] = 0;
      else
        theta_zero[i, j] = theta[IR[i], j];
    }
  }
  beta1[1] = 1;
  beta1[2:M] = beta;
}

model{
  ## additional variables
  row_vector[N] mu;
  
  ## priors
  kappa ~ cauchy(0, 1);
  tau ~ cauchy(0, 1);
  phi ~ cauchy(0, 1);
  Psi ~ lkj_corr_cholesky(1);
  omega ~ inv_gamma(0.01, 0.01);
  eta ~ inv_gamma(0.01, 0.01);
  nu ~ normal(0, 100);
  
  ## model
  beta ~ normal(0, tau);
  for(j in 1:G)
    T[j, ] ~ normal(nu[j], phi);
  a ~ normal(0, omega);
  b ~ normal(0.5, eta);
  lambda ~ normal(0, kappa);
  for(i in 1:R)
    theta[i, ] ~ multi_normal(T[IG[i], IW[i]] * beta1, Sigma);

  ##likelihood
  for(i in 1:N)
    mu[i] = b[IK[i]] * (theta_zero[i, ] * lambda - a[IK[i]]);
  
  Y ~ bernoulli_logit(mu);
}




