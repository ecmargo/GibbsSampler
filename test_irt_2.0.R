library(MCMCpack)
#library(rstan)
#rstan_options(auto_write=T)
#options(mc.cores = parallel::detectCores())

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

sigmoid <- function(x){
  1/(1 + exp(-x))
}


mse <- function(yhat, y, root=T){
  if(root)
    sqrt(mean((yhat - y)^2))
  else
    mean((yhat - y)^2)
}

g = 10
w = 5
k = 10
m = 3
ng = 10
r = 500
n = g * w * k * ng

M = matrix(NA, k, m)
for(j in 1:k){
  M[j, ] = sample(0:1, m, replace=TRUE, prob = c(0.25, 0.75))
}

W = rep(1:w, each=ng*g)
K = rep(1:k, ng*g*w)
R = rep(rep(1:(ng*g*w), each=k))
G = rep(1:g, each=ng*w)

true_kappa = rinvgamma(1, 5, 5)
true_tau = rinvgamma(1, 5, 5)
true_phi = rinvgamma(1, 5, 5)
true_Sigma = Posdef(m)
true_Omega = Posdef(k)
true_beta = rnorm(m, 0, true_tau)
true_nu = rnorm(g, 0, 10)
true_T = matrix(NA, g, w)
for(j in 1:g)
  true_T[j, ] = rnorm(w, true_nu[j], true_phi)
true_a = mvrnorm(1, rep(0,k), true_Omega)
true_b = rinvgamma(k, 5, 5)
true_lambda = rnorm(m, 0, true_kappa)
true_theta = matrix(NA, r, m)
for(i in 1:r){
  true_theta[i, ] = mvrnorm(1, true_T[G[i], W[i]] %*% true_beta, true_Sigma)
  true_theta[i, !M[G[i], ]] = 0
}
Y = rep(NA, n)
for(i in 1:n){
  mu = sigmoid((true_b[K[i]]) * (t(true_theta[R[i], ]) %*% true_lambda  - true_a[K[i]]))
  Y[i] = rbinom(1, 1, mu)
}
write.csv(data.frame(Y, W, K, R, G), '~/Desktop/School/Bayesian/vectors.csv')
write.csv(data.frame(M), '~/Desktop/School/Bayesian/matrix.csv')

#model = stan_model(file = "irt_factor_2.0.stan")
#fit = sampling(model, data=list(N=n, G=g, W=w, K=k, R=r, M=m, IG=G, IW=W, IK=K, IR=R, IM=M, Y=Y), chains=4,
               #iter=2000, verbose=TRUE)

#save.image(test_2.0.RData)
