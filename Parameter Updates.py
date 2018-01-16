import numpy as np
from scipy.stats import invgamma, invwishart, norm, truncnorm
import pandas as pd
import sys
from numba import jit

N = 0   # num individuals
W = 0   # num of waves
G = 0   # num of groups
K = 0   # num of questions
P = 3   # num of factors
N_Iterations = 10    # num iterations

Y = [] #Y[N] answers as 0/1 choices
IG = [] #IG[N] group the observation belongs to
IR = [] #IR[N] unique id of the respondent
IW = [] #IW[N] wave the respondent belongs to
IK = [] #IK[N] question being answered
IM = [] #IM[K, M] factors each question belongs to
NGW = [] #Number of individuals in group G at Wave 0

burnpoint = 5


#parameters
Y_iu = []
Miu =[]
a_k = []
b_k = []
lambda_p = []
theta_i = []
beta_p = []
T_gw = []

omega_sqr = []
tau_sqr = []
sigma = []
zeta_sqr = []
sigma_sqr = []
nu_g = []
phi_sqr = []
kappa_sqr = []

#Hyperparameters
#can expand these later to be specific
hyper_a = 5
hyper_b = 5
h_sqr = 1
s_0 = P
S_0 = np.identity(P)
#positive definite matrix will need to create later


def update_omega(it):  #DONE
    mean = hyper_a + (K/2)
    var = hyper_b + .5 * np.sum(x ** 2 for x in a_k[it])
    return invgamma.rvs(mean, var)

def update_kappa(it): #DONE
    mean  = hyper_a + (P/2)
    var = hyper_b + .5*np.sum(x**2 for x in lambda_p[it])
    return invgamma.rvs(mean,var)

def update_tau(it): #DONE
    mean = hyper_a + (K/2)
    var = hyper_b + .5 * np.sum(x ** 2 for x in b_k[it])
    return invgamma.rvs(mean, var)

def update_zeta(it): #DONE
    mean = hyper_a + (P/ 2)
    var = hyper_b + .5 * np.sum(x ** 2 for x in beta_p[it])
    return invgamma.rvs(mean, var)

def update_sigmasqr(it): #DONE --> NC
    mean = hyper_a + ((N*K) / 2)
    miusum = 0
    for i in range(0,N):
        for u in range(0,K):
            miusum += (Miu[it][i,u]-b_k[it][u]*(np.dot(lambda_p[it],np.squeeze(np.asarray(theta_i[it][i])))-a_k[it][u]))**2
    var = hyper_b + .5*miusum
    return invgamma.rvs(mean, var)

def update_phi(it): #DONE --> NC
    mean = hyper_a + ((G * W) / 2)
    tsum = 0
    for g in range(0,G):
        for w in range(0,W):
            tsum += ((T_gw[it][g,w]-nu_g[it][g]) **2)
    var = hyper_b + .5*(tsum)
    return invgamma.rvs(mean, var)

def update_sigma(it): #DONE
    mean = s_0 + N
    matsum = 0
    for i in range(0,N):
        thetaibeta = np.matrix(theta_i[it][i]-(beta_p[it]*T_gw[it][IG[i]-1,IW[i]-1]))
        matsum += np.matmul(thetaibeta.transpose(), thetaibeta)
    var = np.linalg.inv(np.matrix(S_0 + matsum))
    return invwishart.rvs(mean, var)

def update_nu(it): #DONE
    var = (1 / phi_sqr[it]) + (1 / h_sqr)
    mean = (1/var) * (1/phi_sqr[it])
    return np.array([norm.rvs(mean * np.sum(T_gw[it][g]), 1/var) for g in range(0,G)])

def update_miu(it): #DONE
    update = np.copy(Miu[it])
    for i in range(0,N):
        for u in range(0, K):
            mean = b_k[it][u]*(np.dot(lambda_p[it],theta_i[it][i])-a_k[it][u])
            var = sigma_sqr[it]
            if Y_iu[0][i, u] == 1:
                update[i, u] = truncnorm.rvs(-np.inf, 0, mean, var)
            else:
                update[i,u] = truncnorm.rvs(0,np.inf,mean,var)
    return update

def update_ak(it): #DONE
    update = np.copy(a_k[it])
    for k in range(0,K):
        var = 1/((N*b_k[it][k])**2/sigma_sqr[it] + omega_sqr[it]**-1)
        aksum = 0
        for i in range(0,N):
            aksum += (b_k[it][k]*np.dot(lambda_p[it],np.squeeze(np.asarray(theta_i[it][i])))- Miu[it][i,k])
        mean = (var*b_k[it][k]*aksum)/sigma_sqr[it]
        update[k] = norm.rvs(mean, var)
    return update

def update_bk(it): #DONE --> NC
    update = np.copy(b_k[it])
    for k in range(0, K):
        varsum = 0
        meansum = 0
        for i in range(0, N):
            varsum += (np.dot(lambda_p[it],np.squeeze(np.asarray(theta_i[it][i])))-a_k[it][k])**2
            meansum += Miu[it][i, k]*(np.dot(lambda_p[it],np.squeeze(np.asarray(theta_i[it][i])))-a_k[it][k])
        var = (tau_sqr[it]**-1 + varsum/ sigma_sqr[it]) ** -1
        mean = meansum / sigma_sqr[it]
        update[k] = norm.rvs(mean, var)
    return update

def update_lambda(it): #THERE ARE ISSUES WITH THIS PRODUCING SINGULAR MATRICES
    var = np.zeros((P,P))
    mean = 0
    for i in range(0,N):
        thetai = theta_i[it][i]
        for k in range(0, K):
            var += ((b_k[it][k]**2/sigma_sqr[it])*np.matmul(np.transpose(thetai),thetai)) + kappa_sqr[it]**-1
            mean += (b_k[it][k]/sigma_sqr[it])*theta_i[it][i]*(Miu[it][i, k] - (a_k[it][k]*b_k[it][k]))
    var = np.linalg.inv(var)
    mean = np.squeeze(np.asarray(np.matmul(mean, var)))
    return np.random.multivariate_normal(np.matmul(mean, var),var)

def update_theta(it):
    update = np.copy(theta_i[it])
    lambdatrans = np.outer(lambda_p[it],lambda_p[it])
    invsigma = np.linalg.inv(sigma[it])
    for i in range(0, N):
        mean = 0
        var = invsigma
        for k in range(0, K):
            var += (b_k[it][k]**2/sigma_sqr[it])**lambdatrans
            mean += ((b_k[it][k] / sigma_sqr[it]) * np.transpose(lambda_p[it]) * (Miu[it][i, k] -
                        (a_k[it][k] * b_k[it][k]))*(np.matmul(invsigma*T_gw[it][IG[i]-1,IW[i]-1],beta_p[it])))
        var = np.linalg.inv(var)
        update[i] = np.random.multivariate_normal(np.matmul(var,mean), var)
    return update

def update_beta(it): #will need to check this, dimensions look off
    invsigma = np.linalg.inv(sigma[it])
    varsum = 0
    mean = 0
    for i in range(0,N):
        varsum += (T_gw[it][IG[i]-1,IW[i]-1]**2) * invsigma
        mean += (T_gw[it][IG[i]-1,IW[i]-1]*theta_i[it][i])
    var = np.linalg.inv(varsum + zeta_sqr[it])
    mean = np.squeeze(np.asarray(np.matmul(mean, var)))
    return np.random.multivariate_normal(mean, var)

def update_T(it): #DONE
    update = np.copy(T_gw[it])
    invsigma = np.linalg.inv(sigma[it])
    for g in range(0,G):
        for w in range(0, W):
            var = 1/phi_sqr[it] + NGW[g,w]*np.matmul(np.matmul(beta_p[it].transpose(),invsigma),beta_p[it])
            mean = 0
            for i in range(0,N):
                mean += (np.matmul(np.matmul(beta_p[it],invsigma),theta_i[it][i]))
            mean += (nu_g[it][g]/phi_sqr[it])
            update[g,w] = norm.rvs(mean/var, var**-1)
    return update

def readin(filename1, filename2):
    vectordata = np.delete(pd.read_csv(filename1).as_matrix(),0,1)
    global Y, IW, IK, IR, IG, N, W, G, K, NGW, IM

    Y = vectordata[:,0]
    IW = vectordata[:,1]
    IK = vectordata[:,2]
    IR = vectordata[:,3]
    IG = vectordata[:,4]

    N = len(set(IR))
    W = max(IW)
    G = max(IG)
    K = max(IK)

    IM = np.zeros((K, N))
    NGW = np.zeros((G,W))

    for response in range(0,N):
        NGW[IG[response]-1,IW[response]-1] += 1

    IM = np.delete((pd.read_csv(filename2)).as_matrix(),0,1)


def main(): #Pass in CSV files with data as args
    #readin(sys.argv[1], sys.argv[2])

    readin('vectors.csv','matrix.csv')

    global Y_iu, Miu, a_k, b_k, lambda_p, theta_i, beta_p, T_gw, omega_sqr, tau_sqr, sigma, zeta_sqr, sigma_sqr, nu_g,\
        phi_sqr, kappa_sqr

    #INITIALIZATION
    omega_sqr.append(invgamma.rvs(hyper_a,hyper_b))
    tau_sqr.append(invgamma.rvs(hyper_a,hyper_b))
    sigma.append(invwishart.rvs(s_0,S_0))
    zeta_sqr.append(invgamma.rvs(hyper_a,hyper_b))
    sigma_sqr.append(invgamma.rvs(hyper_a,hyper_b))
    nu_g.append(np.array([norm.rvs(0, h_sqr) for g in range(0,G)]))
    phi_sqr.append(invgamma.rvs(hyper_a,hyper_b))
    kappa_sqr.append(invgamma.rvs(hyper_a,hyper_b))

    a_k.append(np.array([norm.rvs(0, omega_sqr[0]) for u in range(0,K)]))
    b_k.append(np.array([norm.rvs(0, tau_sqr[0]) for u in range(0,K)]))
    lambda_p.append(np.array([norm.rvs(0, kappa_sqr[0]) for p in range(0,P)]))
    beta_p.append(np.array([norm.rvs(0, zeta_sqr[0]) for p in range(0,P)]))
    T_gw.append(np.matrix([[norm.rvs(nu_g[0][g], phi_sqr[0]) for w in range(0, W)] for g in range(0, G)]))
    theta_i.append(np.matrix([np.random.multivariate_normal(beta_p[0]*T_gw[0][IG[i]-1,IW[i]-1], sigma[0]) for i in range(0,N)]))
    Miu.append(np.matrix([[norm.rvs(b_k[0][u]*(np.dot(lambda_p[0],np.squeeze(np.asarray(theta_i[0][i])))-a_k[0][u]), sigma_sqr[0])
                           for u in range(0,K)] for i in range(0,N)]))
    Y_iu.append(np.matrix([[ 0 if Miu[0][i,u] < 0 else 1 for u in range(0,K)] for i in range(0,N)]))


    for i in range(0, burnpoint):
        print(i)
        omega_sqr[0] = update_omega(0)
        tau_sqr[0] = update_tau(0)
        sigma[0] = update_sigma(0)  #SOMETHING IS WRONG! BREAKS PART OF THE TIME??
        zeta_sqr[0] = update_zeta(0)
        sigma_sqr[0] = update_sigmasqr(0)
        nu_g[0]  = update_nu(0)
        phi_sqr[0] = update_phi(0)
        kappa_sqr[0] = update_kappa(0)


        a_k[0] = update_ak(0)
        b_k[0] = update_bk(0)
        lambda_p[0] = update_lambda(0)
        beta_p[0] = update_beta(0)
        theta_i[0] = update_theta(0)
        T_gw[0] = update_T(0)
        Miu[0] = update_miu(0)
        Y_iu[0] = np.matrix([[ 0 if Miu[0][n,u] < 0 else 1 for u in range(0,K)] for n in range(0,N)])


  #  for numit in range(0 , N_Iterations - burnpoint):
   #     omega_sqr.append(update_omega(numit))
    #    tau_sqr.append(update_tau(numit))
     #   kappa_sqr.append(update_kappa(numit))
      #  sigma.append(update_sigma(numit))
      #  zeta_sqr.append(update_zeta(numit))
      #  sigma_sqr.append(update_sigmasqr(numit))
       # nu_g.append(update_nu(numit))
       # phi_sqr.append(update_phi(numit))
      #  a_k.append(update_ak(numit))
   #     b_k.append(update_bk(numit))
    #    lambda_p.append(update_lambda(numit))
    #    beta_p.append(update_beta(numit))
    #    theta_i.append(update_theta(numit))
    #    T_gw.append(update_T(numit))
    #    Y_iu.append(np.matrix([[ 0 if Miu[numit][i,u] < 0 else 1 for u in range(0,K)] for i in range(0,N)]))
    #    Miu.append(update_miu(numit))

if __name__ == '__main__':
    main()