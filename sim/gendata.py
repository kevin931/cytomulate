###################################
#functions were modified based on the code from https://github.com/luiarthur/cytofresearch
#Credits to the original authors: 
#Arthur Lui et al. A Bayesian Feature Allocation Model for Identification of Cell Subpopulations Using Cytometry Data. 2020. arXiv: 2002.08609 [stat.AP].
#For the purpose of comparison
###################################

import math
import numpy as np
import scipy.stats as stats
import random

def logit(p,a,b):
    return math.log(p-a)-math.log(b-p)

def Categorical(p):
    return random.choices(range(1,len(p)+1),p)[0]

def sigmoid(x,a=0,b=1):
    out=0
    if a==0 and b==1:
        out=1/(1+math.exp(-x))
    else:
        ex=math.exp(x)
        out=(b*ex+a)/(1+ex)
    return out

def prob_miss(y,beta):
    n=len(beta)
    np.array(range(n))
    x=sum(y**np.array(range(n))*beta)
    return sigmoid(x)

def eye(n):
    return np.identity(n)

def leftOrder(Z):
    return Z[:,np.argsort(Z[0,:])[::-1]]
    
def genSimpleZ(J,K):
    g = J//K
    assert g==J/K
    return np.outer(eye(K),np.ones(g))

def isValidZ(Z):
    hasRepeatedColumns = np.unique(Z,axis=1).shape[1] < Z.shape[1]
    hasColOfZero = np.any(np.sum(Z,axis=0)==0)
    return ~(hasColOfZero|hasRepeatedColumns)

def genZ(J,K,prob1):
    assert 0 < prob1 < 1
    Z=np.random.randn(J,K)>prob1
    Z=Z[np.argsort(Z[:,0])]
    Z=leftOrder(Z)
    
    if isValidZ(Z):
        return Z
    else:
        return genZ(J,K,prob1)

def z_get(Z,i,n,j,lam):
    return Z[j,lam[i][n]-1]

def mu_get(gam,i,n,j,Z,lam,mus):
    l = gam[i][n,j]
    if l>0:
        z = z_get(Z,i,n,j,lam)
        out = mus[z][1]
    else:
        out = 0.0
    return out

def genData(J,N,K,L,useSimpleZ=True,prob1=0.6,
           sortLambda=False,propMissingScale=0.7):
    I = len(N)
    Z = genSimpleZ(J,K) if useSimpleZ else genZ(J,K,prob1)
    a_W=[float(i) for i in range(1,K+1)]
    a_eta={}
    for z in [0,1]:
        a_eta[z]=range(1,L[z]+1)
    
    genData(J=J,N=N,K=K,L=L,Z=Z,
           beta=[-9.2,-2.3],
           sig2=np.full(I,0.1),
           mus={0:np.linspace(-5,1,L[0]),1:np.linspace(1,5,L[1])},
           a_W=a_W,
           a_eta=a_eta,
           sortLambda=sortLambda,propMissingScale=propMissingScale)
    

def genData(J,N,K,L,Z,beta,sig2,mus,a_W,a_eta,sortLambda,propMissingScale,eps):

    assert len(eps) == len(N)
    assert min(eps)>=0 and max(eps)<=1

    assert Z.shape[1] == K and Z.shape[0] == J
    
    assert np.all(N>0)
    
    I = len(N)
    
    assert (np.array(sig2)>0).all()
    assert len(sig2) == I
    
    print(np.array(mus[0]),np.array(mus[1]))
    assert (np.array(mus[0]) < 0).all() and (np.array(mus[1]) > 0).all()
    assert len(mus[0]) == L[0] and len(mus[1]) == L[1]
    
    assert len(a_W) == K
    assert (np.array(a_W)>0).all()
    
    assert len(a_eta[0]) == L[0]
    assert len(a_eta[1]) == L[1]
    assert (np.array(a_eta[0]) > 0).all()
    assert (np.array(a_eta[1]) > 0).all()
    
    W = np.zeros((I,K))
    for i in range(I):
        random.shuffle(a_W)
        W[i,:] = stats.dirichlet.rvs(alpha=a_W,size=1)[0]

    lam=[]
    for i in range(I):
        p=np.concatenate([[eps[i]],(1-eps[i])*W[i,:]])
        lam.append([Categorical(p)-1 for ni in range(N[i])])
    lam=np.array(lam)
    print(lam.shape)
    
    if sortLambda:
        lam=[np.sort(lami) for lami in lam]
    
    eta={0:np.zeros((I,J,L[0])),1:np.zeros((I,J,L[1]))}
    for i in range(I):
        for j in range(J):
            random.shuffle(a_eta[0])
            eta[0][i,j,:] = stats.dirichlet.rvs(alpha=a_eta[0])
            random.shuffle(a_eta[1])
            eta[1][i,j,:] = stats.dirichlet.rvs(alpha=a_eta[1])
    
    gam = [np.zeros((Ni,J)) for Ni in N]
    for i in range(I):
        for j in range(J):
            for n in range(N[i]):
                k=lam[i][n]
                if k>0:
                    z=Z[j,k-1]
                    gam[i][n,j] = Categorical(eta[z][i,j,:])
                else:
                    gam[i][n,j] = 0
                    
    
    y = [np.zeros((Ni,J)) for Ni in N]
    y_complete = [np.zeros((Ni,J)) for Ni in N]

    for i in range(I):
        for j in range(J):
            for n in range(N[i]):
                sig_i = np.sqrt(sig2[i]) if lam[i][n]>0 else 3.0
                y_complete[i][n,j] = np.random.normal(loc=mu_get(gam,i,n,j,Z,lam,mus),scale=sig_i)
            
            p_miss = np.array([prob_miss(y_complete[i][n,j],beta) for n in range(N[i])])
            prop_missing = np.random.rand()*propMissingScale*sum(W[i,:]*(1-Z[j,:]))
            num_missing = int(round(N[i]*prop_missing))
            idx_missing = np.random.choice(range(N[i]),num_missing,replace=False,p=p_miss/sum(p_miss))
            y[i][:,j] = y_complete[i][:,j]+0
            y[i][idx_missing,j]=None
    
    return {'y':y, 'y_complete':y_complete, 'Z':Z, 'W':W,
           'eta':eta, 'mus':mus, 'sig2':sig2, 'lam':lam, 'gam':gam,
           'beta':beta, 'eps':eps}

def flipbits(x,prob):
    return [1-xi if prob > np.random.rand() else xi for xi in x]

def createSimData(seed,nfac,K=5,L={0:3,1:3},J=20,eps=np.zeros(3),beta=[-9.2,-2.3],
                  sig2=[0.2,0.1,0.3],sortLambda=False,propMissingScale=0.7):
    N = np.array([8,3,2])*nfac
    random.seed(seed)
    I=len(N)
    Z=np.zeros((J,K))
    while ~isValidZ(Z):
        Z = genZ(J,K,0.5)
        Z[:,2] = flipbits(Z[:,1],prob=0.1)
    
    mus = {0:[-1.0,-2.3,-3.5],
          1:[1.0,2.0,3.0]}
    
    a_W = np.random.rand(K)*10
    
    a_eta = {0:np.random.rand(L[0]), 1:np.random.rand(L[1])}
    
    simdat = genData(J=J, N=N, K=K, L=L, Z=Z,
                                beta=beta,
                                sig2=sig2,
                                mus=mus,
                                a_W=a_W,
                                a_eta=a_eta,
                                sortLambda=sortLambda, propMissingScale=propMissingScale,
                                eps=eps)
    return simdat
