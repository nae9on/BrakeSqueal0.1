'''
Function definitions for obtaining the residual

'''

import numpy as np
from numpy import linalg as LA
import scipy.sparse.linalg as norm

def residual_qevp(M,C,K,la,evec):
        n = la.shape[0]
        res = np.zeros(n,np.complex_)   
        for i in range(0, n):
            vec = evec[:,i:i+1]
            A = la[i]*la[i]*M + la[i]*C + K
            res[i] = LA.norm(A.dot(vec))
            #calculate relative residual
            A_norm=norm.onenormest(A, t=3, itmax=5, compute_v=False, compute_w=False)
            res[i] = res[i]/(A_norm*LA.norm(vec))
        return res;

def residual_qevp_nonsparse(M,C,K,la,evec):
        n = la.shape[0]
        res = np.zeros(n,np.complex_)   
        for i in range(0, n):
            vec = evec[:,i:i+1]
            A = la[i]*la[i]*M + la[i]*C + K
            res[i] = LA.norm(A.dot(vec))
            #calculate relative residual
            A_norm=LA.norm(A)
            res[i] = res[i]/(A_norm*LA.norm(vec))
        return res;


def residual_gevp(A,B,la,evec):
        n = la.shape[0]
        res = np.zeros(n,np.complex_)   
        for i in range(0, n):
            vec = evec[:,i:i+1]
            AB = A - la[i]*B 
            res[i] = LA.norm(AB.dot(vec))
            #calculate relative residual
            AB_norm=norm.onenormest(AB, t=3, itmax=5, compute_v=False, compute_w=False)
            res[i] = res[i]/(AB_norm*LA.norm(vec)) 
        return res;

def residual_gevp_nonsparse(A,B,la,evec):
        n = la.shape[0]
        res = np.zeros(n,np.complex_)   
        for i in range(0, n):
            vec = evec[:,i:i+1]
            AB = A - la[i]*B 
            res[i] = LA.norm(AB.dot(vec))
            #calculate relative residual
            AB_norm=LA.norm(AB)
            res[i] = res[i]/(AB_norm*LA.norm(vec)) 
        return res;


