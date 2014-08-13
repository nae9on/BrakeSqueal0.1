'''
Diagonal Scale QEVP

DiagScaleMatrices
Input
1. M
2. C
3. K

Output
1. M
2. C
3. K
4. D1
5. D2

'''

import numpy as np
from scipy.sparse.linalg import eigs
from scipy.sparse import issparse, identity, dia_matrix

maxiter = 1000;
epsilon = 1e-2; # in matlab implementation it is 1e-2

def normalize_cols(A):  #assuming A is real
        n = A.shape[1]
	if(1):
 	  CC = A.copy()
	  if(issparse(CC)):
		 CC.data = (CC.data ** 2)
		 diag = CC.sum(axis=0)
		 diag = diag.A1
	  else:
		 CC = np.square(CC)
	 	 diag = CC.sum(axis=0)
          	 	
	diag = 1/np.sqrt(diag)
	
	D = dia_matrix(([diag], [0]), shape=(n,n))
        return D;

def norm_rc(Y):
	n = Y.shape[0]
        error = float("inf")
	I = identity(n, dtype='float64', format='csr')
	DL = identity(n, dtype='float64', format='csr')
	DR = identity(n, dtype='float64', format='csr')
        i = 1
	while (error > epsilon) and (i<maxiter):
	       # normalize columns	
	       Dcol = normalize_cols(Y);
   	       DR = DR*Dcol;
	       Y = Y*Dcol		
	       # normalize rows
               Drow = normalize_cols(Y.T);    
    	       DL = Drow*DL;
               Y = Drow*Y;
               e1 = abs(Drow-I)
	       e2 = abs(Dcol-I) 
    	       error=max(e1.max(),e2.max());
               i = i + 1;
	#print i, error
        return Y,DL,DR,i
	     
def diag_scale_matrices(obj, M, C, K):
	LOG_LEVEL = obj.log_level
	logger_t = obj.logger_t
	logger_i = obj.logger_i
        
	if(LOG_LEVEL==10): #Debug Mode		
		logger_i.debug("\n"+"\n"+"\n"+'In diagonal scale matrices (scale_matrices.py)')
		logger_i.debug('----------------------------------------------------------------')
		K = K.tocsc()
		la1,v1 = eigs(K, k=1, M=None)
		la2,v2 = eigs(K, k=1, M=None, sigma=0, which='LM')
		logger_i.debug('Condition Number of K before scaling '+str(abs(la1)/abs(la2)))
		
        n = M.shape[0]
        I = identity(n, dtype=complex, format='csr')
        Y=abs(M)+abs(C)+abs(K);
        Y, DL, DR, i = norm_rc(Y)
        M=DL*M*DR;
        C=DL*C*DR;
        K=DL*K*DR;
	
	if(LOG_LEVEL==10): #Debug Mode
		K = K.tocsc()
		la1,v1 = eigs(K, k=1, M=None)
		la2,v2 = eigs(K, k=1, M=None, sigma=0, which='LM')
		logger_i.debug('Condition Number of K after scaling '+str(abs(la1)/abs(la2)))
        	K = K.tocsr()
		logger_i.debug('The condition number of K decreases after diagonal scaling')
		
        return M, C, K, DL, DR;



#previous routine implemented by sarosh
'''
def normalize_cols(A):  #assuming A is real
        n = A.shape[1]
	if(0): #slower
          AA = A.T * A
          diag = AA.diagonal()
        if(1):
 	  CC = A.copy()
	  if(scipy.sparse.issparse(CC)):
		 CC.data = (CC.data ** 2)
		 diag = CC.sum(axis=0)
		 diag = diag.A1
	  else:
		 CC = np.square(CC)
	 	 diag = CC.sum(axis=0)
          	 #diag = diag.A1	
	diag = 1/np.sqrt(diag)
	
	D = dia_matrix(([diag], [0]), shape=(n,n))
        B = A*D
        return B, D;

def norm_rc(Y):
	n = Y.shape[0]
        error = float("inf")
	I = sparse.identity(n, dtype='float64', format='csr')
	D1 = sparse.identity(n, dtype='float64', format='csr')
	D2 = sparse.identity(n, dtype='float64', format='csr')
        i = 1
	while (error > epsilon) and (i<maxiter):
	       Y, Dtemp2 = normalize_cols(Y);
   	       D2 = D2*Dtemp2;
               # normalize rows
               Y, Dtemp1 = normalize_cols(Y.T);    
    	       D1 = D1*Dtemp1;
               Y=Y.T;
               e1 = abs(Dtemp1-I)
	       e2 = abs(Dtemp2-I) 
    	       error=max(e1.max(),e2.max());
               #print i, error
               i = i + 1;
        return Y,D1,D2,i
	     
def DiagScaleMatrices(M,C,K,level=0):
	LOG_LEVEL = level
	logger_t = logging.getLogger('time_logger')
	logger_i = logging.getLogger('info_logger')
        
	if(LOG_LEVEL==10): #Debug Mode		
		logger_i.debug("\n"+"\n"+"\n"+'In diagonal scale matrices (scale_matrices.py)')
		logger_i.debug('----------------------------------------------------------------')
		K = K.tocsc()
		la1,v1 = eigs(K, k=1, M=None)
		la2,v2 = eigs(K, k=1, M=None, sigma=0, which='LM')
		logger_i.debug('Condition Number of K before scaling '+str(abs(la1)/abs(la2)))
		
        n = M.shape[0]
        I = sparse.identity(n, dtype=complex, format='csr')
        MCKmat=abs(M)+abs(C)+abs(K);
        MCKmat, D1, D2, i = norm_rc(MCKmat)
        M=D1*M*D2;
        C=D1*C*D2;
        K=D1*K*D2;
	
	if(LOG_LEVEL==10): #Debug Mode
		K = K.tocsc()
		la1,v1 = eigs(K, k=1, M=None)
		la2,v2 = eigs(K, k=1, M=None, sigma=0, which='LM')
		logger_i.debug('Condition Number of K after scaling '+str(abs(la1)/abs(la2)))
        	K = K.tocsr()
		logger_i.debug('The condition number of K decreases after diagonal scaling')
		
        return M, C, K, D1, D2;
'''
