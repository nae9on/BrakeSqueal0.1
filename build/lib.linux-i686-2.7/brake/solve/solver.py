'''
Function for the generalized eigenvalue problem

Input
A x = lamda B x
evs_per_shift: no of eigenvalues required
kind: largest or smallest in magnitude. parameters 'LM', 'SM' respectively
flag: flag should be passed true when B is positive definite

Output
1. la (Array of evs_per_shift eigenvalues.)
2. v (An array of evs_per_shift eigenvectors. v[:, i] is the eigenvector 
				corresponding to the eigenvalue la[i])

Additional Info
Example of eigs
eigs(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0,
return_eigenvectors=True, Minv=None, OPinv=None, OPpart=None)

The eigs function of PYTHON can calculate the eigenvalues of the generalized eigenvalue
problem A*x=lamda*M*x with the following conditions.
M must represent a real, symmetric matrix if A is real, and must represent a complex, 
hermitian matrix if A is complex.
If sigma is None, M has to be positive definite 
If sigma is specified, M has to be positive semi-definite

When sigma is specified 'say 0' then eigs function will calculate the eigenvalues nearest to sigma. 
The 'LM' clause along with sigma = 0 can be used to calculate the reciprocal eigenvalues of 
Largest Magnitude.

'''


import csv
import array
import time
import timeit
import unittest
import warnings
import logging
from brake.analyze import residual
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import *
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve
from scipy.sparse import vstack, hstack
import scipy
#from scipy.linalg import inv

#Define exceptions
class SOLVERError(Exception): pass
class SOLVER_BadInputError(SOLVERError): pass

#in QEVP: la, evec = solver.qev_optimized_sparse(M,C,K,evs_per_shift);
def qev_optimized_sparse(obj,M,C,K):
    
   LOG_LEVEL = obj.log_level
   logger_t = obj.logger_t
   logger_i = obj.logger_i
   
   evs_per_shift = obj.evs_per_shift

   if(sparse.issparse(M) & sparse.issparse(C) & sparse.issparse(K)):
   	K = K.tocsc()

        #Computing the LU factors of Sparse K
   	try:
   	    begin_LU = timeit.default_timer()
            #LUfactors = splu(K, permc_spec=None, diag_pivot_thresh=None, drop_tol=None, relax=None, panel_size=None, options={})
   	    LUfactors = spilu(K, permc_spec=None, diag_pivot_thresh=None, drop_tol=None, relax=None, panel_size=None, options={})
            #LUfactors = factorized(K) 
   	    end_LU = timeit.default_timer()
   	    logger_t.info("\n"+"\n"+'SOLVER: LU Factorization of sparse K using spilu in time : '+"%.2f" % (end_LU-begin_LU)+' sec')
   	except RuntimeError:
   	    raise SOLVER_BadInputError, 'The matrix K is sparse but singular'
	
   	n = M.shape[0]
   	# Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
   	def mv(x):
   	 y1 = LUfactors.solve(-C.dot(x[0:n])-M.dot(x[n:2*n]))
   	 y2 = x[0:n]
   	 y= np.hstack((y1,y2))	
   	 return y
   	Aoperator = LinearOperator((2*n,2*n),matvec=mv)
   	
   	# computation
   	begin_eigen_solve = timeit.default_timer()  
   	la,v = eigs(Aoperator, k=evs_per_shift, M=None, sigma=None, which='LM')
   	end_eigen_solve = timeit.default_timer()
   	logger_t.info('SOLVER: Call to python eigs command: '+"%.2f" % (end_eigen_solve-begin_eigen_solve)+' sec') 
   	
        #Inverting the obtained eigen values   	      
   	for i in range(0, la.shape[0]):
   	   x1 = la[i].real
   	   x2 = la[i].imag
   	   la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
    
   	 # -----------------------------------------------------------------------------------
   else:
	raise SOLVER_BadInputError, 'The matrix M,C,K are not all in sparse format'
   return la, v;

#IN main:  la, evec = solver.qev_optimized_dense(M,C,K,evs_per_shift);
def qev_optimized_dense(obj,M,C,K):
    
   LOG_LEVEL = obj.log_level
   logger_t = obj.logger_t
   logger_i = obj.logger_i
   
   evs_per_shift = obj.evs_per_shift

   if(~sparse.issparse(M) & ~sparse.issparse(C) & ~sparse.issparse(K)):
   	
        #Computing the LU factors of Dense K
        try:
 	  LUfactors = lu_factor(K, overwrite_a=False, check_finite=True)
	except RuntimeWarning:
	  raise SOLVER_BadInputError, 'The matrix K is dense but singular'

	n = M.shape[0]
        
	# Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
        def mv(x):
          y1 = lu_solve(LUfactors, -C.dot(x[0:n])-M.dot(x[n:2*n]), trans=0, overwrite_b=False, check_finite=True)
          y2 = x[0:n]  
	  y= np.hstack((y1,y2))
          return y
	Aoperator = LinearOperator((2*n,2*n),matvec=mv)

	# computation
        begin_eigen_solve = timeit.default_timer()  
        la,v = eigs(Aoperator, k=evs_per_shift, M=None, sigma=None, which='LM')
        end_eigen_solve = timeit.default_timer()
        
	#logger_t.info('SOLVER: Call to python eigs command: '+"%.2f" % (end_eigen_solve-begin_eigen_solve)+' sec') 
        
	#Inverting the obtained eigen values        
        for i in range(0, la.shape[0]):
          x1 = la[i].real
          x2 = la[i].imag
          la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
        # -----------------------------------------------------------------------------------
   else:
	raise SOLVER_BadInputError, 'The matrix M,C,K are not all dense'
   return la, v;


def qev_sparse(obj,M,C,K):
    
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i
   
    evs_per_shift = obj.evs_per_shift

    if(sparse.issparse(M) & sparse.issparse(C) & sparse.issparse(K)):
    	
	#-----Creating matrix A and B
    	n = M.shape[0]
    	begin_create = timeit.default_timer()
    	I = sparse.identity(n, dtype=complex, format='csr')
    	Z = I-I
    	A = sparse.block_diag((K,I))
	B = vstack([hstack([-1*C,-1*M]),hstack([I,Z])])
    	end_create = timeit.default_timer()
    	logger_t.info("\n"+"\n"+'Creating A and B for GEV: '+"%.2f" % (end_create-begin_create)+' sec')

    	A = A.tocsc()	# warn('splu requires CSC matrix format', SparseEfficiencyWarning)

        #Computing the LU factors of Sparse A
    	try:
          begin_LU = timeit.default_timer() 
          LUfactors = splu(A, permc_spec=None, diag_pivot_thresh=None, drop_tol=None, relax=None, panel_size=None, options={})
          #LUfactors = spilu(A, drop_tol=None, fill_factor=None, drop_rule=None, permc_spec=None, \
          #								diag_pivot_thresh=None, relax=None, panel_size=None, options=None)
          end_LU = timeit.default_timer()
          logger_t.info('SOLVER: LU Factorization of A in time : '+"%.2f" % (end_LU-begin_LU)+' sec')
        except RuntimeError:
          raise SOLVER_BadInputError, 'The matrix A is sparse but singular' 

        # Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
    	def mv(x):
          y = B.dot(x)
          z = LUfactors.solve(y)
          return z	
    	Aoperator = LinearOperator((2*n,2*n),matvec=mv)

    	# computation
    	begin_eigen_solve = timeit.default_timer()  
   	la,v = eigs(Aoperator, k=evs_per_shift, M=None, sigma=None, which='LM')
    	end_eigen_solve = timeit.default_timer()
    	logger_t.info('SOLVER: Call to python eigs command : '+"%.2f" % (end_eigen_solve-begin_eigen_solve)+' sec') 
    	
        #Inverting the obtained eigen values  
    	for i in range(0, la.shape[0]):
      	 x1 = la[i].real
      	 x2 = la[i].imag
      	 la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
    
    	# -----------------------------------------------------------------------------------
    else:
        raise SOLVER_BadInputError, 'The matrix M,C,K are not all in sparse format'	
    return la, v;

def qev_dense(obj,M,C,K):
    
   LOG_LEVEL = obj.log_level
   logger_t = obj.logger_t
   logger_i = obj.logger_i
   
   evs_per_shift = obj.evs_per_shift

   if(~sparse.issparse(M) & ~sparse.issparse(C) & ~sparse.issparse(K)):
   	
	n = M.shape[0]
	I = np.identity(n, dtype=complex)
	Z = I-I
	#Creating matrix A
        AU = np.concatenate((K,Z), axis=1)
	AL = np.concatenate((Z,I), axis=1)
	A = np.concatenate((AU,AL), axis=0)
        #Creating matrix B
	BU = np.concatenate((-1*C,-1*M), axis=1)
	BL = np.concatenate((I,Z), axis=1)
	B = np.concatenate((BU,BL), axis=0)
	

	#Computing the LU factors of Dense A
	warnings.simplefilter("error", RuntimeWarning)		
	try:
 	  LUfactors = lu_factor(A, overwrite_a=False, check_finite=True)
	except RuntimeWarning:
	  raise GEV_BadInputError, 'The matrix A is dense but singular'

        # Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
	def mv(x):
	  y = B.dot(x)
	  y = y.T
	  z = lu_solve(LUfactors, y, trans=0, overwrite_b=False, check_finite=True)
          return z
        Aoperator = LinearOperator((2*n,2*n),matvec=mv)

	# computation
        begin_eigen_solve = timeit.default_timer()  
        la,v = eigs(Aoperator, k=evs_per_shift, M=None, sigma=None, which='LM')
        end_eigen_solve = timeit.default_timer()
        
        #res_gevp = residual_check.residual_gevp_nonsparse(A,B,la,v)
	#print 'Maximum Residual error for the GEVP is ',max(res_gevp)
        
	#logger_t.info('GEV: Call to python eigs command: '+"%.2f" % (end_eigen_solve-begin_eigen_solve)+' sec') 
  
        #Inverting the obtained eigen values 
        for i in range(0, la.shape[0]):
          x1 = la[i].real
          x2 = la[i].imag
          la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
    
   else:
	raise SOLVER_BadInputError, 'M,C,K are not all dense'
   return la, v;
