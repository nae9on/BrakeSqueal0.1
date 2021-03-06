r"""
This module defines the following functions::

  - qev_sparse:
  
    Obtains the eigenvalues(smallest in magnitude) and corresponding eigenvectors 
    for the given Quadratic Eigenvalue Problem(QEVP) with sparse M, C, K matrices.
   
  - qev_dense
  
    Obtains the eigenvalues(smallest in magnitude) and corresponding eigenvectors 
    for the given Quadratic Eigenvalue Problem(QEVP) with dense M, C, K matrices.
      
  - gev_sparse

    Obtains the eigenvalues(smallest in magnitude) and corresponding eigenvectors 
    for the Generalized Eigenvalue Problem(GEVP) with sparse A, B matrices. 
    
  - gev_dense
  
    Obtains the eigenvalues(smallest in magnitude) and corresponding eigenvectors 
    for the Generalized Eigenvalue Problem(GEVP) with dense A, B matrices.
    
Note::

 The eigs function of PYTHON can calculate the eigenvalues of the generalized eigenvalue
 problem A*x=lamda*M*x with the following conditions.
 M must represent a real, symmetric matrix if A is real, and must represent a complex, 
 hermitian matrix if A is complex.
 If sigma is None, M has to be positive definite 
 If sigma is specified, M has to be positive semi-definite

 When sigma is specified 'ex 0' then eigs function will calculate the eigenvalues nearest to sigma. 
 The 'LM' clause along with sigma = 0 can be used to calculate the reciprocal eigenvalues of 
 Largest Magnitude.    

"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import timeit
import unittest
import warnings
import numpy as np
from scipy import sparse
from scipy.sparse import vstack, hstack
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import eigs, spilu, LinearOperator, splu

#----------------------------Application Specific Imports----------------------------------------- 
from brake.analyze import residual

#Define exceptions
class SOLVERError(Exception): pass
class SOLVER_BadInputError(SOLVERError): pass

#in QEVP: la, evec = solver.qev_sparse(obj,M,C,K);
def qev_sparse(obj,M,C,K):
   r"""
   :param obj: object of the class ``BrakeClass``
   :param M: Mass Matrix
   :param C: Damping Matrix
   :param K: Stiffness Matrix
   :return: ``la`` - eigenvalues, ``v`` - eigenvectors
   :raises: ``SOLVER_BadInputError``, When the matrix M,C,K are not all in sparse format
   
   Example::
   
    .
   
   """
    
   LOG_LEVEL = obj.log_level
   logger_t = obj.logger_t
   logger_i = obj.logger_i
   
   evs_per_shift = obj.evs_per_shift
   
   if(sparse.issparse(M) & sparse.issparse(C) & sparse.issparse(K)):
        
        K = K.tocsc()

        #Computing the LU factors of Sparse K
        try:
            begin_LU = timeit.default_timer()
            
            LUfactors = splu(K, permc_spec=None, diag_pivot_thresh=None, drop_tol=None, 
            relax=None, panel_size=None, options={})
            
            '''
            LUfactors = spilu(K, permc_spec=None, diag_pivot_thresh=None, drop_tol=None, 
            relax=None, panel_size=None, options={})
            
            LUfactors = factorized(K) 
            '''
            
            end_LU = timeit.default_timer()
            logger_t.info("\n"+"\n"+'SOLVER: LU Factorization of sparse K using splu in time : '\
            +"%.2f" % (end_LU-begin_LU)+' sec')
            
        except RuntimeError:
            raise SOLVER_BadInputError('The matrix K is sparse but singular')
        
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
        logger_t.info('SOLVER: Call to python eigs command: '+"%.2f" % \
        (end_eigen_solve-begin_eigen_solve)+' sec') 
        
        #Inverting the obtained eigen values          
        for i in range(0, la.shape[0]):
           x1 = la[i].real
           x2 = la[i].imag
           la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
   else:
        raise SOLVER_BadInputError('The matrix M,C,K are not all in sparse format')
   return la, v;

#in main:  la, evec = solver.qev_dense(M,C,K,no_of_evs);
def qev_dense(obj,M,C,K,no_of_evs):
   r"""
   :param obj: object of the class ``BrakeClass``
   :param M: Mass Matrix
   :param C: Damping Matrix
   :param K: Stiffness Matrix
   :param no_of_evs: No of eigenvalues to be computed
   :return: ``la`` - eigenvalues, ``v`` - eigenvectors
   :raises: ``SOLVER_BadInputError``, When the matrix M,C,K are not all dense
   
   Example::
   
    . 
    
   """
   
   LOG_LEVEL = obj.log_level
   logger_t = obj.logger_t
   logger_i = obj.logger_i
   
   if(~sparse.issparse(M) & ~sparse.issparse(C) & ~sparse.issparse(K)):
        
        #Computing the LU factors of Dense K
        try:
          LUfactors = lu_factor(K, overwrite_a=False, check_finite=True)
        except RuntimeWarning:
          raise SOLVER_BadInputError('The matrix K is dense but singular')

        n = M.shape[0]
        
        # Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
        def mv(x):
          y1 = lu_solve(LUfactors, -C.dot(x[0:n])-M.dot(x[n:2*n]), trans=0, overwrite_b=False, 
          check_finite=True)
          y2 = x[0:n]  
          y= np.hstack((y1,y2))
          return y
        Aoperator = LinearOperator((2*n,2*n),matvec=mv)

        # computation
        begin_eigen_solve = timeit.default_timer()  
        la,v = eigs(Aoperator, k=no_of_evs, M=None, sigma=None, which='LM')
        end_eigen_solve = timeit.default_timer()
        
        '''
        logger_t.info('SOLVER: Call to python eigs command: '+"%.2f" % \
        (end_eigen_solve-begin_eigen_solve)+' sec') 
        '''
        
        #Inverting the obtained eigen values        
        for i in range(0, la.shape[0]):
          x1 = la[i].real
          x2 = la[i].imag
          la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
   else:
        raise SOLVER_BadInputError('The matrix M,C,K are not all dense')
   return la, v;


def gev_sparse(obj,A,B):
    r"""
    :param obj: object of the class ``BrakeClass``
    :param A: ``A x = lamda B x``
    :param B: ``A x = lamda B x``
    :return: ``la`` - eigenvalues, ``v`` - eigenvectors
    :raises: ``SOLVER_BadInputError``, When the matrix A, B are not all in sparse format
    
    Example::
   
     .
    
    """    
    
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i
   
    evs_per_shift = obj.evs_per_shift

    if(sparse.issparse(A) & sparse.issparse(B)):
        
        '''
        Example Implementation of the 	QEVP
        #-----Creating matrix A and B
        n = 2*M.shape[0]
        begin_create = timeit.default_timer()
        I = sparse.identity(n, dtype=complex, format='csr')
        Z = I-I
        A = sparse.block_diag((K,I))
        B = vstack([hstack([-1*C,-1*M]),hstack([I,Z])])
        end_create = timeit.default_timer()
        #logger_t.info("\n"+"\n"+'Creating A and B for GEV: '+"%.2f" % \
        #(end_create-begin_create)+' sec')
        '''
        
        A = A.tocsc()   # warn('splu requires CSC matrix format', SparseEfficiencyWarning)
        n = A.shape[0]
        
        #Computing the LU factors of Sparse A
        try:
          begin_LU = timeit.default_timer() 
          LUfactors = splu(A, permc_spec=None, diag_pivot_thresh=None, drop_tol=None, 
          relax=None, panel_size=None, options={})
          '''
          #LUfactors = spilu(A, drop_tol=None, fill_factor=None, drop_rule=None, 
          permc_spec=None, diag_pivot_thresh=None, relax=None, panel_size=None, options=None)
          '''
          end_LU = timeit.default_timer()
          logger_t.info('SOLVER: LU Factorization of A in time : '+"%.2f" % \
          (end_LU-begin_LU)+' sec')
          
        except RuntimeError:
          raise SOLVER_BadInputError('The matrix A is sparse but singular' )

        # Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
        def mv(x):
          y = B.dot(x)
          z = LUfactors.solve(y)
          return z      
        Aoperator = LinearOperator((n,n),matvec=mv)

        # computation
        begin_eigen_solve = timeit.default_timer()  
        la,v = eigs(Aoperator, k=evs_per_shift, M=None, sigma=None, which='LM')
        end_eigen_solve = timeit.default_timer()
        logger_t.info('SOLVER: Call to python eigs command : '+"%.2f" % \
        (end_eigen_solve-begin_eigen_solve)+' sec') 
        
        #Inverting the obtained eigen values  
        for i in range(0, la.shape[0]):
         x1 = la[i].real
         x2 = la[i].imag
         la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
       
    else:
        raise SOLVER_BadInputError('The matrix A, B are not all in sparse format')
    return la, v;

def gev_dense(obj,A,B):
   r"""
   :param obj: object of the class ``BrakeClass``
   :param A: ``A x = lamda B x``
   :param B: ``A x = lamda B x``
   :return: ``la`` - eigenvalues, ``v`` - eigenvectors
   :raises: ``SOLVER_BadInputError``, When the matrix A, B are not all dense
   
   Example::
   
    .
   
   """
   
   LOG_LEVEL = obj.log_level
   logger_t = obj.logger_t
   logger_i = obj.logger_i
   
   evs_per_shift = obj.evs_per_shift

   if(~sparse.issparse(A) & ~sparse.issparse(B)):
        
        '''
        Example Implementation of the 	QEVP
        n = 2*M.shape[0]
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
        '''
        n = A.shape[0]
        
        #Computing the LU factors of Dense A
        warnings.simplefilter("error", RuntimeWarning)          
        try:
          LUfactors = lu_factor(A, overwrite_a=False, check_finite=True)
        except RuntimeWarning:
          raise GEV_BadInputError('The matrix A is dense but singular')

        # Operator defining the matrix vector product A*x 
        # which is passed to eigs(python inbuilt solver) for calculating the eigenvalues
        def mv(x):
          y = B.dot(x)
          y = y.T
          z = lu_solve(LUfactors, y, trans=0, overwrite_b=False, check_finite=True)
          return z
        Aoperator = LinearOperator((n,n),matvec=mv)

        # computation
        begin_eigen_solve = timeit.default_timer()  
        la,v = eigs(Aoperator, k=evs_per_shift, M=None, sigma=None, which='LM')
        end_eigen_solve = timeit.default_timer()
        
        #res_gevp = residual_check.residual_gevp_nonsparse(A,B,la,v)
        #print 'Maximum Residual error for the GEVP is ',max(res_gevp)
        
        #logger_t.info('GEV: Call to python eigs command: '+"%.2f" % \
        #(end_eigen_solve-begin_eigen_solve)+' sec') 
  
        #Inverting the obtained eigen values 
        for i in range(0, la.shape[0]):
          x1 = la[i].real
          x2 = la[i].imag
          la[i] = complex(x1,-x2)/(x1*x1+x2*x2)
    
   else:
        raise SOLVER_BadInputError('A,B are not all dense')
   return la, v;
