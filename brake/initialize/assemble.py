r"""
This module defines the following functions::

  - create_MCK:
  
    Assembles the various component matrices together(for the given angular frequency
    ``omega``) to form the mass(M), damping(C) and stiffness matrix(K).
   
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import math
import scipy.sparse
import scipy.sparse.linalg as norm

#Define exceptions for unittesting
class AssembleError(Exception): pass
class Assemble_BadInputError(AssembleError): pass

def create_MCK(obj, sparse_list, omega):
        r"""
         
        INPUT:
        
        - ``obj`` -- object of the class ``BrakeClass``
        - ``sparse_list`` -- a python list of matrices in Compressed Sparse Column format 
          of type '<type 'numpy.float64'>'
        - ``omega`` -- angular frequency
              
        OUTPUT:
         
        - ``M`` -- Mass Matrix
        - ``C`` -- Damping Matrix
        - ``K`` -- Stiffness Matrix
     
        The ``M`` , ``C`` , ``K`` are assembled as follows:
        
        - ``M`` = m
        - ``C`` = c1+c2*(omega/omegaRef)+c3*(omegaRef/omega)
        - ``K`` = k1+k2+k3*math.pow((omega/omegaRef),2)

        """
        
        #object attributes used in the function
        LOG_LEVEL = obj.log_level
        logger_t = obj.logger_t
        logger_i = obj.logger_i
        omegaRef = obj.omegaRef
        fRef = obj.fRef
        
        #unittesting
        if len(sparse_list) != 8:
         raise Assemble_BadInputError('The sparse list is not of length 8')
        else:
         for i in range(0,len(sparse_list)):
                if(not scipy.sparse.issparse(sparse_list[i])):
                  raise Assemble_BadInputError('The list is not sparse')
                if (sparse_list[i].shape[0]!=sparse_list[i].shape[1]):
                  raise Assemble_BadInputError('The matrix is not square')
                if (sparse_list[i].shape[0]!=sparse_list[0].shape[0]):
                  raise Assemble_BadInputError('The matrix are not of the same size')
                  
        

        if(LOG_LEVEL==10): #Debug Mode
                logger_i.debug("\n"+"\n"+"\n"+'In assemble matrices (assemble.py)')
                logger_i.debug('----------------------------------------------------------------') 
        
        m = sparse_list[0]
        c1 = sparse_list[1]
        c2 = sparse_list[2]
        c3 = sparse_list[3]
        c4 = sparse_list[4]
        k1 = sparse_list[5]
        k2 = sparse_list[6]
        k3 = sparse_list[7]

        #converting csc to csr format
        m = m.tocsr()
        c1 = c1.tocsr()
        c2 = c2.tocsr()
        c3 = c3.tocsr()
        c4 = c4.tocsr()
        k1 = k1.tocsr()
        k2 = k2.tocsr()
        k3 = k3.tocsr()

        if(LOG_LEVEL==10): #Debug Mode
	  
          logger_i.debug('Matrices in CSC format converted to CSR')
          
          #Calculating norm
          m_norm = norm.onenormest(m, t=3, itmax=5, compute_v=False, compute_w=False)
          c1_norm = norm.onenormest(c1, t=3, itmax=5, compute_v=False, compute_w=False) 
          c2_norm = norm.onenormest(c2, t=3, itmax=5, compute_v=False, compute_w=False)
          c3_norm = norm.onenormest(c3, t=3, itmax=5, compute_v=False, compute_w=False)
          c4_norm = norm.onenormest(c4, t=3, itmax=5, compute_v=False, compute_w=False)
          k1_norm = norm.onenormest(k1, t=3, itmax=5, compute_v=False, compute_w=False)
          k2_norm = norm.onenormest(k2, t=3, itmax=5, compute_v=False, compute_w=False) 
          k3_norm = norm.onenormest(k3, t=3, itmax=5, compute_v=False, compute_w=False)

          logger_i.debug("\n"+'Properties of various component matrices')
          logger_i.debug('m '+'symmetric '+'Nonzeros = '+str(m.nnz)+' 1-Norm = '+str(m_norm))
          logger_i.debug('c1 '+'symmetric '+'NNZ = '+str(c1.nnz)+' 1-Norm = '+str(c1_norm))
          logger_i.debug('c2 Dg '+'antisymmetric '+'NNZ = '+str(c2.nnz)+' 1-Norm = '+str(c2_norm))
          logger_i.debug('c3 Dr '+'symmetric '+'NNZ = '+str(c3.nnz)+' 1-Norm = '+str(c3_norm))
          logger_i.debug('c4 '+' '+'NNZ = '+str(c4.nnz)+' 1-Norm = '+str(c4_norm))
          logger_i.debug('k1 '+'symmetric '+'NNZ = '+str(k1.nnz)+' 1-Norm = '+str(k1_norm))
          logger_i.debug('k2 Kr '+'none '+'NNZ = '+str(k2.nnz)+' 1-Norm = '+str(k2_norm))
          logger_i.debug('k3 Kgeo '+'symmetric '+'NNZ = '+str(k3.nnz)+' 1-Norm = '+str(k3_norm))

        
        #Old implementation from sarosh's code
        '''
        M = m
        C = c1+c2*(omega/omegaRef)+c3*((omegaRef/omega)-1)+c4/(2*math.pi*fRef)
        K = k1+k2+k3*(math.pow((omega/omegaRef),2)-1);
        '''
        
        M = m
        C = c1+c2*(omega/omegaRef)+c3*(omegaRef/omega)
        K = k1+k2+k3*math.pow((omega/omegaRef),2)
        
        if(LOG_LEVEL==10): #Debug Mode
          m_norm = norm.onenormest(M, t=3, itmax=5, compute_v=False, compute_w=False)
          c_norm = norm.onenormest(C, t=3, itmax=5, compute_v=False, compute_w=False)
          k_norm = norm.onenormest(K, t=3, itmax=5, compute_v=False, compute_w=False)

          logger_i.debug("\n"+"\n"+'Properties of assembled matrices')
          logger_i.debug('M '+' '+'No of nonzeros = '+str(M.nnz)+' 1-Norm = '+str(m_norm))
          logger_i.debug('C '+' '+'No of nonzeros = '+str(C.nnz)+' 1-Norm = '+str(c_norm))
          logger_i.debug('K '+' '+'No of nonzeros = '+str(K.nnz)+' 1-Norm = '+str(k_norm))

        return M, C, K;
