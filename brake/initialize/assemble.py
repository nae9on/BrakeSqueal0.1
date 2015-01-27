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
         
        :param obj: object of the class ``BrakeClass``
        :param sparse_list: a python list of matrices in Compressed Sparse Column format 
          of type '<type 'numpy.float64'>',
        :param omega: angular frequency
        :return: ``M`` - Mass Matrix, ``C`` - Damping Matrix, ``K`` - Stiffness Matrix
        :raises: Assemble_BadInputError, When a matrix in the list is not sparse
        :raises: Assemble_BadInputError, When a matrix in the list is not square
        :raises: Assemble_BadInputError, When the matrix  are not of the same size
        
        Procedure::
         
         The M , C , K are assembled as follows:
        
         - M = M1
         - C = D1+DR*(omegaRef/omega)+DG*(omega/omegaRef)
         - K = K1+KR+KGeo*math.pow((omega/omegaRef),2)

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

        #converting csc to csr format
        M1 = sparse_list[0].tocsr()
        D1 = sparse_list[1].tocsr()
        DG = sparse_list[2].tocsr()
        DR = sparse_list[3].tocsr()
        D4 = sparse_list[4].tocsr()
        K1 = sparse_list[5].tocsr()
        KR = sparse_list[6].tocsr()
        KGeo = sparse_list[7].tocsr()

        if(LOG_LEVEL==10): #Debug Mode
	  
          logger_i.debug('Matrices in CSC format converted to CSR')
	  
	  for i in range(0,len(sparse_list)):
		componentMatrix = sparse_list[i]
		normMatrix = norm.onenormest(componentMatrix.tocsr(), t=3, itmax=5, compute_v=False, compute_w=False)
		logger_i.debug(obj.data_file_list[i]+' '+obj.data_file_name[i]+' Nonzeros = '+str(componentMatrix.nnz)+' 1-Norm = '+str(normMatrix))


        #Old implementation from sarosh's code
        '''
        M = m
        C = c1+c2*(omega/omegaRef)+c3*((omegaRef/omega)-1)+c4/(2*math.pi*fRef)
        K = k1+k2+k3*(math.pow((omega/omegaRef),2)-1);
        '''

        M = M1
        C = D1+DR*(omegaRef/omega)+DG*(omega/omegaRef)
        K = K1+KR+KGeo*math.pow((omega/omegaRef),2)

        if(LOG_LEVEL==10): #Debug Mode
          m_norm = norm.onenormest(M, t=3, itmax=5, compute_v=False, compute_w=False)
          c_norm = norm.onenormest(C, t=3, itmax=5, compute_v=False, compute_w=False)
          k_norm = norm.onenormest(K, t=3, itmax=5, compute_v=False, compute_w=False)

          logger_i.debug("\n"+"\n"+'Properties of assembled matrices')
          logger_i.debug('M '+' '+'Nonzeros = '+str(M.nnz)+' 1-Norm = '+str(m_norm))
          logger_i.debug('C '+' '+'Nonzeros = '+str(C.nnz)+' 1-Norm = '+str(c_norm))
          logger_i.debug('K '+' '+'Nonzeros = '+str(K.nnz)+' 1-Norm = '+str(k_norm))

        return M, C, K;
