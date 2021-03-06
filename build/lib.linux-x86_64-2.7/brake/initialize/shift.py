r"""
This module defines the following functions::

  - shift_matrices:
  
    Transform the qevp using shift and invert spectral transformations.
   
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import math
import scipy.sparse.linalg as norm


def shift_matrices(obj, m, c, k, tau):
        r"""
         
        :param obj: object of the class ``BrakeClass``
        :param m: Mass Matrix
        :param c: Damping Matrix
        :param k: Stiffness Matrix
        :param tau: Shift
        :return: M, C, K - Shifted Mass, Damping and Stiffness Matrix respectively
        
        Procedure::
        
         The M , C , K are obtained as follows:
        
          M = m
          C = 2 * tau * m + c
          K = tau_squared * m + tau * c + k

        """
    
        #object attributes used in the function
        LOG_LEVEL = obj.log_level
        logger_t = obj.logger_t
        logger_i = obj.logger_i
                
        if(LOG_LEVEL==10): #Debug Mode
                logger_i.debug("\n"+"\n"+"\n"+'In shift matrices (shift_matrices.py)')
                logger_i.debug('----------------------------------------------------------------')
                logger_i.debug('Shifting about the point '+str(tau))
                
        M=m;
        C=2*tau*m+c;
        tau_squared = tau*tau
        K=tau_squared*m+tau*c+k; 
        
        if(LOG_LEVEL==10): #Debug Mode
                m_norm = norm.onenormest(M, t=3, itmax=5, compute_v=False, compute_w=False)
                c_norm = norm.onenormest(C, t=3, itmax=5, compute_v=False, compute_w=False)
                k_norm = norm.onenormest(K, t=3, itmax=5, compute_v=False, compute_w=False)
        
                logger_i.debug('Properties of shifted matrices')
                logger_i.debug("\n"+'M '+' '+'No of nonzeros = '+str(M.nnz)+' 1-Norm = '+str(m_norm))
                logger_i.debug('C '+' '+'No of nonzeros = '+str(C.nnz)+' 1-Norm = '+str(c_norm))
                logger_i.debug('K '+' '+'No of nonzeros = '+str(K.nnz)+' 1-Norm = '+str(k_norm)) 
                
        return M, C, K;
