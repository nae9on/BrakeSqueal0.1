r"""
This module defines the following functions::

  - scale_matrices:
  
    Scales the M, C, K matrices using 2-scalers before linearization.
   
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import math
import scipy.sparse.linalg as norm

def scale_matrices(obj, m , c, k):
        r"""
         
        INPUT:
        
        - ``obj`` -- object of the class ``BrakeClass``
        - ``m`` -- Mass Matrix
        - ``c`` -- Damping Matrix
        - ``k`` -- Stiffness Matrix
                      
        OUTPUT:
         
        - ``M`` -- Scaled Mass Matrix
        - ``C`` -- Scaled Damping Matrix
        - ``K`` -- Scaled Stiffness Matrix
        - ``gamma`` -- scaling parameter
        - ``delta`` -- scaling parameter
     
        The ``M`` , ``C`` , ``K``, ``gamma``, ``delta``  are obtained as follows:
        
        - ``M`` = gamma*gamma*delta*m;
        - ``C`` = gamma*delta*c;
        - ``K`` = delta*k;
        - ``gamma`` = math.sqrt(k_norm/m_norm);
        - ``delta`` = 2/(k_norm+c_norm*gamma);       

        """
        
        #object attributes used in the function
        LOG_LEVEL = obj.log_level
        logger_t = obj.logger_t
        logger_i = obj.logger_i
        
        if(LOG_LEVEL==10): #Debug Mode  
                logger_i.debug("\n"+"\n"+"\n"+'In scale matrices (scale_matrices.py)')
                logger_i.debug('----------------------------------------------------------------')
        
        '''
        Remark  
        scipy.sparse.linalg.onenormest - Computes a lower bound of the 1-norm of a sparse matrix.
        In the disk brake modelling theory they have used spectral norm but since I could not
        find a better(less computational cost) way to obtain this in python, I have used 1-norm
        approximation
        '''
        m_norm=norm.onenormest(m, t=3, itmax=5, compute_v=False, compute_w=False)
        c_norm=norm.onenormest(c, t=3, itmax=5, compute_v=False, compute_w=False)
        k_norm=norm.onenormest(k, t=3, itmax=5, compute_v=False, compute_w=False)
        
        gamma=math.sqrt(k_norm/m_norm);
        delta=2/(k_norm+c_norm*gamma);
        M=gamma*gamma*delta*m;
        C=gamma*delta*c;
        K=delta*k;
        
        if(LOG_LEVEL==10): #Debug Mode
                m_scaled_norm = gamma*gamma*delta*m_norm
                c_scaled_norm = gamma*delta*c_norm
                k_scaled_norm = delta*k_norm
                rho=max(m_scaled_norm,c_scaled_norm,k_scaled_norm)/min(m_scaled_norm,k_scaled_norm)
                logger_i.debug('Scaling factors gamma = '+str(gamma)+' delta = '+str(delta))
                logger_i.debug('The scaled problem has rho = '+str(rho))
                
        return M, C, K, gamma;
