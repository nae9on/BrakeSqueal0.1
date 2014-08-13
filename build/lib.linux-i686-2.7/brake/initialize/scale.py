'''
Scale QEVP
Input
1. M
2. C
3. K
4. tau (shift point)

Output
1. M
2. C
3. K
4. gamma
5. delta

'''

import math
import scipy.sparse.linalg as norm

def scale_matrices(obj, M , C, K):

	LOG_LEVEL = obj.log_level
	logger_t = obj.logger_t
	logger_i = obj.logger_i
        
	if(LOG_LEVEL==10): #Debug Mode	
		logger_i.debug("\n"+"\n"+"\n"+'In scale matrices (scale_matrices.py)')
		logger_i.debug('----------------------------------------------------------------')
	
	'''
	Remark	
	scipy.sparse.linalg.onenormest - Computes a lower bound of the 1-norm of a sparse matrix.
	In the disk brake modelling paper they used spectral norm but since I could not
	find a better(less computational cost) way to obtain this in python, I have used 1-norm
	approximation
	'''
	m_norm=norm.onenormest(M, t=3, itmax=5, compute_v=False, compute_w=False)
        c_norm=norm.onenormest(C, t=3, itmax=5, compute_v=False, compute_w=False)
        k_norm=norm.onenormest(K, t=3, itmax=5, compute_v=False, compute_w=False)
        gamma=math.sqrt(k_norm/m_norm);
        delta=2/(k_norm+c_norm*gamma);
        M=gamma*gamma*delta*M;
        C=gamma*delta*C;
        K=delta*K;
        
	if(LOG_LEVEL==10): #Debug Mode
		m_scaled_norm = gamma*gamma*delta*m_norm
		c_scaled_norm = gamma*delta*c_norm
		k_scaled_norm = delta*k_norm
		rho=max(m_scaled_norm,c_scaled_norm,k_scaled_norm)/min(m_scaled_norm,k_scaled_norm)
		logger_i.debug('Scaling factors gamma = '+str(gamma)+' delta = '+str(delta))
	        logger_i.debug('The scaled problem has rho = '+str(rho))
	        
        return M, C, K, gamma, delta;
