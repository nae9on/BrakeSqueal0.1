r"""
This module defines the following functions::

  - obtain_projection_matrix:
  
   This function forms the Projection Matrix by solving the quadratic eigenvalue problem 
   for each base angular frequency.
   
   - obtain_measurment_matrix
   
   This function forms the Measurment Matrix by solving the quadratic eigenvalue problem 
   for each base angular frequency.
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import timeit
import numpy
import scipy

#----------------------------Application Specific Imports----------------------------------------- 
from brake.solve import qevp

def obtain_projection_matrix(obj,X):
        r"""
        
        :param obj: object of the class ``BrakeClass``
	:param X: measurment matrix
        :return: ``Q`` - projection matrix
        
        Procedure::
        
         The projection matrix is obtained as follows:
        
         - Compute the thin svd of the measurment matrix. X = U * s * V
         - Set Q = truncated(U), where the truncation is done to take only the significant 
           singular values(provided by user in obj.projectionDimension) into account

        """
        
        #object attributes used in the function
        LOG_LEVEL = obj.log_level
        logger_t = obj.logger_t
        logger_i = obj.logger_i
        omega_basis = obj.omega_basis
 
        if(LOG_LEVEL):
          logger_i.info("\n"+'Creating the Measurment matrix X and Projection matrix Q')

        #X = obtain_measurment_matrix(obj)
        '''
        Obtaining the projection matrix Q
        1. Compute the thin svd of the measurment matrix. X = U * s * V
        2. Set Q = truncated(U), where the truncation is done to take only the 
           significant singular values(based on a certain tolerance) into account
        '''
        start_svd = timeit.default_timer()
        U, s, V = scipy.linalg.svd(X, full_matrices=False)
	stop_svd = timeit.default_timer()
	
	if(LOG_LEVEL):
           logger_i.info("\n"+'The shapes of U,s,V of the measurment matrix are as follows '\
                                                +str(U.shape)+' '+str(s.shape)+' '+str(V.shape))
           logger_t.info("\n"+"\n"+'\t\tTime taken to compute svd = '+"%.2f" % (stop_svd-start_svd))    

        
	if obj.projectionDimension > s.shape[0]:
	 print 'Warning: projection dimension exceeds the no of available singular values '
	
	s_truncated = s[0:obj.projectionDimension]
        Q = U[:,0:obj.projectionDimension]
		
	print '\nExtracted '+str(obj.projectionDimension)+' significant singular values from '+str(s.shape[0])+' to obtain projection matrix of dimension '+str(Q.shape)
        
        if(LOG_LEVEL):
           logger_i.info('The no of singular values = '+str(s.shape[0]))
           logger_i.info('The no of singular values after truncation = '+str(s_truncated.shape[0]))
           logger_i.info('The dimensions of the projection matrix = '+str(Q.shape))

        return Q, s

def obtain_measurment_matrix(obj):
        r"""
        
        :param obj: object of the class ``BrakeClass``
        :return: ``X`` - measurment matrix
        
        Procedure::
        
         The measurment matrix is obtained as follows:
        
         - Obtain the measurment matrix X = [X_real X_imag], with X_real as a list of 
           real parts of eigenvectors and X_imag as a list of imaginary parts of 
           eigenvectors, corresponding to each base angular frequency in omega_basis.
        """

        #object attributes used in the function
        LOG_LEVEL = obj.log_level
        logger_t = obj.logger_t
        logger_i = obj.logger_i
        omega_basis = obj.omega_basis
        
        
        ''' Creating the measurment matrix with X_real as a list of real part of eigenvectors 
        and X_imag as a list of imaginary part of eigenvectors corresponding to all the base 
        frequencies 
        '''
        X_real = []
        X_imag = []
        for i in range(0,len(omega_basis)):
                
                begin = timeit.default_timer()
                la, evec = qevp.brake_squeal_qevp(obj, i+1, omega_basis[i])
                end = timeit.default_timer()
                
                if(LOG_LEVEL):
                 logger_t.info("\n"+"\n"+'\t\tTotal time taken by Brake Squeal Qevp = '+"%.2f" % (end-begin)+' sec')
                        
                X_real.append(evec.real)
                X_imag.append(evec.imag)
        
        
        if(LOG_LEVEL):
         logger_i.info("\n"+"\n")
         logger_i.info('------------------ Now creating the Measurment matrix --------')
         logger_i.info('The number of QEVP solved = '+str(len(omega_basis)))

        for i in range(0,len(omega_basis)):
                if(LOG_LEVEL):
                   logger_i.info("\n"+'The number of eigenvectors computed for frequency '\
                                                        +str(i+1)+' = '+str(X_real[i].shape[1]))
        '''
        The Measurment matrix X = [X_real X_imag]
        '''
        #adding the real part of computed eigenvectors to the measurment matrix
        X = X_real[0]
        for i in range(1,len(X_real)):
                X = numpy.concatenate((X,X_real[i]), axis=1)
        #adding the imaginary part of computed eigenvectors to the measurment matrix
        for i in range(0,len(X_imag)):
                X = numpy.concatenate((X,X_imag[i]), axis=1)
        
	return X
