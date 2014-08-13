r"""
This module defines the following functions::

  - unlinearize_matrices:
    
    To obtain the eigenvectors prior linearization of the QEVP from the 
    resulting eigenvectors (after classical companion linearization of the QEVP
    to obtain the generalized eigenvalue propblem).
  
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import numpy
from brake.initialize import diagscale


def unlinearize_matrices(evec):
    r"""
        
    INPUT:
       
    - ``evec`` -- eigenvectors of the generalized eigenvalue propblem
                 
    OUTPUT:
        
    - ``evec_prior`` -- eigenvectors prior linearization of the QEVP
        
    The ``evec_prior`` are obtained as follows:
    
    - check for every vector ``i`` of the GEVP(``evec``) (consider MATLAB convention)
    - if ``norm(evec[1:n,i]) > norm(evec[n+1:2*n,i]) set evec_prior[:,i] = evec[1:n,i]``
    - else ``set evec_prior[:,i] = evec[n+1:2*n,i]``

    """
    
    #dimension of the QEVP
    n = evec.shape[0]/2
    
    #number of eigenvectors
    numev = evec.shape[1]
    
    evec_top = evec[0:n]
    evec_bottom = evec[n:2*n]
    CC = evec.copy()
    CC = numpy.square(CC)
    diag_top = CC[0:n].sum(axis=0)  
    diag_bottom = CC[n:2*n].sum(axis=0)
    diff = diag_top - diag_bottom
    for i in range(0,numev):
       if (diff[i]<0):
          evec_top[:,i:i+1] = evec_bottom[:,i:i+1]
    D = diagscale.normalize_cols(evec_top)
    evec_prior = evec_top*D
    return evec_prior;
