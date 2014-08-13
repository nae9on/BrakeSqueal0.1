'''
Unliniarize GEVP

DiagScaleMatrices
Input
1. evec
2. n
3. numev

Output
1. B

'''

import numpy as np
from brake.initialize import diagscale


def unlinearize_matrices(evec,n,numev):
	evec_top = evec[0:n]
	evec_bottom = evec[n:2*n]
	CC = evec.copy()
	CC = np.square(CC)
	diag_top = CC[0:n].sum(axis=0)	
	diag_bottom = CC[n:2*n].sum(axis=0)
	diff = diag_top - diag_bottom
        for i in range(0,numev):
		if (diff[i]<0):
		  evec_top[:,i:i+1] = evec_bottom[:,i:i+1]
	D = diagscale.normalize_cols(evec_top)
	B = evec_top*D
	return B;
