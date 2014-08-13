'''
Function for loading the .mat data files

Input
1. path (where the data files are located)
2. data_file_list (data file names for the component matrix m,c1,c2,c3,c4,k1,k2,k3)

Output
1. sparse_list (a python list of matrices in Compressed Sparse Column 
		format of type '<type 'numpy.float64'>')

'''

import scipy.io as sio

def load_matrices(obj):
        
        LOG_LEVEL = obj.log_level
	logger_t = obj.logger_t
	logger_i = obj.logger_i
	
	path = obj.input_path
        data_file_list = obj.data_file_list
	
	if(LOG_LEVEL==10): #Debug Mode
		logger_i.debug("\n"+"\n"+"\n"+'In load matrices (load_matrices.py)')
		logger_i.debug('-------------------------------------------------------------------')
		logger_i.debug('loading the following .mat files '+str(data_file_list))
		logger_i.debug('load returns a list of matrices in Compressed Sparse Column format \
		of type numpy.float64')
		
	sparse_list = []
	
	for file in data_file_list:
	    f_list = sio.loadmat(path+file+'00.mat') # ex f_list = sio.loadmat(path+'BMLL00.mat')
            extract = f_list[file]			  # extract = f_list['BMLL']		
            sparse_list.append(extract)
            
	return sparse_list






