r"""
This module defines the following functions::

  - load_matrices:
  
    This function loads the sparse data matrices in .mat format and adds them to a python list
    
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import scipy.io

def load_matrices(obj):
        r"""
        
        :param obj: object of the class ``BrakeClass``
        :return: sparse_list -- a python list of matrices in Compressed Sparse Column format 
          of type '<type 'numpy.float64'>'
        
        Procedure::
        
         The sparse_list is obtained by loading the vaious .mat files present in the
         data_file_list attribute of the BrakeClass and then appending them into
         a python list sparse_list.
        
        """
        
        LOG_LEVEL = obj.log_level
        logger_t = obj.logger_t
        logger_i = obj.logger_i
        
        path = obj.input_path
        data_file_list = obj.data_file_list
        
        if(LOG_LEVEL==10): #Debug Mode
          logger_i.debug("\n"+"\n"+"\n"+'In load matrices (load_matrices.py)')
          logger_i.debug('-------------------------------------------------------------------')
          logger_i.debug('loading the following .mat files '+str(data_file_list))
          logger_i.debug('load returns a list of matrices in Compressed Sparse Column format of type numpy.float64')
                
        sparse_list = []
        
        for file in data_file_list:
            # ex f_list = sio.loadmat(path+'BMLL00.mat')
            f_list = scipy.io.loadmat(path+file+'00.mat')
            # extract = f_list['BMLL']
            extract = f_list[file]                                    
            sparse_list.append(extract)
            
        return sparse_list






