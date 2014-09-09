r"""
This module defines the following functions::

  - BrakeClass:
  
    Python Parent Class for the BrakeSqueal Project
    
  - print_eigs:
    
    Prints all the eigenvalues on the terminal (with two floating points).
    
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import logging
import numpy as np

#----------------------------Application Specific Imports----------------------------------------- 
from initialize import logger


__all__ = [
    'BrakeClass',
    'print_eigs'
    ]

class BrakeClass:
   r"""
   
   ``Member Functions of the BrakeClass:``
   
   - __init__
   - createInfoLogger
   - createTimeLogger
   - displayCount
   - displayParametersConsole
   - displayParametersLog
   
   """
   
   objCount = 0

   def __init__(self, input_path, output_path, info_log_file, time_log_file,
                      log_level):
                      
      self.input_path = input_path
      self.output_path = output_path
      self.info_log_file = info_log_file
      self.time_log_file = time_log_file
      self.log_level = log_level
      
      BrakeClass.createInfoLogger(self)
      BrakeClass.createTimeLogger(self)
     
      BrakeClass.objCount += 1
      
   def createInfoLogger(self):
     self.logger_i = logger.return_info_logger(self)
     
   def createTimeLogger(self):
     self.logger_t = logger.return_time_logger(self) 
     
   def displayCount(self):
     print("Total Object count %d" % BrakeClass.objCount)

   def displayParametersConsole(self):
      attrs = vars(self)
      print("\n",'BrakeClass Attributes and there values')
      for item in attrs.items():
          print(item[0], ' = ', item[1])
          
   def displayParametersLog(self):
      attrs = vars(self)
      self.logger_i.info("\n"+'BrakeClass Attributes and there values')
      for item in attrs.items():
          self.logger_i.info(item[0]+'   '+str(item[1]))       

def print_eigs(obj,arg,which_flag,where_flag):
    r"""
        
    :param obj: object of the class ``BrakeClass``
    :param arg: eigenvalues
    :param which_flag: which eigenvalues are needed('all' or 'target' or 'critical' or 'positive')
    :param where_flag: where to print the eigenvalues('terminal' or 'file')
    :return: prints eigenvalues in the desired format(upto 2 decimal places)
    
    """
    
    n = arg.shape[0]
    arg=sorted(arg,reverse=True)
    logger_i = obj.logger_i
    
    #create a temporary bool array to indicate required eigenvalues
    bool_array = np.zeros(n)
    
    if which_flag == 'all': #no condition
      bool_array = np.ones(n)
      
    if which_flag == 'target':
     for i in range(0, n):
      if (arg[i].real >= obj.target[0] and arg[i].real <= obj.target[1] and \
          arg[i].imag >= obj.target[2] and arg[i].imag <= obj.target[3]):
           bool_array[i] = 1
           
    if which_flag == 'critical':
     for i in range(0, n):
      if (arg[i].real >= obj.target[0] and arg[i].real <= obj.target[1] and \
          arg[i].imag >= obj.target[2] and arg[i].imag <= obj.target[3] and \
          arg[i].real >= 0):
           bool_array[i] = 1

    if which_flag == 'positive':
     for i in range(0, n):
      if (arg[i].real >= 0):
           bool_array[i] = 1
           
    
    for i in range(0, n):
     if bool_array[i] == 1:
       if arg[i].imag < 0:
          if(where_flag == 'file'):
            x = "%04.03e" % arg[i].real,  "%04.03e" % arg[i].imag
            logger_i.info(x)
          if(where_flag == 'terminal'):
            print("%04.03e" % arg[i].real, '-',  "%04.03e" % abs(arg[i].imag), 'I')
       else:
          if(where_flag == 'file'):
            x = "%04.03e" % arg[i].real,  "%04.03e" % arg[i].imag
            logger_i.info(x)
          if(where_flag == 'terminal'):
            print("%04.03e" % arg[i].real, '+',  "%04.03e" % arg[i].imag, 'I')
