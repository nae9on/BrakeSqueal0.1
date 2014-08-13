r"""
This module defines the following functions::

  - BrakeClass:
  
    Python Parent Class for the BrakeSqueal Project
    
  - print_eigs:
    
    Prints all the eigenvalues on the terminal (with two floating points).
    
  - print_target_eigs:
  
    Prints the eigenvalues in the target region on the terminal (when flag is false)
    , in the info file (when flag is true).
  
    
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import logging

#----------------------------Application Specific Imports----------------------------------------- 
from initialize import logger


__all__ = [
    'BrakeClass'
    ]

class BrakeClass:
   r"""
        
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

def print_eigs(arg):
    r"""
        
    INPUT:
       
    - ``arg`` -- eigenvalues
             
    OUTPUT:
        
    - prints eigenvalues on the terminal
    
    """
    
    n = arg.shape[0]
    arg=sorted(arg,reverse=True)
    print("\n")
    for i in range(0, n):
       # print "%.5f" % arg[i].real, '+',  "%.5f" % arg[i].imag, 'I'
       if arg[i].imag < 0:
          print("%04.03e" % arg[i].real, '-',  "%04.03e" % abs(arg[i].imag), 'I')
       else:
          print("%04.03e" % arg[i].real, '+',  "%04.03e" % arg[i].imag, 'I')
    print("\n")

def print_target_eigs(obj,arg,flag):
    r"""
        
    INPUT:
    
    - ``obj`` -- object of the class ``BrakeClass``
    - ``arg`` -- eigenvalues
    - ``flag`` -- 0/1 for output on the terminal/in info file
             
    OUTPUT:
        
    - prints eigenvalues on the terminal when flag = 0 else prints output in the 
      info file.
    
    """
    
    logger_i = obj.logger_i
    n = arg.shape[0]
    arg=sorted(arg,reverse=True)
    
    for i in range(0, n):
      if arg[i].real > 0:
       if arg[i].imag < 0:
          if(flag):
            x = "%04.03e" % arg[i].real,  "%04.03e" % arg[i].imag
            logger_i.info(x)
          else:
            print("%04.03e" % arg[i].real, '-',  "%04.03e" % abs(arg[i].imag), 'I')
       else:
          if(flag):
            x = "%04.03e" % arg[i].real,  "%04.03e" % arg[i].imag
            logger_i.info(x)
          else:
            print("%04.03e" % arg[i].real, '+',  "%04.03e" % arg[i].imag, 'I')
      else:
        break
