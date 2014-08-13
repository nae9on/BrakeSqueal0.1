import logging
from initialize import logger
__all__ = [
    'BrakeClass'
    ]

class BrakeClass:
   'base class for the brake squeal project'
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
     print "Total Object count %d" % BrakeClass.objCount

   def displayParameters(self):
      attrs = vars(self)
      for item in attrs.items():
          print item

def my_print(arg):
    n = arg.shape[0]
    arg=sorted(arg,reverse=True)
    print "\n"
    for i in range(0, n):
       # print "%.5f" % arg[i].real, '+',  "%.5f" % arg[i].imag, 'I'
       if arg[i].imag < 0:
          print "%04.03e" % arg[i].real, '-',  "%04.03e" % abs(arg[i].imag), 'I'
       else:
          print "%04.03e" % arg[i].real, '+',  "%04.03e" % arg[i].imag, 'I'
    print "\n"

def my_print_few_infile(obj,arg,flag):
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
	    print "%04.03e" % arg[i].real, '-',  "%04.03e" % abs(arg[i].imag), 'I'
       else:
          if(flag):
	    x = "%04.03e" % arg[i].real,  "%04.03e" % arg[i].imag
	    logger_i.info(x)
	  else:
	    print "%04.03e" % arg[i].real, '+',  "%04.03e" % arg[i].imag, 'I'		
      else:
        break
