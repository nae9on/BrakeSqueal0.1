r"""
This module defines the following functions::

  - return_info_logger:
    
    Creates and returns a python logger object for information logging
    
  - return_time_logger:
    
    Creates and returns a python logger object for time logging 

Note::

 various logging levels
 LEVELS = {'notset':logging.NOTSET, #0 --> numerical value (for no logging)
          'debug': logging.DEBUG, #10 (to capture detailed debug information)
          'info': logging.INFO, #20 (to capture essential information)
          'warning': logging.WARNING, #30
          'error': logging.ERROR, #40
          'critical': logging.CRITICAL} #50

"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import logging


def return_info_logger(obj):
    r"""
        
    :param obj: object of the class ``BrakeClass``
    :return: ``logger_i`` python logger object for information logging
        
    """ 
    
    INFO_LOG_FILENAME = obj.info_log_file
    TIME_LOG_FILENAME = obj.time_log_file
    LOG_LEVEL = obj.log_level
    
    # creating info logger 
    formatter = logging.Formatter('%(message)s')
    logger_i = logging.getLogger('info_logger')
    hdlr_i = logging.FileHandler(INFO_LOG_FILENAME,mode='w')    
    hdlr_i.setFormatter(formatter)
    logger_i.addHandler(hdlr_i)
    logger_i.setLevel(LOG_LEVEL)
    
    if(LOG_LEVEL == 0):
        logger_i.setLevel(logging.INFO)
        logger_i.info('The current level of logging is '+str(LOG_LEVEL))
        logger_i.info('To enable detailed logging set log level to non zero during object creation:')
	logger_i.info('Set log level = 10(debug mode), log level = 20(info mode)')
        logger_i.info('The Info Data is logged in the file : '+INFO_LOG_FILENAME)
        logger_i.info('The Time Data is logged in the file : '+TIME_LOG_FILENAME)       
        logger_i.setLevel(logging.NOTSET)

    if(LOG_LEVEL):
        logger_i.info('Auto logged Info for the Quadratic Eigen Value Problem(Brake Squeal Project)')
        logger_i.info("\n"+'In Driver file')
    
    return logger_i

def return_time_logger(obj):
    r"""
        
    :param obj: object of the class ``BrakeClass``
    :return: ``logger_t`` python logger object for time logging
        
    """
    
    TIME_LOG_FILENAME = obj.time_log_file
    INFO_LOG_FILENAME = obj.info_log_file
    LOG_LEVEL = obj.log_level

    # creating time logger
    formatter = logging.Formatter('%(message)s') 
    logger_t = logging.getLogger('time_logger')
    hdlr_t = logging.FileHandler(TIME_LOG_FILENAME,mode='w')
    hdlr_t.setFormatter(formatter)
    logger_t.addHandler(hdlr_t)
    logger_t.setLevel(LOG_LEVEL)
    
    if(LOG_LEVEL == 0):
        logger_t.setLevel(logging.INFO)
        logger_t.info('The current level of logging is '+str(LOG_LEVEL))
        logger_t.info('To enable detailed logging set log level to non zero during object creation')
	logger_t.info('Set log level = 10(debug mode), log level = 20(info mode)')
        logger_t.info('The Info Data is logged in the file : '+INFO_LOG_FILENAME)
        logger_t.info('The Time Data is logged in the file : '+TIME_LOG_FILENAME)       
        logger_t.setLevel(logging.NOTSET)

    if(LOG_LEVEL):
        logger_t.info('Auto logged Info for the Quadratic Eigen Value Problem(Brake Squeal Project)')
        logger_t.info("\n"+'Time = 0')
        
    return logger_t

