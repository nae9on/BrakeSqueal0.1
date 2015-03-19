# Creates and returns the Pseudospectrum Object

#----------------------------------Standard Library Imports---------------------------------------
import os
import math
import numpy
import socket
import timeit
import datetime

#----------------------------Application Specific Imports-----------------------------------------
import brake

def returnObject(logLevel):

    begin_program = timeit.default_timer()


  #Logging Parameters
  '''
  'notset':logging.NOTSET, #0 (for no logging - fastest)
  'debug': logging.DEBUG, #10 (to capture detailed debug information - slowest)
  'info': logging.INFO, #20 (to capture essential information)
  '''

  log_level = logLevel
  dt = datetime.date.today()

  #Create output directory if it does not exist
  if not os.path.exists(output_path):
    os.makedirs(output_path)

  '''
  This would create the first object of BrakeClass with certain limited attributes.
  More attributes can be added anywhere in the program using the setattr command.
  These attributes can then be referred later in a module.
  '''
  obj = brake.BrakeClass(input_path, output_path, info_log_file, time_log_file, log_level)

  '''
  Adding more attributes to the BrakeClass object 'obj' created above.

  Problem parameters described below
  M = m
  C = c1+c2*(omega/omegaRef)+c3*(omegaRef/omega)
  K = k1+k2+k3*math.pow((omega/omegaRef),2)
  [m c1 c2 c3 c4 k1 k2 k3]
  ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL']
  [M D1 DG DR D4 K1 KR KGEO]

  omegaRef - reference omega.
  fRef - reference frequency.
  #omega_basis - array of base frequencies that would be used for creating the measurment matrix
  omega_range - array of frequencies for simulating the brake squealing(using projection matrix from omega_basis)
  target - target rectangular region.


  projectionDimension(dimension of the projection matrix) - The number of significant
  singular values to be taken into consideration while creating the projection matrix
  from the measurment matrix.

  evs_per_shift - number of eigenvalues to be calculated per shift
  desired_area_fraction - area fraction of the target region to be covered
  desiredCount - the minimum count of the eigenvalues to be captured for each frequency
  omegaTest - test omega for comparison of the traditional approach with the POD approach
  '''

  setattr(obj, 'data_file_list', ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL'])
  setattr(obj, 'data_file_name', ['M1','D1','DG','DR','D4','K1','KR','KGeo'])
  setattr(obj, 'omegaRef', omegaRef)
  setattr(obj, 'fRef', 1600.0)
  obj.M=lambda M1,omega: M1
  obj.C=lambda D1,DG,DR,D4,omega:  D1+DG*((omega/obj.omegaRef))+DR*((obj.omegaRef/omega)-1)+D4/(2.0*math.pi*obj.fRef)
  obj.K=lambda K1,KR,KGeo,omega: K1+KR+KGeo*(math.pow((omega/obj.omegaRef),2)-1)




  setattr(obj, 'omega_basis', numpy.array([1,4])*2.0*math.pi)
  #setattr(obj, 'omega_basis', numpy.array([1,5,10,15,20])*2*math.pi)
  #setattr(obj, 'omega_basis', numpy.array([1,2.5,5,7.5,10,12.5,15,17.5,20])*2*math.pi)

  setattr(obj, 'omega_range', numpy.linspace(1, 4, num=3)*2.0*math.pi)
  setattr(obj, 'target', numpy.array([-10,1000,-50,42000]))
  setattr(obj, 'projectionDimension', 100)
  setattr(obj, 'evs_per_shift', 50)
  setattr(obj, 'desired_area_fraction', 0.99)
  setattr(obj, 'desiredCount', 0)
  setattr(obj, 'omegaTest', 2*2.0*math.pi)

  #Log all the BrakeClass attribute values in the info log file
  if(obj.log_level):
     obj.displayParametersLog()

  return obj

