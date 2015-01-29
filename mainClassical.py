# Driver File
# Compares Tradional Approach with POD Approach

#----------------------------------Standard Library Imports---------------------------------------
import math
import numpy
import timeit
import socket
import datetime
import matplotlib.pyplot as plt

#----------------------------Application Specific Imports-----------------------------------------
import brake
import createBrakeClassObject
from brake.initialize import load, assemble, scale, diagscale
from brake.solve import projection, classicalProjection, solver, qevp
from brake.analyze import residual, visual


#----------------------------Create Object--------------------------------------------------------
obj = createBrakeClassObject.returnObject()


#Log all the BrakeClass attribute values in the info log file
if(obj.log_level):
  obj.displayParametersLog()
  
'''  
####################################################################
####################################################################
####################### Exact Eigenvalues ##########################
####################################################################
####################################################################
  
#Calculating exact eigenpairs of the QEVP for omega = omegaTest in the target region
#--------------------------------------------------------------------------------------------------
la_exact, evec_exact = qevp.brake_squeal_qevp(obj, 1, obj.omegaTest)

#print '\n\n Eigenvalues in the target region are'
#brake.print_eigs(obj,la_exact,'target','terminal')

#Extracting eigenpairs in the target region from the computed eigenpairs
laTarget, evecTarget = brake.extractEigs(obj,la_exact,evec_exact,'target')

  
####################################################################
####################################################################
####################### POD Approach ###############################
####################################################################
####################################################################
  
#Begin creating the projection matrix
#################################################
if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix'

begin = timeit.default_timer()
Q = projection.obtain_projection_matrix(obj)
end = timeit.default_timer()
#------------------------------------------------done creating the projection matrix

sparse_list = load.load_matrices(obj)
print '\n\n Testing for omega = ',str(obj.omegaTest/(2*math.pi))
M, C, K = assemble.create_MCK(obj, sparse_list, obj.omegaTest)
print 'Projecting the QEVP having dimension '+str(M.shape)
#Projection
QT = Q.T.conjugate()
M =  QT.dot(M.dot(Q))
C =  QT.dot(C.dot(Q))
K =  QT.dot(K.dot(Q))
print 'Onto a smaller subspace of dimension '+str(M.shape)

n = M.shape[0]
no_of_evs = 2*n-2;
la_pod, evec_pod = solver.qev_dense(obj,M,C,K,no_of_evs)

res_qevp = residual.residual_qevp(M,C,K,la_pod,evec_pod[0:n,:])
print 'Maximum Residual error for the QEVP is ',max(res_qevp)

print '\n\n Eigenvalues in the target region computed using projection'
brake.print_eigs(obj,la_pod,'target','terminal')

print '\n\n\nCalculating Delta'
delta = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)
for itr in range(0,len(laTarget)):
  delta[itr] = numpy.min(numpy.absolute(la_pod-laTarget[itr]))/numpy.absolute(laTarget[itr])
  
print delta 
print '\n\n',numpy.max(delta)

'''

####################################################################
####################################################################
####################### Classical Approach #########################
####################################################################
####################################################################
  
#Begin creating the projection matrix
#################################################
if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix'

begin = timeit.default_timer()
Q = classicalProjection.obtain_projection_matrix(obj)
end = timeit.default_timer()
#------------------------------------------------done creating the projection matrix

sparse_list = load.load_matrices(obj)
print '\n\n Testing for omega = ',str(obj.omegaTest/(2*math.pi))
M, C, K = assemble.create_MCK(obj, sparse_list, obj.omegaTest)
print 'Projecting the QEVP having dimension '+str(M.shape)
#Projection
QT = Q.T.conjugate()
M =  QT.dot(M.dot(Q))
C =  QT.dot(C.dot(Q))
K =  QT.dot(K.dot(Q))
print 'Onto a smaller subspace of dimension '+str(M.shape)

n = M.shape[0]
no_of_evs = 2*n-2;
la_pod, evec_pod = solver.qev_dense(obj,M,C,K,no_of_evs)

res_qevp = residual.residual_qevp(M,C,K,la_pod,evec_pod[0:n,:])
print 'Maximum Residual error for the QEVP is ',max(res_qevp)

print '\n\n Eigenvalues in the target region computed using projection'
brake.print_eigs(obj,la_pod,'target','terminal')

print '\n\n\nCalculating Delta'
delta = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)
for itr in range(0,len(laTarget)):
  delta[itr] = numpy.min(numpy.absolute(la_pod-laTarget[itr]))/numpy.absolute(laTarget[itr])
  
print delta 
print '\n\n',numpy.max(delta)


  
  
  
  