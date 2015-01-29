''' 
1. Compares Tradional Approach with POD Approach.
2. Does Convergence anaylysis with increasing dimension of the projection subspace.
3. Compares POD approach for different samples.
'''
#----------------------------------Standard Library Imports---------------------------------------
import os
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

####################################################################
####################################################################
####################### Exact Eigenvalues ##########################
####################################################################
####################################################################

#Calculating exact eigenpairs of the QEVP for omega = omegaTest in the target region
#--------------------------------------------------------------------------------------------------
print "\n"+"\n"+'Calculating exact eigenpairs of the QEVP for omega = '+str(obj.omegaTest/(2*math.pi))+' in the target region'
laExact, evecExact = qevp.brake_squeal_qevp(obj, 1, obj.omegaTest)
#brake.printEigs(obj,la_exact,'target','terminal')
#Extracting eigenpairs in the target region from the computed eigenpairs
laTarget, evecTarget = brake.extractEigs(obj, laExact, evecExact, 'target')

####################################################################
####################################################################
####################### POD Approach ###############################
####################################################################
####################################################################

#Begin creating the Measurment matrix
#################################################
if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Setup Phase: Creating the Measurment Matrix')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Setup Phase: Creating the Measurment Matrix'
begin = timeit.default_timer()

X = []

obj.omega_basis = numpy.array([1,20])*2*math.pi
X_wset1 = projection.obtain_measurment_matrix(obj)
X.append(X_wset1)

obj.omega_basis = numpy.array([1,5,10,15,20])*2*math.pi
X_wset2 = projection.obtain_measurment_matrix(obj)
X.append(X_wset2)

obj.omega_basis = numpy.array([1,2.5,5,7.5,10,12.5,15,17.5,20])*2*math.pi
X_wset3 = projection.obtain_measurment_matrix(obj)
X.append(X_wset3)


end = timeit.default_timer()
#------------------------------------------------done creating the Measurment matrix

#Obtaining the original M,C,K matrices
sparse_list = load.load_matrices(obj)
print '\n\n Testing for omega = ',str(obj.omegaTest/(2*math.pi))
M_orig, C_orig, K_orig = assemble.create_MCK(obj, sparse_list, obj.omegaTest)

#Relative error plot with increasing dimension of the projection matrix
xDim = numpy.arange(10,310,10)

yDeltaList = []

for setItr in range(0,len(X)):
	
	yDelta = numpy.zeros(xDim.shape[0],dtype=numpy.float64)
	
	for i1 in range(0,len(xDim)):
		
		obj.projectionDimension = xDim[i1]
		Q = projection.obtain_projection_matrix(obj,X[setItr])
	        
		print 'wset '+str(setItr)+' Projecting the QEVP having dimension '+str(M_orig.shape)
		
		#Projection
		QT = Q.T.conjugate()
		M =  QT.dot(M_orig.dot(Q))
		C =  QT.dot(C_orig.dot(Q))
		K =  QT.dot(K_orig.dot(Q))
		print 'Onto a smaller subspace of dimension '+str(M.shape)
	
		n = M.shape[0]
		no_of_evs = 2*n-2;
		la_pod, evec_pod = solver.qev_dense(obj,M,C,K,no_of_evs)
		
		'''
		res_qevp = residual.residual_qevp(M,C,K,la_pod,evec_pod[0:n,:])
		print 'Maximum Residual error for the QEVP is ',max(res_qevp)
		
		print '\n\n Eigenvalues in the target region computed using projection'
		brake.printEigs(obj,la_pod,'target','terminal')
		'''
		
		#print '\n\n\nCalculating Delta'
		delta = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)
		for itr in range(0,len(laTarget)):
		    delta[itr] = numpy.min(numpy.absolute(la_pod-laTarget[itr]))
		
		yDelta[i1] = numpy.max(delta)/numpy.max(numpy.absolute(laTarget))
		
        yDeltaList.append(yDelta)

print '\n\n The delta values for the three sets are as follows \n\n'
for setItr in range(0,len(X)):
	print yDeltaList[setItr]

plt.figure(1)
plt.xlabel('dimension')
plt.ylabel('relative error')
plt.title('logarithmic decay of relative error with dimension')
plt.plot(xDim, yDeltaList[0], '-ro')
plt.plot(xDim, yDeltaList[1], '-bD')
plt.plot(xDim, yDeltaList[2], '-gs')
plt.yscale('log')
plt.grid(True)
brake.save(obj.output_path+'errorDecay', ext="png", close=True, verbose=False)
