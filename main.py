''' 
1. Compares Tradional Approach with POD Approach.
2. Convergence anaylysis test with increasing dimension of the projection matrix.
3. Compares POD approach for different size of the samples.
'''
#----------------------------------Standard Library Imports---------------------------------------
import os
import math
import numpy
import timeit
import socket
import datetime
from pylab import *
import numpy.linalg as LA
from numpy.linalg import norm
import matplotlib.pyplot as plt

#----------------------------Application Specific Imports-----------------------------------------
import brake
import createBrakeClassObject
from brake.initialize import load, assemble, scale, diagscale, unlinearize
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
print '\n\n----------------------------------------------------------------------------------------'
print "\n"+"\n"+'Calculating exact eigenpairs of the QEVP for omega = '+str(obj.omegaTest/(2*math.pi))+' in the target region \n'
laExact, evecExact = qevp.brake_squeal_qevp(obj, 1, obj.omegaTest)
#brake.printEigs(obj,la_exact,'target','terminal')
#Extracting eigenpairs in the target region from the computed eigenpairs
laTarget, evecTarget = brake.extractEigs(obj, laExact, evecExact, 'target')
  
  
####################################################################
####################################################################
####################### Classical Projection #######################
####################################################################
####################################################################
print '\n\n----------------------------------------------------------------------------------------'
laClassical, evecClassical = classicalProjection.Obtain_eigs(obj,150)
XClassical = numpy.concatenate((evecClassical.real,evecClassical.imag), axis=1)

#X is list of measurment matrices form which projection matrices of varying dimension will
#be later extracted in the loop using singular value decomposition
X = []
X.append(XClassical)

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

print '\n\nBeginning Setup Phase: Creating the Measurment Matrix'

print '\n\nCreating the Measurment Matrix for wset1'
obj.evs_per_shift = 75
obj.omega_basis = numpy.array([1,20])*2*math.pi
X_wset1 = projection.obtain_measurment_matrix(obj)
X.append(X_wset1)

print '\n\nCreating the Measurment Matrix for wset2'
obj.evs_per_shift = 30
obj.omega_basis = numpy.array([1,5,10,15,20])*2*math.pi
X_wset2 = projection.obtain_measurment_matrix(obj)
X.append(X_wset2)

print '\n\nCreating the Measurment Matrix for wset3'
obj.evs_per_shift = 17
obj.omega_basis = numpy.array([1,2.5,5,7.5,10,12.5,15,17.5,20])*2*math.pi
X_wset3 = projection.obtain_measurment_matrix(obj)
X.append(X_wset3)

#------------------------------------------------done creating the Measurment matrix

#Obtaining the original M,C,K matrices
sparse_list = load.load_matrices(obj)
print '\n\nTesting for omega = '+str(obj.omegaTest/(2*math.pi))+'\n\n'
M_orig, C_orig, K_orig = assemble.create_MCK(obj, sparse_list, obj.omegaTest)

#Relative error plot with increasing dimension of the projection matrix
xDim = numpy.arange(10,310,10)

yDeltaList = []
yTheta1List = []
yTheta2List = []
ySingularList = []

for setItr in range(0,len(X)):
	
	yDelta = numpy.zeros(xDim.shape[0],dtype=numpy.float64)
	yTheta1 = numpy.zeros(xDim.shape[0],dtype=numpy.float64)
	yTheta2 = numpy.zeros(xDim.shape[0],dtype=numpy.float64)
	ySingular = numpy.zeros(xDim.shape[0],dtype=numpy.float64)
	
	if setItr == 0:
	 print '\n\n Classical Projection with size of measurment matrix = '+str(X[setItr].shape)
	else:
         print '\n\n POD Approach with wset'+str(setItr)+' size of measurment matrix = '+str(X[setItr].shape)
	
	print '\nProjecting the QEVP having dimension '+str(M_orig.shape)+' onto a smaller subspace of dimension'
	
	for i1 in range(0,len(xDim)):
		
		
		obj.projectionDimension = xDim[i1]
		
		#obtaining the projection matrix
		Q, singularValues = projection.obtain_projection_matrix(obj,X[setItr])
	        
		#Projection
		QT = Q.T.conjugate()
		M =  QT.dot(M_orig.dot(Q))
		C =  QT.dot(C_orig.dot(Q))
		K =  QT.dot(K_orig.dot(Q))
		
		n = M.shape[0]
		no_of_evs = 2*n-2;
		la_pod, evec_pod = solver.qev_dense(obj,M,C,K,no_of_evs)
		
		evec_pod = unlinearize.unlinearize_matrices(evec_pod)
		
		#print '\n\n\nCalculating Delta'
		delta = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)
		Theta1 = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)
		Theta2 = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)
		for itr in range(0,len(laTarget)):
		    		    
		    #computing cosine of the angle between evecTarget[:,itr] and the subspace spanned by columns of Q
		    #QR factorization
		    qFactor, rFactor = numpy.linalg.qr(Q)	    
		    cosineAngle1 = numpy.zeros(qFactor.shape[1],dtype=numpy.float64)
		    for jtr in range(0,len(cosineAngle1)):
	    		#v1 = qFactor[:,jtr]
			v2 = evecTarget[:,itr]
			svd = norm(numpy.dot(numpy.transpose(v2), qFactor))
			cosineAngle1[jtr] = numpy.arccos(min(svd,1.0))
			
			#cosineAngle1[jtr] = numpy.absolute(numpy.dot(v1,v2))/(norm(v1)*norm(v2))
			#cosineAngle1[jtr] = 1.0 if cosineAngle1[jtr] > 1.0 else cosineAngle1[jtr]
			#svd = LA.svd(numpy.dot(numpy.transpose(v2), qFactor))		    
			#cosineAngle1[jtr] = numpy.arccos(min(svd[1].min(), 1.0))
			
							    
		    cosineAngle2 = numpy.zeros(evec_pod.shape[1],dtype=numpy.float64)
		    for jtr in range(0,len(cosineAngle2)):
	    		v1 = Q.dot(evec_pod[:,jtr])
			v2 = evecTarget[:,itr]
			cosineAngle2[jtr] = numpy.absolute(numpy.dot(v1,v2))/(norm(v1)*norm(v2))
			cosineAngle2[jtr] = numpy.arccos(min(cosineAngle2[jtr], 1.0))
			
		    delta[itr] = numpy.min(numpy.absolute(la_pod-laTarget[itr]))
		    Theta1[itr] = numpy.min(cosineAngle1)
		    Theta2[itr] = numpy.min(cosineAngle2)
		    
		
		yDelta[i1] = numpy.max(delta)/numpy.max(numpy.absolute(laTarget))
		yTheta1[i1] = max(numpy.max(Theta1),math.pow(0.1,30)) #setting 10^-30 if 0 for log plot
		yTheta2[i1] = max(numpy.max(Theta2),math.pow(0.1,30)) #setting 10^-30 if 0 for log plot
		ySingular[i1] = singularValues[obj.projectionDimension-1] #min(singularValues)
		
		
        yDeltaList.append(yDelta)
	yTheta1List.append(yTheta1)
	yTheta2List.append(yTheta2)
	ySingularList.append(ySingular)

print '\n\n The delta values for the three sets are as follows \n\n'
for setItr in range(0,len(X)):
	print yTheta1List[setItr]

fig = plt.figure('logarithmic decay of different parameters with dimension')

plt.subplot(221)
plt.plot(xDim, yDeltaList[0], '-ro')
plt.plot(xDim, yDeltaList[1], '-kD')
plt.plot(xDim, yDeltaList[2], '-bs')
plt.plot(xDim, yDeltaList[3], '-go')
plt.yscale('log')
plt.xlabel('dimension')
plt.ylabel('delta')
#plt.title('relative error between computed and exact eigenvalues')
plt.grid(True)

plt.subplot(222)
plt.plot(xDim, yTheta1List[0], '-ro')
plt.plot(xDim, yTheta1List[1], '-kD')
plt.plot(xDim, yTheta1List[2], '-bs')
plt.plot(xDim, yTheta1List[3], '-go')
plt.yscale('log')
plt.xlabel('dimension')
plt.ylabel('theta1')
#plt.title('maximum angle between exact eigenvector and the projection space')
plt.grid(True)

plt.subplot(223)
plt.plot(xDim, yTheta2List[0], '-ro')
plt.plot(xDim, yTheta2List[1], '-kD')
plt.plot(xDim, yTheta2List[2], '-bs')
plt.plot(xDim, yTheta2List[3], '-go')
plt.yscale('log')
plt.xlabel('dimension')
plt.ylabel('theta2')
#plt.title('maximum angle between exact eigenvector and the computed eigenvector')
plt.grid(True)

plt.subplot(224)
#plt.plot(xDim, ySingularList[0], '-ro')
plt.plot(xDim, ySingularList[1], '-kD')
plt.plot(xDim, ySingularList[2], '-bs')
plt.plot(xDim, ySingularList[3], '-go')
plt.yscale('log')
plt.xlabel('dimension')
plt.ylabel('singular values')
plt.title('decay of singular values')
plt.grid(True)

plt.subplot_tool()
plt.show()
brake.save(obj.output_path+'errorDecay', ext="png", close=True, verbose=False)