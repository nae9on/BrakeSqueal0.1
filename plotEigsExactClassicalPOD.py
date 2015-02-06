''' 
1. Compares the exact eigenvalues in the target region with the eigenvalues 
obtained from the Tradional and POD Approach.
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

#----------------------------Create Object with no logging ---------------------------------------
logLevel = 0
obj = createBrakeClassObject.returnObject(logLevel)

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
#Extracting the eigenvalues and the corresponding eigenvectors that lie in the target region
laTarget, evecTarget = brake.extractEigs(obj, laExact, evecExact, 'target')
deltaTarget = numpy.zeros(laTarget.shape[0],dtype=numpy.float64)


#Obtaining the original M,C,K matrices for the test omega
sparse_list = load.load_matrices(obj)
M_orig, C_orig, K_orig = assemble.create_MCK(obj, sparse_list, obj.omegaTest)

####################################################################
####################################################################
####################### Classical Projection #######################
####################################################################
####################################################################
print '\n\nObtaining the classical projection matrix'
la, evec = classicalProjection.Obtain_eigs(obj,150)
XClassical = numpy.concatenate((evec.real,evec.imag), axis=1)
print '\nClassical Projection with size of measurment matrix = '+str(XClassical.shape)
#Setting Dimension of the projection matrix
obj.projectionDimension = 300
Q, singularValues = projection.obtain_projection_matrix(obj,XClassical)
QT = Q.T.conjugate()
M =  QT.dot(M_orig.dot(Q))
C =  QT.dot(C_orig.dot(Q))
K =  QT.dot(K_orig.dot(Q))
print 'Projected the QEVP having dimension '+str(M_orig.shape)+' onto a smaller subspace of dimension '+str(M.shape[1])
n = M.shape[0]
no_of_evs = 2*n-2;
laClassical, evecClassical = solver.qev_dense(obj,M,C,K,no_of_evs)
evecClassical = unlinearize.unlinearize_matrices(evecClassical)

#print '\n\n\nCalculating Delta'
deltaClassical = numpy.zeros(laClassical.shape[0],dtype=numpy.float64)
for itr in range(0,len(laClassical)):
	deltaClassical[itr] = numpy.min(numpy.absolute(laClassical[itr]-laTarget))/numpy.max(numpy.absolute(laTarget))

####################################################################
####################################################################
####################### POD Approach ###############################
####################################################################
####################################################################
print '\n\nObtaining the POD projection matrix'
print '\nCreating the Measurment Matrix for wset3'
obj.evs_per_shift = 17
obj.omega_basis = numpy.array([1,2.5,5,7.5,10,12.5,15,17.5,20])*2*math.pi
X_wset3 = projection.obtain_measurment_matrix(obj)
print '\nPOD Projection with size of measurment matrix = '+str(X_wset3.shape)
#Setting Dimension of the projection matrix
obj.projectionDimension = 300
Q, singularValues = projection.obtain_projection_matrix(obj,X_wset3)
QT = Q.T.conjugate()
M =  QT.dot(M_orig.dot(Q))
C =  QT.dot(C_orig.dot(Q))
K =  QT.dot(K_orig.dot(Q))
print 'Projected the QEVP having dimension '+str(M_orig.shape)+' onto a smaller subspace of dimension '+str(M.shape[1])
n = M.shape[0]
no_of_evs = 2*n-2;
laPOD, evecPOD = solver.qev_dense(obj,M,C,K,no_of_evs)
evecPOD = unlinearize.unlinearize_matrices(evecPOD)

#print '\n\n\nCalculating Delta'
deltaPOD = numpy.zeros(laPOD.shape[0],dtype=numpy.float64)
for itr in range(0,len(laPOD)):
	deltaPOD[itr] = numpy.min(numpy.absolute(laPOD[itr]-laTarget))/numpy.max(numpy.absolute(laTarget))
	
print numpy.sort(deltaClassical)
print numpy.sort(deltaPOD)
####################################################################
####################################################################
'''
x1 = obj.target[0]
x2 = obj.target[1]
y1 = obj.target[2]
y2 = obj.target[3]
# calculate center of the obj.target rectangular region
tau = complex((x1+x2)/2,(y1+y2)/2)
fig = plt.figure(figsize=(16.0, 10.0))
react = plt.Rectangle((x1,y1),abs(x2-x1),abs(y2-y1),color='black',fill=False)


plt.plot(laTarget.real, laTarget.imag, 'ko', label='exact')
plt.plot(laClassical.real, laClassical.imag, 'r+', markersize=10, fillstyle='none', label='traditional')
plt.plot(laPOD.real, laPOD.imag, 'go', markersize=10, fillstyle='none', label='POD')
plt.legend(loc='upper right', shadow=True)
brake.save(obj.output_path+'eigsExactClassicalPOD'+str(obj.omegaTest/(2*math.pi)), ext="png", close=False, verbose=True)
plt.show()
'''

fig = plt.figure()
ax = plt.gca()
minLim = min(deltaPOD)
maxLim = pow(10,-5)
hdl = plt.scatter(laTarget.real, laTarget.imag, s=10, c=deltaTarget, marker='o', vmin=minLim, vmax=maxLim, label='exact')
hdl = plt.scatter(laClassical.real, laClassical.imag, s=100, c=deltaClassical, marker='+', facecolors='none', linewidth='3', vmin=minLim, vmax=maxLim, label='traditional')
hdl = plt.scatter(laPOD.real, laPOD.imag, s=100, c=deltaPOD, marker='o', facecolors='none', linewidth='3', vmin=minLim, vmax=maxLim, label='POD')
plt.legend(loc='upper left', shadow=True)
xWidth = max(laTarget.real)-min(laTarget.real)
yWidth = max(laTarget.imag)-min(laTarget.imag)
ax.set_xlim((min(laTarget.real)-0.1*xWidth,max(laTarget.real)+0.1*xWidth))
ax.set_ylim((min(laTarget.imag)-0.1*yWidth,max(laTarget.imag)+0.1*yWidth))
plt.axhline(0, color='blue')
plt.axvline(0, color='blue')
clevs = [-15, -10, -5]
cb1 = plt.colorbar(hdl, orientation='vertical')#, ticks=clevs)
cb1.ax.set_yticklabels(np.arange(-15,-4))# vertically oriented colorbar
brake.save(obj.output_path+'eigsExactClassicalPOD'+str(obj.omegaTest/(2*math.pi)), ext="png", close=False, verbose=True)
plt.show()
####################################################################
####################################################################
