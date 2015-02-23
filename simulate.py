'''
Simulates the transition of eigenvalues from the left half plane into the
target rectangular region. The plots are obtained in plotEigsCover.png 
and plotEigsTransition.png which gets updated on fly thus showing the transition.
'''

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
from brake.solve import projection, solver
from brake.analyze import residual, visual

begin_program = timeit.default_timer()


#----------------------------Create Object with no logging ---------------------------------------
obj = createBrakeClassObject.returnObject(0)

#Log all the BrakeClass attribute values in the info log file
if(obj.log_level):
  obj.displayParametersLog()
  

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

begin = timeit.default_timer()
X_wset = projection.obtain_measurment_matrix(obj)
Q, singularValues = projection.obtain_projection_matrix(obj,X_wset)
end = timeit.default_timer()
#------------------------------------------------done creating the projection matrix


'''
Using the projection matrix(Q) constructed above for the base frequencies, we can now obtain results 
for several intermediate frequencies by projecting each of them onto a smaller subspace 
'''
#-------------------------------------------Simulation-------------------------------------------
print "\n","\n","\n",'Beginning Simulation Phase'

if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Simulation Phase')
  obj.logger_i.info('-------------------------------------------------------------------------')
  obj.logger_i.info('Performing Proper Orthogonal Decomposition for the following frequencies')
	
sparse_list = load.load_matrices(obj)

obj.logger_t.info("\n"+"\n"+'Logging Simulation Times')
for i in range(0,len(obj.omega_range)):
        begin_sim = timeit.default_timer()
	print "\n"
	print '---------------------------------- QEVP '+str(i+1),' --------------------------'
	omega = obj.omega_range[i]
	print 'omega = ',str(omega/(2*math.pi))
	if(obj.log_level):
          obj.logger_i.info("\n"+"\n"+'POD for omega = '+str(omega/(2*math.pi)))
	M, C, K = assemble.create_MCK(obj, sparse_list, omega)
	print 'Projecting the QEVP having dimension '+str(M.shape)
	#Projection
	QT = Q.T.conjugate()
	M =  QT.dot(M.dot(Q))
	C =  QT.dot(C.dot(Q))
	K =  QT.dot(K.dot(Q))
	print 'Onto a smaller subspace of dimension '+str(M.shape)
	
	n = M.shape[0]
	no_of_evs = n-5;
	begin_solver = timeit.default_timer()	
        la, evec = solver.qev_dense(obj,M,C,K,no_of_evs);
        end_solver = timeit.default_timer()
	res_qevp = residual.residual_qevp(M,C,K,la,evec[0:n,:])
	print 'Maximum Residual error for the QEVP is ',max(res_qevp)
	
        brake.printEigs(obj,la,'target','terminal')
	radius = visual.plot_eigs_cover(obj,la)
        radius = visual.plot_eigs_transition(obj,la)

        end_sim = timeit.default_timer() 
        obj.logger_t.info('Simulation time for the '+str(i+1)+'th simulation: '+"%.2f" % \
        (end_sim-begin_sim)+' sec')
        
end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'
