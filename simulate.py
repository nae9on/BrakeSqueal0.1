# Driver file for the brakesqueal project

#----------------------------------Standard Library Imports---------------------------------------
import math
import numpy
import timeit
import socket
import datetime
 
#----------------------------Application Specific Imports-----------------------------------------
import brake
from brake.initialize import load, assemble, scale, diagscale
from brake.solve import projection, solver
from brake.analyze import residual, visual

begin_program = timeit.default_timer()

#Initializing Parameters
#################################################        

#Set Input/Output Path
#Get the host system name to define the path
host=socket.gethostname()
print "Working on :", host, "computer"

if host == 'laplace':
   input_path = '/homes/numerik/quraishi/work/paperprojects/brake_squeel_project/ProblemData/\
   Matrix/5kOmegaRef%(row)d/' % {'row': omegaRef}
elif host == 'aif-server':
     input_path = 'd:/eigenwerte/ProblemData/Matrix/800kOmegaRef%(row)d/' % {'row': omegaRef}
elif host == 'frenet':
     input_path = '/homes/extern/kadar/Desktop/project/python_source/data/5koref1/'
     output_path = '/homes/extern/kadar/Desktop/BrakeSqueal0.1/output/'
elif host == 'ali-Inspiron-1525':
     input_path = './data/5koref1/'
     output_path = './output/'  
else:
     input_path = './data/5koref1/'
     output_path = './output/'

#Logging Parameters
'''
'notset':logging.NOTSET, #0 (for no logging - fastest)
'debug': logging.DEBUG, #10 (to capture detailed debug information - slowest)
'info': logging.INFO, #20 (to capture essential information)
'''
log_level = 20
dt = datetime.date.today()
info_log_file = output_path+'info_'+dt.strftime("%d%b")+'.log' #example info_02Aug.log
time_log_file = output_path+'time_'+dt.strftime("%d%b")+'.log' #example time_02Aug.log


'''
This would create the first object of BrakeClass with certain limited attributes. 
More attributes can be added anywhere in the program using the setattr command.
These attributes can then be referred later in a module.
'''
obj = brake.BrakeClass(input_path, output_path, info_log_file, time_log_file, log_level)


'''
Adding more attributes to the BrakeClass object 'obj' created above.

enable - flag 0/1 for comparing POD results with actual results, very 'time intensive.

Problem parameters described below
M = m
C = c1+c2*(omega/omegaRef)+c3*(omegaRef/omega)
K = k1+k2+k3*math.pow((omega/omegaRef),2)
[m c1 c2 c3 c4 k1 k2 k3]
['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL']
[M D1 DG DR ~ K1 KR KGEO]

omegaRef - reference omega.
fRef - reference frequency.
omega_basis - base angular velocity's to be used for creating the projection matrix.
omega_range - range of angular velocity's for simulation.
target - target rectangular region.


cutoff - cutoff for truncating the less sgnificant singular values
as comapred to the most significant s[0]
ex s[i] < s[0]*cutoff

evs_per_shift - number of eigenvalues to be calculated per shift
desired_area_fraction - area fraction of the target region to be covered

'''

setattr(obj, 'enable', 0) 
setattr(obj, 'data_file_list', ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL'])
setattr(obj, 'omegaRef', 1)
setattr(obj, 'fRef', 1600)
setattr(obj, 'omega_basis', numpy.array([1,20])*2*math.pi)  
setattr(obj, 'omega_range', numpy.linspace(1, 20, num=2)*2*math.pi)
setattr(obj, 'target', numpy.array([-10,1000,-50,12000]))
setattr(obj, 'cutoff', 0.0001) 
setattr(obj, 'evs_per_shift', 30)
setattr(obj, 'desired_area_fraction', 0.99)




#------------------------------------------------done initialization

#Log all the BrakeClass attribute values in the info log file
if(obj.log_level):
  obj.displayParametersLog()

#Begin creating the projection matrix
#################################################
if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Setup Phase')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix'

begin = timeit.default_timer()
Q = projection.obtain_projection_matrix(obj)
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
	
        brake.print_eigs(obj,la,'target','terminal')
	radius = visual.plot_eigs_cover(obj,la)
        radius = visual.plot_eigs_transition(obj,la)

	#set enable = 1 for verification
	if(obj.enable):
	  obj.logger_i.info("\n"+'The eigenvalues(in the target region) approximated using \
										POD are :'+"\n") 
	  brake.print_eigs(obj,la,'target','terminal')
	  la, evec = brake_squeal_qevp.BrakeSquealQevp(i,path,data_file_list,omega,
								megaRef, fRef,target, evs_per_shift)
	  obj.logger_i.info("\n"+'The eigenvalues(in the target region) obtained after solving \
								the original QEVP are :'+"\n") 
	  brake.print_eigs(obj,la,'target','terminal')
	  obj.logger_i.info("\n"+'----------------------------------------------------------'+"\n")
        
        end_sim = timeit.default_timer() 
        obj.logger_t.info('Simulation time for the '+str(i+1)+'th simulation: '+"%.2f" % \
        (end_sim-begin_sim)+' sec')
        
end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'
