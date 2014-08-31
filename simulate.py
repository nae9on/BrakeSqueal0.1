#Driver file for the brakesqueal project

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
from brake.analyze import residual, py_plot

begin_program = timeit.default_timer()

#Initializing Parameters
#################################################        

#Set Input Path
#Get the host system name to define the input path
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
elif host == 'ubuntu':
     input_path = '/home/ali/Desktop/project/python_source/data/5koref1/'
     output_path = '/home/ali/Desktop/BrakeSqueal0.1/output/'
else:
     input_path = '/homes/extern/kadar/Desktop/project/python_source/data/5koref1/'
     output_path = '/homes/extern/kadar/Desktop/BrakeSqueal0.1/output/'

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
More attributes can be added using the setattr command.
These attributes can then be referred later in a module
'''
obj = brake.BrakeClass(input_path, output_path, info_log_file, time_log_file, log_level)

#Add more attributes to the BrakeClass object 'obj'

#flags
setattr(obj, 'enable', 0) #flag 0/1 for comparing results, very 'time intensive'

#Problem Parameters
'''
M = m
C = c1+c2*(omega/omegaRef)+c3*(omegaRef/omega)
K = k1+k2+k3*math.pow((omega/omegaRef),2)
[m c1 c2 c3 c4 k1 k2 k3]
['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL']
[M D1 DG DR ~ K1 KR KGEO]
'''

setattr(obj, 'data_file_list', ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL'])
setattr(obj, 'omegaRef', 1) #Reference omega
setattr(obj, 'fRef', 1600) #Reference omega


#set base angular velocity's for creating the projection matrix
setattr(obj, 'omega_basis', numpy.array([1,20])*2*math.pi)  

#set range of angular velocity's for simulation
setattr(obj, 'omega_range', numpy.linspace(1, 20, num=2*2*math.pi))

setattr(obj, 'target', numpy.array([-10,1000,-50,12000])) #target rectangular region

'''
cutoff for truncating the less sgnificant singular values
as comapred to the most significant s[0]
s[i] < s[0]*cutoff
'''
setattr(obj, 'cutoff', 0.0001) 

#number of eigenvalues to be calculated per shift
setattr(obj, 'evs_per_shift', 30)

#The area fraction of the target region to be covered
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
	
        brake.print_target_eigs(obj,la,0)
	radius = py_plot.plot_eigs_cover(obj,la)
        radius = py_plot.plot_eigs_transition(obj,la)

	#set enable = 1 for verification
	if(obj.enable):
	  obj.logger_i.info("\n"+'The eigenvalues(in the target region) approximated using \
										POD are :'+"\n") 
	  brake.print_target_eigs(obj,la,1)
	  la, evec = brake_squeal_qevp.BrakeSquealQevp(i,path,data_file_list,omega,
								megaRef, fRef,target, evs_per_shift)
	  obj.logger_i.info("\n"+'The eigenvalues(in the target region) obtained after solving \
								the original QEVP are :'+"\n") 
	  brake.print_target_eigs(obj,la,1)
	  obj.logger_i.info("\n"+'----------------------------------------------------------'+"\n")
        
        end_sim = timeit.default_timer() 
        obj.logger_t.info('Simulation time for the '+str(i+1)+'th simulation: '+"%.2f" % \
        (end_sim-begin_sim)+' sec')
        
end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'
