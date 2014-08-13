'''
Driver file for the brakesqueal project

__author__ = "Ali Kadar (kadar@math.tu-berlin.de)"
__version__ = ""
__date__ = ""
__copyright__ = ""
__license__ = "Python"

Input
1. path (where the data files are located)
2. data_file_list (data file names for the component matrix m,c1,c2,c3,c4,k1,k2,k3)
3. omega_basis (array of base frequencies for creating the projection matrix, 
								ex array([1,10,20])*2*math.pi)
4. omega_range (range of frequencies for which the simulation is to be done, 
							ex np.linspace(1, 20, num=1)*2*math.pi)
4. omegaRef (reference omega value, ex 1)
5. fRef (reference frequency value, ex 1600)
6. target (target region for the shift points, ex array([-10,1000,-50,12000]))
7. cutoff (limit to take into account only the significant singular values, ex 0.0001)
8. enable (flag which should be set true for comparing results, 'time intensive', ex 0/1)
9. evs_per_shift (number of eigenvalues per shift, ex 10)

Output

'''





#----------------------------------------------System Imports-------------------------------------- 
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import sys
import math
import time
import timeit
import socket
import logging
import numpy as np
from numpy import array
from datetime import date
from scipy.sparse import vstack, hstack 





#----------------------------------------User Defined Imports--------------------------------------- 
# Please ensure that the following files are placed in the working directory
import solver
import py_plot
import assemble
import def_list
import projection
import load_matrices
import residual_check
import brake_squeal_qevp


begin_program = timeit.default_timer()


#-----------------------------------------------Logging---------------------------------------------
# will create a log file, ex record09Feb.log
d = date.today()
INFO_LOG_FILENAME = 'info_'+d.strftime("%d%b")+'.log'
TIME_LOG_FILENAME = 'time_'+d.strftime("%d%b")+'.log'

#define logging levels
LEVELS = {'notset':logging.NOTSET, #0 --> numerical value
	  'debug': logging.DEBUG, #10
          'info': logging.INFO, #20
          'warning': logging.WARNING, #30
          'error': logging.ERROR, #40
          'critical': logging.CRITICAL} #50


'''
sets logging parameter passed by user(as sys argument)
'DEBUG' to capture detailed debug information
'INFO' to capture essential information
'NOTSET'(=0) if no parameter is passed by user
'''

if len(sys.argv) > 1:
    level_name = sys.argv[1]
    level = LEVELS.get(level_name, logging.NOTSET)
else:
    level=logging.NOTSET 


formatter = logging.Formatter('%(message)s')

# creating info logger
logger_i = logging.getLogger('info_logger')
hdlr_i = logging.FileHandler(INFO_LOG_FILENAME,mode='w')    
hdlr_i.setFormatter(formatter)
logger_i.addHandler(hdlr_i)
logger_i.setLevel(level)

# creating time logger
logger_t = logging.getLogger('time_logger')
hdlr_t = logging.FileHandler(TIME_LOG_FILENAME,mode='w')
hdlr_t.setFormatter(formatter)
logger_t.addHandler(hdlr_t)
logger_t.setLevel(level)

LOG_LEVEL = level

if(LOG_LEVEL == 0):
	logger_i.setLevel(logging.INFO)
	logger_i.info('The current level of logging is '+str(LOG_LEVEL))
	logger_i.info('To enable detailed logging provide one of the following parameters:')
	logger_i.info('info,debug')
	logger_i.info('The Info Data is logged in the file : '+INFO_LOG_FILENAME)
	logger_i.info('The Time Data is logged in the file : '+TIME_LOG_FILENAME)	
	logger_i.info('ex execute following in shell: python main.py info')
	logger_i.setLevel(logging.NOTSET)

if(LOG_LEVEL):
	logger_i.info('Auto logged Info for the Quadratic Eigen Value Problem(Brake Squeal Project)')
	logger_i.info("\n"+'In Driver file')





# ------------------------------------------Initialization------------------------------------------ 
# Specify the system path where the data files are located
host=socket.gethostname()
print "Working on :", host, "computer"
if host == 'laplace':
    path = '/homes/numerik/quraishi/work/paperprojects/brake_squeel_project/\
					ProblemData/Matrix/5kOmegaRef%(row)d/' % {'row': omegaRef}
elif host == 'aif-server':
    path = 'd:/eigenwerte/ProblemData/Matrix/800kOmegaRef%(row)d/' % {'row': omegaRef}
elif host == 'frenet':
    path = '/homes/extern/kadar/Desktop/project/python_source/data/5koref1/'
elif host == 'ubuntu':
    path = '/home/ali/Desktop/project/python_source/data/5koref1/'	
else:
    path = '/homes/extern/kadar/Desktop/project/python_source/data/5koref1/'

#path = '/homes/extern/kadar/store/800koref5/'

data_file_list = ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL']

# array of base frequencies 
omega_basis = array([1,20])*2*math.pi
omega_range = np.linspace(1, 20, num=2)*2*math.pi
omegaRef = 1
fRef = 1600
target = array([-10,1000,-50,12000])
cutoff = 0.0001; # s[i] < s[0]*cutoff
enable = 0
evs_per_shift = 30

# calculate center of the target rectangular region
tau = complex((target[0]+target[1])/2,(target[2]+target[3])/2)

if(LOG_LEVEL):
    logger_i.info("\n"+'Parameters: '+"\n"+'omega = '+str(omegaRef)+"\n"+'fRef = '+str(fRef))
    logger_i.info('Target rectangle = '+str(target)+"\n"+'tau = '+str(tau))
    logger_i.info("\n"+'Reading data from '+path)





# ---------------------------------------------Projection------------------------------------------
if(LOG_LEVEL):
	logger_i.info("\n"+"\n"+'Beginning Setup Phase')
        logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix'


begin = timeit.default_timer()
Q = projection.obtain_projection_matrix(path, data_file_list, omega_basis, omegaRef, fRef, target, 
								     cutoff, evs_per_shift, level)
end = timeit.default_timer()


if(LOG_LEVEL):
	logger_i.info('---------------------------------------------------------------'+"\n"+"\n")
	logger_t.info("\n"+"\n"+"\n"+'\t\t\tProjection Matrix Created in time: '\
								+"%.2f" % (end-begin)+' sec')
	logger_i.info("\n"+'Projection matrix created for the following base frequencies '\
								+str(omega_basis/(2*math.pi))+"\n")
	
'''
Using the projection matrix(Q) constructed above for the base frequencies, we can now obtain results 
for several intermediate frequencies by projecting each of them onto a smaller subspace 
'''
#-------------------------------------------Simulation----------------------------------------------
print "\n","\n","\n",'Beginning Simulation Phase'

if(LOG_LEVEL):
	logger_i.info("\n"+"\n"+'Beginning Simulation Phase')
	logger_i.info('-------------------------------------------------------------------------')
	logger_i.info('Performing Proper Orthogonal Decomposition for the following frequencies')
	
sparse_list = load_matrices.load(path,data_file_list, level)

logger_t.info("\n"+"\n"+'Logging Simulation Times')
for i in range(0,len(omega_range)):
        begin_sim = timeit.default_timer()
	print "\n"
	print '---------------------------------- QEVP '+str(i+1),' --------------------------'
	omega = omega_range[i]
	print 'omega = ',str(omega/(2*math.pi))
	if(LOG_LEVEL):
          logger_i.info("\n"+"\n"+'POD for omega = '+str(omega/(2*math.pi)))
	M, C, K = assemble.create_MCK(sparse_list,omega,omegaRef,fRef, level)
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
        la, evec = solver.qev_optimized_dense(M,C,K,no_of_evs);
        end_solver = timeit.default_timer()
	res_qevp = residual_check.residual_qevp_nonsparse(M,C,K,la,evec[0:n,:])
	print 'Maximum Residual error for the QEVP is ',max(res_qevp)
	
        def_list.my_print_few_infile(la,0)
	radius = py_plot.plot_eigs_cover(la,tau,target)
        radius = py_plot.plot_eigs_transition(la,tau,target)

	#set enable = 1 for verification
	if(enable):
		logger_i.info("\n"+'The eigenvalues(in the target region) approximated using \
										POD are :'+"\n") 
		def_list.my_print_few_infile(la,1)
		la, evec = brake_squeal_qevp.BrakeSquealQevp(i,path,data_file_list,omega,
								megaRef, fRef,target, evs_per_shift)
		logger_i.info("\n"+'The eigenvalues(in the target region) obtained after solving \
								the original QEVP are :'+"\n") 
		def_list.my_print_few_infile(la,1)
	 	logger_i.info("\n"+'----------------------------------------------------------'+"\n")
        end_sim = timeit.default_timer() 
        logger_t.info('Simulation time for the '+str(i+1)+'th simulation: '+"%.2f" % (end_sim-begin_sim)+' sec')
end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'
