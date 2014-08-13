#Driver file for the brakesqueal project

#Standard Library Imports
import math
import timeit
import numpy
from numpy import array
from datetime import date
from socket import gethostname

#Application Specific Imports
import brake
from brake.initialize import load, assemble, scale, diagscale
from brake.solve import projection

begin_program = timeit.default_timer()

#Initializing Parameters
#---------------------------------------------------------------------------------------------------        

# Specify the system path where the data files are located
host=gethostname()
print "Working on :", host, "computer"

if host == 'laplace':
   path = '/homes/numerik/quraishi/work/paperprojects/brake_squeel_project/ProblemData/\
   Matrix/5kOmegaRef%(row)d/' % {'row': omegaRef}
elif host == 'aif-server':
     path = 'd:/eigenwerte/ProblemData/Matrix/800kOmegaRef%(row)d/' % {'row': omegaRef}
elif host == 'frenet':
     path = '/homes/extern/kadar/Desktop/project/python_source/data/5koref1/'
elif host == 'ubuntu':
     path = '/home/ali/Desktop/project/python_source/data/5koref1/'	
else:
     path = '/homes/extern/kadar/Desktop/project/python_source/data/5koref1/'

#Input/Output Path
input_path = path
output_path = '/home/ali/Desktop/project/python_source/brake_package/output/'

#Logging Parameters
log_level = 0
dt = date.today()
info_log_file = output_path+'info_'+dt.strftime("%d%b")+'.log'
time_log_file = output_path+'time_'+dt.strftime("%d%b")+'.log'


'''
This would create the first object of BrakeClass class
with certain limited attributes. More attributes can
be added using the setattr command
'''
obj = brake.BrakeClass(input_path, output_path, info_log_file, time_log_file, log_level)

#flags
setattr(obj, 'enable', 0)

#Problem Parameters
setattr(obj, 'fRef', 1600)
setattr(obj, 'omegaRef', 1)
setattr(obj, 'omega_basis', array([1,20])*2*math.pi)
setattr(obj, 'omega_range', numpy.linspace(1, 20, num=2*2*math.pi))
setattr(obj, 'target', array([-10,1000,-50,12000]))
setattr(obj, 'tau', complex((obj.target[0]+obj.target[1])/2,(obj.target[2]+obj.target[3])/2))
setattr(obj, 'cutoff', 0.0001) # s[i] < s[0]*cutoff
setattr(obj, 'evs_per_shift', 30)
setattr(obj, 'data_file_list', ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL'])


if(obj.log_level):
  obj.logger_i.info("\n"+'Parameters: '+"\n"+'omega = '+str(obj.omegaRef)+"\n"+'fRef = '+str(obj.fRef))      
  obj.logger_i.info('Target rectangle = '+str(obj.target)+"\n"+'tau = '+str(obj.tau))
  obj.logger_i.info("\n"+'Reading data from '+obj.input_path)


# ---------------------------------------------Projection------------------------------------------
if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Setup Phase')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrix'


begin = timeit.default_timer()
Q = projection.obtain_projection_matrix(obj)
end = timeit.default_timer()

'''
if(obj.log_level):
  obj.logger_i.info('---------------------------------------------------------------'+"\n"+"\n")
  obj.logger_t.info("\n"+"\n"+"\n"+'\t\t\tProjection Matrix Created in time: '\
								+"%.2f" % (end-begin)+' sec')
  obj.logger_i.info("\n"+'Projection matrix created for the following base frequencies '\
								+str(omega_basis/(2*math.pi))+"\n")
								
'''
'''
obj.displayParameters()

print "Total object count %d" % brake.BrakeClass.objCount

sparse_list = load.load_matrices(obj)

M, C, K = assemble.create_MCK(obj, sparse_list, obj.omega_basis[1])

M, C, K, tau, gamma, delta = scale.scale_matrices(obj, M , C, K, tau)

M, C, K, DL, DR = diagscale.diag_scale_matrices(obj,M,C,K)
'''

