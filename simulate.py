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
from brake.solve import projection

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