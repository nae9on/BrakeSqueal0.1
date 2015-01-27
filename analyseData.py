 
# This file reads the component matrices and displays there sparsity pattern
# and other properties

#----------------------------------Standard Library Imports---------------------------------------
import math
import numpy
import timeit
import socket
import datetime
import matplotlib.pylab as plt
import scipy.sparse.linalg as norm
 
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
log_level = 10
dt = datetime.date.today()
info_log_file = output_path+'dataProperties'+'.log' #example info_02Aug.log
time_log_file = output_path+'dataAnalysisTime'+'.log' #example time_02Aug.log


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
setattr(obj, 'data_file_name', ['M1','D1','DG','DR','D4','K1','KR','KGeo'])
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
  obj.logger_i.info("\n"+"\n"+'Beginning Data Analysis')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Data Analysis'

sparse_list = load.load_matrices(obj)

obj.logger_i.debug("\n"+'Properties of various component matrices')

obj.logger_i.debug('Matrices in CSC format converted to CSR')

fig = plt.gcf()

for i in range(0,len(sparse_list)):
	componentMatrix = sparse_list[i]
	normMatrix = norm.onenormest(componentMatrix.tocsr(), t=3, itmax=5, compute_v=False, compute_w=False)
	print obj.data_file_list[i], obj.data_file_name[i], ' NonZeros = ', componentMatrix.nnz, ' 1-Norm = ', normMatrix
	obj.logger_i.debug(obj.data_file_list[i]+' '+obj.data_file_name[i]+' Nonzeros = '+str(componentMatrix.nnz)+' 1-Norm = '+str(normMatrix))
	fig.clf()
	fig.gca().add_artist(plt.spy(componentMatrix))
        fig.savefig(obj.output_path+'sparsity_'+obj.data_file_name[i]+'.png')


end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'


