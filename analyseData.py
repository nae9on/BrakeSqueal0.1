 
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


log_level = 10
dt = datetime.date.today()
info_log_file = output_path+'dataAnalysis/dataProperties'+'.log' 
time_log_file = output_path+'dataAnalysis/dataAnalysisTime'+'.log' 

obj = brake.BrakeClass(input_path, output_path, info_log_file, time_log_file, log_level)

setattr(obj, 'data_file_list', ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL'])
setattr(obj, 'data_file_name', ['M1','D1','DG','DR','D4','K1','KR','KGeo'])


if(obj.log_level):
  obj.logger_i.info("\n"+"\n"+'Beginning Data Analysis')
  obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Data Analysis'

sparse_list = load.load_matrices(obj)

obj.logger_i.debug("\n"+'Matrices in CSC format converted to CSR')

obj.logger_i.debug("\n\n"+'Properties of various component matrices'+"\n")

for i in range(0,len(sparse_list)):
	componentMatrix = sparse_list[i]
	normMatrix = norm.onenormest(componentMatrix.tocsr(), t=3, itmax=5, compute_v=False, compute_w=False)
	print obj.data_file_list[i], obj.data_file_name[i], ' NonZeros = ', componentMatrix.nnz, ' 1-Norm = ', normMatrix
	obj.logger_i.debug(obj.data_file_list[i]+' '+obj.data_file_name[i]+' Nonzeros = '+str(componentMatrix.nnz)+' 1-Norm = '+str(normMatrix))
	
	#saving the sparsity pattern
	fig = plt.gcf()
	fig.clf()
	fig.gca().add_artist(plt.spy(componentMatrix))
        fig.savefig(obj.output_path+'dataAnalysis/'+obj.data_file_name[i]+'.png')


end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'

if(obj.log_level):
    obj.logger_i.info("\n"+"\n"+'Finished Data Analysis')

    #----------------Logging Time Complexity-----
    obj.logger_t.info('Time Taken for data analysis: '+"%.2f" % (end_program-begin_program)+' sec')

