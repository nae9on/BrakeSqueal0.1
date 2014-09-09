#Experiment demonstrating that the POD approach is better than 
#the classcial model reduction approach

#----------------------------------Standard Library Imports---------------------------------------
import math
import numpy
import timeit
import socket
import datetime
 
#----------------------------Application Specific Imports-----------------------------------------
import brake
from brake.initialize import load, assemble, scale, diagscale
from brake.solve import projection, solver, qevp
from brake.analyze import residual, visual

begin_program = timeit.default_timer()

#Initializing Parameters

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

log_level = 0
dt = datetime.date.today()
info_log_file = output_path+'info_'+dt.strftime("%d%b")+'.log' #example info_02Aug.log
time_log_file = output_path+'time_'+dt.strftime("%d%b")+'.log' #example time_02Aug.log

## creating object
obj = brake.BrakeClass(input_path, output_path, info_log_file, time_log_file, log_level)
setattr(obj, 'data_file_list', ['BMLL','BDLL','BYLL','BDIWLL','BHLL','BKLL','BKQLL','BWLL'])
setattr(obj, 'omegaRef', 1) #Reference omega
setattr(obj, 'omega_range', 17*2*math.pi)
setattr(obj, 'fRef', 1600) #Reference omega
setattr(obj, 'target', numpy.array([-10,1000,-50,12000])) #target rectangular region
setattr(obj, 'cutoff', 0.00001) 
setattr(obj, 'evs_per_shift', 30)
setattr(obj, 'desired_area_fraction', 0.99)

#Calculating exact eigenpairs of the QEVP
#--------------------------------------------------------------------------------------------------
setattr(obj, 'omega_basis', 11*2*math.pi)
la_exact, evec_exact = qevp.brake_squeal_qevp(obj, 1, obj.omega_basis)


#POD Approach
#--------------------------------------------------------------------------------------------------

wset = []
wset.append(numpy.array([1,20])*2*math.pi)
wset.append(numpy.array([1,10,20])*2*math.pi)
wset.append(numpy.array([1,5,10,15,20])*2*math.pi)
wset.append(numpy.array([1,2.5,5,7.5,10,12.5,15,17.5,20])*2*math.pi)

la_pod_set = []
evec_pod_set = []

for ctr in range(0,len(wset)):
  #Begin creating the projection matrix
  print "\n"+"\n"+'Beginning Setup Phase: Creating the Projection Matrices'

  setattr(obj, 'omega_basis', wset[ctr])  
  Q = projection.obtain_projection_matrix(obj)

  print "\n","\n","\n",'Beginning POD', Q.shape

  sparse_list = load.load_matrices(obj)

  omega = obj.omega_range
  print 'omega = ',str(omega/(2*math.pi))
  M, C, K = assemble.create_MCK(obj, sparse_list, omega)
  print 'Projecting the QEVP having dimension '+str(M.shape)
  #Projection
  QT = Q.T.conjugate()
  M =  QT.dot(M.dot(Q))
  C =  QT.dot(C.dot(Q))
  K =  QT.dot(K.dot(Q))
  print 'Onto a smaller subspace of dimension '+str(M.shape)

  n = M.shape[0]
  no_of_evs = 2*n-2;
  la_pod, evec_pod = solver.qev_dense(obj,M,C,K,no_of_evs)
  la_pod_set.append(la_pod)
  evec_pod_set.append(evec_pod)
  res_qevp = residual.residual_qevp(M,C,K,la_pod,evec_pod[0:n,:])
  print 'Maximum Residual error for the QEVP is ',max(res_qevp)


print "\n","\n",'The exact eigenvalues computed for frequency ',  obj.omega_range
brake.print_eigs(obj,la_exact,'target','terminal')

for ctr in range(0,len(wset)):
 print "\n","\n",'The projected eigenvalues for frequency ', obj.omega_range,\
      ' computed using base frequency ', wset[ctr]  
 brake.print_eigs(obj,la_pod_set[ctr],'target','terminal')

end_program = timeit.default_timer()
print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'


