'''
Function Definition
BrakeSquealQevp: Implementation of the Quadratic Eigenvalue Problem for the brakesqueal project

BrakeSquealQevp
Input
1. freq_i (the index of the base frequency)
2. path (where the data files are located)
3. data_file_list (data file names for the component matrix m,c1,c2,c3,c4,k1,k2,k3)
4. omega (omega_basis[freq_i])
5. omegaRef (reference omega value)
6. fRef (reference frequency value)
7. target (target region for the shift points)
8. evs_per_shift (number of eigenvalues per shift)

Output
1. assembled_la (eigenvalues for many shifts)
2. assembled_evec (eigenvectors for many shifts)

'''





#---------------------------------------- System Imports -----------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import sys
import math
import cmath
import timeit
import numpy as np
import scipy.io as sio
from scipy import sparse
import matplotlib.pyplot as plt
from scipy.sparse import vstack, hstack





#------------------------------------User defined imports-----------------------------------------
# Please ensure that the following files are placed in the directory where the program is run


import brake

from brake.initialize import load, assemble, shift, scale, diagscale, unlinearize
from brake.solve import cover, solver
from brake.analyze import residual



def BrakeSquealQevp(obj, freq_i, omega):
    
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i
        
    path = obj.input_path
    data_file_list = obj.data_file_list
    fRef = obj.fRef
    omegaRef = obj.omegaRef
    target = obj.target
    evs_per_shift = obj.evs_per_shift
    
    if(LOG_LEVEL):
	logger_i.info("\n"+"\n"+'-----------------------------------------------------------------')
	logger_i.info('Taking eigenvector snapshots of brakesqueal problem for frequence '\
							+str(freq_i+1)+' = '+str(omega/(2*math.pi)))
    x1 = target[0]
    x2 = target[1]
    y1 = target[2]
    y2 = target[3]

    total_area = math.fabs((x2-x1)*(y2-y1))
    desired_area_fraction = 0.99
    area_fraction_covered = 0 
    assembled_la = []
    assembled_evec = []	 
    previous_shifts = []
    previous_radius = []
    qevp_j= 1	
    while (area_fraction_covered < desired_area_fraction):
	   next_shift = cover.next_shift(target,previous_shifts,previous_radius)

	   begin = timeit.default_timer()
	   la, evec = Obtain_eigs(obj, freq_i, qevp_j, omega, next_shift)
	   end = timeit.default_timer()
	   if(LOG_LEVEL):
	     logger_t.info('\tTotal time taken by Obtain_eigs = '+"%.2f" % (end-begin)+' sec')

	   if(qevp_j==1):
	      assembled_la = la
	      assembled_evec = evec
	   else:		
	      assembled_la = np.concatenate((assembled_la,la), axis=1)
	      assembled_evec = np.concatenate((assembled_evec,evec), axis=1)
	   
	   farthest_eval_dist = 0
	   for itr in range(0,len(la)):
	      rad = cmath.polar(la[itr]-next_shift)[0] 
	      if (rad > farthest_eval_dist):
	         farthest_eval_dist = rad
           next_radius = farthest_eval_dist
	
	   previous_shifts.append(next_shift)
	   previous_radius.append(next_radius)
	   area_fraction_covered = cover.calculate_area_fraction(target, previous_shifts, \
										previous_radius)
	   if(LOG_LEVEL):
	      logger_i.info("\n"+'for Shift '+str(qevp_j)+' = '+str(next_shift)+\
' approximate area fraction covered = '+str(area_fraction_covered)) 	
	   cover.draw_circles(target,next_shift,next_radius)
	   qevp_j = qevp_j + 1
 
    ax = plt.gca()	
    ax.cla()
    return assembled_la, assembled_evec


def Obtain_eigs(obj, freq_i, qevp_j, omega, next_shift):
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i
    path = obj.input_path
    data_file_list = obj.data_file_list
    omegaRef = obj.omegaRef
    fRef = obj.fRef
    evs_per_shift = obj.evs_per_shift

    #------------------------------Initialization---------------------------------------
    tau = next_shift
    check_flag = 0
       
    #------------------------------Reading Data Files-----------------------------------
    begin_read = timeit.default_timer()
    sparse_list = load.load_matrices(obj)
    end_read = timeit.default_timer()
    
    #------------------------------Preprocessing---------------------------------------
    #---Assembling---------------
    begin_assemble = timeit.default_timer()
    M, C, K = assemble.create_MCK(obj, sparse_list, omega)
    end_assemble = timeit.default_timer()

    n = M.shape[0]
    
    M_orig = M
    C_orig = C
    K_orig = K
    
    #---Shifting---------------
    begin_shift = timeit.default_timer()
    M, C, K = shift.shift_matrices(obj, M, C, K, tau)
    end_shift = timeit.default_timer()

    #---Scaling----------------
    begin_scale = timeit.default_timer()
    M, C, K, gamma, delta = scale.scale_matrices(obj, M, C, K)
    end_scale = timeit.default_timer()\

    #---Diagonal Scaling----------------
    begin_diagscale = timeit.default_timer()
    M, C, K, D1, D2 = diagscale.diag_scale_matrices(obj,M,C,K)
    end_diagscale = timeit.default_timer()

    #----------------------------------Solver-------------------------------------------
    begin_solver = timeit.default_timer()	
    la, evec = solver.qev_optimized_sparse(obj,M,C,K);
    end_solver = timeit.default_timer()

    #Residual check for the quadratic eigenvalue problem
    if(LOG_LEVEL):
        begin_rescheck1 = timeit.default_timer()	
        res_qevp_prior = residual.residual_qevp(M,C,K,la,evec[0:n])
        end_rescheck1 = timeit.default_timer()	
        
    #---------------------------Unfolding-----------------------------------------------
    begin_unfolding = timeit.default_timer()	
    la = la*gamma
    la = tau+la
    evec = unlinearize.unlinearize_matrices(evec,n,evs_per_shift)
    evec = D2*evec;
    DD = diagscale.normalize_cols(evec)
    evec = evec * DD	
    end_unfolding = timeit.default_timer()	
    
    #---------------------------Post Processing-----------------------------------------
    #---------------------------Error Analysis------------------------------------------
    if(LOG_LEVEL):
    	begin_rescheck2 = timeit.default_timer()	
    	res_qevp_post = residual.residual_qevp(M_orig,C_orig,K_orig,la,evec)
    	end_rescheck2 = timeit.default_timer()
	
    '''
    #print 'eigen values sorted by real part'
    #def_list.my_print(la)
    '''

    print 'Eigenvalues in the positive plane for shift',qevp_j,'=',next_shift,\
    		' and frequency',freq_i+1,'=',"%.2f" % omega,'are'	
    brake.my_print_few_infile(obj,la,0)
    
    if(LOG_LEVEL):
	logger_i.info('Eigenvalues in the positive plane: ')
   	brake.my_print_few_infile(obj,la,1)

    if(check_flag):
    	sio.savemat('evec_py.mat', mdict={'data': evec})
	sio.savemat('eval_py.mat', mdict={'data': la})

           
    #---------------Logging-----------------
    if(LOG_LEVEL):
    	#----------------Logging Result--------------
    	logger_i.info('Maximum Residual error(solver) for the QEVP is '+str(max(res_qevp_prior)))
    	logger_i.info('Maximum Residual error = '+str(max(res_qevp_post)))
    
	#----------------Logging Time Complexity-----
	logger_t.info("\n"+"\n"+'------------------------------------')
	logger_t.info('Reading: '+"%.2f" % (end_read-begin_read)+' sec')
	logger_t.info('Assembling: '+"%.2f" % (end_assemble-begin_assemble)+' sec')
	logger_t.info('Shifting: '+"%.2f" % (end_shift-begin_shift)+' sec')
	logger_t.info('Scaling: '+"%.2f" % (end_scale-begin_scale)+' sec')
	logger_t.info('Diagonal Scaling: '+"%.2f" % (end_diagscale-begin_diagscale)+' sec')
	logger_t.info('Total Eigenvalue Computation Time: '+"%.2f" % (end_solver-begin_solver)+' sec')
	logger_t.info('Unfolding: '+"%.2f" % (end_unfolding-begin_unfolding)+' sec')
        logger_t.info('Error Analysis1(standard QEVP): '+"%.2f" % (end_rescheck1-begin_rescheck1)+' sec')
	logger_t.info('Error Analysis2(after unfolding): '+"%.2f" % (end_rescheck2-begin_rescheck2)+' sec')
	logger_t.info('------------------------------------')    
    return la, evec 	

