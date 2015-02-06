r"""
This module defines the following functions::

  - brake_squeal_qevp:
  
    For a particuar base angular frequency this function assembles the eigenvalues
    and eigenvectors for different shift points in the target region.
    
  - Obtain_eigs:
    
    For a particuar base angular frequency and for a particular shift point this 
    function obtains obj.evs_per_shift eigenvalues and the corresponding eigenvectors
  
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import math
import cmath
import numpy
import timeit
import scipy.io
import matplotlib.pyplot

#----------------------------Application Specific Imports----------------------------------------- 
import brake
from brake.initialize import load, assemble, shift, scale, diagscale, unlinearize
from brake.solve import cover, solver
from brake.analyze import residual



def brake_squeal_qevp(obj, freq_i, omega):
    r"""
        
    :param obj: object of the class ``BrakeClass``
    :param freq_i: the index of the base angular freq
    :param omega: ith base angular freq
    :return: ``assembled_la`` - assembled eigenvalues, ``assembled_evec`` -- assembled eigenvectors
    
    Procedure::
    
     The assembled_la and assembled_evec are obtained as follows:
       
     - calculate the next shift point in the target region
     - obtain eigenvalues and eigenvectors for that particular shift point
     - add the eigenvalues and eigenvectors to assembled_la and assembled_evec respectively
     - check if the required area fraction of the target region has been covered. If yes return
       assembled_la and assembled_evec else calculate the next shift point in the target
       region and repeat

    """
    
    #object attributes used in the function
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i
    target = obj.target
    desired_area_fraction = obj.desired_area_fraction
    desiredCount = obj.desiredCount
        
    if(LOG_LEVEL):
        logger_i.info("\n"+"\n"+'-----------------------------------------------------------------')
        logger_i.info('Taking eigenvector snapshots of brakesqueal problem for frequence '\
                                                        +str(freq_i)+' = '+str(omega/(2*math.pi)))
    x1 = target[0]
    x2 = target[1]
    y1 = target[2]
    y2 = target[3]

    total_area = math.fabs((x2-x1)*(y2-y1))
    area_fraction_covered = 0
    eigCount = 0 
    assembled_la = []
    assembled_evec = []  
    previous_shifts = []
    previous_radius = []
    qevp_j= 1   
    while ((area_fraction_covered < desired_area_fraction) | (eigCount < desiredCount)):
           next_shift = cover.next_shift(obj, previous_shifts, previous_radius)

           begin = timeit.default_timer()
           la, evec = Obtain_eigs(obj, freq_i, qevp_j, omega, next_shift)
           end = timeit.default_timer()
           
           if(LOG_LEVEL):
             logger_t.info('\tTotal time taken by Obtain_eigs = '+"%.2f" % (end-begin)+' sec')

            #obtain radius = distance of the farthest eigenvalue from next_shift
           farthest_eval_dist = 0
           for itr in range(0,len(la)):
              rad = cmath.polar(la[itr]-next_shift)[0] 
              if (rad > farthest_eval_dist):
                 farthest_eval_dist = rad
           next_radius = farthest_eval_dist
           
           #Update Phase
           if(qevp_j==1):
              assembled_la = la
              assembled_evec = evec
           else:                
              assembled_la = numpy.concatenate((assembled_la,la), axis=1)
              assembled_evec = numpy.concatenate((assembled_evec,evec), axis=1)
           
           previous_shifts.append(next_shift)
           previous_radius.append(next_radius)
           
           area_fraction_covered = cover.calculate_area_fraction(obj, previous_shifts, \
                                                                                previous_radius)
           
	   eigCount = eigCount + len(la)
	   if(LOG_LEVEL):
              logger_i.info("\n"+'for Shift '+str(qevp_j)+' = '+str(next_shift)+\
                             ' approximate area fraction covered = '+str(area_fraction_covered))
              
           cover.draw_circles(obj, next_shift, next_radius)
           qevp_j = qevp_j + 1
 
    ax = matplotlib.pyplot.gca()        
    ax.cla()
    return assembled_la, assembled_evec


def Obtain_eigs(obj, freq_i, qevp_j, omega, next_shift):
    r"""
        
    :param obj: object of the class ``BrakeClass``
    :param freq_i: the index of the base angular freq
    :param qevp_j: the index of the shift point
    :param omega: ith base angular freq
    :param next_shift: jth shift point in the target region
    :return: ``la`` - eigenvalues, ``evec`` - eigenvectors
    
    Procedure::
    
     The la and evec are obtained as follows:
       
     - load the various component matrices
     - assemble the various component matrices together(for the given angular frequency
       omega) to form the mass(M), damping(C) and stiffness matrix(K).
     - because we are interested in inner eigenvalues around certain shift points
       next_shift, so we transform the qevp using shift and invert spectral
       transformations.
     
    """
  
    #object attributes used in the function
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i
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
    
    #save the matrices before shift
    M_orig = M
    C_orig = C
    K_orig = K
    
    #---Shifting---------------
    begin_shift = timeit.default_timer()
    M, C, K = shift.shift_matrices(obj, M, C, K, tau)
    end_shift = timeit.default_timer()

    #---Scaling----------------
    begin_scale = timeit.default_timer()
    M, C, K, gamma = scale.scale_matrices(obj, M, C, K)
    end_scale = timeit.default_timer()\

    #---Diagonal Scaling----------------
    begin_diagscale = timeit.default_timer()
    M, C, K, D2 = diagscale.diag_scale_matrices(obj,M,C,K)
    end_diagscale = timeit.default_timer()
    
    #----------------------------------Solver-------------------------------------------
    begin_solver = timeit.default_timer()       
    la, evec = solver.qev_sparse(obj,M,C,K);
    end_solver = timeit.default_timer()
    
    
    #Residual check for the quadratic eigenvalue problem
    if(LOG_LEVEL):
        begin_rescheck1 = timeit.default_timer()        
        res_qevp_prior = residual.residual_qevp(M,C,K,la,evec[0:n])
        end_rescheck1 = timeit.default_timer()  
    
    
    #---------------------------Unfolding-----------------------------------------------
    begin_unfolding = timeit.default_timer()
    
    #obtaining the original eigenvalues
    #undo scaling
    la = la*gamma
    #undo shifting
    la = tau+la
    
    #obtaining the original eigenvectors
    evec = unlinearize.unlinearize_matrices(evec) #evec = evec_after_diagscale
    evec = D2*evec; #D2 = DR
    
    #normalizing the eigenvector(why is it needed ?)
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
    
    print('Eigenvalues for shift',qevp_j,'=',next_shift,\
                ' and frequency',freq_i,'=',"%.2f" % omega,'are')      
    brake.printEigs(obj,la,'target','terminal')
    
    if(LOG_LEVEL):
        logger_i.info('Eigenvalues : ')
        brake.printEigs(obj,la,'target','file')

    if(check_flag):
        scipy.io.savemat('evec_py.mat', mdict={'data': evec})
        scipy.io.savemat('eval_py.mat', mdict={'data': la})

           
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

