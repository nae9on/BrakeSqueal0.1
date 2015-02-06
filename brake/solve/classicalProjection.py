r"""
This module defines the following functions::

  - obtain_projection_matrix:
  
   This function forms the Projection Matrix by solving the quadratic eigenvalue problem 
   for each base angular frequency.
   
   - obtain_measurment_matrix
   
   This function forms the Measurment Matrix by solving the quadratic eigenvalue problem 
   for each base angular frequency.
"""

#----------------------------------Standard Library Imports---------------------------------------
# Please ensure that the following libraries are installed on the system prior 
# running the program.
import timeit
import numpy
import scipy

#----------------------------Application Specific Imports----------------------------------------- 
import brake
from brake.initialize import load, assemble, shift, scale, diagscale, unlinearize
from brake.solve import solver
from brake.analyze import residual

def Obtain_eigs(obj,no_of_evs):

    #object attributes used in the function
    LOG_LEVEL = obj.log_level
    logger_t = obj.logger_t
    logger_i = obj.logger_i

    #------------------------------Initialization---------------------------------------
    #Making the shift point as the center of the target rectangular region
    tau = complex((obj.target[0]+obj.target[1])/2,(obj.target[2]+obj.target[3])/2)
    check_flag = 0

    #------------------------------Reading Data Files-----------------------------------
    begin_read = timeit.default_timer()
    sparse_list = load.load_matrices(obj)
    end_read = timeit.default_timer()

    #------------------------------Preprocessing---------------------------------------
    #---Assembling---------------
    begin_assemble = timeit.default_timer()
    M = sparse_list[0].tocsr()
    C = 0*M
    K = sparse_list[5].tocsr() + sparse_list[7].tocsr()
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
    saveEVS = obj.evs_per_shift
    obj.evs_per_shift = no_of_evs
    begin_solver = timeit.default_timer()
    la, evec = solver.qev_sparse(obj,M,C,K);

    end_solver = timeit.default_timer()
    obj.evs_per_shift = saveEVS

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

    print('Calculating eigenpairs of the QEVP corressponding to the Classical Projection')
    brake.printEigs(obj,la,'target','terminal')
    print('largest real part of obtained eigenvalue = '),numpy.max(la.real)

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