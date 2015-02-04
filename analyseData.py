# This file reads the component matrices and displays there sparsity pattern
# and other properties

#----------------------------------Standard Library Imports---------------------------------------
import timeit
import matplotlib.pylab as plt
import scipy.sparse.linalg as norm

#----------------------------Application Specific Imports-----------------------------------------
import brake
import createBrakeClassObject
from brake.initialize import load


begin_program = timeit.default_timer()


#----------------------------Create Object--------------------------------------------------------
obj = createBrakeClassObject.returnObject(20)

obj.logger_i.info("\n"+"\n"+'Beginning Data Analysis')
obj.logger_i.info('------------------------------------------------------------------------')

print "\n"+"\n"+'Beginning Data Analysis'

sparse_list = load.load_matrices(obj)

obj.logger_i.info("\n"+'Matrices in CSC format converted to CSR')

obj.logger_i.info("\n\n"+'Properties of various component matrices'+"\n")

for i in range(0,len(sparse_list)):
	componentMatrix = sparse_list[i]
	normMatrix = norm.onenormest(componentMatrix.tocsr(), t=3, itmax=5, compute_v=False, compute_w=False)
	print obj.data_file_list[i], obj.data_file_name[i], ' NonZeros = ', componentMatrix.nnz, ' 1-Norm = ', normMatrix
	obj.logger_i.info(obj.data_file_list[i]+' '+obj.data_file_name[i]+' Nonzeros = '+str(componentMatrix.nnz)+' 1-Norm = '+str(normMatrix))
	
	#saving the sparsity pattern
	fig = plt.figure(figsize=(24.0, 15.0))
	fig.clf()
	fig.gca().add_artist(plt.spy(componentMatrix))
	brake.save(obj.output_path+'dataAnalysis/'+obj.data_file_name[i], ext="png", close=True, verbose=False)
        

end_program = timeit.default_timer()

print "\n","\n","\n",'Total Run Time = : '+"%.2f" % (end_program-begin_program)+' sec'

if(obj.log_level):
    obj.logger_i.info("\n"+"\n"+'Finished Data Analysis')

    #----------------Logging Time Complexity-----
    obj.logger_t.info('Time Taken for data analysis: '+"%.2f" % (end_program-begin_program)+' sec')
