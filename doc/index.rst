.. BrakeSqueal documentation master file, created by
   sphinx-quickstart on Mon Aug  4 16:04:03 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BrakeSqueal's documentation!
=======================================

Contents:

.. toctree::
   :maxdepth: 2
   
brake
=======================================

__init__
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/__init__.py?at=master>`__

.. automodule:: brake.__init__
    :members:
    
.. seealso:: `logger`_  brake.initialize.logger
   
initialize
=======================================

logger
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/logger.py?at=master>`__

.. automodule:: brake.initialize.logger
    :members:
    
.. seealso:: `Python Logging <https://docs.python.org/2/library/logging.html>`__  
        
load
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/load.py?at=master>`__

.. automodule:: brake.initialize.load
    :members:
    
.. seealso:: `scipy.io.loadmat <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.io.loadmat.html>`__

assemble
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/assemble.py?at=master>`__

.. automodule:: brake.initialize.assemble
    :members:
    
.. seealso:: `scipy.sparse.linalg.onenormest <http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.sparse.linalg.onenormest.html>`__

shift
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/shift.py?at=master>`__

.. automodule:: brake.initialize.shift
    :members:
    
.. seealso:: `scipy.sparse.linalg.onenormest <http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.sparse.linalg.onenormest.html>`__
      
scale
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/scale.py?at=master>`__

.. automodule:: brake.initialize.scale
    :members:
    
.. seealso:: `scipy.sparse.linalg.onenormest <http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.sparse.linalg.onenormest.html>`__
      
diagscale
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/diagscale.py?at=master>`__

.. automodule:: brake.initialize.diagscale
    :members:
    
.. seealso::

   `scipy.sparse.linalg.eigs <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.eigs.html>`__  
      Documentation of the Python eigs command
      
unlinearize
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/initialize/unlinearize.py?at=master>`__

.. automodule:: brake.initialize.unlinearize
    :members:   

.. seealso:: `diagscale`_  brake.initialize.diagscale
      
solve
=======================================

projection
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/solve/projection.py?at=master>`__

.. automodule:: brake.solve.projection
    :members:

.. seealso:: `scipy.linalg.svd <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.linalg.svd.html>`_

qevp
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/solve/qevp.py?at=master>`__

.. automodule:: brake.solve.qevp
    :members:

solver
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/solve/solver.py?at=master>`__

.. automodule:: brake.solve.solver
    :members:
    
.. seealso::

   `scipy.sparse.linalg.eigs <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.eigs.html>`__  
      Documentation of the Python eigs command
      
   `scipy.sparse.linalg.splu <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.splu.html>`__
      Documentation of the Python splu command       
   
   `scipy.sparse.linalg.spilu <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.spilu.html>`__
      Documentation of the Python spilu command 
      
   `scipy.linalg.lu_factor <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.linalg.lu_factor.html>`__
      Documentation of the Python lu_factor command 
      
   `scipy.linalg.lu_solve <http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.linalg.lu_solve.html>`__
      Documentation of the Python lu_solve command       
      
cover
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/solve/cover.py?at=master>`__

.. automodule:: brake.solve.cover
    :members:
    
.. seealso::

   `matplotlib.pyplot <http://matplotlib.org/api/pyplot_api.html>`__  
      Documentation of the Python Matplotlib library
 
   `randon.py <https://docs.python.org/2/library/random.html>`__
      Documentation of the Python random number generator module
      
analyze
=======================================

residual
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/analyze/residual.py?at=master>`__

.. automodule:: brake.analyze.residual
    :members:
    
.. seealso::

   `scipy.sparse.linalg.onenormest <http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.sparse.linalg.onenormest.html>`__

   `numpy.linalg.norm <http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.norm.html>`__
     
visual
---------------------------------------
`Source <https://bitbucket.org/akadar/brakesqueal0.1/src/master/brake/analyze/visual.py?at=master>`__

.. automodule:: brake.analyze.visual
    :members:
    
.. seealso::

   `matplotlib.pyplot <http://matplotlib.org/api/pyplot_api.html>`__  
      Documentation of the Python Matplotlib library
      
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

