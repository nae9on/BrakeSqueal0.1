BrakeSqueal: Python Implementation of Parametric Model Reduction in Disc Brake Modelling
========================================

Requirements
------------
Standard Python Packages (docutils, matplotlib, numpy, scipy, sphinx, vtk)


Getting Started
---------------


Examples
---------------

Installation
---------------
Install using `python setup.py install`



Documentation
---------------
http://sphinx-doc.org/tutorial.html

To update the current documentation perform the following steps.

1. When documenting Python code, it is common to put a lot of documentation in the source files, 
in documentation strings. Sphinx supports the inclusion of docstrings from the modules with an 
extension called “autodoc”. The autodoc feature has already been activated in conf.py. 
To add suitable docstrings in the current python code follow
http://sphinx-doc.org/ext/autodoc.html#module-sphinx.ext.autodoc 

2. When adding a new module with docstrings update doc/index.rst

3. Update the git code repository.

4. Update documentation using `BrakeSqueal0.1/doc/ make html` (or make latexpdf, for more options check the inculded Makefile).

html documenation available in `BrakeSqueal0.1/doc/_build/html/index.html`
latex documenation available in `BrakeSqueal0.1/doc/_build/latex/BrakeSqueal.pdf`