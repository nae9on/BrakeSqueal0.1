#!/usr/bin/env python

from distutils.core import setup

setup(
    name='BrakeSqueal',
    version='0.1',
    description='Implementation of Parametric Model Reduction in Disc Brake Modelling',
    author='Ali H. Kadar',
    author_email='kadar@math.tu-berlin.de',
    url='https://bitbucket.org/akadar/brakesqueal',
    packages=['brake','brake.initialize','brake.solve','brake.analyze'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
)
