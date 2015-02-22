r"""
This module defines the following functions::

  - next_shift
  
    Implementation of the MonteCarlo Algorithm for choosing the next shift point 
    in the target region.
  
  - calculate_area_fraction
  
    Calculates the area fraction covered(of the target rectangle) by the chosen shift points.
    
  - draw_circles
  
    plots a circle corresponding to the next shift point and appends it to the existing plot.
    Thus creates a simulation showing how the target region is covered by the various shift points
    chosen on fly.
  
     
"""

import math
import cmath
import random
import numpy as np
import matplotlib.pyplot as plt

import brake
#Define exceptions
class CoverError(Exception): pass
class Cover_BadInputError(CoverError): pass

def next_shift(obj, previous_shifts = [], previous_radius = []):
        r"""
        :param obj: object of the class ``BrakeClass``
        :param previous_shifts: python list for the previous shift points already calculated
        :param previous_radius: corresponding radius of the previous shift points
        :return: ``next_shift`` - next shift point in the target region
        :raises: ``Cover_BadInputError``, When the provided input is not as expected
        
        """
        
        target = obj.target
        
        # Exceptions
        #------------------------------------------------------------
        if len(target) == 0:
         raise Cover_BadInputError('The target is not specified')
        elif len(target) < 4:
         raise Cover_BadInputError('target region does not define a reactangle')
        elif target[1] <= target[0] or target[3] <= target[2]:
         raise Cover_BadInputError('the target reactangle is not defined')
          
        if (len(previous_shifts)==0) and  (len(previous_radius)==0):
                first_shift = complex((target[0]+target[1])/2,(target[2]+target[3])/2)
                next_shift = first_shift
                
        else:
                
                x1 = target[0]
                x2 = target[1]
                y1 = target[2]
                y2 = target[3]
                
                if len(previous_shifts) != len(previous_radius):
                    raise Cover_BadInputError('The shift and radius vector are not of the same length')
                
                #check if previous shifts lie in the target region
                for i in range(0,len(previous_shifts)):
                  if previous_shifts[i].real < x1 \
                        or previous_shifts[i].real > x2 \
                        or previous_shifts[i].imag < y1 \
                        or previous_shifts[i].imag > y2:
                    raise Cover_BadInputError('shift point not in target region')
                          
                if np.amin(previous_radius) <= 0:
                 raise Cover_BadInputError('The radius is nonpositive')
        
                #generate 1000 random sampling points and radius in the target region               
                sampling_points = []
                sampling_radius = []
                for i in range(1000):
                  x = random.randint(x1,x2)
                  y = random.randint(y1,y2)
                  sampling_points.append(complex(x,y))
                  sampling_radius.append(random.uniform(0,1))

                mean_radius = np.mean(sampling_radius)
                mean_previous_radius = np.mean(previous_radius)
                sampling_radius[:] = [x*mean_previous_radius/mean_radius for x in sampling_radius]  

                #scanning for the best shift
                area_increment = 0
                for i in range(0,len(sampling_points)):
                  increment = 10000000
                  for j in range(0,len(previous_shifts)):
                      delta = cmath.polar(sampling_points[i]-previous_shifts[j])[0] - previous_radius[j]                
                      if delta < 0:
                        # Point i is inside Circle j   
                        increment = 0                   
                        break
                      else:
                        if delta < increment:
                          increment =  delta
                  max_increment = [increment]
                  max_increment.append(np.fabs(sampling_points[i].real-x1))
                  max_increment.append(np.fabs(sampling_points[i].real-x2))
                  max_increment.append(np.fabs(sampling_points[i].imag-y1))
                  max_increment.append(np.fabs(sampling_points[i].imag-y2))
                  increment = min(max_increment)        
                  if(increment > area_increment):
                     area_increment = increment
                     next_shift = sampling_points[i]
                
        return next_shift

def calculate_area_fraction(obj, previous_shifts, previous_radius):
        r"""
        :param obj: object of the class ``BrakeClass``
        :param previous_shifts: python list for the previous shift points already calculated
        :param previous_radius: corresponding radius of the previous shift points
        :return: ``area_fraction_covered`` - the total area fraction covered with the chosen shift points
        :raises: ``Cover_BadInputError``, When the provided input is not as expected
                
        """
        
        target = obj.target
        
        # Exceptions
        #------------------------------------------------------------
        if len(target)==0:
         raise Cover_BadInputError('The target is not specified')
        elif len(target) < 4:
         raise Cover_BadInputError('target region does not define a reactangle')
        elif target[1] <= target[0] or target[3] <= target[2]:
         raise Cover_BadInputError('the target reactangle is not defined')

        if len(previous_shifts)==0:
         raise Cover_BadInputError('The previous_shifts is not specified')
        
        if len(previous_radius)==0:
         raise Cover_BadInputError('The previous_radius is not specified')
        
        if len(previous_shifts) != len(previous_radius):
            raise Cover_BadInputError('The shift and radius vector are not of the same length')

        x1 = target[0]
        x2 = target[1]
        y1 = target[2]
        y2 = target[3]
                
        #check if previous shifts lie in the target region
        for i in range(0,len(previous_shifts)):
          if previous_shifts[i].real < x1 \
                or previous_shifts[i].real > x2 \
                or previous_shifts[i].imag < y1 \
                or previous_shifts[i].imag > y2:
            raise Cover_BadInputError('shift point not in target region')
                          
        if np.amin(previous_radius) <= 0:
          raise Cover_BadInputError('The radius is nonpositive')
        #------------------------------------------------------------
                  
        inn = 0.0
        pool = 10000
        for i in range(0,pool):
                  x = random.randint(x1,x2)
                  y = random.randint(y1,y2)
                  for j in range(0,len(previous_shifts)):
                     r = cmath.polar(complex(x,y)-previous_shifts[j])[0]
                     if r<previous_radius[j]:
                      inn = inn+1
                      break
        area_fraction_covered = inn/pool
        return area_fraction_covered                       

def draw_circles(obj, next_shift, next_radius):
        r"""
        :param obj: object of the class ``BrakeClass``
        :param next_shift: the shift point to be shown on the plot
        :param next_radius: the radius of the shift point to be shown
        :return: plots a circle corresponding to the next shift point and appends it to the existing plot.
        :raises: ``Cover_BadInputError``, When the provided input is not as expected
        """
        
        target = obj.target
        # Exceptions
        #------------------------------------------------------------
        if len(target)==0:
         raise Cover_BadInputError('The target is not specified')
        elif len(target) < 4:
         raise Cover_BadInputError('target region does not define a reactangle')
        elif target[1] <= target[0] or target[3] <= target[2]:
         raise Cover_BadInputError('the target reactangle is not defined')
        
        if not next_radius:
         raise Cover_BadInputError('The previous_radius is not specified')
        elif next_radius <= 0:
         raise Cover_BadInputError('The radius is nonpositive')
                                                                                
        x1 = target[0]
        x2 = target[1]
        y1 = target[2]
        y2 = target[3]

        if not next_shift:
         raise Cover_BadInputError('The next_shift array is not length 1')
        elif next_shift.real < x1 \
                or next_shift.real > x2 \
                or next_shift.imag < y1 \
                or next_shift.imag > y2:
         raise Cover_BadInputError('shift point not in target region')
        #------------------------------------------------------------

        x1 = target[0]
        x2 = target[1]
        y1 = target[2]
        y2 = target[3]

        fig = plt.gcf()
        circle=plt.Circle((next_shift.real,next_shift.imag),next_radius,color='b',fill=False)
        react = plt.Rectangle((x1,y1),abs(x2-x1),abs(y2-y1),color='black',fill=False)
        ax = plt.gca()
        #ax.cla() # clear things for fresh plot
        ax.set_xlim((y1-1000,y2+1000))
        ax.set_ylim((y1-1000,y2+1000))
        plt.axhline(0, color='blue')
        plt.axvline(0, color='blue')
        fig.gca().add_artist(circle)
        fig.gca().add_artist(react)
        brake.save(obj.output_path+'drawCircles', ext="png", close=True, verbose=False)
        return
