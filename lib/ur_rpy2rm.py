
import numpy as np

import cmath
import math
from math import cos 
from math import sin 
from math import atan2 
from math import acos
from math import asin
from math import sqrt
from math import pi

# Z, Y, X
def theta2rm(alfa, beta, gamma):
  rm = np.zeros(3,3)
  rm[0, 0] = cos(alfa) * cos(beta)
  rm[0, 1] = sin(alfa) * cos(beta)
  rm[0, 2] = -sin(beta)

  rm[1, 0] = cos(alfa) * sin(beta) * sin(gamma) - sin(alfa) * cos(gamma)
  rm[1, 1] = sin(alfa) * sin(beta) * sin(gamma) + cos(alfa) * cos(gamma)
  rm[1, 2] = cos(beta) * sin(gamma)
  
  rm[2, 0] = cos(alfa) * sin(beta) * cos(gamma) + sin(alfa) * sin(gamma)
  rm[2, 1] = sin(alfa) * sin(beta) * cos(gamma) - cos(alfa) * sin(gamma)
  rm[2, 2] = cos(beta) * cos(gamma)

  return rm
