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

def rm2theta(rm):
  beta = atan2(-rm[2, 0], math.sqrt(rm[0, 0]*rm[0, 0] + rm[1, 0]*rm[1, 0])) # y
  alfa = atan2(rm[1, 0]/cos(beta), rm[0, 0]/cos(beta)) # z
  gamma = atan2(rm[2, 1]/cos(beta), rm[2, 2]/cos(beta)) # x
  print("beta, alfa, gamma:" + str(beta)  + str(alfa)  + str(gamma))
  return np.array([gamma, beta, alfa])
