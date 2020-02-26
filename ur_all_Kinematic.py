#!/usr/bin/python2
__author__ = "Lai-Hsuan Liu"
__UR_version__ = "3.11.0"
__maintainer__ = "Lai-Hsuan Liu""
__email__ = "carlosstarw@gmail.com"
__status__ = "Prototype"

import numpy as np
from numpy import linalg

import cmath
import math
from math import cos 
from math import sin 
from math import atan2 
from math import acos
from math import asin
from math import sqrt
from math import pi


# ****** D-H Table values ******

global mat
global d1, a2, a3, d4, d5, d6
global d, a, alph

mat=np.matrix
a1 = a4 = a5 = a6 =0
d2 = d3 = 0

# ur5, ur10
# d1 =  0.1273
# d4 =  0.163941
# d5 =  0.1157
# d6 =  0.0922

# a2 = -0.612
# a3 = -0.5723

# ur3
d1 = 0.1519 
d4 = 0.11235
d5 = 0.08535
d6 = 0.0819

a2 = -0.24365
a3 = -0.21325

d = mat([0.1519 , 0, 0, 0.11235, 0.08535, 0.0819]) #ur3
a = mat([0 ,-0.24365 ,-0.21325 ,0 ,0 ,0]) #ur3
alph = mat([math.pi/2, 0, 0, math.pi/2, -math.pi/2, 0 ])  #ur3

# d = mat([0.089159, 0, 0, 0.10915, 0.09465, 0.0823]) #ur5
# a =mat([0 ,-0.425 ,-0.39225 ,0 ,0 ,0]) #ur5
# alph = mat([math.pi/2, 0, 0, math.pi/2, -math.pi/2, 0 ])  #ur5

# d = mat([0.1273, 0, 0, 0.163941, 0.1157, 0.0922])#ur10 mm
# a =mat([0 ,-0.612 ,-0.5723 ,0 ,0 ,0])#ur10 mm
# alph = mat([pi/2, 0, 0, pi/2, -pi/2, 0 ]) # ur10

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

def rm2theta(rm):
  beta = atan2(-rm[2, 0], math.sqrt(rm[0, 0]*rm[0, 0] + rm[1, 0]*rm[1, 0])) # y
  alfa = atan2(rm[1, 0]/cos(beta), rm[0, 0]/cos(beta)) # z
  gamma = atan2(rm[2, 1]/cos(beta), rm[2, 2]/cos(beta)) # x
  print("beta, alfa, gamma:" + str(beta)  + str(alfa)  + str(gamma))
  return np.array([gamma, beta, alfa])


# ************************************************** FORWARD KINEMATICS
# homogeneous coordinate
def HTM(n ,th ,c):
  T_a = mat(np.identity(4), copy=False)
  T_a[0,3] = a[0,n-1]
  T_d = mat(np.identity(4), copy=False)
  T_d[2,3] = d[0,n-1]

  Rzt = mat([[cos(th[n-1,c]), -sin(th[n-1,c]), 0 ,0],
	         [sin(th[n-1,c]),  cos(th[n-1,c]), 0, 0],
	         [0,               0,              1, 0],
	         [0,               0,              0, 1]],copy=False)
      

  Rxa = mat([[1, 0,                 0,                  0],
			 [0, cos(alph[0,n-1]), -sin(alph[0,n-1]),   0],
			 [0, sin(alph[0,n-1]),  cos(alph[0,n-1]),   0],
			 [0, 0,                 0,                  1]],copy=False)

  A_i = T_d * Rzt * T_a * Rxa
  return A_i

def HTrans(th ,c):  
  A_1=HTM( 1,th,c  )
  A_2=HTM( 2,th,c  )
  A_3=HTM( 3,th,c  )
  A_4=HTM( 4,th,c  )
  A_5=HTM( 5,th,c  )
  A_6=HTM( 6,th,c  )
      
  T_06=A_1*A_2*A_3*A_4*A_5*A_6

  return T_06

# ************************************************** INVERSE KINEMATICS 
def invKine(desired_pos):# Get 8 solution from InvKine
  th = mat(np.zeros((6, 8)))
  P_05 = (desired_pos * mat([0,0, -d6, 1]).T-mat([0,0,0,1 ]).T)
  
  # **** theta1 ****
  
  psi = atan2(P_05[2-1,0], P_05[1-1,0])
  phi = acos(d4 /sqrt(P_05[2-1,0]*P_05[2-1,0] + P_05[1-1,0]*P_05[1-1,0]))
  #The two solutions for theta1 correspond to the shoulder
  #being either left or right
  th[0, 0:4] = pi/2 + psi + phi
  th[0, 4:8] = pi/2 + psi - phi
  th = th.real
  
  # **** theta5 ****
  
  cl = [0, 4]# wrist up or down
  for i in range(0,len(cl)):
	      c = cl[i]
	      T_10 = linalg.inv(HTM(1,th,c))
	      T_16 = T_10 * desired_pos
	      th[4, c:c+2] = + acos((T_16[2,3]-d4)/d6)
	      th[4, c+2:c+4] = - acos((T_16[2,3]-d4)/d6)

  th = th.real
  
  # **** theta6 ****
  # theta6 is not well-defined when sin(theta5) = 0 or when T16(1,3), T16(2,3) = 0.

  cl = [0, 2, 4, 6]
  for i in range(0,len(cl)):
	      c = cl[i]
	      T_10 = linalg.inv(HTM(1,th,c))
	      T_16 = linalg.inv( T_10 * desired_pos )
	      th[5, c:c+2] = atan2((-T_16[1,2]/sin(th[4, c])),(T_16[0,2]/sin(th[4, c])))
		  
  th = th.real

  # **** theta3 ****
  cl = [0, 2, 4, 6]
  for i in range(0,len(cl)):
	      c = cl[i]
	      T_10 = linalg.inv(HTM(1,th,c))
	      T_65 = HTM( 6,th,c)
	      T_54 = HTM( 5,th,c)
	      T_14 = ( T_10 * desired_pos) * linalg.inv(T_54 * T_65)
	      P_13 = T_14 * mat([0, -d4, 0, 1]).T - mat([0,0,0,1]).T
	      t3 = cmath.acos((linalg.norm(P_13)**2 - a2**2 - a3**2 )/(2 * a2 * a3)) # norm ?
	      th[2, c] = t3.real
	      th[2, c+1] = -t3.real

  # **** theta2 and theta 4 ****

  cl = [0, 1, 2, 3, 4, 5, 6, 7]
  for i in range(0,len(cl)):
	      c = cl[i]
	      T_10 = linalg.inv(HTM( 1,th,c ))
	      T_65 = linalg.inv(HTM( 6,th,c))
	      T_54 = linalg.inv(HTM( 5,th,c))
	      T_14 = (T_10 * desired_pos) * T_65 * T_54
	      P_13 = T_14 * mat([0, -d4, 0, 1]).T - mat([0,0,0,1]).T
	      
	      # theta 2
	      th[1, c] = -atan2(P_13[1], -P_13[0]) + asin(a3* sin(th[2,c])/linalg.norm(P_13))
	      # theta 4
	      T_32 = linalg.inv(HTM( 3,th,c))
	      T_21 = linalg.inv(HTM( 2,th,c))
	      T_34 = T_32 * T_21 * T_14
	      th[3, c] = atan2(T_34[1,0], T_34[0,0])
  th = th.real

  return th

# Test """Joint Angle (in degrees) reading from encoders."""
theta1 = np.radians(45.0)   
theta2 = np.radians(-60.0)  
theta3 = np.radians(90.0)
theta4 = np.radians(0.0)
theta5 = np.radians(0.0)
theta6 = np.radians(0.0)

th = np.matrix([[theta1], [theta2], [theta3], [theta4], [theta5], [theta6]])
c = [0]
location = HTrans(th,c )

print(location)
print("rm2theta:" + str(rm2theta(location)))


# desired_pos = location

# th = invKine(desired_pos)
# print(th)

desired_pos = np.matrix(
[[ 3.06161700e-17, 8.66025404e-01, -5.00000000e-01, 3.63537488e-01],
[-1.00000000e+00, 0.00000000e+00, -5.55111512e-17, -1.09150000e-01],
[-5.55111512e-17, 5.00000000e-01, 8.66025404e-01, 4.25598256e-01],
[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00]]
)

th = invKine(desired_pos)
print(th)


