# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 13:14:07 2018

@author: Vegard Broen
"""

import numpy as np
import matplotlib.pyplot as plt

h = 0.0001 #step size
N = 50000 #number of steps


g = 9.81 #Gravitational acceleration
m = 0.01 #Rass
k = 0.01 #Friction constant (air resistance, vibrations etc.)
r = 0.02 #Radius of object
I_0 = 2*(m*r**2)/5 #Moment of inertia solid sphere


t_0 = 0
v_0 = 0
x_0 = 10

t = np.zeros(N+1)
v = np.zeros(N+1)
x = np.zeros(N+1)
E = np.zeros(N+1)
a = np.zeros(N+1)


t[0] = t_0
v[0] = v_0
x[0] = x_0
E[0] = m*g*(x_0-5)**2


t_old = t_0
v_old = v_0
x_old = x_0

#return angle for a position x
def alpha(x):
    return -np.arctan(2*x-10)

#time step for position and velocity
def newx(x_old, v_old):
    return (x_old + v_old*h/np.sqrt(1 + (alpha(x_old))**2))

def newv(x_old, v_old):
    return (v_old + h*(g*np.sin(alpha(x_old)) - k*v_old/m)/(1+I_0/(m*r**2)))

for n in range(N): #Euler's method
    
    x_new = newx(x_old, v_old)
    v_new = newv(x_old, v_old)
    
    t[n+1] = t_old + h
    v[n+1] = v_new
    x[n+1] = x_new
    E[n+1] = m*(g*(x_new-5)**2 + 0.5*(v_new)**2*(1+I_0/(r**2)))
    a[n] = (v[n+1]-v[n])/h 
    
    t_old = t[n+1]
    v_old = v_new
    x_old = x_new






plt.figure()
plt.plot(t, x, 'r')
plt.show()
    
plt.figure()
plt.plot(t, v, 'b')
plt.show()

plt.figure()
plt.plot(t, a, 'y')
plt.show()

plt.figure()
plt.plot(t, E, 'g')
plt.show()


    
    
    
    