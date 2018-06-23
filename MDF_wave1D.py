#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
d2u    1    d2u
--- = --- X ---   Solving Wave equation 1D by Finite Differences
dx2    c2   dt2


Code Written by Felipe Timoteo | felipetimoteo@id.uff.br
                 Last update: May 23th, 2018 

Copyright (C) 2016 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
                   Departamento de Geologia e Geofísica
                   Universidade Federal Fluminense
"""

from numpy import arange,size,pi,sqrt,exp,zeros,copy
from matplotlib.pylab import figure,plot,show,draw,pause,clf
from time import sleep
import matplotlib.animation as animation

# Parameters
dx  = 0.5
dt  = 0.0002
T   = 5
L   = 1000 


x  = arange(0,L+dx,dx)
t  = arange(0,T+dt,dt)
Nx = size(x)
Nt = size(t)

fc = 30.
c  = 1500.
r = (c*dt/dx)*(c*dt/dx)


# Source
td  = t - 2*sqrt(pi)/fc
fcd = fc/(sqrt(pi)*3)
source = (1 - 2*pi*(pi*fcd*td)*(pi*fcd*td))*exp(-pi*(pi*fcd*td)*(pi*fcd*td))


#figure(1)
#plot(source)
#draw()

# Initial Conditions

Uc = zeros(Nx)
Up = zeros(Nx)
Uf = zeros(Nx)

# Solving Differential Equation
for k in arange(1,Nt+1):
    
    if k < size(source)-1:
        Uc[int(round(Nx/2))] = Uc[int(round(Nx/2))] + source[k]
          
    for i in arange(1,Nx-1):
        Uf[i] = r*(Uc[i-1] + Uc[i+1]) + 2*(1-r)*Uc[i] -Up[i]
    
    Up=copy(Uc)
    Uc=copy(Uf)
    
    if k%200==1:
        #print 'iteration', k
        clf()
        plot(Up)
        draw()
        pause(1e-17)
        sleep(0.1)
        
show()      
