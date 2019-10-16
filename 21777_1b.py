#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 14:31:43 2019

@author: alexfassone1
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants

#define all the constants to be used in all the equations
M = 1.99*10**30
G = scipy.constants.G
#put the constants into an array
cnsts = [M,G]

#define the order of the variables that will be used as if in an array
#Var = [x,A,y,B,t]

#define the differential equations using the array of constants and the array of variables
def dxdt(Cn,Var):
    return Var[1]

def dAdt(Cn,Var):
    return -1*Cn[1]*Cn[0]*(Var[0]/(Var[0]**2+Var[2]**2)**(3/2))

def dydt(Cn,Var):
    return Var[3]

def dBdt(Cn,Var):
    return -1*Cn[1]*Cn[0]*(Var[2]/(Var[0]**2+Var[2]**2)**(3/2))

#create an array of the functions 
drdt = [dxdt,dAdt,dydt,dBdt]

#create an array of the initial conditions (keeping the order the same as defined above)
#x0, dxdt0, y0, dydt0, t0
iniConds = [5.2*10**12,0,0,880,0]

#arguments(array of functions, array of initial conditions, array of constants, time to integrate to, number of points)
def Rk4(dRdt,Ics,Cn,t_f,N):
    
    n=len(dRdt) #find the number of equations 
    dt=(t_f-Ics[n])/N #calculate step size 
    dt_min=dt/(100) #use the intial stepsize to calculate a minimum acceptable stepsize
    
    #two 2D arrays to store the values for all the equations and time from 0 to t_f 
    Str1 = np.zeros((n+1, 10))
    Str2 = np.zeros((n+1, 10))
    
    #populate the first value of each of the equations with the initial conditions for the first array
    for i in range(0,n+1):
        Str1[i][0] = Ics[i]
    
    #define a new function to calculate the value of every functions at a time h later
    #arguments(integer position, stepsize, array of equations, number of equations)
    def step(i,h,Str,n):
        
        #create arrays to store all 4 k values for each of the equations 
        k1_= np.zeros(n)
        k2_= np.zeros(n)
        k3_= np.zeros(n)
        k4_= np.zeros(n)
        
        #create an array to store the input values for k1 function
        Var1 = np.zeros(n+1)
        
        #populate the array with stored values from the functions 
        for k in range(0,n+1):
            Var1[k] = Str[k][i-1]
        
        #calculate the k1 constants for each function
        for j in range(0,n):
            k1_[j]=h*dRdt[j](Cn,Var1)
            
        #repeat as above but for the k2 constants
        Var2 = np.zeros(n+1)
            
        for k in range(0,n):
            Var2[k] = Str[k][i-1] + k1_[k]/2
            
        Var2[n] = Str[n][i-1] + h/2
        
        for j in range(0,n):
            k2_[j]=h*dRdt[j](Cn,Var2)
        
        #repeat as above but for the k3 constants
        Var3 = np.zeros(n+1)
            
        for k in range(0,n):
            Var3[k] = Str[k][i-1] + k2_[k]/2
            
        Var3[n] = Str[n][i-1] + h/2
        
        for j in range(0,n):
            k3_[j]=h*dRdt[j](Cn,Var3)
            
        #repeat as above but for the k4 constants
        Var4 = np.zeros(n+1)
            
        for k in range(0,n):
            Var4[k] = Str[k][i-1] + k3_[k]
            
        Var4[n] = Str[n][i-1] + h
        
        for j in range(0,n):
            k4_[j]=h*dRdt[j](Cn,Var4)
        
        
        #use the k constants to calculate the next values in each of the functions 
        for j in range(0,n):
            Str[j][i] = Str[j][i-1] + 1/6*(k1_[j]+2*k2_[j]+2*k3_[j]+k4_[j])
        
        #step forward one stepseize in time
        Str[n][i] = Str[n][i-1]+h
        
        #return the updated equations array
        return Str
    
    c = 1 #set counter to 1
    time = dt #set time to one stepsize (first time after the initial conditions)
    tol = 500000 #set the error tolerance
    
    #create an empty array used to prevent function getting trapped at the minimum stepsize
    t_check = np.zeros(3) 
    
    #continue to iterate unitl the time has reached the final time set by the user
    while time < t_f:
        
        #update the check array so that it containes the 3 most recent times 
        t_check[2] = t_check[1]
        t_check[1] = t_check[0]
        t_check[0] = time
        
        #if the counter is up to the limit of the size storage array then increase the size of the array
        if c == Str1.shape[1] - 1:
            
            Str3 = np.zeros((n+1, Str1.shape[1] + 10)) 
            Str3[:,0:Str1.shape[1]] = Str1[:,:]
            Str1 = Str3.copy()
            
        
        #take two steps in time each of size dt 
        Str1 = step(c,dt,Str1,n) 
        Str1 = step(c+1,dt,Str1,n)
        
        #calculate the distance from the sun
        r1 = (Str1[0][c+1]**2 + Str1[2][c+1]**2)**0.5
        
        Str2 = Str1.copy() #copy over the data to a temporary array
        
        #take one step in time of size 2*dt
        Str2 = step(c,2*dt,Str2,n)    
        #print(Str2[0][c])
        
        #calculate the distance from the sun
        r2 = (Str2[0][c]**2 + Str2[2][c]**2)**0.5
        
        err = abs(r1 - r2)/30 #use the two distances to calculate the error 
        
        
        #if the error is bigger than the set tolerance level 
        if err > tol:
            
            #if the stepsize is at its minimum possible value and the time hasnt changed in 3 iterations
            #this shows that the function is trapped at a point
            if dt == dt_min and t_check[2] == t_check[0]:
                    
                c = c+1 #increment the counter by 1
                time = time + dt #increase the time by dt
            
            #if half dt is less than the minimum value then set dt to the minimum value otherwise half dt
            elif dt/2 < dt_min: 
                dt=dt_min
            else:
                dt = dt/2 
        
        else:
                
            c = c+1 #increment the counter by 1
            time = time + dt #increase the time by dt
            dt = dt*2 #set dt to double its value
            
            
    #use the calculated values to plot as desired               
    plt.scatter(Str1[0]/(1.5*10**11),Str1[2]/(1.5*10**11), c = Str1[n]/(3.154*10**7), cmap='Blues')
    
    clb = plt.colorbar()
    clb.set_label('time (years)', labelpad=-25, y=1.08, rotation=0)
    
    plt.grid()
    plt.xlabel('AU')
    plt.ylabel('AU')
    plt.plot(0,0, 'ro')
    plt.show()
    
Rk4(drdt,iniConds,cnsts,2.4*10**9,10000)