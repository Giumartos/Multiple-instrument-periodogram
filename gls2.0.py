#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:19:04 2020

@author: giumartos
"""

from sympy import Matrix, linsolve
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

#---------------------------------------------------------------------------------------------------------------------------------
def savetofile(x, y, outfile_name):
    
    '''
    Saves vectors of interested into a file.
    '''
    file = open(outfile_name,'w')
    for i in range(len(x)):
        file.write('{0:.5f} {1:.5f}\n'.format(x[i],y[i]))
    file.close()
    
#---------------------------------------------------------------------------------------------------------------------------------
def get_data (instruments: list, freq, filename):
    
    '''
    Extracts the times and rv data from a file and creates a dictionary containing 
    the data from each instrument, and returns a list with these dictionaries.
    '''
    
    file = open(filename,'r')
    lines = file.readlines()
    
    dicts = []
    for k in range(len(instruments)):
        instruments[k] = {'name': instruments[k], 'time':[], 'rv':[], 'rv_err':[]}
        for l in lines:
            l = l.split()
            if l[-1] == instruments[k]['name']:
                instruments[k]['time'].append(float(l[0]))
                instruments[k]['rv'].append(float(l[1]))
                instruments[k]['rv_err'].append(float(l[2]))
                instruments[k]['freq'] = np.array(freq)
        dicts.append(instruments[k])
    file.close()      

    return dicts

#---------------------------------------------------------------------------------------------------------------------------------

def gls (instruments: list, freq, filename):

    '''
    Calculates and plots the multiple-intrument periodogram.
    '''
    
    power = []
    data_dicts = get_data(instruments, freq, filename)


    # Calculating the parameters
    
    for f in freq:

        X2 = X2_0 = p = 0
        omega = 2 * np.pi * f
        W, Ĉs, Ŝs, Ŷ, w = [], [], [], [], []
        time, rvs, rvs_err = [],[],[]
        ĈĈ = ŜŜ = ĈŜ = ŶĈ = ŶŜ = 0
        values = {} #dictionary that will contain the values used for solve the system
        
        
        for k in range(len(instruments)): #for each instrument
            
            
            t = np.array(data_dicts[k]['time'])
            time.append(t)
            rv = np.array(data_dicts[k]['rv'])
            rvs.append(rv)
            rv_err = np.array(data_dicts[k]['rv_err'])
            rvs_err.append(rv_err)
           

            W_k = np.sum(1/(rv_err**2)) 
            W.append(W_k)
            w_k = 1/ (W_k * rv_err**2)
            w.append(w_k)
            Ŷ_k = np.sum(w_k * rv)
            Ŷ.append(Ŷ_k)
            
            Ĉ_k = np.sum(w_k * np.cos(omega * t))
            Ĉs.append(Ĉ_k)
            Ŝ_k = np.sum(w_k * np.sin(omega * t))
            Ŝs.append(Ŝ_k)
         
            ĈĈ += W_k * np.sum(w_k * (np.cos(omega * t))**2)
            ŜŜ += W_k * np.sum(w_k * (np.sin(omega * t))**2)
            ĈŜ += W_k * np.sum(w_k * np.cos(omega * t) * np.sin(omega * t))
            ŶĈ += W_k * np.sum(w_k * rv * np.cos(omega * t))
            ŶŜ += W_k * np.sum(w_k * rv * np.sin(omega * t))
            
            values['Ĉ'+str(k+1)], values['Ŝ'+str(k+1)],values['W'+str(k+1)], values['Ŷ'+str(k+1)] = Ĉs[k], Ŝs[k], W[k], Ŷ[k]
            values['W'+str(k+1)+'Ĉ'+str(k+1)], values['W'+str(k+1)+'Ŝ'+str(k+1)] = W[k]*Ĉs[k], W[k]*Ŝs[k]
        
        values.update({'ĈĈ':ĈĈ, 'ŜŜ':ŜŜ, 'ĈŜ': ĈŜ, 'ŶĈ':ŶĈ, 'ŶŜ':ŶŜ})

        
        # Constructing and solving the system of equations for the X**2 considering the given instruments.
        
        Is = [] # lines of the system
        matrix = [] # coefficients matrix
        I = [values['ĈĈ'], values['ĈŜ']]
        II = [values['ĈŜ'], values['ŜŜ']]
        res = [values['ŶĈ'], values['ŶŜ']]
        variables = ['a','b']
        
        
        # Creating the system, considering all the instruments
        
        for i in range(1,len(instruments)+1):
            equation = [] # equation that will be added to the system 
            variables.append('c'+str(i))
            I.append(values['W'+str(i)+'Ĉ'+str(i)])
            II.append(values['W'+str(i)+'Ŝ'+str(i)])
            equation.extend((values['Ĉ'+str(i)], values['Ŝ'+str(i)]))
            
            for pos in range(len(equation), len(instruments)+2):
                equation.append(0)
            equation[i+1] = 1   
            Is.append(equation)
            res.append(values['Ŷ'+str(i)])  
        matrix.extend((I, II))
        
        for eq in Is:
            matrix.append(eq)
 

    
        # Solving the system
    
        results = linsolve((Matrix(matrix), Matrix(res)), variables)

        coefs = [] # coefficients
        for v in range(len(variables)):
            coefs.append(results.args[0][v])
        coefs = np.array(coefs).astype(float)
                
        offsets = coefs[2:]     
        a, b = coefs[0], coefs[1]


        # Calculating the X**2 and the power
 
        for k in range (len(instruments)):
            cos = np.cos(omega*time[k])
            sin = np.sin(omega*time[k])
            rv = rvs[k]

            X2 += W[k] * np.sum(w[k] * (rv - a*cos - b*sin  - offsets[k])**2)
            X2_0 += W[k] * np.sum(w[k] * (rv - Ŷ[k])**2)
     
        p = 1 - (X2/X2_0)
        power.append(p)
        
    power = np.array(power)
    savetofile(1/freq, power, outfile )
    

#---------------------------------------------------------------------------------------------------------------
    # plotting the rv points and the periodogram 
 
    fig = plt.figure(figsize=(12,8))
    
    fig.subplots_adjust(hspace=0.4)
    
    ax = fig.add_subplot(2,1,1)

    for i in range(len(instruments)):
        ax.errorbar(time[i], rvs[i], rvs_err[i], fmt='.', label=instruments[i]['name'])
    

    ax.set_xlabel('Time (days)') 
    ax.set_ylabel('RVs')
    ax.legend()
    
    ax1 = fig.add_subplot(2,1,2)
    ax1.set_xlabel('Period (days)') # Or frequency depending on the plot
    ax1.set_ylabel('Power')
    ax1.minorticks_on()
    ax1.grid('minor', linestyle=':', axis ='both')
    ax1.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labeltop=True)
    # Limits of the plot
    ax1.set_xlim(0.1, 6)
    ax1.set_ylim(0.,1.)
    ax1.set_xticks(np.arange(0, 6, 1.0))
    
    # Plot 1/freq to have a period x-axis or freq to have frequency x-axis
    ax1.plot(1/freq, power, '-', color='k', linewidth = 0.7)
    plt.savefig(imagename)
    plt.show()

#---------------------------------------------------------------------------------------------------------------
instruments = [] # List with the name of the instruments
filename = # Name of the rv file
freq = np.linspace(1/15, 10, 1e4) # Frequency to evaluate the data
imagename = # Name to save the periodogram image plot
outfile = # Name of the file to save the frequency and the power 

gls(instruments, freq, filename)
