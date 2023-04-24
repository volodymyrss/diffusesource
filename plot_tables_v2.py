import healpy as hp
import numpy as np
from numpy import sqrt, sin, cos, pi
import scipy as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from mhealpy import HealpixMap
from astropy.io import fits
import csv
import glob
from scipy.interpolate import interp1d

import pandas as pd
import os
import math
import sys


files = sorted(glob.glob("*.csv"))

tables = []

for i in files:
    test = pd.read_csv('/Users/sokolenko/Dropbox/KICP/Dan/Elena_tab_2/alpha_15/'+i)
    tables.append(test)


#Read the first table
main_table = pd.read_csv('dphi_theta0.0.csv')

#Adding all other tables to the first one
iterator = 0.1
for i in tables:
    main_table[str(iterator)] = i['dphi [1/TeV/cm^2/s]']
    iterator = round(iterator+0.1,1)
  
main_table.rename(columns={'dphi [1/TeV/cm^2/s]' : '0.0'}, inplace = True)

# Tables to the 2D numpy matrix: 1st [] - row; 2nd [] column
numpy_table = main_table.to_numpy()

# This example gives 1st row and all the columns but the 1st one (the first one is the energy), 
# i.e. all values for the flux 
#print(numpy_table[energy_bin][2:])


Radius=10 #degrees


def f(x):
    return 2*pi*sin(x)


energy_bin = 0
k=0
flux=0
while k < Radius:
    res = integrate.quad(f, pi/180 * (k), pi/180 * (k+0.1))[0]
    if k<Radius-1:
        flux = flux + res * (numpy_table[energy_bin][2+int(10*k)]+numpy_table[energy_bin][2+int(10*(k+1))])/2
    else:
        flux = flux + res * numpy_table[energy_bin][2+int(10*k)]
    k = k+0.1
#flux=flux/(3.14*Radius**2) 


Tabs_E_F=[numpy_table[energy_bin][0],flux]


for i in range(1,len(numpy_table)):
    energy_bin = i
    k=0
    flux=0
    while k < Radius:
        res = integrate.quad(f, pi/180 * (k), pi/180 * (k+0.1))[0]
        if k<Radius-1:
            flux = flux + res * (numpy_table[energy_bin][2+int(10*k)]+numpy_table[energy_bin][2+int(10*(k+1))])/2
        else:
            flux = flux + res * numpy_table[energy_bin][2+int(10*k)]
        k = k+0.1
    #flux=flux/(3.14*Radius**2) 
    Tabs_E_F = np.vstack([Tabs_E_F, [numpy_table[energy_bin][0],flux]])



# print('Radius',Radius)
# print('E,[GeV]    F[1/TeV/cm^2/s]')
# for i in range(0,len(Tabs_E_F[::,0])):
#     print(str(Tabs_E_F[i,0]), " ", str(Tabs_E_F[i,1]))

# sys.exit()

#Load data from Moskalenko et al 2019

x = []
y = []
with open('/Users/sokolenko/Dropbox/KICP/Dan/Elena_tab_2/Moskalenko_IC_iso_gamma1_2p0_rz_rt_Inf.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

x1 = []
y1 = []
with open('/Users/sokolenko/Dropbox/KICP/Dan/Elena_tab_2/Moskalenko_IC_iso_gamma1_2p0.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[1]))  

x_HAWC_u = []
y_HAWC_u = []
with open('/Users/sokolenko/Dropbox/KICP/Dan/Elena_tab_2/HAWC_upper.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_HAWC_u.append(float(row[0]))
        y_HAWC_u.append(float(row[1])) 

x_HAWC_l = []
y_HAWC_l = []
with open('/Users/sokolenko/Dropbox/KICP/Dan/Elena_tab_2/HAWC_lower.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_HAWC_l.append(float(row[0]))
        y_HAWC_l.append(float(row[1]))

x_Magic = []
y_Magic = []
with open('/Users/sokolenko/Dropbox/KICP/Dan/Elena_tab_2/Magic.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_Magic.append(float(row[0]))
        y_Magic.append(float(row[1]))
                
Int_HAWC_u=interp1d(x_HAWC_u, y_HAWC_u,fill_value="extrapolate")
Int_HAWC_l=interp1d(x_HAWC_l, y_HAWC_l,fill_value="extrapolate")


plt.figure()
plt.title("$r=10^\\circ$")
plt.plot(x,y, label='Paper, $\\gamma_1=2.0, r_z, r_t =$ Infinity',color = 'magenta')
plt.plot(x1,y1, label='Paper, $\\gamma_1=2.0, r_z, r_t =(30,50)$ ',color = 'gray')
plt.plot(Tabs_E_F[::,0], 10**6 * Tabs_E_F[::,1]*Tabs_E_F[::,0]**2, label='Elena', linestyle="--", color = 'b')
plt.fill_between(x_HAWC_l, Int_HAWC_u(x_HAWC_l),Int_HAWC_l(x_HAWC_l), color='green',alpha=.4,label='HAWC')
plt.scatter(x_Magic,y_Magic,color = 'purple',label='Magic')
plt.xscale("log")
plt.yscale("log")
plt.xlim(10, 4e4)
plt.ylim(0.1, 400)
plt.legend(frameon=False, loc='lower left')
plt.xlabel("$E, [GeV]$")
plt.ylabel("Flux, [eV cm$^{-2} \cdot s^{-1}$]")

plt.show()