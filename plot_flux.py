# !pip install pandas
import healpy as hp
import numpy as np
from numpy import sqrt, sin, cos, pi
import scipy as sp
import scipy.integrate as integrate
import matplotlib.pyplot as plt
# from mhealpy import HealpixMap
from astropy.io import fits
import csv
import glob
from scipy.interpolate import interp1d

import pandas as pd
import os
import math
import sys


files = sorted(glob.glob("dphi_*.csv"))

Norm = 2
name = "alpha15_model_Flux.dat"

tables = []

for i in files:
    test = pd.read_csv(i)
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


Radius=5 #degrees


def f(x):
    return 2*pi*sin(x)


energy_bin = 0
k = 0
flux = 0



Tabs_E_F = None


for i in range(len(numpy_table)):
    energy_bin = i
    radius_deg = 0
    flux=0
    while radius_deg < Radius:
        # res = integrate.quad(f, pi/180 * (k), pi/180 * (k+0.1))[0]
        # precise form for area, better than an approximation
        ring_area = 2*pi * (cos(pi/180 * (radius_deg)) - cos(pi/180 * (radius_deg+0.1)))
        
        if radius_deg<Radius-1:
            flux += Norm * ring_area * (numpy_table[energy_bin][2+int(10*radius_deg)])
            # flux += Norm * ring_area * (numpy_table[energy_bin][2+int(10*radius_deg)]+numpy_table[energy_bin][2+int(10*(radius_deg)+1)])/2
            # flux = flux + Norm * ring_area * (numpy_table[energy_bin][2+int(10*radius_deg)]+numpy_table[energy_bin][2+int(10*(radius_deg+1))])/2
            # flux = flux + Norm * ring_area * (numpy_table[energy_bin][2+int(10*radius_deg)])
        else:
            flux += Norm * ring_area * numpy_table[energy_bin][2+int(10*radius_deg)]
        
        radius_deg = radius_deg + 0.1
        
    #flux=flux/(3.14*Radius**2) 
    
    row = [numpy_table[energy_bin][0],flux]
    
    if Tabs_E_F is not None:
        Tabs_E_F = np.vstack([Tabs_E_F, row])
    else:
        Tabs_E_F = row

    if numpy_table[energy_bin][0] > 133 and numpy_table[energy_bin][0] < 134:
        print("energy", numpy_table[energy_bin][0], flux)
        

# print('Radius',Radius)
# print('E,[GeV]    F[1/TeV/cm^2/s]')
# for i in range(0,len(Tabs_E_F[::,0])):
#     print(str(Tabs_E_F[i,0]), " ", str(Tabs_E_F[i,1]))

# sys.exit()

#Load data from Moskalenko et al 2019

# x = []
# y = []
# with open('Moskalenko_IC_iso_gamma1_2p0_rz_rt_Inf.csv','r') as csvfile:
#     plots = csv.reader(csvfile, delimiter=',')
#     for row in plots:
#         x.append(float(row[0]))
#         y.append(float(row[1]))

# x1 = []
# y1 = []
# with open('Moskalenko_IC_iso_gamma1_2p0.csv','r') as csvfile:
#     plots = csv.reader(csvfile, delimiter=',')
#     for row in plots:
#         x1.append(float(row[0]))
#         y1.append(float(row[1]))  

x_HAWC_u = []
y_HAWC_u = []
with open('HAWC_upper.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_HAWC_u.append(float(row[0]))
        y_HAWC_u.append(float(row[1])) 

x_HAWC_l = []
y_HAWC_l = []
with open('HAWC_lower.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_HAWC_l.append(float(row[0]))
        y_HAWC_l.append(float(row[1]))

x_Magic = []
y_Magic = []
with open('Magic.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x_Magic.append(float(row[0]))
        y_Magic.append(float(row[1]))
                
Int_HAWC_u=interp1d(x_HAWC_u, y_HAWC_u,fill_value="extrapolate")
Int_HAWC_l=interp1d(x_HAWC_l, y_HAWC_l,fill_value="extrapolate")


plt.figure()
plt.title("$r=5^\\circ$")
#plt.plot(x,y, label='Paper, $\\gamma_1=2.0, r_z, r_t =$ Infinity',color = 'magenta')
#plt.plot(x1,y1, label='Paper, $\\gamma_1=2.0, r_z, r_t =(30,50)$ ',color = 'gray')

plt.plot(Tabs_E_F[::,0], 10**6 * Tabs_E_F[::,1]*Tabs_E_F[::,0]**2, label='Elena', linestyle="--", color = 'b')

with open(name,'w') as outfile:
    print('E,[GeV]    F[1 eV/cm^2/s]')
    for i in range(0,len(Tabs_E_F[::,0])):
        print(str(Tabs_E_F[i,0]), " ", str(10**6 * Tabs_E_F[i,1]*Tabs_E_F[i,0]**2))
        outfile.write(str(Tabs_E_F[i,0]))
        outfile.write('    ')
        outfile.write(str(10**6 * Tabs_E_F[i,1]*Tabs_E_F[i,0]**2))
        outfile.write('\n')
outfile.close()

##

# flux_from_2map = fits.open("allenergies.fits")[0].data.sum(1).sum(1)
# plt.plot(Tabs_E_F[::,0], 
#          10**6 * flux_from_2map*Tabs_E_F[::,0]**2, lw=5,alpha=0.5,
#          label="from 3d map")
##

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

# Spectrum of IC emission averaged over a 10° wide region