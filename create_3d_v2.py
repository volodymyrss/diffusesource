import astropy.wcs as pywcs
import astropy.io.fits as pyfits
import numpy
from numpy import sqrt, sin, cos, pi
import scipy as sp
import scipy.integrate as integrate
import os
import csv
import pandas as pd
import glob

import sys

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

#coordinates of the image center for Geminga (FK4)
ra0 = 97.746636
dec0 = +17.809522
print("start")
binn = 0.05 #deg/pix
size = 10 #deg -- image size
#image size and central pixels
npix=int(size/binn)
center_x=float(npix)/2 #center x
center_y=float(npix)/2 #center y
    
# Read csv file
data = pd.read_csv("dphi_theta0.0.csv")
Energies=data["E [GeV]"]
Energies=Energies*1000 # E in MeV
tmpFlux=data["dphi [1/TeV/cm^2/s]"]
tmpFlux=tmpFlux/1e6 # Flux in 1/MeV/cm^2/s

#sys.exit()

# Creat a text file with energies
outfile=open('Energies_Elena.txt','w')
for ii in range(1,len(Energies)+1):
    outfile.write(str(Energies[ii-1])+'\n')

outfile.close()

# Get CSV files list from a folder
csv_files = sorted(glob.glob("*.csv"))

# Read each CSV file into DataFrame. This creates a list of dataframes
dfs = list()
for filename in csv_files:    
    df = pd.read_csv(filename)    
    dfs.append(df)


fout = 'allenergies.fits'#'E={:.2f}GeV.fits'.format(Energies[k])    
try:
    os.remove(fout)
except:
    pass


cube = numpy.zeros((len(Energies), npix, npix))


def f(x):
    return 2*pi*sin(x)
k=0
while k < len(Energies):

    #defining WCS header
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [float(npix)/2, float(npix)/2]
    wcs.wcs.cdelt = numpy.array([-binn, binn])
    wcs.wcs.crval = [ra0, dec0]
    wcs.wcs.ctype = ["RA---CAR", "DEC--CAR"]

    Fluxes = [0, dfs[0]["dphi [1/TeV/cm^2/s]"][k]/1e6]
    
    #generate image data
    dat=numpy.zeros((npix,npix)) #image data will be stored here ; just zeros at the moment


    #go pixel-by-pixel and fill the image data with the required value
    for xx in range(npix):
        for yy in range(npix):
            r=(xx-center_x+0.5)**2+(yy-center_y+0.5)**2 #distance to the image center in pixels
            r=r**0.5
            if r*binn<9.9:
                #area = 3.14 * ( ((int(r*binn*10)+1)/10)**2-(int(r*binn*10)/10)**2 ) * 3.046*10**-4 # sr
                res = integrate.quad(f, pi/180 * int(r*binn*10), pi/180 * (int(r*binn*10)+1))[0]
                Narea = int(3.14 * ( (((int(r*binn*10)+1)/10)/binn)**2 - ((int(r*binn*10)/10)/binn)**2 )) # pixels
                dat[xx,yy]= res*(dfs[int(r*binn*10)+1]["dphi [1/TeV/cm^2/s]"][k]+dfs[int(r*binn*10)]["dphi [1/TeV/cm^2/s]"][k])/2/1e6/Narea #fill the data
            else:
                dat[xx,yy]=0


    #writing down the image
    header = wcs.to_header()
    cube[k, :, :] = dat

    print('Done with E =', Energies[k]/1000, 'GeV')

    k = k+1

cube = numpy.array(cube)
hdu = pyfits.PrimaryHDU(data=cube,header=header)
hdu.writeto(fout)

# hdu2 = pyfits.ImageHDU(data=Energies,header=header)
# hdulist = pyfits.HDUList([hdu, hdu2])
# hdulist.writeto(fout)

# sys.exit()   

fname = 'allenergies.fits'
spec_name = 'spec.txt'
fluxes = []

#energies=numpy.loadtxt('Energies_Elena.txt',unpack=True)

with pyfits.open(fname,'update') as f:
    data = f[0].data
    for ii in range(data.shape[0]):
        data_slice = data[ii]
        tot_sum = numpy.ma.sum(data_slice)
        fluxes.append(tot_sum)
        if( tot_sum > 0 ):
            f[0].data[ii] = data_slice/tot_sum



with open(spec_name,'w') as f:
    for ii in range(len(Energies)):
        ss = str(Energies[ii]) +' '+str(fluxes[ii])
        f.write(ss+'\n')

