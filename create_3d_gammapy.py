# !pip install astropy scipy numpy pandas matplotlib

import astropy.wcs as pywcs
import astropy.io.fits as pyfits
from astropy.io import fits
from astropy.io.fits import getheader
import numpy
import numpy as np
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

#coordinates of the image center for Geminga [1711.06223]
ra0 = 98.48
dec0 = 17.77
print("start")
binn = 0.05 #deg/pix
size = 10 #deg -- image size
#image size and central pixels
npix=int(size/binn)
center_x=float(npix)/2 #center x
center_y=float(npix)/2 #center y

Norm = 2
    
# Read csv file
data = pd.read_csv("dphi_theta0.0.csv")
Energies=data["E [GeV]"]
Energies=Energies*1000 # E in MeV
tmpFlux=Norm*data["dphi [1/TeV/cm^2/s]"]

#sys.exit()

# Create a text file with energies
outfile=open('Energies_Elena.txt','w')
for ii in range(1,len(Energies)+1):
    outfile.write(str(Energies[ii-1])+'\n')

outfile.close()

# Get CSV files list from a folder
csv_files = sorted(glob.glob("dphi*.csv"))

# Read each CSV file into DataFrame. This creates a list of dataframes
dfs = list()
for filename in csv_files:    
    df = pd.read_csv(filename)    
    dfs.append(df)
    
dfs_2d = np.array([np.array(r["dphi [1/TeV/cm^2/s]"]) for r in dfs])    
    
def psf_value_by_index(theta_index, energy_index):
    #print(theta_index.shape)
    m = theta_index < dfs_2d.shape[0]
    r = np.zeros_like(theta_index, dtype=float)
    r[m] = dfs_2d[theta_index[m], energy_index]    
    return r
    

fout = 'allenergies_TeV.fits'#'E={:.2f}GeV.fits'.format(Energies[k])    
try:
    os.remove(fout)
except:
    pass


cube = numpy.zeros((len(Energies), npix, npix))


def f(x):
    return 2*pi*sin(x)

k=0
# k = 10

while k < len(Energies):
# while k < 11:

    #defining WCS header
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [float(npix)/2, float(npix)/2]
    wcs.wcs.cdelt = numpy.array([-binn, binn])
    wcs.wcs.crval = [ra0, dec0]
    wcs.wcs.ctype = ["RA---CAR", "DEC--CAR"]

    #Fluxes = [0, Norm*dfs[0]["dphi [1/TeV/cm^2/s]"][k]/1e6]
    
    #generate image data
    dat=numpy.zeros((npix,npix)) #image data will be stored here ; just zeros at the moment

    X, Y = np.mgrid[:npix, :npix]
    R = ((X - center_x)**2 + (Y - center_y)**2)**0.5
    #Dat = psf_value_by_index(np.round(R*binn*10).astype(int), k)

    #go pixel-by-pixel and fill the image data with the required value
    for xx in range(npix):
        for yy in range(npix):
            r=((xx-center_x+0.5)**2+(yy-center_y+0.5)**2)**0.5 #distance to the image center in pixels
            if r*binn<=9.9:
                #area = 3.14 * ( ((int(r*binn*10)+1)/10)**2-(int(r*binn*10)/10)**2 ) * 3.046*10**-4 # sr
                # VS: this "f" function is analytical, quad is taking resource unnecessarily
                # res = integrate.quad(f, pi/180 * int(r*binn*10)/10, pi/180 * (int(r*binn*10)+1)/10)[0]
                #res = integrate.quad(f, pi/180 * r*binn, pi/180 * (r*binn+0.1))[0]
                #print(res)
                
                flux_sky_density_rad2 = dfs[int(r*binn*10)]["dphi [1/TeV/cm^2/s]"][k]
                
                # dat is flux integrated over pixel area, i.e. just multiplied by it's area
                pixel_area_deg2 = binn**2
                pixel_area_rad2 = pixel_area_deg2 * (np.pi/180)**2
                dat[xx,yy] = Norm * flux_sky_density_rad2 * pixel_area_rad2 
                #dat[xx,yy] = dat[xx,yy] / 10**6 # 10**6 gives the correct flux inside 10 deg
                
                # Narea = int(3.14 * ( (((int(r*binn*10)+1)/10)/binn)**2 - ((int(r*binn*10)/10)/binn)**2 )) # pixels
                # dat[xx,yy] = Norm*res*(dfs[int(r*binn*10)+1]["dphi [1/TeV/cm^2/s]"][k]+dfs[int(r*binn*10)]["dphi [1/TeV/cm^2/s]"][k])/2/1e6/Narea #fill the data
            else:
                dat[xx,yy] = 0
                
                     
    dfs_theta_deg = np.linspace(0, 10-0.1, 100)
    dfs_theta_rad = dfs_theta_deg*(1./180*np.pi)
    annulus_area_rad2 = dfs_theta_rad*np.pi*(0.1/180*np.pi)*2    
    
    # if Energies[k]/1000 > 133 and Energies[k]/1000 < 134:
    #     print("input interpretted:", np.sum([df["dphi [1/TeV/cm^2/s]"][k]*area for df, area in zip(dfs, annulus_area_rad2)]))
    #     print("sum over pixels", dat.sum())


    #writing down the image
    header = wcs.to_header()
    cube[k, :, :] = dat / (7.615 * 1e-7) #(in 1 / cm2 TeV s sr; 7.615 * 1e-7 = 0.05*0.05 degrees^2 in sr)

    print('Done with E =', Energies[k], 'MeV')

    k = k + 1
    
cube = numpy.array(cube)
hdu = pyfits.PrimaryHDU(data=cube,header=header)
hdu.writeto(fout)

# hdu2 = pyfits.ImageHDU(data=Energies,header=header)
# hdulist = pyfits.HDUList([hdu, hdu2])
# hdulist.writeto(fout)

# sys.exit()   

fname = 'allenergies_TeV.fits'
spec_name = 'spec.txt'
# fluxes = []


# #energies=numpy.loadtxt('Energies_Elena.txt',unpack=True)

# with pyfits.open(fname,'update') as f:
#     data = f[0].data
#     for ii in range(data.shape[0]):
#         data_slice = data[ii]
#         tot_sum = numpy.ma.sum(data_slice) * pixel_area_rad2
#         fluxes.append(tot_sum)
#         if( tot_sum > 0 ):
#             f[0].data[ii] = data_slice/tot_sum / pixel_area_rad2
#             if(ii>10):
#                 f[0].data[ii] = 0


# with open(spec_name,'w') as f:
#     for ii in range(len(Energies)):
#         ss = str(Energies[ii]) +' '+str(fluxes[ii])
#         f.write(ss+'\n')


with open(spec_name,'w') as f:
    for ii in range(len(Energies)):
        ss = str(Energies[ii]) +' '+str(1)
        f.write(ss+'\n')


# Open the original FITS file
with fits.open('allenergies_TeV.fits', mode='update') as hdul:

    # Create a new table HDU
    col1 = fits.Column(name='Energy', format='1D', array=Energies) #np.array([1, 2, 3]))
    table_hdu = fits.BinTableHDU.from_columns([col1],name='ENERGIES')



    # Add the table HDU to the original file
    hdul.append(table_hdu)



    # Save the changes
    hdul.flush()

