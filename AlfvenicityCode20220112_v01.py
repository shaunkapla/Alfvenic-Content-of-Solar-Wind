#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:19:43 2022

"""

import matplotlib.pyplot as plt
import numpy as np

## import scripts from pytplot

from pytplot import cdf_to_tplot
from pytplot import get_data

## import scripts from PySpedas

import pyspedas as ps
from pyspedas import time_double
from pyspedas import time_string

import scipy
from scipy import signal


# --------------------------------------------------------------------------------

"""Functions"""

def magnitude(x,y,z):
    """This function accepts three one-dimensional lists as inputs.
    These inputs correspond to the directional components of a single vector

    The function computes and returns the magnitude of the vector
    as a one-dimensional list.

    Parameters
    ----------
    x : a 1D list corresponding to the first vector component (RTN, XYZ, etc.).

    y : a 1D list corresponding to the second vector component.

    z : a 1D list corresponding to the third vector component.

    Returns
    -------
    The magnitude of the vector at all indices as a 1D list.
    """

    return np.sqrt(x**2+y**2+z**2)

def mag3d(a):
    """This function accepts one two-dimensional (Nx3) array as an input.
    It then computes and returns the magnitude of the array at each row.

    Parameters
    ----------
    a : a 2D array (Nx3: N rows, 3 columns).

    Returns
    -------
    The magnitude of the 2D array at each N-index.
    """

    return(np.sqrt(a[:,0]**2 + a[:,1]**2 + a[:,2]**2))

def dot_arrays(times, a,b):
    """This function accepts two two-dimmendional arrays (Nx3) and a list of
    times each index corresponds to. It then computes product of the two arrays
    at each time index, and stores it in a separate array, dotArr.

    This function was not used in the final algorithm.

    Parameters
    ----------
    times : a 1D list of times, len(times) == len(a) == len(b).

    a : an array (Nx3) representing a vector, len(a) == len(b) == len(times).

    b : an array (Nx3) representing a vector, len(b) == len(a) == len(times).

    Returns
    -------
    dotArr : a 1D list, the computed product of the two arrays at each index.
    """

    dotArr = []
    for i in range(len(times)):
        dotArr.append(a[i]*b[i])
    return dotArr

def nanFunky(a, b, c, epoch):
    """
    This function accepts three 1D lists, corresponding to each component of a
    vector, and a list of times each index corresponds to. It then indexes each
    invalid 'nan' value, and stores all valid (non-nan) values in a separate array.

    It then returns a new array populated by each unique value for every index.
    This algorithm functionally removes all 'nan' values from downloaded PSP data,
    and recombines the valid indices into a useable array for further computation
    and analysis.

    Parameters
    ----------


    """

    Arrnan=[]
    aNan=np.isnan(a)
    for i in range(len(aNan)):
        if aNan[i]==False:
            Arrnan.append(i)

    bNan=np.isnan(b)
    for i in range(len(bNan)):
        if bNan[i]==False:
            Arrnan.append(i)

    cNan=np.isnan(c)
    for i in range(len(cNan)):
        if cNan[i]==False:
            Arrnan.append(i)

    nan1=np.unique(Arrnan)

    return nan1

def dotproduct(a,b):
    return a[:,0]*b[:,0]+ a[:,1]*b[:,1]+a[:,2]*b[:,2]
#----------------------------------------------------------------------------------

path = '/Users/shaunkapla/Downloads/Scientific Notebook/projects/SolarProbe/Data/'

"""B-FIELD DATA"""

file = 'psp_fld_l2_mag_RTN_4_Sa_per_Cyc_20210916_v02.cdf'

# Open the file

tvars = cdf_to_tplot(path+file, get_support_data=1)

# View contents

print('')
for ind in range(len(tvars)):
    print(ind, '', tvars[ind])
print('\n')

atemp = get_data('psp_fld_l2_mag_RTN_4_Sa_per_Cyc')

### B-Field Data RTN ### ------------------------------------------------------
B_epoch = atemp[0]
B_rtn = atemp[1] # radial, tangential, and normal B field data

B_r = B_rtn[:,0] # isolate the radial B field data
B_t = B_rtn[:,1] # isolate the tangential B-field data
B_n = B_rtn[:,2] # isolate the normal B-field data

### Redifining B-Field RTN ### ------------------------------------------------

NAN1=nanFunky(B_r, B_t, B_n, B_epoch)

B_r= B_r[NAN1]
B_t=B_t[NAN1]
B_n=B_n[NAN1]
B_epoch=B_epoch[NAN1]

B_rtn_mag = magnitude(B_r,B_t,B_n) # compute the total RTN magnitude of B-field data
#print(B_epoch[0])
T1=time_string(B_epoch[0])
T1=T1[0:19]
T_end=time_string(B_epoch[-1])
T_end=T_end[0:19]
T_0 = time_double(T1) # start time June 16th, 2020
time_local = (B_epoch - T_0) / 3600 # in hours
T0_string = time_string(T_0)[0:10]

# plot the raw RTN B-field data

plt.style.use('default')
fig, (ax1, ax2, ax3) = plt.subplots(3)

ax1.plot(time_local, B_rtn_mag)
ax1.set_xlim(np.min(time_local), np.max(time_local))
ax1.set_title('Solar Wind RTN B-Field Magnitude')
ax1.set_xlabel('Hours after ' + T0_string)
ax1.set_ylabel('B-Field Magnitude (nT)')
ax1.set_xticklabels(time_local, fontsize = 1)
ax1.set_yticklabels(B_rtn_mag, fontsize = 10)

"""VELOCITY DATA"""

### Velocity Data RTN ### -----------------------------------------------------

file1='psp_swp_spc_l3i_20210916_v02.cdf'

# Open the file

tvars = cdf_to_tplot(path+file1, get_support_data=1)

# View contents

print('')
for ind in range(len(tvars)):
    print(ind, '', tvars[ind])
print('\n')

### Velocity Data RTN ### -----------------------------------------------------

atemp1        = get_data("vp_moment_RTN")
epoch_vel_rtn = atemp1[0]
data_vel_rtn  = atemp1[1]

rvel=data_vel_rtn[:,0] # isolate the radial component of velocity
tvel=data_vel_rtn[:,1] # isolate the tangential component of the velocity
nvel=data_vel_rtn[:,2] # isolate the normal component of the velocity

###Redifining RTN Velocity ###-------------------------------------------------
NAN3=nanFunky(rvel, tvel, nvel, epoch_vel_rtn)

rvel = rvel[NAN3]
tvel = tvel[NAN3]
nvel = nvel[NAN3]
epoch_vel_rtn = epoch_vel_rtn[NAN3]

#------------------------------------------------------------------------------

time_rtn_vel_local = (epoch_vel_rtn - T_0) /3600
T = time_string(T_0)

vel_rtn_mag =magnitude(rvel,tvel,nvel) # compute the magnitude of the full RTN velocity vector

### Plot the RTN velocity data to stack plot

ax2.plot(time_rtn_vel_local, vel_rtn_mag, label = 'Calculated Mean')
ax2.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax2.set_title('Solar Wind RTN Velocity Magnitude (km/s)')
ax2.set_xlabel('Hours after ' + T)
ax2.set_ylabel('Veloctiy (km/s)')

plt.tight_layout()

### Spacecraft Position Data ### ----------------------------------------------------------

atemp2 = get_data("sc_pos_HCI")
epoch_sc_pos = atemp2[0]
sc_position = atemp2[1]

# Use the following line only for SPC
#sc_pos = mag3d(sc_position)

s_rad = 1.4374e-6 # km to solar radii conversion factor

time_sc_pos_local = (epoch_sc_pos - T_0) /3600
T1 = time_string(T_0)

# Use the following line only for SPI
sc_pos = sc_position[0:len(epoch_sc_pos)]

ax3.plot(time_sc_pos_local, sc_pos * s_rad, label = 'Spacecraft Position from Sun')
ax3.set_title('Position of PSP Relative to Sun')
ax3.set_xlabel('Hours after ' + T)
ax3.set_ylabel('Position from Sun (R${_s}$)')

plt.tight_layout()

plt.show()

av_dist = (np.mean(sc_pos * s_rad))

### --------------------------------------------------------------------------------------

### INTERPOLATION ### --------------------------------------------------------
"""
This portion of the code interpolates the B-field data onto the velocity data,
separated into the individual RTN components. Each returned array has 386
populated indices.
"""

B_rad_interp = np.interp(time_rtn_vel_local, time_local, B_r)
r_samples = len(B_rad_interp)

B_tan_interp = np.interp(time_rtn_vel_local, time_local, B_t)
t_samples = len(B_tan_interp)

B_nor_interp = np.interp(time_rtn_vel_local, time_local, B_n)
n_samples = len(B_nor_interp)
# interpolate the RTN B-field data onto the RTN velocity data

print("Interpolated Radial B-Field Data Samples: ",r_samples)
print("Interpolated Tangential B-Field Data Samples: ",t_samples)
print("Interpolated Normal B-Field Data Samples: ",n_samples)
print("RTN Velocity Samples",len(rvel))

### Plot Separated B-Field Data (R, T, N)

plt.plot(time_rtn_vel_local, B_rad_interp, label = 'Radial')
plt.plot(time_rtn_vel_local, B_tan_interp, label = 'Tangential')
plt.plot(time_rtn_vel_local, B_nor_interp, label = 'Normal')
plt.xlabel('Hours after ' + T, size = 13)
plt.ylabel('|B| (nT)', size = 13)
plt.title('Interpolated RTN B-Field Components', size = 13)
plt.legend(bbox_to_anchor=(1.0, 1.0))
#plt.savefig("20210428_sepB.png", dpi = 250)
plt.show()

### Plot Separated Velocity Data (R, T, N)

plt.plot(time_rtn_vel_local, rvel, label = 'Radial')
plt.plot(time_rtn_vel_local, tvel, label = 'Tangential')
plt.plot(time_rtn_vel_local, nvel, label = 'Normal')
plt.xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
plt.title('Separated RTN Solar Wind Velocity Components', size = 13)
plt.xlabel('Hours after ' + T, size = 13)
plt.ylabel('|V| (km/s)', size = 13)
plt.legend(bbox_to_anchor=(1.0, 1.0))
#plt.savefig("20210428_sepV", dpi=250)
plt.show()


### Interpolated RTN Plots ### ------------------------------------------------------
# Plot the interpolated data as stack plots

fig, (ax4, ax5) = plt.subplots(2)

ax4.plot(time_rtn_vel_local, B_rad_interp)
ax4.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax4.set_title('Interpolated Radial B-Field Magnitude (n='+str(r_samples)+')', size = 14)
ax4.set_xlabel('Hours after ' + T, size = 10)
ax4.set_ylabel('|Br| (nT)', size = 14)

ax5.plot(time_rtn_vel_local, rvel)
ax5.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax5.set_title('Radial-Component Solar Wind Velocity', size = 14)
ax5.set_xlabel('Hours after '+ T, size = 10)
ax5.set_ylabel('Vr (km/s)', size = 14)

plt.tight_layout()
#plt.savefig("20210428_Rad_B.png", dpi = 250)

plt.show()

fig, (ax6, ax7) = plt.subplots(2)

ax6.plot(time_rtn_vel_local, B_tan_interp)
ax6.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax6.set_title('Interpolated Tangential B-Field Magnitude (n='+str(t_samples)+')', size = 14)
ax6.set_xlabel('Hours after ' + T, size = 10)
ax6.set_ylabel('|Br| (nT)', size = 14)

ax7.plot(time_rtn_vel_local, tvel)
ax7.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax7.set_title('Tangent-Component Solar Wind Velocity', size = 14)
ax7.set_xlabel('Hours after '+ T, size = 10)
ax7.set_ylabel('Vr (km/s)', size = 14)

plt.tight_layout()
#plt.savefig("20210428_Tan_B.png", dpi=250)

plt.show()

fig, (ax8, ax9) = plt.subplots(2)

ax8.plot(time_rtn_vel_local, B_nor_interp)
ax8.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax8.set_title('Interpolated Normal B-Field Magnitude (n='+str(n_samples)+')', size = 14)
ax8.set_xlabel('Hours after ' + T, size = 10)
ax8.set_ylabel('|Br| (nT)', size = 14)

ax9.plot(time_rtn_vel_local, nvel)
ax9.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax9.set_title('Normal-Component Solar Wind Velocity ', size = 14)
ax9.set_xlabel('Hours after '+ T, size = 10)
ax9.set_ylabel('Vr (km/s)', size = 14)

plt.tight_layout()
#plt.savefig("20210428_Nor_B.png", dpi=250)

plt.show()

# ------------------------------------------------------------------------------------------

### Correlation Coefficient ### -------------------------------------------

# caclulate magnetic field in velocity units
# plug into corellation coefficient formula from
# https://www.aanda.org/articles/aa/full_html/2020/01/aa37064-19/aa37064-19.html

# Get Velocity-Density Data
atemp3 = get_data("np_moment")
epoch_dens = atemp3[0]
dens = atemp3[1]

### Redifining Velocity Density ### --------------------------------------------
dens=dens[NAN3]
epoch_dens=epoch_dens[NAN3]

# Recombine Interpolated B-Field Data into useable 2D Array
B_rtn_interp = np.zeros(shape = (len(B_rad_interp), 3))
B_rtn_interp[:,0] = B_rad_interp
B_rtn_interp[:,1] = B_tan_interp
B_rtn_interp[:,2] = B_nor_interp


# Calculate the magnitude of the interpolated 2d Array, returns a 1d array
B_rtn_interp_mag = np.sqrt(B_rtn_interp[:,0]**2 + B_rtn_interp[:,1]**2 + B_rtn_interp[:,2]**2)

# Constants
mu_0 = (4*(np.pi)) * (10**(-7)) # H / m
proton_mass = 1.67 * (10**(-27))

#Interpolating Density to B#
#dens_interp=np.interp(time_rtn_vel_local,epoch_dens,dens)


# Alfven Velocity
Va = B_rtn_interp_mag * 1e-9 / np.sqrt(mu_0 * proton_mass * dens * 1e6)
Va = Va / 1e3

# Magnetic Field in Velocity Units
B_r_vel = Va * (B_rad_interp / B_rtn_interp_mag)
B_r_vel=np.abs(B_r_vel)

B_t_vel = Va * (B_tan_interp / B_rtn_interp_mag)
B_t_vel=np.abs(B_r_vel)

B_n_vel = Va * (B_nor_interp / B_rtn_interp_mag)
B_n_vel=np.abs(B_r_vel)

#### Cube Time ###--------------------------------------------------------------
B_rtn_vel_cube=np.zeros(shape=(len(B_r_vel),3))
B_rtn_vel_cube[:,0] = B_r_vel
B_rtn_vel_cube[:,1] = B_t_vel
B_rtn_vel_cube[:,2] = B_n_vel




vel_RTN_cube=np.zeros(shape=(len(rvel),3))
vel_RTN_cube[:,0] = rvel
vel_RTN_cube[:,1] = tvel
vel_RTN_cube[:,2] = nvel



vel_RTN_cube_T=vel_RTN_cube.reshape(3,len(rvel))
###Magnitude ### --------------------------------------------------------------

velmagcube=mag3d(vel_RTN_cube)
Bmagcube=mag3d(B_rtn_vel_cube)


VelDotB=dotproduct(vel_RTN_cube, B_rtn_vel_cube)




# Calculation of Correlation Coefficient
# This is a 1D list of calculated correlation coefficients at each data point
# The correlation is between radial solar wind velocity and the magnetic field in velocity units
sigma = 2 * VelDotB / ((Bmagcube**2 + velmagcube**2))

NaN = np.where(np.isnan(sigma)==False)[0]

sigmaNaN = sigma[NaN]

AvergeSigmaValue= np.mean(np.abs(sigmaNaN))
print('')
print('Mean Distance from Sun (Solar Radii):', av_dist)
print('Mean Sigma Value for',T1[0:10], 'is:', AvergeSigmaValue)
print('Lower Sigma Quartile: ', np.quantile(sigmaNaN, .25))
print('Upper Sigma Quartile: ', np.quantile(sigmaNaN, .75))
print('Samples: ', r_samples)
print('')

plt.scatter(time_rtn_vel_local, sigma, s = 12)
plt.xlabel('Hours after '+T, size = 14)
plt.ylabel('\u03C3', size = 17)
plt.title('B-Field-Velocity Correlation Coefficient', size = 14)
#plt.savefig("20210428_sigma.png", dpi=250)

plt.show()

# -----------------------------------------------------------------------------------------
# Generate Stack plot of all Key Figure Data
# Radial Interpolated B-Field Magnitude
# Radial Velocity Magnitude
# Proton Density Data (used for Corellation)
# Correlation Coefficient
# Radial PSP Distance

fig, (ax10, ax11, ax12, ax13, ax14) = plt.subplots(5, figsize = (10,15), sharex=True)

# B-Field
ax10.plot(time_rtn_vel_local, B_rad_interp)
ax10.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
ax10.set_title('2021-04-28 Parker Solar Probe (SPI) (n='+str(r_samples)+')', size = 23)
#ax10.set_xlabel('Hours after ' + T)
ax10.set_ylabel('Br (nT)', size = 23)

# Velocity

ax11.plot(time_rtn_vel_local, rvel)
ax11.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
#ax11.set_title('Radial-Component Solar Wind Velocity ')
#ax11.set_xlabel('Hours after '+ T)
ax11.set_ylabel('Vr (km/s)', size = 23)

# Proton Density

ax12.plot(time_rtn_vel_local, dens)
ax12.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
#ax12.set_title('Proton Density (1/cm^3)')
#ax12.set_xlabel('Hours after '+ T)
ax12.set_ylabel('Np (1/cm$^3$)', size = 23)

# Corellation Coefficient

ax13.scatter(time_rtn_vel_local, sigma, s = 12)
ax13.set_xlim(np.min(time_rtn_vel_local), np.max(time_rtn_vel_local))
#ax13.set_xlabel('Hours after '+T)
ax13.set_ylabel('\u03C3', size = 23)
#ax13.set_title('B-Field-Velocity Correlation Coefficient')

# PSP Distance

ax14.plot(time_sc_pos_local, sc_pos * s_rad, label = 'Spacecraft Position from Sun')
#ax14.set_title('Position of PSP Relative to Sun')
ax14.set_xlabel('Hours after ' + T, size = 23)
ax14.set_ylabel('R (R${_s}$)', size = 23)

plt.tight_layout()
#plt.savefig("20210428_combined.png", dpi=300)

plt.show()
