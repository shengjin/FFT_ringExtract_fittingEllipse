#!/usr/bin/env python2
#import sys
import numpy as np
import math
#import os

import matplotlib as mpl
mpl.use('Agg')   # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

debug = False

filename = "JyPixel_Ra_Dec.out"

vis_r_d = np.genfromtxt(filename, skip_header=0, dtype=float)

center_ra = (vis_r_d[:,1].max() + vis_r_d[:,1].min())*0.5
center_de = (vis_r_d[:,2].max() + vis_r_d[:,2].min())*0.5
ref_de_top = vis_r_d[:,2].max()
print "ra.max,ra.min", vis_r_d[:,1].max(), vis_r_d[:,1].min()
print "de.max,de.min", vis_r_d[:,2].max(), vis_r_d[:,2].min()
print "center_ra,center_de", center_ra, center_de
print "ref_de_top", ref_de_top, "\n"

n = vis_r_d.shape[0]
vis_r_d_phi = np.zeros((n,4), dtype=np.float64)
vis_r_d_phi[:,:-1] = vis_r_d

# define the angle calculated function
def angle_between_points( p0, p1, p2 ):
    a = (p1[0]-p0[0])**2 + (p1[1]-p0[1])**2
    b = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2
    c = (p2[0]-p0[0])**2 + (p2[1]-p0[1])**2
    #print a, b, c
    angle = math.acos( (a+b-c) / math.sqrt(4*a*b) ) * 180/math.pi
    return angle

def once_between_points( p0, p1, p2 ):
    a = (p1[0]-p0[0])**2 + (p1[1]-p0[1])**2
    b = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2
    c = (p2[0]-p0[0])**2 + (p2[1]-p0[1])**2
    #print a, b, c
    value = (a+b-c) / math.sqrt(4*a*b) 
    if (abs(value) < 1.0):
        once = True
    else:
        once = False
        if debug:
            print a, b, c
    return once

p_ori = [center_ra,center_de]
p_ref = [center_ra,ref_de_top]
#p_one = [vis_r_d[:,1].min(),center_de]
#p_one = [vis_r_d[:,1].max(),center_de]
for i in range(n_fin):
    p_one = [vis_r_d_phi[i,1],vis_r_d_phi[i,2]]
    once = once_between_points(p_ref,p_ori,p_one)
    if once:
        angle = angle_between_points(p_ref,p_ori,p_one) 
        if p_one[0] < center_ra:
            angle = 360 - angle
        vis_r_d_phi[i,3] = angle


#sort the array according to azimuthal angle
vis_r_d_phi = vis_r_d_phi[vis_r_d_phi[:,3].argsort()] 

#save vis_r_d_phi
np.savetxt('JyPixel_Ra_Dec_ang.out', np.transpose([vis_r_d_phi[:,0],vis_r_d_phi[:,1],vis_r_d_phi[:,2],vis_r_d_phi[:,3]]))


n_phi = 720
vis_r_d_phi_Fin = np.zeros((n_phi,4), dtype=np.float64)
vis_r_d_phi_Fin[:,3] = np.linspace(0,360,n_phi, endpoint=False)


def LinIntP(x,x1,y1,x2,y2):
    """
    makes linear interpolation for f(x) at the point x, if (x1,f(x1)=y1) and (x2,f(x2)=y2)
    using Lagrange's formula
    """
    return ((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2



#for i in range(2):
for i in range(n_phi):

    ## print a progress bar
    if(i%max(1,int(0.1*n_phi))==0):
        print "%s%s" % (100*i/n_phi, "% done ...")

    phi_tmp = vis_r_d_phi_Fin[i,3]
    ##################
    # find phi_up,phi_low
    for ii in range(n):
        if vis_r_d_phi[ii,3] >= phi_tmp:
            phi_up = ii
            break
        elif (ii == n-1) and vis_r_d_phi[ii,3] < phi_tmp:
            phi_up = 0
            print"(ii = n_phi-1) and vis_r_d_phi[ii,3] < phi_tmp:phi_up", phi_up
    if phi_up == 0:
        phi_low = n - 1
        print "phi_up==0:phi_low", phi_low
    else:
        phi_low = phi_up - 1
    if debug:
        print "(i,phi_tmp),(vis_r_d_phi[phi_up,3],phi_up),(vis_r_d_phi[phi_low,3],phi_low)"
        print i,phi_tmp,vis_r_d_phi[phi_up,3],phi_up,vis_r_d_phi[phi_low,3],phi_low
    ################
    # do interpolation
    vis_r_d_phi_Fin[i,1] = LinIntP(phi_tmp, vis_r_d_phi[phi_low,3], vis_r_d_phi[phi_low,1], vis_r_d_phi[phi_up,3], vis_r_d_phi[phi_up,1])
    vis_r_d_phi_Fin[i,2] = LinIntP(phi_tmp, vis_r_d_phi[phi_low,3], vis_r_d_phi[phi_low,2], vis_r_d_phi[phi_up,3], vis_r_d_phi[phi_up,2])
    vis_r_d_phi_Fin[i,0] = LinIntP(phi_tmp, vis_r_d_phi[phi_low,3], vis_r_d_phi[phi_low,0], vis_r_d_phi[phi_up,3], vis_r_d_phi[phi_up,0])


np.savetxt('JyPixel_Ra_Dec_ang_Fin.out', np.transpose([vis_r_d_phi_Fin[:,0],vis_r_d_phi_Fin[:,1],vis_r_d_phi_Fin[:,2],vis_r_d_phi_Fin[:,3]]))

# plot azimuthal flux density
plt.xlabel('Azimuthal angle')
plt.ylabel('Flux density [Jy/pixel]')
plt.plot(vis_r_d_phi_Fin[:,3], vis_r_d_phi_Fin[:,0])
plt.savefig('flux_phi.png')
plt.clf



#### FFT
dt = vis_r_d_phi_Fin[1,3] - vis_r_d_phi_Fin[0,3]
print "in FFT, dt=: ", dt

vis = vis_r_d_phi_Fin[:,0]
FFT_x = np.fft.fft(vis)

Freq = np.fft.fftfreq(FFT_x.shape[0], dt)

FFT_x = np.fft.fftshift(FFT_x)   # Shift zero freq to center
Freq = np.fft.fftshift(Freq)   # Shift zero freq to center

np.savetxt('FFT.out', np.transpose([abs(FFT_x), FFT_x.real, FFT_x.imag, Freq]))

# plot FFT
# plot
plt.xlabel('Frequency')
plt.ylabel('abs(FFT)')
plt.xlim(-0.1,0.1)
plt.plot(Freq, abs(FFT_x))
plt.savefig('fft_abs.png')

quit()


