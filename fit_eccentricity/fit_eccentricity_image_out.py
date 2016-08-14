#!/usr/bin/env python2.7

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from readimage import *
from natconst import *
from matplotlib import cm
import math
from math import sin, cos


#debug = True
debug = False

#### THE Spectral Index got from ALMA
spectral_index = 2.72256367113

## 
# frequencies
mu_mod = cc*10/1.0
mu_obs = 343500000000.0
#av4    343500000000.0
#spw0	337244575000.0
#spw1	339182075000.0
#spw2	347744575000.0
#spw3	349744575000.0

# wavelength in micron
u_wav = 870.0
#av4 870.0
#0 888.947
#1 883.869
#2 862.105
#3 857.175
v_wav = u_wav

######### radius wanted
gap_au = 95.0 # AU
gap_width = 12.0 # AU

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

image = readImage()
dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out
ImageJyppix_scaled = image.imageJyppix[:,:,n_image]/dpc/dpc
print "Total flux mJy (original):"
print np.sum(ImageJyppix_scaled, axis=None)*1000.0, "\n"

# scale the flux using Fv2/Fv1 = (v2/v1)**spectral_index
ImageJyppix_scaled = ImageJyppix_scaled*(mu_obs/mu_mod)**spectral_index 
print "Total flux scaled mJy (Fv2/Fv2 = (mu2/mu1)**spectral_index):"
print np.sum(ImageJyppix_scaled, axis=None)*1000.0
ImageJyTotal = np.sum(ImageJyppix_scaled, axis=None)

m,n = ImageJyppix_scaled.shape
print "dimension check (m,n):", m, n, "\n"


############################################ 
# write the flux in Jy/pixel for debug
############################################ 
if debug:
    image_output = np.zeros(m*n)
    for i in range(m):
        for j in range(n):
            image_output[n*i+j]=ImageJyppix_scaled[i,j]
    np.savetxt('ImageJyppix_scaled', np.transpose(image_output))



##########################################
# set the grids in pixel
########################################## 

#################################
# NOTE: in this case, nx=ny=m=n

# assuming the same number of grids in x,y axes
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
fft_dx = pix_x[1] - pix_x[0]
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)
fft_dy = pix_y[1] - pix_y[0]


 
############################################ 
############ info printting
############################################ 

print "pixel-size (x) in cm :"
print "     ", image.sizepix_x
print "pixel-size (y) in cm :"
print "     ", image.sizepix_y
print "distance from the source (cm), PC x dpc(nPCs) :"
print "     ", dpc*3.0856775e18
print "radian in x :"
print "     ", image.sizepix_x/dpc/3.0856775e18
print "radian in y :"
print "     ", image.sizepix_y/dpc/3.0856775e18, "\n"


############################################ 
####### set the u_max,v_max for FFT
############################################ 

u_max = u_wav/1e6/(image.sizepix_x/dpc/3.0856775e18)
v_max = v_wav/1e6/(image.sizepix_y/dpc/3.0856775e18)
print "u_Max (lamda_obs / radian_x) :"
print "     ", u_wav/1e6/(image.sizepix_x/dpc/3.0856775e18)
print "y_Max (lamda_obs / radian_y) :"
print "     ", v_wav/1e6/(image.sizepix_y/dpc/3.0856775e18), "\n"


############################################ 
######### Shift the image to solve the [-1,-1] 
#########  centering problem in FFT
############################################ 

i_in = (gap_au-0.5*gap_width)*au/image.sizepix_x
print i_in, "i_in"
i_ou = (gap_au+0.5*gap_width)*au/image.sizepix_x
print i_ou, "i_ou"
#

# screen out the point in the ring at [i_in, i_out]
fin0_x = []
fin0_y = []
fin0_vis = []

i_2_au = image.sizepix_x/au
for i in range(1,m):
    for j in range(1,n):
        if (((i-501)**2+(j-501)**2)**0.5 >= i_in) and (((i-501)**2+(j-501)**2)**0.5 < i_ou):
            fin0_x.append((i-501)*i_2_au)
            fin0_y.append((j-501)*i_2_au)
            fin0_vis.append(ImageJyppix_scaled[i,j])


fin_x = np.asarray(fin0_x)
fin_y = np.asarray(fin0_y)
fin_vis = np.asarray(fin0_vis)

np.savetxt('x_y_vis', np.transpose([fin_x,fin_y,fin_vis]))


n_fin = fin_x.shape[0]
vis_r_d_phi = np.zeros((n_fin,4), dtype=np.float64)
vis_r_d_phi[:,0] = fin_vis
vis_r_d_phi[:,1] = fin_x
vis_r_d_phi[:,2] = fin_y

# define the angle calculated function
def angle_between_points( p0, p1, p2 ):
    a = (p1[0]-p0[0])**2 + (p1[1]-p0[1])**2
    b = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2
    c = (p2[0]-p0[0])**2 + (p2[1]-p0[1])**2
    #print a, b, c
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


center_ra = 0.0
center_de = 0.0
p_ori = [center_ra,center_de]
p_ref = [center_ra,fin_y.max()]
#p_one = [vis_r_d[:,1].min(),center_de]
#p_one = [vis_r_d[:,1].max(),center_de]
ang_x = []
ang_y = []
ang_vis = []
ang_ang = []

for i in range(n_fin):
    p_one = [vis_r_d_phi[i,1],vis_r_d_phi[i,2]]
    once = once_between_points(p_ref,p_ori,p_one) 
    if once:
        angle = angle_between_points(p_ref,p_ori,p_one) 
        if p_one[0] < center_ra:
            angle = 360 - angle
        ang_x.append(vis_r_d_phi[i,1])
        ang_y.append(vis_r_d_phi[i,2])
        ang_vis.append(vis_r_d_phi[i,0])
        ang_ang.append(angle)

#sort the array according to azimuthal angle
ang_x = np.asarray(ang_x)
ang_y = np.asarray(ang_y)
ang_vis = np.asarray(ang_vis)
ang_ang = np.asarray(ang_ang)
n_ang = ang_x.shape[0]
vis_r_d_phi_Fin = np.zeros((n_ang,4), dtype=np.float64)

vis_r_d_phi_Fin[:,0] = ang_vis
vis_r_d_phi_Fin[:,1] = ang_x
vis_r_d_phi_Fin[:,2] = ang_y
vis_r_d_phi_Fin[:,3] = ang_ang
vis_r_d_phi_Fin = vis_r_d_phi_Fin[vis_r_d_phi_Fin[:,3].argsort()] 

#save vis_r_d_phi
np.savetxt('JyPixel_Ra_Dec_ang.out', np.transpose([vis_r_d_phi_Fin[:,0],vis_r_d_phi_Fin[:,1],vis_r_d_phi_Fin[:,2],vis_r_d_phi_Fin[:,3]]))


############# Find the Minimium value at a phi_angle

# defien phi grid
n_phi = 721

# creat grid      
phi_min = 0.0
phi_max = 360.0
dphi  = (phi_max-phi_min)/(n_phi-1)
phi = zeros(n_phi)
for i in range(n_phi):
    phi[i] = phi_min+i*dphi


vis_r_d_phi_Nphi = np.zeros((n_phi-1,4), dtype=np.float64)


for i in range(n_phi-1):

    ## print a progress bar
    if(i%max(1,int(0.1*(n_phi-1)))==0):
        print "%s%s" % (100*i/(n_phi-1), "% done ...")

    j_l = 0
    for j in range(n_ang):
        if (vis_r_d_phi_Fin[j,3]>=phi[i]):
            j_l = j
            break
    j_u = 0
    for j in range(n_ang): 
        if (vis_r_d_phi_Fin[j,3]>=phi[i+1]):
            j_u = j-1
            break
        elif (j>=(n_ang-1)):
            j_u = j
            break
    if debug:
        print j_l, j_u

    vis_min = vis_r_d_phi_Fin[j_l,0]
    n_vis_min = j_l
    for j in range(j_l,j_u+1):
        if (vis_r_d_phi_Fin[j,0] > vis_min):
            vis_min = vis_r_d_phi_Fin[j,0]
            n_vis_min = j

    vis_r_d_phi_Nphi[i,0] = vis_min
    vis_r_d_phi_Nphi[i,1] = vis_r_d_phi_Fin[n_vis_min,1]
    vis_r_d_phi_Nphi[i,2] = vis_r_d_phi_Fin[n_vis_min,2]
    vis_r_d_phi_Nphi[i,3] = vis_r_d_phi_Fin[n_vis_min,3]

np.savetxt('vis_r_d_phi.Min', np.transpose([vis_r_d_phi_Nphi[:,0],vis_r_d_phi_Nphi[:,1],vis_r_d_phi_Nphi[:,2],vis_r_d_phi_Nphi[:,3]]))
