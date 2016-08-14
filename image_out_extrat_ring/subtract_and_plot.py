#!/usr/bin/env python2

import os
import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # to avoid window

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm

from readimage import *
from natconst import *
from matplotlib import cm
import math
from math import sin, cos


#debug = True
debug = False

colorlog = False

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

image = readImage('./withline_free_PD/image.out')
dpc = 140.0  # distance in pc

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

im_x_au = (pix_x)*image.sizepix_x/au
im_y_au = (pix_y)*image.sizepix_y/au

######################
print "INFO:","\n"
print "image_dust.nwav"
print image.nwav
print "image.wav"
print image.wav


#######################################################
#####  Remove the dust continuum at each wav for gas image
#########################################


nwav = image.nwav
for i in range(image.nwav):
    imagedust = "%s%s" % ("image.out_lambda_", i)
    image_dust = readImage(imagedust)
    image.image[:,:,i] = image.image[:,:,i] - image_dust.image[:,:,0]
    image.imageJyppix[:,:,i] = image.imageJyppix[:,:,i] - image_dust.imageJyppix[:,:,0]



#######################################################
#####  Create average image
#########################################

image_ave = readImage('./image.out_lambda_0')

image_ave.image[:,:,0] = image_ave.image[:,:,0]*0.0
image_ave.imageJyppix[:,:,0] = image_ave.imageJyppix[:,:,0]*0.0
for i in range(image.nwav):
    image_ave.image[:,:,0] = image_ave.image[:,:,0] + image.image[:,:,i]
    image_ave.imageJyppix[:,:,0] = image_ave.imageJyppix[:,:,0] + image.imageJyppix[:,:,i]


#######################################################
#####  Plotting

plt.pcolormesh(im_x_au, im_y_au, image_ave.imageJyppix[:,:,0]/dpc/dpc, cmap='RdBu')
plt.xlabel('AU')
plt.ylabel('AU')
plt.ylim(-180,180)
plt.xlim(-180,180)
cbar1=plt.colorbar()
cbar1.set_label("Janskey/pixel")
plt.title("Flux density")
plt.savefig('flux_density_ave.png')
plt.clf()


#######################################################
#####  Make some Plot 
#########################################
#plt, axes = plt.subplots(nrows=1, ncols=3, figsize=(13, 6), dpi=80, sharex=True, sharey=True)
plt, axes = plt.subplots(nrows=5, ncols=8, figsize=(13, 8), dpi=80, sharex=True, sharey=True)


axisadd = 10
v_min   = 0.0
v_max   = image.imageJyppix[:,:,:].max()/dpc/dpc


n_row = 5
n_col = image.nwav/n_row
for i in range(n_row):
    for j in range(n_col):
        print i,j
        if colorlog:
            im=axes[i,j].pcolormesh(im_x_au, im_y_au, image.imageJyppix[:,:,n_col*i+j]/dpc/dpc, cmap='jet', norm=SymLogNorm(linthresh=0.03, linscale=0.03, vmin=v_min, vmax=v_max))
        else:
            im=axes[i,j].pcolormesh(im_x_au, im_y_au, image.imageJyppix[:,:,n_col*i+j]/dpc/dpc, cmap='jet', vmin=v_min, vmax=v_max)
        axes[i,j].axis([im_x_au.min()-axisadd, im_x_au.max()+axisadd, im_y_au.min()-axisadd, im_y_au.max()+axisadd])

#plt.tight_layout()
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, hspace=0.0, wspace=0.0)
cax = plt.add_axes([0.8, 0.1, 0.02, 0.8])
plt.colorbar(im, cax=cax, label=' Janskey/pixel ')
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
#plt.colorbar(cax=cax)

#plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('gas_lines_radmc.png')

plt.clf()





