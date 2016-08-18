#!/usr/bin/env python2

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


dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out

gaus_ma = 6.323864890469E-05    # maj axis of guassian in deg
gaus_mi = 4.604819748137E-05    # min
gaus_pa = 2.383443450928E+01    # PA of maj axis of Guass measure East from North

gaus_ma_arcsec = gaus_ma*3600
gaus_mi_arcsec = gaus_mi*3600
gaus_pa_ctclk  = gaus_pa  # PA in imConv

fwhm = [gaus_ma_arcsec, gaus_mi_arcsec]

#debug = True
debug = False


AU = 1.496e13

#sizeAU = math.sin(pixel_in_rad)*dpc*pc/AU

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

image = readImage()

image_conv = image.imConv(fwhm, gaus_pa_ctclk, dpc)

ImagePlot = image_conv.imageJyppix[:,:,n_image]
print "Total flux mJy (original):"
print np.sum(ImagePlot, axis=None)*1000.0, "\n"


############################################ 
############ info printting
############################################ 

############################################ 
# write the flux in Jy/pixel for debug
############################################ 
if debug:
    image_output = np.zeros(m*n)
    for i in range(m):
        for j in range(n):
            image_output[n*i+j]=ImagePlot[i,j]
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

#####################
# set the grids in AU
# change
im_x_au = (pix_x+0.5)*image.sizepix_x/au
im_y_au = (pix_y+0.5)*image.sizepix_y/au




#######################################################
#####  Make some Plot 
#########################################

plt.xlabel('AU')
plt.ylabel('AU')
plt.ylim(-180,180)
plt.xlim(-180,180)
plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000,cmap='RdBu')
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
cbar1=plt.colorbar()
cbar1.set_label("m Janskey / beam")
plt.title("Flux density")
plt.savefig('flux_density.png')
plt.clf()






