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


# inclination of the disk
#incli = 51.0
incli = 224.0

# position angle
# use it to scale the center Radius VS T_B
#posang = 61.0
posang = 222.0

# rotation angle used in the center ellipse to
rotang = posang
#rotang = 90-posang

# radius_wanted
n_rad = 15
inneravoid = 3
gap_width = 0.9 # AU

dpc = 122.0  # distance in pc

# LkCa 15
#gaus_ma = 6.323864890469E-05    # maj axis of guassian in deg
#gaus_mi = 4.604819748137E-05    # min
#gaus_pa = 2.383443450928E+01    # PA of maj axis of Guass measure East from North

# HL Tau
gaus_ma = 8.367259158856E-06*3.5
gaus_mi = 5.335457002123E-06*3.5
gaus_pa = -87.3

gaus_ma_arcsec = 0.65
gaus_mi_arcsec = 0.42
if gaus_pa > 0:
    gaus_pa_ctclk  = 180.0-gaus_pa  # PA in imConv
else:
    gaus_pa_ctclk  = -gaus_pa  # PA in imConv

fwhm = [gaus_ma_arcsec, gaus_mi_arcsec]

#debug = True
debug = False

AU = 1.496e13

# the n_th image for azimuthal extraction
n_image = 0

#########################
#### print the gaussian
print "gaus_pa_ctclk"
print gaus_pa_ctclk
print "gaus_ma_arcsec, gaus_mi_arcsec"
print gaus_ma_arcsec, gaus_mi_arcsec
gaus_ma_rad = gaus_ma*math.pi/180.0
gaus_mi_rad = gaus_mi*math.pi/180.0
gaus_ma_AU = math.sin(gaus_ma_rad)*dpc*pc/AU
gaus_mi_AU = math.sin(gaus_mi_rad)*dpc*pc/AU
print "gaus_ma_AU, gaus_mi_AU"
print gaus_ma_AU, gaus_mi_AU


#debug = True
debug = False

colorlog = False

###################### ###################### 
########### Read image and scale the flux
###################### ###################### 

#image = readImage('./image.out')
image = readImage('./image.out')
dpc = 122.0  # distance in pc
print AU/pc/dpc
au2arcsec = AU/pc/dpc*206265
print au2arcsec

#################################
# NOTE: in this case, nx=ny=m=n
# assuming the same number of grids in x,y axes
nx = image.nx
NN = nx
pix_x = np.linspace(-NN/2,NN/2,nx, endpoint=False)
ny = image.ny
pix_y = np.linspace(-NN/2,NN/2,ny, endpoint=False)

im_x_au = (pix_x)*image.sizepix_x/au
im_y_au = (pix_y)*image.sizepix_y/au

arcsec_x = im_x_au*au2arcsec
arcsec_y = im_y_au*au2arcsec

np.savetxt('arcsec_x', np.transpose(arcsec_x))

######################
print "INFO:","\n"
print "image_dust.nwav"
print image.nwav
print "image.wav"
print image.wav




#######################################################
#####  Create average image
#########################################

image_ave = readImage('./image_dust.out')

image_ave.image[:,:,0] = image_ave.image[:,:,0]*0.0
image_ave.imageJyppix[:,:,0] = image_ave.imageJyppix[:,:,0]*0.0
for i in range(image.nwav):
    image_ave.image[:,:,0] = image_ave.image[:,:,0] + image.image[:,:,i]
    image_ave.imageJyppix[:,:,0] = image_ave.imageJyppix[:,:,0] + image.imageJyppix[:,:,i]



#######################################################
#####  Plotting

print " Plotting... ", "\n"

plt.pcolormesh(im_x_au, im_y_au, image_ave.imageJyppix[:,:,0]/dpc/dpc, cmap='RdBu')
plt.xlabel('AU')
plt.ylabel('AU')
#plt.ylim(-180,180)
#plt.xlim(-180,180)
cbar1=plt.colorbar()
cbar1.set_label("Jansky/pixel")
plt.title("Flux density")
plt.savefig('flux_density_ave.png')
plt.clf()


#######################################################
#####  Make some Plot 
#########################################
n_row = 2
n_col = image.nwav/n_row

#plt, axes = plt.subplots(nrows=1, ncols=3, figsize=(13, 6), dpi=80, sharex=True, sharey=True)
plt, axes = plt.subplots(n_row, n_col, figsize=(16, 4), dpi=80, sharex=True, sharey=True)


axisadd = 0
v_min   = 0.0
v_max   = image.imageJyppix[:,:,:].max()/dpc/dpc*1000


for i in range(n_row):
    for j in range(n_col):
        if colorlog:
            im=axes[i,j].pcolormesh(arcsec_x, arcsec_y, image.imageJyppix[:,:,n_col*i+j]/dpc/dpc*1000, cmap='jet', norm=SymLogNorm(linthresh=0.03, linscale=0.03, vmin=v_min, vmax=v_max))
        else:
            im=axes[i,j].pcolormesh(arcsec_x, arcsec_y, image.imageJyppix[:,:,n_col*i+j]/dpc/dpc*1000, cmap='jet', vmin=v_min, vmax=v_max)
        axes[i,j].axis([arcsec_x.min()-axisadd, arcsec_x.max()+axisadd, arcsec_y.min()-axisadd, arcsec_y.max()+axisadd])
        axes[i,j].tick_params(axis='both', which='major', labelsize=8)
        axes[i,j].set_yticks( [-3, 0, 3] )
        axes[i,j].set_xticks( [-3, 0, 3] )
        #axes[i,j].set_xticks( [-150, -100, -50, 0, 50, 100, 150] )
        axes[i,j].set_xticklabels(['-3', '0', '3'],rotation=45)
        axes[i,j].set_yticklabels(['-3', '0', '3'])
        if j==0:
            axes[i,j].set_ylabel('arcsec')
        axes[i,j].set_xlabel('arcsec')

#plt.tight_layout()
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, hspace=0.0, wspace=0.0)
cax = plt.add_axes([0.82, 0.1, 0.02, 0.8])
plt.colorbar(im, cax=cax, label='m Jansky/pixel')
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
#plt.colorbar(cax=cax)

#plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('gas_lines_radmc.png')

plt.clf()


#######################################################
#####  Convolve image with gaussian
#########################################

print "Convolve the image wiht Guassian"
image_conv = image.imConv(fwhm, gaus_pa_ctclk, dpc)
image_ave_conv = image_ave.imConv(fwhm, gaus_pa_ctclk, dpc)



############################################ 
# write the flux in Jy/pixel for debug
############################################ 
if debug:
    image_output = np.zeros(nx*nx)
    for i in range(nx):
        for j in range(nx):
            image_output[nx*i+j]=image_ave.image[i,j,0]
    np.savetxt('ImageJyppix_scaled', np.transpose(image_output))


##########################################
# define azimuthal extract function
#  could be ellipse or circle
##########################################

def azimuthal_Jy_avg(image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix):
    # input: image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix, nx
    # output: fin_vis
    ##################################
    # screen out the point in the ring at [i_in, i_out]

    au = 1.496e13

    fin0_x = []
    fin0_y = []
    fin0_vis = []
    
    # parameters for the center ellipse
    gap_min = gap_au - gap_width*0.5
    gap_max = gap_au + gap_width*0.5
    
    # assuming sizepix_x = sizepix_y
    e_a_grid_max = gap_max*au/sizepix  # long semi-axis
    e_a_grid_min = gap_min*au/sizepix  # long semi-axis
    
    inclination = math.radians(incli)
    e_b_grid_max = e_a_grid_max * cos(inclination)     # short semi-axis
    e_b_grid_min = e_a_grid_min * cos(inclination)     # short semi-axis
    rotang = math.radians(rotang)
    m = image.shape[0]
    # convert integer to float in order to make sure we find
    #    the center of the image.
    # image.out: 1) 20 grids
    #               python array: 0, 1, ..., 19
    #               center is at point 10
    #            2) 19 grids
    #               python array: 0, 1, ..., 18
    #               center is at point 9.5
    i_2_au = sizepix/au
    for ii in range(m):
        i=float(ii)
        for jj in range(m):
            j=float(jj)
            if ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_max**2 + ((i-e_h)*sin(rotang)-( j-e_k)*cos(rotang))**2/e_b_grid_max**2 <= 1.0**2) and ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_min**2 + ((i-e_h)*sin(rotang)-(j-e_k)*cos(rotang))**2/e_b_grid_min**2 > 1.0**2) :
                fin0_x.append(i)
                fin0_y.append(j)
                fin0_vis.append(image[ii,jj])
                #print nhalfpix, i, j
    fin_x = np.asarray(fin0_x)
    fin_y = np.asarray(fin0_y)
    fin_vis = np.asarray(fin0_vis)
    n_fin = fin_x.shape[0]
    if n_fin > 0: 
        np.savetxt('x_y_vis', np.transpose([fin_x,fin_y,fin_vis]))
    #total=np.sum(fin_vis), "\n"
    #avg = total/n_fin
    return fin_vis


##########################################
##########################################

#### define the ImagePlot nx*nx array to extract and plot azimuthal ring
ImagePlot = image_ave_conv.imageJyppix[:,:,n_image]


sizepix = image_ave_conv.sizepix_x
print "image_ave_conv.sizepix_x/au,image.sizepix_x/au"
print image_ave_conv.sizepix_x/au,image.sizepix_x/au
print nx, image_ave_conv.nx

nfloat = float(nx)
nhalfpix = nfloat/2
e_h = nhalfpix
e_k = nhalfpix

r_Jy = np.zeros([n_rad,2], dtype=float64)

for i in range(n_rad):
    gap_au = float(i)+inneravoid
    r_Jy[i,0] = gap_au
    avg = azimuthal_Jy_avg(ImagePlot, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix)
    n_num = avg.shape[0]
    if n_num == 0:
        r_Jy[i,1] = 0
        print "WARNNING: no points found around ", r_Jy[i,0], " AU!!"
    else:
        f_n_num = float(n_num)
        #print avg
        total = np.sum(avg)
        #print total
        r_Jy[i,1] = total/f_n_num
        print r_Jy[i,0], r_Jy[i,1]
    
    
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6)) 
#mpl.rcParams['xtick.labelsize'] = 22
#mpl.rcParams['ytick.labelsize'] = 24
plt.tick_params(labelsize=20)
plt.xlabel('AU',fontsize=20)
plt.ylabel("m Jansky / beam",fontsize=20)
#plt.ylim(-180,180)
plt.xlim(3,160)
plt.plot(r_Jy[:,0], r_Jy[:,1]*1000)
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
plt.savefig('azimuthalavg_fluxJy_convolved.png')
plt.clf()



#######################################################
#####  Make some Plot 
#########################################

plt.figure()
plt.xlabel('AU')
plt.ylabel('AU')
plt.ylim(-180,180)
plt.xlim(-180,180)
plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000,cmap='RdBu')
#plt.pcolormesh(im_x_au, im_y_au, ImagePlot*1000,cmap='RdBu')
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
cbar1=plt.colorbar()
cbar1.set_label("m Jansky / beam")
plt.title("Flux density")
#plt.tick_params(axis='both', which='major', labelsize=16)
plt.savefig('flux_density_convolved.png')
plt.clf()


#######################################################
#########################################
n_row = 2
n_col = image.nwav/n_row

#plt, axes = plt.subplots(nrows=1, ncols=3, figsize=(13, 6), dpi=80, sharex=True, sharey=True)
plt, axes = plt.subplots(n_row, n_col, figsize=(18, 4), dpi=80, sharex=True, sharey=True)


axisadd = 0
v_min   = 0.0
v_max   = image_conv.imageJyppix[:,:,:].max()*1


for i in range(n_row):
    for j in range(n_col):
        if colorlog:
            im=axes[i,j].pcolormesh(arcsec_x, arcsec_x, image_conv.imageJyppix[:,:,n_col*i+j]*1.0, cmap='jet', norm=SymLogNorm(linthresh=0.03, linscale=0.03, vmin=v_min, vmax=v_max))
        else:
            im=axes[i,j].pcolormesh(arcsec_x, arcsec_y, image_conv.imageJyppix[:,:,n_col*i+j]*1.0, cmap='jet', vmin=v_min, vmax=v_max)
        axes[i,j].tick_params(axis='both', which='major', labelsize=8)
        #tick.label.set_fontsize(14) 
        axes[i,j].axis([arcsec_x.min()-axisadd, arcsec_x.max()+axisadd, arcsec_y.min()-axisadd, arcsec_y.max()+axisadd])
        axes[i,j].set_yticks( [-3, 0, 3] )
        axes[i,j].set_xticks( [-3, 0, 3] )
        #axes[i,j].set_xticks( [-150, -100, -50, 0, 50, 100, 150] )
        axes[i,j].set_xticklabels(['-3', '0', '3'],rotation=45)
        axes[i,j].set_yticklabels(['-3', '0', '3'])
        if j==0:
            axes[i,j].set_ylabel('arcsec')
        axes[i,j].set_xlabel('arcsec')
        #for tick in axes[i,j].get_xticklabels():
        #    tick.set_rotation(45)

#plt.tight_layout()
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, hspace=0.0, wspace=0.0)
cax = plt.add_axes([0.82, 0.1, 0.02, 0.8])
plt.colorbar(im, cax=cax, label='Jansky/beam')
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
#plt.colorbar(cax=cax)

#plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('gas_lines_radmc_convolved.png')

plt.clf()


